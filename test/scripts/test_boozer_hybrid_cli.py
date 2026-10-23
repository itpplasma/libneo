"""Check HYBRID CLI selection and iotas tangent/adjoint serialization."""

import importlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import numpy as np
from netCDF4 import Dataset


SURFACES = np.array([0.15, 0.5, 0.85])
SPECTRUM = (
    "bmnc_b", "bmns_b", "numns_b", "numnc_b", "iota", "G", "I",
    "gradsmnc_b", "gradsmns_b", "psi_edge", "R",
)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def read_output(path):
    with Dataset(path) as data:
        data.set_auto_maskandscale(False)
        arrays = {name: np.array(var[:]) for name, var in data.variables.items()}
        construction = int(data.getncattr("construction"))
    for name, value in arrays.items():
        require(np.all(np.isfinite(value)), f"nonfinite output {name}")
    require(arrays["xm_b"].ndim == arrays["xn_b"].ndim == 1, "wrong mode shape")
    require(arrays["xm_b"].shape == arrays["xn_b"].shape, "mode shapes differ")
    require(np.array_equal(arrays["s"], SURFACES), "wrong output surfaces")
    for name in SPECTRUM:
        expected = (len(SURFACES), arrays["xm_b"].size) if name.endswith("_b") else (
            len(SURFACES),
        )
        require(arrays[name].shape == expected, f"wrong {name} shape")
    return arrays, construction


def write_input(path, construction=None, iotas=None, seed=None):
    with Dataset(path, "w", format="NETCDF3_CLASSIC") as data:
        data.createDimension("s", len(SURFACES))
        data.createVariable("s", "f8", ("s",))[:] = SURFACES
        for name, value in (("mboz", 8), ("nboz", 6), ("grid_refinement", 3)):
            data.createVariable(name, "i4")[...] = value
        if construction is not None:
            data.createVariable("construction", "i4")[...] = construction
        if iotas is not None:
            data.createDimension("radius", len(iotas))
            data.createVariable("iotas", "f8", ("radius",))[:] = iotas
        if seed is not None:
            data.createVariable("iota", "f8", ("s",))[:] = seed


def run_cli(binary, root, label, mode, wout, **inputs):
    source = root / f"{label}-input.nc"
    output = root / f"{label}-output.nc"
    write_input(source, **inputs)
    result = subprocess.run(
        [str(binary), mode, str(wout), str(source), str(output)],
        cwd=root, capture_output=True, text=True, check=False,
    )
    require(
        result.returncode == 0,
        f"{label} exited {result.returncode}: "
        f"{(result.stdout + result.stderr)[-4000:]}",
    )
    return read_output(output)


def copy_profile(original, destination, iotas=None, iotaf=None):
    shutil.copyfile(original, destination)
    with Dataset(destination, "r+") as data:
        data.set_auto_maskandscale(False)
        if iotas is not None:
            data.variables["iotas"][:] = iotas
        if iotaf is not None:
            data.variables["iotaf"][:] = iotaf


def check_selector(binary, original, root, ns):
    # Deliberately inconsistent profiles distinguish input dispatch, not equilibria.
    wout = root / "selector-wout.nc"
    copy_profile(original, wout, np.full(ns, 0.375), np.full(ns, 0.625))
    outputs = []
    errors = {}
    for label, selection, expected in (
        ("default", None, 0.375), ("hybrid", 1, 0.375), ("legacy", 2, 0.625),
    ):
        arrays, construction = run_cli(
            binary, root, label, "primal", wout, construction=selection,
        )
        require(construction == (selection or 1), f"wrong {label} construction")
        error = float(np.max(np.abs(arrays["iota"] - expected)))
        require(error <= 1e-11, f"{label} iota error {error:.3e}")
        errors[label] = error
        outputs.append(arrays)
    require(outputs[0].keys() == outputs[1].keys(), "default output fields differ")
    for name in outputs[0]:
        require(
            np.array_equal(outputs[0][name], outputs[1][name]),
            f"default and explicit HYBRID differ in {name}",
        )
    return errors


def check_derivatives(binary, original, root, profile):
    ns = len(profile)
    radius = np.linspace(0.0, 1.0, ns)
    tangent = 0.03 * (1.0 + radius + 0.25 * np.sin(2.0 * np.pi * radius))
    require(tangent.shape == (ns,), "wrong tangent shape")
    require(np.all(np.isfinite(tangent)) and np.all(tangent != 0.0), "bad tangent")
    require(np.ptp(tangent) > 0.0, "tangent must be nonuniform")
    wout = root / "derivative-wout.nc"
    copy_profile(original, wout)
    jvp, construction = run_cli(binary, root, "jvp", "jvp", wout, iotas=tangent)
    require(construction == 1, "omitted JVP construction is not HYBRID")
    derivatives = {name: jvp[f"d_{name}"] for name in SPECTRUM}
    for name in SPECTRUM:
        require(derivatives[name].shape == jvp[name].shape, f"wrong d_{name} shape")
    norm2 = float(sum(np.sum(value**2) for value in derivatives.values()))
    iota_norm2 = float(np.dot(derivatives["iota"], derivatives["iota"]))
    require(np.isfinite(norm2) and norm2 > 0.0, "zero/nonfinite JVP norm")
    require(np.isfinite(iota_norm2) and iota_norm2 > 0.0, "zero/nonfinite d_iota")

    curve = []
    for k in range(2, 10):
        step = 10.0**(-k)
        primal = []
        for sign, label in ((1.0, "plus"), (-1.0, "minus")):
            perturbed = root / f"{label}-{k}-wout.nc"
            copy_profile(original, perturbed, iotas=profile + sign * step * tangent)
            arrays, construction = run_cli(
                binary, root, f"{label}-{k}", "primal", perturbed,
            )
            require(construction == 1, "omitted FD construction is not HYBRID")
            primal.append(arrays)
        residuals = {
            name: (primal[0][name] - primal[1][name]) / (2.0 * step)
            - derivatives[name] for name in SPECTRUM
        }
        error = float(np.sqrt(
            sum(np.sum(value**2) for value in residuals.values()) / norm2
        ))
        iota_error = float(np.linalg.norm(residuals["iota"]) / np.sqrt(iota_norm2))
        require(np.isfinite(error) and np.isfinite(iota_error), "nonfinite FD error")
        curve.append({"step": step, "relative_error": error,
                      "iota_relative_error": iota_error})
    best = min(row["relative_error"] for row in curve)
    best_iota = min(row["iota_relative_error"] for row in curve)
    require(best <= 1e-8, f"JVP V-curve minimum {best:.3e}")
    require(best_iota <= 1e-8, f"d_iota V-curve minimum {best_iota:.3e}")

    vjp, construction = run_cli(
        binary, root, "vjp", "vjp", wout, seed=derivatives["iota"],
    )
    require(construction == 1, "omitted VJP construction is not HYBRID")
    for name in ("iotas", "iotaf"):
        require(vjp[name].shape == (ns,), f"wrong serialized {name} cotangent shape")
    cotangent_norm = float(np.linalg.norm(vjp["iotas"]))
    require(
        np.isfinite(cotangent_norm) and cotangent_norm > 0.0,
        "zero/nonfinite serialized iotas cotangent norm",
    )
    require(np.all(vjp["iotaf"] == 0.0), "HYBRID returned nonzero iotaf cotangent")
    lhs = iota_norm2
    rhs = float(np.dot(vjp["iotas"], tangent))
    duality = abs(lhs - rhs) / max(abs(lhs), abs(rhs))
    require(np.isfinite(duality) and duality <= 1e-13, f"duality error {duality:.3e}")
    return {"fd_curve": curve, "jvp_norm": float(np.sqrt(norm2)),
            "d_iota_norm": float(np.sqrt(iota_norm2)),
            "iotas_cotangent_norm": cotangent_norm, "duality_lhs": lhs,
            "duality_rhs": rhs, "duality_relative_error": duality}


def import_origins():
    origins = {}
    for name in ("numpy", "netCDF4", "cftime", "netCDF4._netCDF4", "cftime._cftime"):
        module = importlib.import_module(name)
        origins[name] = str(Path(module.__file__).resolve())
    try:
        module = importlib.import_module("numpy._core._multiarray_umath")
    except ModuleNotFoundError:
        module = importlib.import_module("numpy.core._multiarray_umath")
    origins[module.__name__] = str(Path(module.__file__).resolve())
    return {"executable": sys.executable, "import_origins": origins}


def main():
    result_path = Path.cwd() / "boozer-hybrid-cli-results.json"
    result_path.unlink(missing_ok=True)
    require(len(sys.argv) == 3, "usage: test_boozer_hybrid_cli.py <CLI> <wout.nc>")
    binary, original = (Path(arg).resolve() for arg in sys.argv[1:])
    with Dataset(original) as data:
        data.set_auto_maskandscale(False)
        ns = int(np.array(data.variables["ns"][...]).item())
        profile = np.array(data.variables["iotas"][:], dtype=np.float64)
        full_profile = np.array(data.variables["iotaf"][:])
    require(ns > 1 and profile.shape == full_profile.shape == (ns,), "wrong ns shape")
    require(np.all(np.isfinite(profile)), "nonfinite original iotas")
    with tempfile.TemporaryDirectory(
        prefix="boozer-hybrid-cli-", dir=Path.cwd(),
    ) as tmp:
        root = Path(tmp)
        result = {"selector_iota_errors": check_selector(binary, original, root, ns)}
        result.update(check_derivatives(binary, original, root, profile))
    result["python"] = import_origins()
    result_path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print("HYBRID CLI selector, iotas JVP finite differences and VJP duality: PASS")


if __name__ == "__main__":
    main()
