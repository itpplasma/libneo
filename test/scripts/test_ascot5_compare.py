#!/usr/bin/env python
"""Cross-check libneo Cerfon-Freidberg solver against ASCOT5's B_GS implementation."""
from __future__ import annotations

import argparse
import ctypes
import math
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


R0_DEFAULT = 6.2
Z0_DEFAULT = 0.0
BPHI0_DEFAULT = 5.3
PSI_MULT_DEFAULT = 200.0
PSI0_DEFAULT = 0.0
PSI1_DEFAULT = 1.0
EPSILON_DEFAULT = 0.32
A_PARAM_DEFAULT = -0.142


class AscotGS(ctypes.Structure):
    _fields_ = [
        ("R0", ctypes.c_double),
        ("z0", ctypes.c_double),
        ("raxis", ctypes.c_double),
        ("zaxis", ctypes.c_double),
        ("B_phi0", ctypes.c_double),
        ("psi0", ctypes.c_double),
        ("psi1", ctypes.c_double),
        ("psi_mult", ctypes.c_double),
        ("psi_coeff", ctypes.c_double * 14),
        ("Nripple", ctypes.c_int),
        ("a0", ctypes.c_double),
        ("alpha0", ctypes.c_double),
        ("delta0", ctypes.c_double),
    ]


@dataclass(frozen=True)
class CaseConfig:
    name: str
    kappa: float
    delta: float
    csv: str
    png: str


CASES = (
    CaseConfig("circular", 1.0, 0.0, "flux_circular.csv", "flux_circular.png"),
    CaseConfig("shaped", 1.7, 0.33, "flux_shaped.csv", "flux_shaped.png"),
)


def run(cmd: list[str], cwd: Path | None = None, **kwargs) -> None:
    subprocess.run(cmd, cwd=cwd, check=True, **kwargs)


def compile_ascot_library(clone_dir: Path) -> Path:
    compiler = shutil.which("gcc-14") or shutil.which("gcc")
    if compiler is None:
        raise RuntimeError("gcc or gcc-14 not found; required to build ascot5 reference library")

    flags: list[str] = []
    link_flags: list[str] = ["-fopenmp", "-lhdf5_hl", "-lhdf5"]

    # Check for Homebrew (macOS)
    homebrew_root = Path("/opt/homebrew")
    if homebrew_root.exists():
        hb_include = homebrew_root / "include"
        if hb_include.exists():
            flags.append(f"-I{hb_include}")
        omp_include = homebrew_root / "opt" / "libomp" / "include"
        if omp_include.exists():
            flags.append(f"-I{omp_include}")
        hb_lib = homebrew_root / "lib"
        if hb_lib.exists():
            link_flags.append(f"-L{hb_lib}")
        omp_lib = homebrew_root / "opt" / "libomp" / "lib"
        if omp_lib.exists():
            link_flags.append(f"-L{omp_lib}")
            link_flags.append("-lomp")

    # Check for standard Linux HDF5 locations
    for hdf5_include in [Path("/usr/include/hdf5/serial"), Path("/usr/include/hdf5"), Path("/usr/include")]:
        if hdf5_include.exists():
            flags.append(f"-I{hdf5_include}")
            break

    # Architecture-agnostic library paths (x86_64, aarch64, etc.)
    import platform
    machine = platform.machine()
    multiarch = None
    if machine == "x86_64":
        multiarch = "x86_64-linux-gnu"
    elif machine in ("aarch64", "arm64"):
        multiarch = "aarch64-linux-gnu"

    lib_search_paths = []
    if multiarch:
        lib_search_paths.extend([
            Path(f"/usr/lib/{multiarch}/hdf5/serial"),
            Path(f"/usr/lib/{multiarch}")
        ])
    lib_search_paths.extend([Path("/usr/lib/hdf5/serial"), Path("/usr/lib")])

    for hdf5_lib in lib_search_paths:
        if hdf5_lib.exists() and (hdf5_lib / "libhdf5.so").exists():
            link_flags.append(f"-L{hdf5_lib}")
            break

    run(["make", "clean"], cwd=clone_dir)
    run([
        "make",
        f"CC={compiler}",
        f"FLAGS={' '.join(flags)}",
        f"LFLAGS={' '.join(link_flags)}",
        "libascot",
        "-j2",
    ], cwd=clone_dir)

    for candidate in (clone_dir / "build" / "libascot.so", clone_dir / "build" / "libascot.dylib"):
        if candidate.exists():
            return candidate
    raise FileNotFoundError("ASCOT5 build did not produce libascot shared library")


def load_ascot_functions(lib_path: Path) -> tuple[ctypes.CDLL, AscotGS]:
    lib = ctypes.CDLL(str(lib_path))
    lib.B_GS_init.argtypes = [
        ctypes.POINTER(AscotGS),
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_double,
    ]
    lib.B_GS_init.restype = ctypes.c_int
    lib.B_GS_eval_psi.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.POINTER(AscotGS),
    ]
    lib.B_GS_eval_psi.restype = ctypes.c_int
    lib.B_GS_get_axis_rz.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(AscotGS)]
    lib.B_GS_get_axis_rz.restype = ctypes.c_int
    lib.B_GS_free.argtypes = [ctypes.POINTER(AscotGS)]
    lib.B_GS_free.restype = None
    return lib, AscotGS()


def evaluate_ascot_flux(lib: ctypes.CDLL, data: AscotGS, points: np.ndarray) -> np.ndarray:
    flux = np.empty(points.shape[0], dtype=np.float64)
    psi_val = ctypes.c_double()
    for idx, (r, z) in enumerate(points):
        err = lib.B_GS_eval_psi(ctypes.byref(psi_val), float(r), 0.0, float(z), ctypes.byref(data))
        if err != 0:
            raise RuntimeError(f"B_GS_eval_psi failed at point index {idx}")
        flux[idx] = psi_val.value
    return flux


def read_flux_csv(csv_path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    data = np.loadtxt(csv_path, delimiter=",", skiprows=1)
    R = data[:, 0]
    Z = data[:, 1]
    psi = data[:, 2]
    return R, Z, psi


def reshape_grid(R: np.ndarray, Z: np.ndarray, values: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    r_unique = np.unique(R)
    z_unique = np.unique(Z)
    if r_unique.size * z_unique.size != values.size:
        raise ValueError("Flux data cannot be reshaped to a rectangular grid.")
    order = np.argsort(R * 1e6 + Z)
    values_sorted = values[order]
    R_sorted = R[order]
    Z_sorted = Z[order]
    V_grid = values_sorted.reshape((r_unique.size, z_unique.size))
    R_grid = R_sorted.reshape((r_unique.size, z_unique.size))
    Z_grid = Z_sorted.reshape((r_unique.size, z_unique.size))
    return R_grid, Z_grid, V_grid


def compute_axis_shift(R_grid: np.ndarray, Z_grid: np.ndarray, psi_grid: np.ndarray, R0: float) -> float:
    z_axis_idx = np.argmin(np.abs(Z_grid[0, :]))
    psi_slice = psi_grid[:, z_axis_idx]
    axis_idx = np.argmin(psi_slice)
    return float(R_grid[axis_idx, z_axis_idx] - R0)


def build_coeff_array(analytic_gs, kappa: float, delta: float) -> np.ndarray:
    coeffs = analytic_gs(
        A=A_PARAM_DEFAULT,
        R0=R0_DEFAULT,
        epsilon=EPSILON_DEFAULT,
        kappa=kappa,
        delta=delta,
        Xpointx=0.0,
        Xpointy=0.0,
        sym=True,
    )
    coeff_array = np.zeros(14, dtype=np.float64)
    coeff_array[: coeffs.size] = coeffs
    coeff_array[12] = A_PARAM_DEFAULT
    return coeff_array


def initialise_case(lib: ctypes.CDLL, coeff_array: np.ndarray) -> AscotGS:
    coeff_ptr = (ctypes.c_double * len(coeff_array))(*coeff_array)
    data = AscotGS()
    init_code = lib.B_GS_init(
        ctypes.byref(data),
        R0_DEFAULT,
        Z0_DEFAULT,
        R0_DEFAULT,
        Z0_DEFAULT,
        BPHI0_DEFAULT,
        PSI0_DEFAULT,
        PSI1_DEFAULT,
        PSI_MULT_DEFAULT,
        coeff_ptr,
        0,
        1.0,
        0.0,
        0.0,
    )
    if init_code != 0:
        raise RuntimeError("B_GS_init failed")
    return data


def overlay_plot(R_grid: np.ndarray, Z_grid: np.ndarray, psi_lib: np.ndarray, psi_asc: np.ndarray, png_path: Path, title: str) -> None:
    psi_lib_aligned = psi_lib - psi_lib.min()
    psi_asc_aligned = psi_asc - psi_asc.min()
    levels = np.linspace(
        min(psi_lib_aligned.min(), psi_asc_aligned.min()),
        max(psi_lib_aligned.max(), psi_asc_aligned.max()),
        15,
    )
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.contour(R_grid, Z_grid, psi_lib_aligned, levels=levels, colors="C0", linestyles="solid", linewidths=1.1)
    ax.contour(R_grid, Z_grid, psi_asc_aligned, levels=levels, colors="C1", linestyles="dashed", linewidths=1.1)
    ax.contour(R_grid, Z_grid, psi_lib_aligned, levels=[0.0], colors="red", linewidths=2.0)
    ax.set_xlabel("R [m]")
    ax.set_ylabel("Z [m]")
    ax.set_title(title)
    ax.set_aspect("equal")
    handles = [
        plt.Line2D([], [], color="C0", linestyle="-", label="libneo"),
        plt.Line2D([], [], color="C1", linestyle="--", label="ascot5"),
    ]
    ax.legend(handles=handles, loc="best")
    fig.tight_layout()
    fig.savefig(png_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path, default=None, help="Path to libneo build directory")
    parser.add_argument("--output-dir", type=Path, default=None, help="Directory for artifacts")
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[2]
    build_dir = (args.build_dir or repo_root / "build").resolve()
    if not build_dir.exists():
        raise FileNotFoundError(f"Build directory {build_dir} does not exist")

    fortran_exe = build_dir / "test" / "test_analytical_circular.x"
    if not fortran_exe.exists():
        raise FileNotFoundError(f"Required executable not found: {fortran_exe}")

    artifact_root = (args.output_dir or build_dir / "test" / "ascot5_compare").resolve()
    artifact_root.mkdir(parents=True, exist_ok=True)

    libneo_dir = artifact_root / "libneo_flux"
    libneo_dir.mkdir(exist_ok=True)

    run([str(fortran_exe), str(libneo_dir)])

    ascot_clone = artifact_root / "ascot5"
    if ascot_clone.exists():
        shutil.rmtree(ascot_clone)
    run(["git", "clone", "--depth", "1", "https://github.com/ascot4fusion/ascot5.git", str(ascot_clone)])

    lib_path = compile_ascot_library(ascot_clone)
    lib, _ = load_ascot_functions(lib_path)

    # Check for required Python dependencies
    try:
        import unyt
    except ImportError:
        raise ImportError(
            "The 'unyt' module is required for ASCOT5 comparison.\n"
            "Please install it on the system:\n"
            "  sudo apt install python3-unyt\n"
            "or install manually if not on Debian/Ubuntu."
        )

    # Add ASCOT5 repo root to path so 'a5py' package can be imported
    ascot_root = str(ascot_clone)
    if ascot_root in sys.path:
        sys.path.remove(ascot_root)
    sys.path.insert(0, ascot_root)

    # Clear any previously imported a5py modules to force reload from correct path
    for module_name in list(sys.modules.keys()):
        if module_name.startswith('a5py') or module_name.startswith('physlib'):
            del sys.modules[module_name]

    from a5py.physlib.analyticequilibrium import analyticGS  # type: ignore

    rms_tol = 3.0e-2
    shift_tol = 2.0e-3

    results = []
    for case in CASES:
        csv_path = libneo_dir / case.csv
        if not csv_path.exists():
            raise FileNotFoundError(f"Expected flux output missing: {csv_path}")

        R, Z, psi_lib = read_flux_csv(csv_path)
        R_grid, Z_grid, psi_lib_grid = reshape_grid(R, Z, psi_lib)

        coeff_array = build_coeff_array(analyticGS, kappa=case.kappa, delta=case.delta)
        data = initialise_case(lib, coeff_array)
        psi_asc = evaluate_ascot_flux(lib, data, np.column_stack((R, Z)))
        psi_asc_grid = psi_asc.reshape(R_grid.shape)

        psi_lib_aligned = psi_lib_grid - psi_lib_grid.min()
        psi_asc_aligned = psi_asc_grid - psi_asc_grid.min()
        diff = psi_asc_aligned - psi_lib_aligned
        rms = float(np.sqrt(np.mean(diff**2)))

        shift_lib = compute_axis_shift(R_grid, Z_grid, psi_lib_aligned, R0_DEFAULT)
        shift_asc = compute_axis_shift(R_grid, Z_grid, psi_asc_aligned, R0_DEFAULT)
        shift_diff = abs(shift_lib - shift_asc)

        if not math.isfinite(rms) or rms > rms_tol:
            raise AssertionError(f"[{case.name}] Flux RMS difference {rms:.3e} exceeds tolerance {rms_tol:.3e}")
        if not math.isfinite(shift_diff) or shift_diff > shift_tol:
            raise AssertionError(f"[{case.name}] Axis shift mismatch {shift_diff:.3e} exceeds tolerance {shift_tol:.3e}")

        png_path = libneo_dir / case.png
        overlay_plot(
            R_grid,
            Z_grid,
            psi_lib_grid,
            psi_asc_grid,
            png_path,
            f"{case.name.capitalize()} Tokamak Flux Surfaces (Cerfon-Freidberg vs ASCOT5)",
        )

        lib.B_GS_free(ctypes.byref(data))
        results.append((case.name, rms, shift_diff, png_path))

    extra_plot = artifact_root / "ascot5_flux_comparison.png"
    if extra_plot.exists():
        extra_plot.unlink()

    for name, rms, shift_diff, png_path in results:
        print(f"[{name}] Flux RMS difference: {rms:.6e} Wb")
        print(f"[{name}] Axis shift |Δ|:     {shift_diff:.6e} m")
        print(f"[{name}] Saved overlay to {png_path}")


if __name__ == "__main__":
    main()
