#!/usr/bin/env python3
"""Require VMEC ``signgs`` relabeling to preserve Cartesian geometry and B."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
from pathlib import Path

import netCDF4
import numpy as np
from scipy.interpolate import CubicSpline

from manufacture_vmec_theta_relabel import relabel_wout


def clone_without_variable(source: Path, output: Path, omitted: str) -> None:
    """Clone a NetCDF file while omitting one scalar metadata variable."""

    with netCDF4.Dataset(source) as src, netCDF4.Dataset(output, "w") as dst:
        for name, dimension in src.dimensions.items():
            dst.createDimension(
                name, None if dimension.isunlimited() else len(dimension)
            )
        dst.setncatts({name: src.getncattr(name) for name in src.ncattrs()})
        for name, variable in src.variables.items():
            if name == omitted:
                continue
            kwargs = {}
            if "_FillValue" in variable.ncattrs():
                kwargs["fill_value"] = variable.getncattr("_FillValue")
            copied = dst.createVariable(
                name, variable.datatype, variable.dimensions, **kwargs
            )
            copied.setncatts(
                {
                    attr: variable.getncattr(attr)
                    for attr in variable.ncattrs()
                    if attr != "_FillValue"
                }
            )
            copied[:] = variable[:]


def read_chartmap(path: Path) -> dict[str, np.ndarray | float | int]:
    with netCDF4.Dataset(path) as dataset:
        result: dict[str, np.ndarray | float | int] = {
            name: np.asarray(dataset.variables[name][:], dtype=float)
            for name in (
                "rho",
                "s",
                "theta",
                "zeta",
                "x",
                "y",
                "z",
                "A_phi",
                "B_theta",
                "B_phi",
                "Bmod",
            )
        }
        result["torflux"] = float(dataset.torflux)
        result["signgs"] = int(dataset.signgs)
    return result


def mapped(array: np.ndarray, theta_map: np.ndarray) -> np.ndarray:
    return array[:, theta_map, :]


def relative(first: np.ndarray, second: np.ndarray) -> float:
    denominator = np.linalg.norm(first)
    if denominator == 0.0:
        return float(np.linalg.norm(second))
    return float(np.linalg.norm(first - second) / denominator)


def cartesian_field(
    chart: dict[str, np.ndarray | float | int],
    izeta: int,
    itheta: int,
    irho: int,
) -> np.ndarray:
    theta = np.asarray(chart["theta"])
    zeta = np.asarray(chart["zeta"])
    s = np.asarray(chart["s"])
    ntheta = theta.size
    nzeta = zeta.size
    dtheta = 2.0 * np.pi / ntheta
    nfp = round(2.0 * np.pi / (zeta[-1] + (zeta[1] - zeta[0])))
    dzeta = 2.0 * np.pi / (nfp * nzeta)
    position = np.stack((chart["x"], chart["y"], chart["z"]))
    e_theta = (
        position[:, izeta, (itheta + 1) % ntheta, irho]
        - position[:, izeta, (itheta - 1) % ntheta, irho]
    ) / (2.0 * dtheta)
    e_zeta = (
        position[:, (izeta + 1) % nzeta, itheta, irho]
        - position[:, (izeta - 1) % nzeta, itheta, irho]
    ) / (2.0 * dzeta)

    torflux = float(chart["torflux"])
    d_aphi_ds = CubicSpline(s, np.asarray(chart["A_phi"]))(s[irho], 1)
    iota = -float(d_aphi_ds) / torflux
    btheta = float(np.asarray(chart["B_theta"])[irho])
    bzeta = float(np.asarray(chart["B_phi"])[irho])
    bmod = float(np.asarray(chart["Bmod"])[izeta, itheta, irho])
    bcontra_zeta = bmod**2 / (bzeta + iota * btheta)
    return bcontra_zeta * (iota * e_theta + e_zeta)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--wout", type=Path, required=True)
    parser.add_argument("--exporter", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    relabeled_wout = output_dir / "wout_theta_relabel.nc"
    missing_signgs_wout = output_dir / "wout_missing_signgs.nc"
    invalid_signgs_wout = output_dir / "wout_invalid_signgs.nc"
    original_chart = output_dir / "signgs_original.chartmap.nc"
    relabeled_chart = output_dir / "signgs_relabel.chartmap.nc"
    missing_signgs_chart = output_dir / "signgs_missing.chartmap.nc"
    relabel_wout(args.wout.resolve(), relabeled_wout)
    clone_without_variable(args.wout.resolve(), missing_signgs_wout, "signgs")
    shutil.copyfile(args.wout.resolve(), invalid_signgs_wout)
    with netCDF4.Dataset(invalid_signgs_wout, "a") as dataset:
        dataset.variables["signgs"].assignValue(0)

    for source, output in (
        (args.wout.resolve(), original_chart),
        (relabeled_wout, relabeled_chart),
    ):
        subprocess.run([args.exporter.resolve(), source, output], check=True)
    missing_completed = subprocess.run(
        [args.exporter.resolve(), missing_signgs_wout, missing_signgs_chart],
        check=True,
        text=True,
        capture_output=True,
    )
    invalid_completed = subprocess.run(
        [args.exporter.resolve(), invalid_signgs_wout, output_dir / "invalid.nc"],
        check=False,
        text=True,
        capture_output=True,
    )

    original = read_chartmap(original_chart)
    relabeled = read_chartmap(relabeled_chart)
    missing = read_chartmap(missing_signgs_chart)
    ntheta = np.asarray(original["theta"]).size
    theta_map = (-np.arange(ntheta)) % ntheta

    geometry_errors = {
        name: float(
            np.max(
                np.abs(
                    np.asarray(original[name])
                    - mapped(np.asarray(relabeled[name]), theta_map)
                )
            )
        )
        for name in ("x", "y", "z")
    }
    bmod_error = relative(
        np.asarray(original["Bmod"]),
        mapped(np.asarray(relabeled["Bmod"]), theta_map),
    )

    same_errors: list[float] = []
    reversed_errors: list[float] = []
    samples = ((5, 3, 10), (13, 7, 20), (27, 11, 30), (41, 17, 40))
    for izeta, itheta, irho in samples:
        mapped_theta = int(theta_map[itheta])
        first = cartesian_field(original, izeta, itheta, irho)
        second = cartesian_field(relabeled, izeta, mapped_theta, irho)
        norm = float(np.linalg.norm(first))
        same_errors.append(float(np.linalg.norm(second - first) / norm))
        reversed_errors.append(float(np.linalg.norm(second + first) / norm))

    with netCDF4.Dataset(args.wout) as dataset:
        phi_edge = float(dataset.variables["phi"][-1])
        input_signgs = int(dataset.variables["signgs"][:])
    expected_torflux = input_signgs * phi_edge / (2.0 * np.pi) * 1.0e8

    metrics = {
        "maximum_geometry_error_cm": max(geometry_errors.values()),
        "relative_l2_Bmod_after_theta_map": bmod_error,
        "relative_Azeta_error_against_same": relative(
            np.asarray(original["A_phi"]), np.asarray(relabeled["A_phi"])
        ),
        "relative_Btheta_error_against_opposite": relative(
            np.asarray(original["B_theta"]), -np.asarray(relabeled["B_theta"])
        ),
        "relative_Bzeta_error_against_same": relative(
            np.asarray(original["B_phi"]), np.asarray(relabeled["B_phi"])
        ),
        "maximum_relative_cartesian_B_error": max(same_errors),
        "minimum_wrong_reversed_B_error": min(reversed_errors),
        "original_torflux_relative_error": abs(
            float(original["torflux"]) - expected_torflux
        )
        / abs(expected_torflux),
        "relabel_torflux_opposite_relative_error": abs(
            float(relabeled["torflux"]) + float(original["torflux"])
        )
        / abs(float(original["torflux"])),
        "missing_signgs_torflux_relative_error": abs(
            float(missing["torflux"]) - float(original["torflux"])
        )
        / abs(float(original["torflux"])),
    }
    gates = {
        "input_signgs_is_minus_one": input_signgs == -1,
        "chartmap_signgs_transforms": (
            int(original["signgs"]),
            int(relabeled["signgs"]),
        )
        == (-1, 1),
        "geometry_is_same": metrics["maximum_geometry_error_cm"] < 2.0e-7,
        "Bmod_is_same": bmod_error < 2.0e-10,
        "vector_potential_covariant_transforms": metrics[
            "relative_Azeta_error_against_same"
        ]
        < 2.0e-9,
        "B_covariants_transform": max(
            metrics["relative_Btheta_error_against_opposite"],
            metrics["relative_Bzeta_error_against_same"],
        )
        < 2.0e-7,
        "Cartesian_B_is_same": metrics["maximum_relative_cartesian_B_error"] < 2.0e-3,
        "wrong_B_reversal_is_rejected": metrics["minimum_wrong_reversed_B_error"] > 1.9,
        "standard_path_is_numerically_unchanged": metrics[
            "original_torflux_relative_error"
        ]
        < 1.0e-14,
        "relabel_torflux_changes_sign": metrics[
            "relabel_torflux_opposite_relative_error"
        ]
        < 1.0e-14,
        "missing_signgs_uses_documented_minus_one_default": (
            int(missing["signgs"]) == -1
            and metrics["missing_signgs_torflux_relative_error"] < 1.0e-14
            and "assuming the historical signgs=-1" in missing_completed.stdout
        ),
        "invalid_signgs_is_rejected": (
            invalid_completed.returncode != 0
            and "signgs" in (invalid_completed.stdout + invalid_completed.stderr)
        ),
    }
    evidence = {"gates": gates, "metrics": metrics}
    (output_dir / "vmec_signgs_relabel.json").write_text(
        json.dumps(evidence, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(evidence, indent=2, sort_keys=True))
    if not all(gates.values()):
        raise RuntimeError(f"VMEC signgs relabel gate failed: {gates}")


if __name__ == "__main__":
    main()
