#!/usr/bin/env python3
"""Validate automatic rho_lcfs estimation for VMEC-theta chartmaps.

This test checks that when rho_lcfs is not provided, libneo derives a sensible
value from geometry by comparing axis-to-LCFS and axis-to-wall distances.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np
from netCDF4 import Dataset


def _import_libneo_from_repo() -> None:
    repo_root = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(repo_root / "python"))


def _median_ratio_rho_lcfs(
    R_axis: np.ndarray,
    Z_axis: np.ndarray,
    R_lcfs: np.ndarray,
    Z_lcfs: np.ndarray,
    R_wall: np.ndarray,
    Z_wall: np.ndarray,
) -> float:
    d_lcfs = np.sqrt((R_lcfs - R_axis) ** 2 + (Z_lcfs - Z_axis) ** 2)
    d_wall = np.sqrt((R_wall - R_axis) ** 2 + (Z_wall - Z_axis) ** 2)
    mask = (d_lcfs > 0.0) & (d_wall > 0.0)
    ratio = d_lcfs[mask] / d_wall[mask]
    ratio = ratio[np.isfinite(ratio)]
    if ratio.size == 0:
        raise RuntimeError("failed to compute rho_lcfs ratio")
    return float(np.median(ratio))


def _snap_to_grid(rho: np.ndarray, val: float) -> float:
    ir = int(np.searchsorted(rho, val, side="right") - 1)
    ir = max(1, min(rho.size - 2, ir))
    return float(rho[ir])


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="validate_vmec_rho_lcfs_auto")
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args(argv)

    out_dir = Path(args.output_dir).resolve()
    wout = out_dir / "wout.nc"
    if not wout.exists() or wout.stat().st_size == 0:
        raise RuntimeError(f"missing wout fixture: {wout}")

    _import_libneo_from_repo()
    from libneo.chartmap import write_chartmap_from_vmec_extended, write_chartmap_from_vmec_to_wall
    from libneo.chartmap_io import build_chartmap_grid
    from libneo.vmec import VMECGeometry

    with Dataset(wout, "r") as ds:
        nfp = int(np.array(ds.variables["nfp"][...]))

    geom = VMECGeometry.from_file(str(wout))

    nrho = 33
    ntheta = 65
    nzeta = 17
    period = 2.0 * np.pi / float(nfp)
    zeta = np.linspace(0.0, period, nzeta, endpoint=False, dtype=float)
    rho = build_chartmap_grid(nrho=nrho, ntheta=ntheta, zeta=zeta, num_field_periods=nfp).rho
    theta_grid = np.linspace(0.0, 2.0 * np.pi, ntheta, endpoint=False, dtype=float)

    # --- extended: synthetic outer boundary defined by VMEC radial derivative
    outfile_ext = out_dir / "wout_vmec_extended_auto_rho_lcfs.chartmap.nc"
    boundary_offset = 0.1
    write_chartmap_from_vmec_extended(
        wout,
        outfile_ext,
        nrho=nrho,
        ntheta=ntheta,
        nzeta=nzeta,
        boundary_offset=boundary_offset,
    )

    with Dataset(outfile_ext, "r") as ds:
        if "rho_lcfs" not in ds.ncattrs():
            raise RuntimeError("missing rho_lcfs global attribute on extended chartmap")
        rho_lcfs_ext = float(ds.getncattr("rho_lcfs"))

    ratios = []
    for phi in zeta:
        R_axis, Z_axis, _ = geom.coords_s(0.0, np.array([0.0]), float(phi), use_asym=True)
        R_axis = float(R_axis[0])
        Z_axis = float(Z_axis[0])
        R_lcfs, Z_lcfs, dR_ds, dZ_ds = geom.coords_s_with_deriv(
            1.0,
            theta_grid,
            float(phi),
            use_asym=True,
        )
        norm = np.sqrt(dR_ds**2 + dZ_ds**2)
        norm = np.where(norm == 0.0, 1.0, norm)
        R_outer = R_lcfs + boundary_offset * (dR_ds / norm)
        Z_outer = Z_lcfs + boundary_offset * (dZ_ds / norm)
        ratios.append(
            _median_ratio_rho_lcfs(
                R_axis=np.full_like(R_lcfs, R_axis),
                Z_axis=np.full_like(Z_lcfs, Z_axis),
                R_lcfs=R_lcfs,
                Z_lcfs=Z_lcfs,
                R_wall=R_outer,
                Z_wall=Z_outer,
            )
        )

    expected_ext = _snap_to_grid(rho, float(np.median(np.array(ratios))))
    if abs(rho_lcfs_ext - expected_ext) > 1.0e-12:
        raise RuntimeError(f"extended rho_lcfs mismatch: got {rho_lcfs_ext}, expected {expected_ext}")

    # --- to_wall: build a synthetic wall by offsetting the LCFS along poloidal-plane normal
    outfile_wall = out_dir / "wout_vmec_to_wall_auto_rho_lcfs.chartmap.nc"
    offset = 0.125

    wall_rz = []
    for phi in zeta:
        dtheta = float(theta_grid[1] - theta_grid[0])

        R_lcfs, Z_lcfs, dR_ds, dZ_ds = geom.coords_s_with_deriv(
            1.0,
            theta_grid,
            float(phi),
            use_asym=True,
        )
        dR_dth = (np.roll(R_lcfs, -1) - np.roll(R_lcfs, 1)) / (2.0 * dtheta)
        dZ_dth = (np.roll(Z_lcfs, -1) - np.roll(Z_lcfs, 1)) / (2.0 * dtheta)

        nR = dZ_dth
        nZ = -dR_dth
        nmag = np.sqrt(nR**2 + nZ**2)
        nmag = np.where(nmag == 0.0, 1.0, nmag)
        nR = nR / nmag
        nZ = nZ / nmag

        outward = nR * dR_ds + nZ * dZ_ds
        flip = outward < 0.0
        nR = np.where(flip, -nR, nR)
        nZ = np.where(flip, -nZ, nZ)

        R_wall = R_lcfs + offset * nR
        Z_wall = Z_lcfs + offset * nZ
        wall_rz.append((R_wall, Z_wall))

    write_chartmap_from_vmec_to_wall(
        wout,
        outfile_wall,
        wall_rz=wall_rz,
        wall_zeta=zeta,
        nrho=nrho,
        ntheta=ntheta,
        rho_lcfs=None,
        num_field_periods=nfp,
        wall_match="poloidal_normal",
    )

    with Dataset(outfile_wall, "r") as ds:
        if "rho_lcfs" not in ds.ncattrs():
            raise RuntimeError("missing rho_lcfs global attribute on to_wall chartmap")
        rho_lcfs_wall = float(ds.getncattr("rho_lcfs"))

    ratios = []
    for phi, (R_wall, Z_wall) in zip(zeta, wall_rz, strict=True):
        R_axis, Z_axis, _ = geom.coords_s(0.0, np.array([0.0]), float(phi), use_asym=True)
        R_axis = float(R_axis[0])
        Z_axis = float(Z_axis[0])
        R_lcfs, Z_lcfs, _ = geom.coords_s(1.0, theta_grid, float(phi), use_asym=True)
        ratios.append(
            _median_ratio_rho_lcfs(
                R_axis=np.full_like(R_lcfs, R_axis),
                Z_axis=np.full_like(Z_lcfs, Z_axis),
                R_lcfs=R_lcfs,
                Z_lcfs=Z_lcfs,
                R_wall=R_wall,
                Z_wall=Z_wall,
            )
        )

    expected_wall = _snap_to_grid(rho, float(np.median(np.array(ratios))))
    if abs(rho_lcfs_wall - expected_wall) > 1.0e-12:
        raise RuntimeError(f"to_wall rho_lcfs mismatch: got {rho_lcfs_wall}, expected {expected_wall}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

