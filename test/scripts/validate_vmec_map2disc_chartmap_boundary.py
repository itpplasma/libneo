#!/usr/bin/env python3
"""
Validate that a map2disc-generated VMEC-boundary chartmap reproduces the VMEC boundary.

For selected zeta slices and all theta grid points at rho=1:
  - Compute R_chart = sqrt(x^2 + y^2)/100, Z_chart = z/100 (convert cm->m)
  - Compare against VMEC boundary evaluation at s=1 (last surface index)

Fails if max abs error in R or Z exceeds tolerances.
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


def _load_chartmap_boundary(chartmap_path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with Dataset(chartmap_path, "r") as ds:
        rho = np.array(ds.variables["rho"][:], dtype=float)
        theta = np.array(ds.variables["theta"][:], dtype=float)
        zeta = np.array(ds.variables["zeta"][:], dtype=float)
        x = np.array(ds.variables["x"][:], dtype=float)
        y = np.array(ds.variables["y"][:], dtype=float)
        z = np.array(ds.variables["z"][:], dtype=float)

        if ds.getncattr("zeta_convention") != "cyl":
            raise RuntimeError("expected zeta_convention=cyl")
        if ds.getncattr("rho_convention") != "unknown":
            raise RuntimeError("expected rho_convention=unknown")

    if rho.ndim != 1 or theta.ndim != 1 or zeta.ndim != 1:
        raise RuntimeError("expected rho/theta/zeta to be 1D")
    if rho.size < 2 or theta.size < 2 or zeta.size < 2:
        raise RuntimeError("expected rho/theta/zeta lengths >= 2")

    # File order is (zeta, theta, rho)
    if x.shape != (zeta.size, theta.size, rho.size):
        raise RuntimeError(f"unexpected x shape {x.shape}")

    ir = int(rho.size - 1)
    R_chart = np.sqrt(x[:, :, ir] ** 2 + y[:, :, ir] ** 2) / 100.0
    Z_chart = z[:, :, ir] / 100.0
    return theta, zeta, np.stack([R_chart, Z_chart], axis=-1)


def _select_zeta_indices(nzeta: int) -> list[int]:
    if nzeta < 2:
        return [0]
    candidates = [0, 1, nzeta // 2, nzeta - 2, nzeta - 1]
    out: list[int] = []
    for i in candidates:
        ii = int(max(0, min(nzeta - 1, i)))
        if ii not in out:
            out.append(ii)
    return out


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="validate_vmec_map2disc_chartmap_boundary")
    p.add_argument("--wout", type=Path, required=True)
    p.add_argument("--chartmap", type=Path, required=True)
    p.add_argument("--s-boundary", type=float, default=1.0)
    p.add_argument("--boundary-scale", type=float, default=1.0)
    p.add_argument("--boundary-padding", type=float, default=0.0)
    p.add_argument("--tol-R", type=float, default=1.0e-10)
    p.add_argument("--tol-Z", type=float, default=1.0e-10)
    args = p.parse_args(argv)

    _import_libneo_from_repo()
    from libneo.vmec import VMECGeometry

    wout = args.wout.resolve()
    chartmap = args.chartmap.resolve()
    if not wout.exists() or wout.stat().st_size == 0:
        raise RuntimeError(f"missing wout file: {wout}")
    if not chartmap.exists() or chartmap.stat().st_size == 0:
        raise RuntimeError(f"missing chartmap file: {chartmap}")

    theta, zeta, rz_chart = _load_chartmap_boundary(chartmap)

    geom = VMECGeometry.from_file(str(wout))
    s_boundary = float(args.s_boundary)
    if not (0.0 < s_boundary <= 1.0):
        raise ValueError("s_boundary must be in (0, 1]")
    boundary_scale = float(args.boundary_scale)
    if boundary_scale <= 0.0:
        raise ValueError("boundary_scale must be > 0")
    boundary_padding = float(args.boundary_padding)
    if boundary_padding < 0.0:
        raise ValueError("boundary_padding must be >= 0")

    iz_list = _select_zeta_indices(int(zeta.size))
    max_dR = 0.0
    max_dZ = 0.0
    for iz in iz_list:
        zeta_val = float(zeta[iz])
        R_vmec, Z_vmec, _ = geom.coords_s(s_boundary, theta, zeta_val, use_asym=True)
        if boundary_scale != 1.0 or boundary_padding != 0.0:
            R_axis, Z_axis, _ = geom.coords_s(0.0, np.array([0.0]), zeta_val, use_asym=True)
            R0 = float(R_axis[0])
            Z0 = float(Z_axis[0])
            dR = R_vmec - R0
            dZ = Z_vmec - Z0
            norm = np.sqrt(dR * dR + dZ * dZ)
            if np.any(norm == 0.0):
                raise ValueError("boundary curve coincides with axis at some theta")
            R_vmec = R0 + boundary_scale * dR + boundary_padding * (dR / norm)
            Z_vmec = Z0 + boundary_scale * dZ + boundary_padding * (dZ / norm)
        R_chart = rz_chart[iz, :, 0]
        Z_chart = rz_chart[iz, :, 1]

        dR = np.max(np.abs(R_chart - R_vmec))
        dZ = np.max(np.abs(Z_chart - Z_vmec))
        max_dR = max(max_dR, float(dR))
        max_dZ = max(max_dZ, float(dZ))

    print(
        f"PASS: {chartmap.name} rho=1 vs VMEC s={s_boundary:.6g} "
        f"scale={boundary_scale:.6g} pad={boundary_padding:.6g} m "
        f"max|dR|={max_dR:.3e} m max|dZ|={max_dZ:.3e} m"
    )

    if max_dR > float(args.tol_R) or max_dZ > float(args.tol_Z):
        raise SystemExit(
            f"FAIL: boundary mismatch exceeds tolerance: max|dR|={max_dR:.3e} "
            f"(tol {args.tol_R:.3e}), max|dZ|={max_dZ:.3e} (tol {args.tol_Z:.3e})"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
