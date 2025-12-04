#!/usr/bin/env python3
"""
Generate a Babin (map2disc-based) circular test volume for libneo.

The output NetCDF file contains a single toroidal slice (zeta=0) of the
map2disc mapping for a circular boundary. It is placed in the specified
output directory as ``babin_circular.nc`` with the following variables:

- rho(nrho)
- theta(ntheta)
- zeta(nzeta)
- R(nrho, ntheta, nzeta)
- Z(nrho, ntheta, nzeta)
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from map2disc import map as m2d
from netCDF4 import Dataset


def generate_babin_circular(outfile: Path) -> None:
    R0 = 1.7
    a = 0.35
    nrho = 63
    ntheta = 64
    nzeta = 65

    def curve(t: np.ndarray) -> np.ndarray:
        return np.array([R0 + a * np.cos(t), a * np.sin(t)])

    bcm = m2d.BoundaryConformingMapping(
        curve=curve,
        M=16,
        Nt=256,
        Ng=(256, 256),
    )
    bcm.solve_domain2disk()
    bcm.solve_disk2domain()

    rho = np.linspace(0.0, 1.0, nrho)
    theta = np.linspace(0.0, 2.0 * np.pi, ntheta, endpoint=False)
    zeta = np.linspace(0.0, 2.0 * np.pi, nzeta, endpoint=False)

    # Use polar interface: rho in [0,1], theta in [0,2π]
    # eval_rt_1d returns xy with shape (2, nrho, ntheta)
    xy = bcm.eval_rt_1d(rho, theta)
    R_2d = xy[0]  # shape (nrho, ntheta)
    Z_2d = xy[1]  # shape (nrho, ntheta)

    # Expand to 3D with zeta axis: shape (nrho, ntheta, nzeta)
    # For axisymmetric case, replicate the same 2D slice for all zeta
    R_3d = np.repeat(R_2d[:, :, np.newaxis], nzeta, axis=2)
    Z_3d = np.repeat(Z_2d[:, :, np.newaxis], nzeta, axis=2)

    outfile.parent.mkdir(parents=True, exist_ok=True)
    with Dataset(outfile, "w") as ds:
        ds.createDimension("rho", nrho)
        ds.createDimension("theta", ntheta)
        ds.createDimension("zeta", nzeta)

        v_rho = ds.createVariable("rho", "f8", ("rho",))
        v_theta = ds.createVariable("theta", "f8", ("theta",))
        v_zeta = ds.createVariable("zeta", "f8", ("zeta",))

        # NetCDF (C-based) uses row-major order; Fortran uses column-major.
        # For Fortran array(nrho, ntheta, nzeta), Python/NetCDF dims are reversed.
        v_R = ds.createVariable("R", "f8", ("zeta", "theta", "rho"))
        v_Z = ds.createVariable("Z", "f8", ("zeta", "theta", "rho"))

        v_rho[:] = rho
        v_theta[:] = theta
        v_zeta[:] = zeta

        # Transpose from (nrho, ntheta, nzeta) to (nzeta, ntheta, nrho) for NetCDF
        v_R[:, :, :] = np.transpose(R_3d, (2, 1, 0))
        v_Z[:, :, :] = np.transpose(Z_3d, (2, 1, 0))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="setup_babin_volume")
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory in which to place babin_circular.nc",
    )
    args = parser.parse_args(argv)

    out_dir = Path(args.output_dir).resolve()
    outfile = out_dir / "babin_circular.nc"

    if outfile.exists() and outfile.stat().st_size > 0:
        return 0

    generate_babin_circular(outfile)
    if not outfile.exists() or outfile.stat().st_size == 0:
        raise RuntimeError(f"failed to create {outfile}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
