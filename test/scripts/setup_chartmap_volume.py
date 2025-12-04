#!/usr/bin/env python3
"""
Generate a chartmap test volume for libneo using map2disc.

The output NetCDF file contains a toroidal volume mapping for a circular
boundary. It is placed in the specified output directory as ``chartmap.nc``
with the following variables:

- rho(nrho)
- theta(ntheta)
- zeta(nzeta)
- pos(xyz, nrho, ntheta, nzeta) - Cartesian X,Y,Z positions
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from map2disc import map as m2d
from netCDF4 import Dataset


def generate_chartmap_circular(outfile: Path) -> None:
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

    # Use polar interface: rho in [0,1], theta in [0,2pi]
    # eval_rt_1d returns xy with shape (2, nrho, ntheta)
    xy = bcm.eval_rt_1d(rho, theta)
    R_2d = xy[0]  # shape (nrho, ntheta)
    Z_2d = xy[1]  # shape (nrho, ntheta)

    # Build full Cartesian pos(3, nrho, ntheta, nzeta)
    pos = np.zeros((3, nrho, ntheta, nzeta))
    for iz, phi in enumerate(zeta):
        pos[0, :, :, iz] = R_2d * np.cos(phi)  # X
        pos[1, :, :, iz] = R_2d * np.sin(phi)  # Y
        pos[2, :, :, iz] = Z_2d  # Z (same for all zeta in axisymmetric case)

    outfile.parent.mkdir(parents=True, exist_ok=True)
    with Dataset(outfile, "w") as ds:
        ds.createDimension("xyz", 3)
        ds.createDimension("rho", nrho)
        ds.createDimension("theta", ntheta)
        ds.createDimension("zeta", nzeta)

        v_rho = ds.createVariable("rho", "f8", ("rho",))
        v_theta = ds.createVariable("theta", "f8", ("theta",))
        v_zeta = ds.createVariable("zeta", "f8", ("zeta",))

        # NetCDF (C-based) uses row-major order; Fortran uses column-major.
        # For Fortran array pos(3, nrho, ntheta, nzeta), NetCDF dims are reversed.
        v_pos = ds.createVariable("pos", "f8", ("zeta", "theta", "rho", "xyz"))

        v_rho[:] = rho
        v_theta[:] = theta
        v_zeta[:] = zeta

        # Transpose from (3, nrho, ntheta, nzeta) to (nzeta, ntheta, nrho, 3)
        v_pos[:, :, :, :] = np.transpose(pos, (3, 2, 1, 0))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="setup_chartmap_volume")
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory in which to place chartmap.nc",
    )
    args = parser.parse_args(argv)

    out_dir = Path(args.output_dir).resolve()
    outfile = out_dir / "chartmap.nc"

    if outfile.exists() and outfile.stat().st_size > 0:
        return 0

    generate_chartmap_circular(outfile)
    if not outfile.exists() or outfile.stat().st_size == 0:
        raise RuntimeError(f"failed to create {outfile}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
