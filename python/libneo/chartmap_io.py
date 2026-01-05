from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class ChartmapGrid:
    rho: np.ndarray
    theta: np.ndarray
    zeta: np.ndarray
    num_field_periods: int


def build_chartmap_grid(
    *,
    nrho: int,
    ntheta: int,
    zeta: np.ndarray,
    num_field_periods: int,
) -> ChartmapGrid:
    if nrho < 2 or ntheta < 2:
        raise ValueError("nrho and ntheta must be >= 2")
    if num_field_periods < 1:
        raise ValueError("num_field_periods must be >= 1")

    z = np.asarray(zeta, dtype=float)
    if z.ndim != 1 or z.size < 2:
        raise ValueError("zeta must be a 1D array with at least 2 values")

    rho = np.linspace(0.0, 1.0, int(nrho), dtype=float)
    theta = np.linspace(0.0, 2.0 * np.pi, int(ntheta), endpoint=False, dtype=float)
    return ChartmapGrid(rho=rho, theta=theta, zeta=z, num_field_periods=int(num_field_periods))


def write_chartmap_netcdf(
    out_path: Path,
    *,
    grid: ChartmapGrid,
    x_rtz: np.ndarray,
    y_rtz: np.ndarray,
    z_rtz: np.ndarray,
    zeta_convention: str,
    rho_convention: str,
    attrs: dict[str, object] | None = None,
) -> None:
    from netCDF4 import Dataset

    if x_rtz.shape != y_rtz.shape or x_rtz.shape != z_rtz.shape:
        raise ValueError("x, y, z arrays must have the same shape")
    if x_rtz.shape != (grid.rho.size, grid.theta.size, grid.zeta.size):
        raise ValueError("x, y, z must have shape (nrho, ntheta, nzeta)")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with Dataset(out_path, "w", format="NETCDF4") as ds:
        ds.setncattr("zeta_convention", str(zeta_convention))
        ds.setncattr("rho_convention", str(rho_convention))
        if attrs:
            for k, v in attrs.items():
                ds.setncattr(str(k), v)

        ds.createDimension("rho", grid.rho.size)
        ds.createDimension("theta", grid.theta.size)
        ds.createDimension("zeta", grid.zeta.size)

        v_rho = ds.createVariable("rho", "f8", ("rho",))
        v_theta = ds.createVariable("theta", "f8", ("theta",))
        v_zeta = ds.createVariable("zeta", "f8", ("zeta",))

        v_x = ds.createVariable("x", "f8", ("zeta", "theta", "rho"))
        v_y = ds.createVariable("y", "f8", ("zeta", "theta", "rho"))
        v_z = ds.createVariable("z", "f8", ("zeta", "theta", "rho"))
        v_nfp = ds.createVariable("num_field_periods", "i4")

        v_x.units = "cm"
        v_y.units = "cm"
        v_z.units = "cm"

        v_rho[:] = grid.rho
        v_theta[:] = grid.theta
        v_zeta[:] = grid.zeta

        v_x[:, :, :] = np.transpose(x_rtz, (2, 1, 0))
        v_y[:, :, :] = np.transpose(y_rtz, (2, 1, 0))
        v_z[:, :, :] = np.transpose(z_rtz, (2, 1, 0))

        v_nfp.assignValue(int(grid.num_field_periods))
