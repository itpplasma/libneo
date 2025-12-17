from __future__ import annotations

from pathlib import Path
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class ChartmapGrid:
    rho: np.ndarray
    theta: np.ndarray
    zeta: np.ndarray
    num_field_periods: int


def _as_path(path: str | Path) -> Path:
    p = Path(path)
    if not p:
        raise ValueError("path must be non-empty")
    return p


def _zeta_period(num_field_periods: int) -> float:
    if num_field_periods < 1:
        raise ValueError("num_field_periods must be >= 1")
    return 2.0 * float(np.pi) / float(num_field_periods)


def _build_grid(
    *,
    nrho: int,
    ntheta: int,
    nzeta: int,
    num_field_periods: int,
) -> ChartmapGrid:
    if nrho < 2 or ntheta < 2 or nzeta < 2:
        raise ValueError("nrho, ntheta, nzeta must be >= 2")

    rho = np.linspace(0.0, 1.0, int(nrho), dtype=float)
    theta = np.linspace(0.0, 2.0 * np.pi, int(ntheta), endpoint=False, dtype=float)
    period = _zeta_period(int(num_field_periods))
    zeta = np.linspace(0.0, period, int(nzeta), endpoint=False, dtype=float)
    return ChartmapGrid(rho=rho, theta=theta, zeta=zeta, num_field_periods=int(num_field_periods))


def _write_chartmap_netcdf(
    out_path: Path,
    *,
    grid: ChartmapGrid,
    x_rtz: np.ndarray,
    y_rtz: np.ndarray,
    z_rtz: np.ndarray,
    zeta_convention: str,
    rho_convention: str,
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


def write_chartmap_from_vmec_boundary(
    wout_path: str | Path,
    out_path: str | Path,
    *,
    nrho: int = 33,
    ntheta: int = 65,
    nzeta: int = 33,
    s_boundary: float = 1.0,
    num_field_periods: int | None = None,
    M: int = 16,
    Nt: int = 256,
    Ng: tuple[int, int] = (256, 256),
    use_asym: bool = True,
) -> None:
    """
    Generate a chartmap NetCDF using map2disc from the VMEC last closed surface.

    The output matches libneo_coordinates chartmap conventions:
    - Dimensions: rho, theta, zeta
    - Variables: rho, theta, zeta, x, y, z, num_field_periods
    - x/y/z stored with file dims (zeta, theta, rho)
    - Units: cm
    """
    from map2disc import map as m2d
    from netCDF4 import Dataset

    from .vmec import VMECGeometry

    wout = _as_path(wout_path)
    out = _as_path(out_path)

    if not (0.0 < float(s_boundary) <= 1.0):
        raise ValueError("s_boundary must be in (0, 1]")

    geom = VMECGeometry.from_file(str(wout))
    s_index = geom.rmnc.shape[1] - 1

    if num_field_periods is None:
        with Dataset(wout, "r") as ds:
            if "nfp" not in ds.variables:
                raise ValueError("wout file missing nfp variable")
            nfp = int(np.array(ds.variables["nfp"][...]))
    else:
        nfp = int(num_field_periods)

    grid = _build_grid(nrho=nrho, ntheta=ntheta, nzeta=nzeta, num_field_periods=nfp)

    x = np.zeros((grid.rho.size, grid.theta.size, grid.zeta.size), dtype=np.float64)
    y = np.zeros_like(x)
    z = np.zeros_like(x)

    for iz, phi in enumerate(grid.zeta):
        phi_val = float(phi)

        def curve(t: np.ndarray) -> np.ndarray:
            th = np.asarray(t, dtype=float)
            R, Zc, _ = geom.coords(s_index, th, phi_val, use_asym=use_asym)
            return np.array([R, Zc])

        bcm = m2d.BoundaryConformingMapping(curve=curve, M=int(M), Nt=int(Nt), Ng=Ng)
        bcm.solve_domain2disk()
        bcm.solve_disk2domain()

        rz = bcm.eval_rt_1d(grid.rho, grid.theta)
        R_grid = rz[0]
        Z_grid = rz[1]

        x[:, :, iz] = R_grid * np.cos(phi_val) * 100.0
        y[:, :, iz] = R_grid * np.sin(phi_val) * 100.0
        z[:, :, iz] = Z_grid * 100.0

    _write_chartmap_netcdf(
        out,
        grid=grid,
        x_rtz=x,
        y_rtz=y,
        z_rtz=z,
        zeta_convention="cyl",
        rho_convention="unknown",
    )


def write_chartmap_from_stl(
    stl_path: str | Path,
    out_path: str | Path,
    *,
    nrho: int = 33,
    ntheta: int = 65,
    nphi: int = 32,
    axis_xy: tuple[float, float] | None = None,
    num_field_periods: int = 1,
    n_boundary_points: int = 512,
    stitch_tol: float = 1.0e-6,
    M: int = 16,
    Nt: int = 256,
    Ng: tuple[int, int] = (256, 256),
) -> None:
    """
    Generate a chartmap NetCDF using map2disc from an STL boundary.
    """
    from map2disc import map as m2d

    from .stl_boundary import extract_boundary_slices, _curve_from_contour

    stl = _as_path(stl_path)
    out = _as_path(out_path)

    if nphi < 2:
        raise ValueError("nphi must be >= 2")
    if num_field_periods < 1:
        raise ValueError("num_field_periods must be >= 1")

    period = _zeta_period(int(num_field_periods))
    phi_vals = np.linspace(0.0, period, int(nphi), endpoint=False, dtype=float)

    slices = extract_boundary_slices(
        stl,
        n_phi=int(nphi),
        n_boundary_points=int(n_boundary_points),
        axis_xy=axis_xy,
        stitch_tol=float(stitch_tol),
        phi_vals=phi_vals,
    )

    grid = _build_grid(
        nrho=nrho, ntheta=ntheta, nzeta=len(slices), num_field_periods=int(num_field_periods)
    )

    x = np.zeros((grid.rho.size, grid.theta.size, grid.zeta.size), dtype=np.float64)
    y = np.zeros_like(x)
    z = np.zeros_like(x)

    for iz, s in enumerate(slices):
        curve = _curve_from_contour(s.outer_filled)
        bcm = m2d.BoundaryConformingMapping(curve=curve, M=int(M), Nt=int(Nt), Ng=Ng)
        bcm.solve_domain2disk()
        bcm.solve_disk2domain()

        rz = bcm.eval_rt_1d(grid.rho, grid.theta)
        R_grid = rz[0]
        Z_grid = rz[1]

        phi = float(s.phi)
        ax_x, ax_y = s.axis_xy
        x[:, :, iz] = (ax_x + R_grid * np.cos(phi)) * 100.0
        y[:, :, iz] = (ax_y + R_grid * np.sin(phi)) * 100.0
        z[:, :, iz] = Z_grid * 100.0

    _write_chartmap_netcdf(
        out,
        grid=grid,
        x_rtz=x,
        y_rtz=y,
        z_rtz=z,
        zeta_convention="cyl",
        rho_convention="unknown",
    )
