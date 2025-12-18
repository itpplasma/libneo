from __future__ import annotations

from pathlib import Path

import numpy as np

from .chartmap_io import ChartmapGrid, build_chartmap_grid, write_chartmap_netcdf


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
    if nzeta < 2:
        raise ValueError("nzeta must be >= 2")
    period = _zeta_period(int(num_field_periods))
    zeta = np.linspace(0.0, period, int(nzeta), endpoint=False, dtype=float)
    return build_chartmap_grid(
        nrho=int(nrho),
        ntheta=int(ntheta),
        zeta=zeta,
        num_field_periods=int(num_field_periods),
    )


def write_chartmap_from_vmec_boundary(
    wout_path: str | Path,
    out_path: str | Path,
    *,
    nrho: int = 33,
    ntheta: int = 65,
    nzeta: int = 33,
    s_boundary: float = 1.0,
    boundary_scale: float = 1.0,
    boundary_padding: float = 0.0,
    num_field_periods: int | None = None,
    M: int = 16,
    Nt: int = 256,
    Ng: tuple[int, int] = (256, 256),
    use_asym: bool = True,
) -> None:
    """
    Generate a chartmap NetCDF using map2disc from a VMEC flux-surface boundary.

    - `s_boundary` selects the VMEC surface by normalized toroidal flux in (0, 1].
      (`s_boundary=1.0` corresponds to the LCFS.)
    - `boundary_scale` and `boundary_padding` optionally enlarge the boundary curve
      outward relative to the magnetic axis at each zeta slice:
        - `boundary_scale` scales the axis-to-boundary vector (dimensionless).
        - `boundary_padding` adds a fixed outward offset in meters.

    The output matches libneo_coordinates chartmap conventions:
    - Dimensions: rho, theta, zeta
    - Variables: rho, theta, zeta, x, y, z, num_field_periods
    - x/y/z stored with file dims (zeta, theta, rho)
    - Units: cm
    - Conventions: zeta_convention=cyl, rho_convention=unknown (geometric disk radius)
    """
    from map2disc import map as m2d
    from netCDF4 import Dataset

    from .vmec import VMECGeometry

    wout = _as_path(wout_path)
    out = _as_path(out_path)

    boundary_scale = float(boundary_scale)
    boundary_padding = float(boundary_padding)

    if not (0.0 < float(s_boundary) <= 1.0):
        raise ValueError("s_boundary must be in (0, 1]")
    if boundary_scale <= 0.0:
        raise ValueError("boundary_scale must be > 0")
    if boundary_padding < 0.0:
        raise ValueError("boundary_padding must be >= 0")

    geom = VMECGeometry.from_file(str(wout))

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
        if boundary_scale != 1.0 or boundary_padding != 0.0:
            R_axis, Z_axis, _ = geom.coords_s(0.0, np.array([0.0]), phi_val, use_asym=use_asym)
            axis_rz = (float(R_axis[0]), float(Z_axis[0]))
        else:
            axis_rz = None

        def curve(t: np.ndarray) -> np.ndarray:
            th = np.asarray(t, dtype=float)
            R, Zc, _ = geom.boundary_rz(
                float(s_boundary),
                th,
                phi_val,
                boundary_scale=boundary_scale,
                boundary_padding=boundary_padding,
                axis_rz=axis_rz,
                use_asym=use_asym,
            )
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

    write_chartmap_netcdf(
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
    from .stl_boundary import extract_boundary_slices, write_chartmap_from_stl as _write_from_slices

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
    _write_from_slices(
        slices,
        out,
        nrho=int(nrho),
        ntheta=int(ntheta),
        num_field_periods=int(num_field_periods),
        M=int(M),
        Nt=int(Nt),
        Ng=Ng,
    )
