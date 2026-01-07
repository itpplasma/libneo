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
    boundary_offset: float = 0.0,
    boundary_param: str = "arc",
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
    - `boundary_offset` optionally displaces the boundary curve outward (meters)
      along the local curve normal in each zeta slice before map2disc.
    - `boundary_param` selects how the boundary is parameterized for map2disc:
      "arc" uses arclength along the offset boundary; "theta" preserves VMEC
      theta labeling and offsets points along the local normal.

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

    boundary_offset = float(boundary_offset)
    boundary_param = str(boundary_param).strip().lower()

    if not (0.0 < float(s_boundary) <= 1.0):
        raise ValueError("s_boundary must be in (0, 1]")
    if boundary_offset < 0.0:
        raise ValueError("boundary_offset must be >= 0")
    if boundary_param not in ("arc", "theta"):
        raise ValueError("boundary_param must be 'arc' or 'theta'")

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

        if boundary_param == "theta":
            if boundary_offset != 0.0:
                theta_base = np.linspace(0.0, 2.0 * np.pi, 4096, endpoint=False, dtype=float)
                R_base, Z_base, _ = geom.coords_s(float(s_boundary), theta_base, phi_val, use_asym=use_asym)
                dR = np.gradient(R_base, theta_base, edge_order=2)
                dZ = np.gradient(Z_base, theta_base, edge_order=2)
                nR = dZ
                nZ = -dR
                norm = np.sqrt(nR**2 + nZ**2)
                norm = np.where(norm == 0.0, 1.0, norm)
                nR = nR / norm
                nZ = nZ / norm
                Rc = float(np.mean(R_base))
                Zc = float(np.mean(Z_base))
                outward = (R_base - Rc) * nR + (Z_base - Zc) * nZ
                flip = np.where(outward < 0.0, -1.0, 1.0)
                nR = nR * flip
                nZ = nZ * flip
                R_off = R_base + boundary_offset * nR
                Z_off = Z_base + boundary_offset * nZ

                theta_ext = np.concatenate([theta_base, [2.0 * np.pi]])
                R_ext = np.concatenate([R_off, [R_off[0]]])
                Z_ext = np.concatenate([Z_off, [Z_off[0]]])

                def curve(t: np.ndarray) -> np.ndarray:
                    tt = np.asarray(t, dtype=float)
                    th = np.mod(tt, 2.0 * np.pi)
                    R = np.interp(th, theta_ext, R_ext)
                    Zc = np.interp(th, theta_ext, Z_ext)
                    return np.array([R, Zc])
            else:
                def curve(t: np.ndarray) -> np.ndarray:
                    th = np.asarray(t, dtype=float)
                    R, Zc, _ = geom.coords_s(float(s_boundary), th, phi_val, use_asym=use_asym)
                    return np.array([R, Zc])
        else:
            theta_base = np.linspace(0.0, 2.0 * np.pi, 4096, endpoint=False, dtype=float)
            try:
                from shapely.geometry import Polygon
            except Exception as exc:  # pragma: no cover
                raise ImportError(
                    "boundary_offset requires shapely; install libneo with the 'chartmap' extra"
                ) from exc

            R_base, Z_base, _ = geom.coords_s(float(s_boundary), theta_base, phi_val, use_asym=use_asym)
            poly = Polygon(np.column_stack([R_base, Z_base]))
            if not poly.is_valid:
                poly = poly.buffer(0.0)
            if poly.is_empty:
                raise ValueError("invalid VMEC boundary polygon for buffering")

            if boundary_offset == 0.0:
                off = poly
            else:
                off = poly.buffer(boundary_offset, join_style=1)
                if off.is_empty:
                    raise ValueError("boundary_offset produced empty geometry")
                if off.geom_type != "Polygon":
                    off = max(list(off.geoms), key=lambda g: g.area)

            coords = np.asarray(off.exterior.coords, dtype=float)
            coords_open = coords[:-1, :]
            i0 = int(np.argmax(coords_open[:, 0]))
            coords_open = np.vstack([coords_open[i0:, :], coords_open[:i0, :]])
            if coords_open.shape[0] >= 2 and coords_open[1, 1] - coords_open[0, 1] < 0.0:
                coords_open = coords_open[::-1, :]
            coords = np.vstack([coords_open, coords_open[0, :]])
            seg = np.sqrt(np.sum((coords[1:] - coords[:-1]) ** 2, axis=1))
            s_coords = np.concatenate(([0.0], np.cumsum(seg)))
            length = float(s_coords[-1])
            if length == 0.0:
                raise ValueError("buffered boundary has zero length")

            def curve(t: np.ndarray) -> np.ndarray:
                tt = np.asarray(t, dtype=float)
                u = (tt % (2.0 * np.pi)) / (2.0 * np.pi)
                s_query = u * length
                R = np.interp(s_query, s_coords, coords[:, 0])
                Zc = np.interp(s_query, s_coords, coords[:, 1])
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
