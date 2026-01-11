from __future__ import annotations

from pathlib import Path

import numpy as np

from .chartmap_io import ChartmapGrid, build_chartmap_grid, write_chartmap_netcdf


def _normalize_boundary_param(value: str) -> str:
    param = str(value).strip().lower()
    if param not in ("arc", "theta"):
        raise ValueError("boundary_param must be 'arc' or 'theta'")
    return param


def _vmec_theta_curve(
    geom,
    s_boundary: float,
    phi_val: float,
    *,
    boundary_offset: float,
    use_asym: bool,
):
    if boundary_offset == 0.0:
        def curve(t: np.ndarray) -> np.ndarray:
            th = np.asarray(t, dtype=float)
            R, Zc, _ = geom.coords_s(float(s_boundary), th, phi_val, use_asym=use_asym)
            return np.array([R, Zc])
        return curve

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

    return curve


def _vmec_arc_curve(
    geom,
    s_boundary: float,
    phi_val: float,
    *,
    boundary_offset: float,
    use_asym: bool,
):
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

    return curve


def _vmec_boundary_curve(
    geom,
    s_boundary: float,
    phi_val: float,
    *,
    boundary_offset: float,
    boundary_param: str,
    use_asym: bool,
):
    if boundary_param == "theta":
        return _vmec_theta_curve(
            geom,
            s_boundary,
            phi_val,
            boundary_offset=boundary_offset,
            use_asym=use_asym,
        )
    return _vmec_arc_curve(
        geom,
        s_boundary,
        phi_val,
        boundary_offset=boundary_offset,
        use_asym=use_asym,
    )


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
    boundary_param = _normalize_boundary_param(boundary_param)

    if not (0.0 < float(s_boundary) <= 1.0):
        raise ValueError("s_boundary must be in (0, 1]")
    if boundary_offset < 0.0:
        raise ValueError("boundary_offset must be >= 0")

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
        curve = _vmec_boundary_curve(
            geom,
            float(s_boundary),
            phi_val,
            boundary_offset=boundary_offset,
            boundary_param=boundary_param,
            use_asym=use_asym,
        )

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


def _cubic_hermite(t: np.ndarray, p0: np.ndarray, m0: np.ndarray, p1: np.ndarray, m1: np.ndarray) -> np.ndarray:
    """
    Cubic Hermite interpolation.

    Given endpoints p0, p1 and tangents m0, m1 (scaled by interval width),
    evaluate the cubic Hermite spline at parameter t in [0, 1].

    H(t) = (2t³ - 3t² + 1)·p0 + (t³ - 2t² + t)·m0 + (-2t³ + 3t²)·p1 + (t³ - t²)·m1
    """
    t = np.asarray(t, dtype=float)
    t2 = t * t
    t3 = t2 * t
    h00 = 2.0 * t3 - 3.0 * t2 + 1.0
    h10 = t3 - 2.0 * t2 + t
    h01 = -2.0 * t3 + 3.0 * t2
    h11 = t3 - t2
    return h00 * p0 + h10 * m0 + h01 * p1 + h11 * m1


def write_chartmap_from_vmec_extended(
    wout_path: str | Path,
    out_path: str | Path,
    *,
    nrho: int = 33,
    ntheta: int = 65,
    nzeta: int = 33,
    nrho_vmec: int | None = None,
    rho_lcfs: float | None = None,
    boundary_offset: float = 0.1,
    num_field_periods: int | None = None,
    use_asym: bool = True,
) -> None:
    """
    Generate a chartmap NetCDF using VMEC coordinates with cubic Hermite extension.

    This mode uses exact VMEC flux coordinates inside the LCFS and extends beyond
    using cubic Hermite interpolation along coordinate line tangents. This provides
    C1 continuity at the LCFS boundary, suitable for symplectic integrators.

    Parameters
    ----------
    wout_path : path to VMEC wout NetCDF file
    out_path : output chartmap NetCDF path
    nrho : total number of radial grid points
    ntheta : number of poloidal grid points
    nzeta : number of toroidal grid points
    nrho_vmec : number of radial points inside LCFS (alternative to rho_lcfs)
    rho_lcfs : radial location of LCFS in [0, 1] (alternative to nrho_vmec)
    boundary_offset : distance beyond LCFS in meters for outer boundary
    num_field_periods : override nfp from wout file
    use_asym : include asymmetric Fourier terms if available

    Grid structure
    --------------
    - rho in [0, 1] with nrho points
    - For rho <= rho_lcfs: VMEC coordinates with s = (rho / rho_lcfs)^2
    - For rho > rho_lcfs: cubic Hermite extrapolation from LCFS tangent

    The cubic Hermite spline matches position and first derivative at the LCFS,
    providing C1 continuity. At the outer boundary, the tangent points radially
    outward to ensure smooth coordinate lines.
    """
    from netCDF4 import Dataset

    from .vmec import VMECGeometry

    wout = _as_path(wout_path)
    out = _as_path(out_path)

    boundary_offset = float(boundary_offset)
    if boundary_offset <= 0.0:
        raise ValueError("boundary_offset must be > 0 for extended mode")

    if nrho_vmec is not None and rho_lcfs is not None:
        raise ValueError("provide only one of nrho_vmec or rho_lcfs")
    if nrho_vmec is None and rho_lcfs is None:
        nrho_vmec = max(2, int(0.8 * nrho))

    geom = VMECGeometry.from_file(str(wout))

    if num_field_periods is None:
        with Dataset(wout, "r") as ds:
            if "nfp" not in ds.variables:
                raise ValueError("wout file missing nfp variable")
            nfp = int(np.array(ds.variables["nfp"][...]))
    else:
        nfp = int(num_field_periods)

    grid = _build_grid(nrho=nrho, ntheta=ntheta, nzeta=nzeta, num_field_periods=nfp)

    if nrho_vmec is not None:
        nrho_vmec = int(nrho_vmec)
        if nrho_vmec < 2 or nrho_vmec >= grid.rho.size:
            raise ValueError("nrho_vmec must be in [2, nrho-1]")
        ir_lcfs = nrho_vmec - 1
        rho_lcfs_val = float(grid.rho[ir_lcfs])
    else:
        assert rho_lcfs is not None
        rho_lcfs_val = float(rho_lcfs)
        if not (0.0 < rho_lcfs_val < 1.0):
            raise ValueError("rho_lcfs must be in (0, 1)")
        ir_lcfs = int(np.searchsorted(grid.rho, rho_lcfs_val, side="right") - 1)
        if ir_lcfs < 1 or ir_lcfs >= grid.rho.size - 1:
            raise ValueError("rho_lcfs must leave room for extended region")
        rho_lcfs_val = float(grid.rho[ir_lcfs])

    x = np.zeros((grid.rho.size, grid.theta.size, grid.zeta.size), dtype=np.float64)
    y = np.zeros_like(x)
    z = np.zeros_like(x)

    for iz, phi in enumerate(grid.zeta):
        phi_val = float(phi)

        R_lcfs, Z_lcfs, dR_ds, dZ_ds = geom.coords_s_with_deriv(
            1.0, grid.theta, phi_val, use_asym=use_asym
        )

        outward_norm = np.sqrt(dR_ds**2 + dZ_ds**2)
        outward_norm = np.where(outward_norm == 0.0, 1.0, outward_norm)
        outward_R = dR_ds / outward_norm
        outward_Z = dZ_ds / outward_norm

        R_outer = R_lcfs + boundary_offset * outward_R
        Z_outer = Z_lcfs + boundary_offset * outward_Z

        ds_drho_at_lcfs = 2.0 / rho_lcfs_val
        drho_dt = 1.0 - rho_lcfs_val
        m0_R = dR_ds * ds_drho_at_lcfs * drho_dt
        m0_Z = dZ_ds * ds_drho_at_lcfs * drho_dt

        m1_R = outward_R * boundary_offset
        m1_Z = outward_Z * boundary_offset

        for ir in range(ir_lcfs + 1):
            rho_val = grid.rho[ir]
            s_val = (float(rho_val) / rho_lcfs_val) ** 2
            R_in, Z_in, _ = geom.coords_s(s_val, grid.theta, phi_val, use_asym=use_asym)
            x[ir, :, iz] = R_in * np.cos(phi_val) * 100.0
            y[ir, :, iz] = R_in * np.sin(phi_val) * 100.0
            z[ir, :, iz] = Z_in * 100.0

        n_outer = grid.rho.size - ir_lcfs - 1
        if n_outer > 0:
            rho_outer = grid.rho[ir_lcfs + 1:]
            t_param = (rho_outer - rho_lcfs_val) / (1.0 - rho_lcfs_val)

            for ith in range(grid.theta.size):
                R_interp = _cubic_hermite(
                    t_param,
                    R_lcfs[ith],
                    m0_R[ith],
                    R_outer[ith],
                    m1_R[ith],
                )
                Z_interp = _cubic_hermite(
                    t_param,
                    Z_lcfs[ith],
                    m0_Z[ith],
                    Z_outer[ith],
                    m1_Z[ith],
                )
                x[ir_lcfs + 1:, ith, iz] = R_interp * np.cos(phi_val) * 100.0
                y[ir_lcfs + 1:, ith, iz] = R_interp * np.sin(phi_val) * 100.0
                z[ir_lcfs + 1:, ith, iz] = Z_interp * 100.0

    write_chartmap_netcdf(
        out,
        grid=grid,
        x_rtz=x,
        y_rtz=y,
        z_rtz=z,
        zeta_convention="cyl",
        rho_convention="vmec_extended",
    )


def _match_wall_to_lcfs_by_fraction(
    R_lcfs: np.ndarray,
    Z_lcfs: np.ndarray,
    R_wall: np.ndarray,
    Z_wall: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Match wall points to LCFS points by fractional position along curve.

    Both curves are aligned to start at their maximum-R point (outboard midplane),
    then wall points are interpolated at the same fractional positions as LCFS.
    """
    i_lcfs_start = int(np.argmax(R_lcfs))
    n_lcfs = len(R_lcfs)
    lcfs_order = np.r_[np.arange(i_lcfs_start, n_lcfs), np.arange(0, i_lcfs_start)]
    R_lcfs_ordered = R_lcfs[lcfs_order]
    Z_lcfs_ordered = Z_lcfs[lcfs_order]

    i_wall_start = int(np.argmax(R_wall))
    n_wall = len(R_wall)
    wall_order = np.r_[np.arange(i_wall_start, n_wall), np.arange(0, i_wall_start)]
    R_wall_ordered = R_wall[wall_order]
    Z_wall_ordered = Z_wall[wall_order]

    if Z_lcfs_ordered[1] < Z_lcfs_ordered[0]:
        R_lcfs_ordered = R_lcfs_ordered[::-1]
        Z_lcfs_ordered = Z_lcfs_ordered[::-1]
    if Z_wall_ordered[1] < Z_wall_ordered[0]:
        R_wall_ordered = R_wall_ordered[::-1]
        Z_wall_ordered = Z_wall_ordered[::-1]

    t_lcfs = np.arange(n_lcfs) / n_lcfs
    t_wall = np.arange(n_wall) / n_wall

    t_wall_ext = np.r_[t_wall - 1.0, t_wall, t_wall + 1.0]
    R_wall_ext = np.r_[R_wall_ordered, R_wall_ordered, R_wall_ordered]
    Z_wall_ext = np.r_[Z_wall_ordered, Z_wall_ordered, Z_wall_ordered]

    R_matched_ordered = np.interp(t_lcfs, t_wall_ext, R_wall_ext)
    Z_matched_ordered = np.interp(t_lcfs, t_wall_ext, Z_wall_ext)

    reverse_order = np.argsort(lcfs_order)
    R_matched = R_matched_ordered[reverse_order]
    Z_matched = Z_matched_ordered[reverse_order]

    return R_matched, Z_matched


def write_chartmap_from_vmec_to_wall(
    wout_path: str | Path,
    out_path: str | Path,
    wall_rz: list[tuple[np.ndarray, np.ndarray]],
    wall_zeta: np.ndarray,
    *,
    nrho: int = 33,
    ntheta: int = 65,
    rho_lcfs: float = 0.8,
    num_field_periods: int | None = None,
    use_asym: bool = True,
) -> None:
    """
    Generate a chartmap from VMEC coordinates extended to wall boundaries.

    This creates a boundary-fitted coordinate system where:
    - rho in [0, rho_lcfs]: VMEC coordinates with rho = sqrt(s) * rho_lcfs
    - rho in [rho_lcfs, 1]: Cubic Hermite interpolation from LCFS to wall
    - rho = 1: Exactly on the wall boundary

    The coordinate system has C1 continuity at the LCFS.

    Parameters
    ----------
    wout_path : path to VMEC wout NetCDF file
    out_path : output chartmap NetCDF path
    wall_rz : list of (R, Z) arrays defining wall boundary at each zeta
              Each R, Z should be a closed curve (first point != last point)
    wall_zeta : array of zeta values corresponding to wall_rz curves
    nrho : number of radial grid points
    ntheta : number of poloidal grid points
    rho_lcfs : radial location of LCFS in (0, 1)
    num_field_periods : override nfp from wout file
    use_asym : include asymmetric Fourier terms if available

    Notes
    -----
    The wall curves are matched to VMEC theta using arc-length parameterization.
    Both the LCFS and wall curves are parameterized by fractional arc length,
    and corresponding points are connected by cubic Hermite interpolation.
    """
    from netCDF4 import Dataset

    from .vmec import VMECGeometry

    wout = _as_path(wout_path)
    out = _as_path(out_path)

    if not (0.0 < rho_lcfs < 1.0):
        raise ValueError("rho_lcfs must be in (0, 1)")

    wall_zeta = np.asarray(wall_zeta, dtype=float)
    if len(wall_rz) != len(wall_zeta):
        raise ValueError("wall_rz and wall_zeta must have same length")

    geom = VMECGeometry.from_file(str(wout))

    if num_field_periods is None:
        with Dataset(wout, "r") as ds:
            if "nfp" not in ds.variables:
                raise ValueError("wout file missing nfp variable")
            nfp = int(np.array(ds.variables["nfp"][...]))
    else:
        nfp = int(num_field_periods)

    grid = ChartmapGrid(
        rho=np.linspace(0.0, 1.0, nrho, dtype=float),
        theta=np.linspace(0.0, 2.0 * np.pi, ntheta, endpoint=False, dtype=float),
        zeta=wall_zeta,
        num_field_periods=nfp,
    )

    ir_lcfs = int(np.searchsorted(grid.rho, rho_lcfs, side="right") - 1)
    ir_lcfs = max(1, min(nrho - 2, ir_lcfs))
    rho_lcfs_actual = float(grid.rho[ir_lcfs])

    x = np.zeros((grid.rho.size, grid.theta.size, grid.zeta.size), dtype=np.float64)
    y = np.zeros_like(x)
    z = np.zeros_like(x)

    for iz, phi_val in enumerate(grid.zeta):
        phi_val = float(phi_val)

        R_wall, Z_wall = wall_rz[iz]
        R_wall = np.asarray(R_wall, dtype=float)
        Z_wall = np.asarray(Z_wall, dtype=float)

        R_lcfs, Z_lcfs, dR_ds, dZ_ds = geom.coords_s_with_deriv(
            1.0, grid.theta, phi_val, use_asym=use_asym
        )

        R_wall_matched, Z_wall_matched = _match_wall_to_lcfs_by_fraction(
            R_lcfs, Z_lcfs, R_wall, Z_wall
        )

        ds_drho_at_lcfs = 2.0 / rho_lcfs_actual
        drho_interval = 1.0 - rho_lcfs_actual
        m0_R = dR_ds * ds_drho_at_lcfs * drho_interval
        m0_Z = dZ_ds * ds_drho_at_lcfs * drho_interval

        chord_R = R_wall_matched - R_lcfs
        chord_Z = Z_wall_matched - Z_lcfs
        m1_R = chord_R
        m1_Z = chord_Z

        for ir in range(ir_lcfs + 1):
            rho_val = grid.rho[ir]
            s_val = (float(rho_val) / rho_lcfs_actual) ** 2
            R_in, Z_in, _ = geom.coords_s(s_val, grid.theta, phi_val, use_asym=use_asym)
            x[ir, :, iz] = R_in * np.cos(phi_val) * 100.0
            y[ir, :, iz] = R_in * np.sin(phi_val) * 100.0
            z[ir, :, iz] = Z_in * 100.0

        n_outer = grid.rho.size - ir_lcfs - 1
        if n_outer > 0:
            rho_outer = grid.rho[ir_lcfs + 1:]
            t_param = (rho_outer - rho_lcfs_actual) / (1.0 - rho_lcfs_actual)

            for ith in range(grid.theta.size):
                R_interp = _cubic_hermite(
                    t_param,
                    R_lcfs[ith],
                    m0_R[ith],
                    R_wall_matched[ith],
                    m1_R[ith],
                )
                Z_interp = _cubic_hermite(
                    t_param,
                    Z_lcfs[ith],
                    m0_Z[ith],
                    Z_wall_matched[ith],
                    m1_Z[ith],
                )
                x[ir_lcfs + 1:, ith, iz] = R_interp * np.cos(phi_val) * 100.0
                y[ir_lcfs + 1:, ith, iz] = R_interp * np.sin(phi_val) * 100.0
                z[ir_lcfs + 1:, ith, iz] = Z_interp * 100.0

    write_chartmap_netcdf(
        out,
        grid=grid,
        x_rtz=x,
        y_rtz=y,
        z_rtz=z,
        zeta_convention="cyl",
        rho_convention="vmec_to_wall",
    )


def write_chartmap_from_vmec_and_stl(
    wout_path: str | Path,
    stl_path: str | Path,
    out_path: str | Path,
    *,
    nrho: int = 33,
    ntheta: int = 65,
    nzeta: int = 33,
    rho_lcfs: float = 0.8,
    n_boundary_points: int = 512,
    stitch_tol: float = 1.0e-6,
    num_field_periods: int | None = None,
    use_asym: bool = True,
) -> None:
    """
    Generate a chartmap from VMEC extended to STL wall boundaries.

    Convenience function that combines STL boundary extraction with
    write_chartmap_from_vmec_to_wall.

    Parameters
    ----------
    wout_path : path to VMEC wout NetCDF file
    stl_path : path to STL wall geometry file
    out_path : output chartmap NetCDF path
    nrho : number of radial grid points
    ntheta : number of poloidal grid points
    nzeta : number of toroidal slices
    rho_lcfs : radial location of LCFS in (0, 1)
    n_boundary_points : points for STL boundary extraction
    stitch_tol : tolerance for stitching open contours
    num_field_periods : override nfp from wout file
    use_asym : include asymmetric Fourier terms if available
    """
    from netCDF4 import Dataset

    from .stl_boundary import extract_boundary_slices

    wout = _as_path(wout_path)
    stl = _as_path(stl_path)

    if num_field_periods is None:
        with Dataset(wout, "r") as ds:
            if "nfp" not in ds.variables:
                raise ValueError("wout file missing nfp variable")
            nfp = int(np.array(ds.variables["nfp"][...]))
    else:
        nfp = int(num_field_periods)

    period = _zeta_period(nfp)
    phi_vals = np.linspace(0.0, period, nzeta, endpoint=False, dtype=float)

    slices = extract_boundary_slices(
        stl,
        n_phi=nzeta,
        n_boundary_points=n_boundary_points,
        stitch_tol=stitch_tol,
        phi_vals=phi_vals,
    )

    wall_rz = []
    wall_zeta = []
    for s in slices:
        R = s.outer_filled[:-1, 0]
        Z = s.outer_filled[:-1, 1]
        wall_rz.append((R, Z))
        wall_zeta.append(s.phi)

    write_chartmap_from_vmec_to_wall(
        wout,
        out_path,
        wall_rz,
        np.array(wall_zeta),
        nrho=nrho,
        ntheta=ntheta,
        rho_lcfs=rho_lcfs,
        num_field_periods=nfp,
        use_asym=use_asym,
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
