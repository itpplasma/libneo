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
    - Conventions: zeta_convention=cyl, rho_convention=unknown (geometric disk radius)
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


def write_chartmap_from_coils_offset_surface(
    coils_path: str | Path,
    out_path: str | Path,
    *,
    offset_cm: float,
    nrho: int = 33,
    ntheta: int = 65,
    nphi: int = 32,
    num_field_periods: int = 1,
    grid_shape: tuple[int, int, int] = (72, 72, 72),
    padding_cm: float = 20.0,
    sample_step_cm: float | None = None,
    axis_xy: tuple[float, float] | None = None,
    seed_rz: tuple[float, float] | None = None,
    window_r_quantiles: tuple[float, float] = (0.0, 1.0),
    window_z_quantiles: tuple[float, float] = (0.0, 1.0),
    smooth_window: int = 11,
    n_boundary_points: int = 512,
    stitch_tol: float = 1.0e-6,
    M: int = 16,
    Nt: int = 256,
    Ng: tuple[int, int] = (256, 256),
) -> None:
    """
    Generate a chartmap NetCDF from a coils file by building an offset surface.

    The coils file must be in SIMPLE format (see neo_biotsavart).
    """
    import trimesh

    from .coils_simple import read_simple_coils, coils_to_segments
    from .coils_offset_surface import build_inner_offset_surface_from_segments
    from .stl_boundary import extract_boundary_slices_from_mesh, write_chartmap_from_stl as _write_from_slices

    coils_file = _as_path(coils_path)
    out = _as_path(out_path)

    if offset_cm <= 0.0:
        raise ValueError("offset_cm must be > 0")
    if num_field_periods < 1:
        raise ValueError("num_field_periods must be >= 1")
    if nphi < 2:
        raise ValueError("nphi must be >= 2")

    coils = read_simple_coils(coils_file)
    a, b = coils_to_segments(coils)

    if axis_xy is None:
        pts = np.vstack(coils.coils)
        axis_xy = (float(np.mean(pts[:, 0])), float(np.mean(pts[:, 1])))

    a = a.copy()
    b = b.copy()
    a[:, 0] -= float(axis_xy[0])
    a[:, 1] -= float(axis_xy[1])
    b[:, 0] -= float(axis_xy[0])
    b[:, 1] -= float(axis_xy[1])

    if seed_rz is None:
        pts0 = np.vstack(coils.coils)
        rr = np.sqrt((pts0[:, 0] - float(axis_xy[0])) ** 2 + (pts0[:, 1] - float(axis_xy[1])) ** 2)
        R0 = float(np.mean(rr))
        z0 = float(np.mean(pts0[:, 2]))
        seed_rz = (R0, z0)

    seed = np.array([float(seed_rz[0]), 0.0, float(seed_rz[1])], dtype=float)

    offset_m = float(offset_cm) / 100.0
    padding_m = float(padding_cm) / 100.0
    sample_step_m = None if sample_step_cm is None else float(sample_step_cm) / 100.0

    surf = build_inner_offset_surface_from_segments(
        a,
        b,
        offset_m=offset_m,
        seed=seed,
        grid_shape=grid_shape,
        padding_m=padding_m,
        sample_step_m=sample_step_m,
        window_r_quantiles=window_r_quantiles,
        window_z_quantiles=window_z_quantiles,
    )

    mesh = trimesh.Trimesh(vertices=surf.vertices, faces=surf.faces, process=True)
    mesh.merge_vertices()
    mesh.remove_unreferenced_vertices()
    trimesh.repair.fix_normals(mesh)
    trimesh.repair.fill_holes(mesh)
    trimesh.smoothing.filter_taubin(mesh, lamb=0.5, nu=-0.53, iterations=12)

    if not mesh.is_watertight:
        raise RuntimeError("offset surface mesh is not watertight; cannot form closed surface")

    period = _zeta_period(int(num_field_periods))
    phi_vals = np.linspace(0.0, period, int(nphi), endpoint=False, dtype=float)

    slices = extract_boundary_slices_from_mesh(
        mesh,
        n_phi=int(nphi),
        n_boundary_points=int(n_boundary_points),
        axis_xy=(0.0, 0.0),
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
        smooth_window=int(smooth_window),
    )
