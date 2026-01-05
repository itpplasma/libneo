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

        if boundary_offset != 0.0:
            theta_base = np.linspace(0.0, 2.0 * np.pi, 4096, endpoint=False, dtype=float)
            R_base, Z_base, _ = geom.coords_s(float(s_boundary), theta_base, phi_val, use_asym=use_asym)

            dtheta = float(theta_base[1] - theta_base[0])
            dR = (np.roll(R_base, -1) - np.roll(R_base, 1)) / (2.0 * dtheta)
            dZ = (np.roll(Z_base, -1) - np.roll(Z_base, 1)) / (2.0 * dtheta)
            area = 0.5 * float(np.sum(R_base * np.roll(Z_base, -1) - np.roll(R_base, -1) * Z_base))

            if area > 0.0:
                nR = dZ
                nZ = -dR
            else:
                nR = -dZ
                nZ = dR

            nrm = np.sqrt(nR**2 + nZ**2)
            nrm = np.where(nrm == 0.0, 1.0, nrm)
            nR = nR / nrm
            nZ = nZ / nrm

            R_off = R_base + boundary_offset * nR
            Z_off = Z_base + boundary_offset * nZ

            cR = np.fft.rfft(R_off) / float(theta_base.size)
            cZ = np.fft.rfft(Z_off) / float(theta_base.size)
            kmax = min(int(M), int(cR.size) - 1)

            def eval_series(th: np.ndarray, c: np.ndarray) -> np.ndarray:
                theta_eval = np.asarray(th, dtype=float)
                out = np.full(theta_eval.shape, float(np.real(c[0])), dtype=float)
                for k in range(1, kmax + 1):
                    out = out + 2.0 * np.real(c[k] * np.exp(1j * float(k) * theta_eval))
                return out

            def curve(t: np.ndarray) -> np.ndarray:
                tt = np.asarray(t, dtype=float)
                return np.array([eval_series(tt, cR), eval_series(tt, cZ)])
        else:

            def curve(t: np.ndarray) -> np.ndarray:
                th = np.asarray(t, dtype=float)
                R, Zc, _ = geom.coords_s(float(s_boundary), th, phi_val, use_asym=use_asym)
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
        attrs={
            "chartmap_source": "vmec_boundary_map2disc",
            "s_boundary": float(s_boundary),
            "boundary_offset_m": float(boundary_offset),
            "map2disc_M": int(M),
            "map2disc_Nt": int(Nt),
            "map2disc_Ng_r": int(Ng[0]),
            "map2disc_Ng_t": int(Ng[1]),
        },
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
