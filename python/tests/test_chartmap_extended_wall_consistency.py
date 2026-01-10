"""Integration tests verifying extended and wall chartmaps are consistent inside VMEC region.

Key finding: Extended (cubic Hermite) and Wall (quintic Hermite) chartmaps should give
IDENTICAL (R, Z) coordinates inside the VMEC region (rho < rho_lcfs). They only differ
in the extrapolation region beyond LCFS.

This test verifies that the coordinate conversion is correct and that any differences
in particle tracing results are due to dynamics in the extrapolation region, not bugs
in the coordinate system.
"""
from __future__ import annotations

from pathlib import Path
from urllib.request import urlopen

import numpy as np
import pytest


def _download(url: str, out_path: Path) -> None:
    with urlopen(url, timeout=60) as resp, open(out_path, "wb") as f:
        f.write(resp.read())


def _figure_path(name: str) -> Path:
    repo_root = Path(__file__).resolve().parents[2]
    target = repo_root / "build" / "test" / "python"
    target.mkdir(parents=True, exist_ok=True)
    return target / name


def _create_synthetic_wall(
    geom, zeta: float, offset: float, ntheta: int = 512
) -> tuple[np.ndarray, np.ndarray]:
    """Create a synthetic wall boundary using shapely buffer."""
    from shapely.geometry import Polygon

    theta = np.linspace(0.0, 2.0 * np.pi, 4096, endpoint=False)
    R_lcfs, Z_lcfs, _ = geom.coords_s(1.0, theta, zeta, use_asym=True)

    poly = Polygon(np.column_stack([R_lcfs, Z_lcfs]))
    if not poly.is_valid:
        poly = poly.buffer(0.0)

    off = poly.buffer(offset, join_style="round")
    if off.geom_type != "Polygon":
        off = max(list(off.geoms), key=lambda g: g.area)  # type: ignore[union-attr]

    coords = np.asarray(off.exterior.coords, dtype=float)
    coords_open = coords[:-1, :]

    i0 = int(np.argmax(coords_open[:, 0]))
    coords_open = np.vstack([coords_open[i0:, :], coords_open[:i0, :]])
    if coords_open.shape[0] >= 2 and coords_open[1, 1] - coords_open[0, 1] < 0.0:
        coords_open = coords_open[::-1, :]

    seg = np.sqrt(np.sum(np.diff(coords_open, axis=0) ** 2, axis=1))
    s_coords = np.concatenate(([0.0], np.cumsum(seg)))
    length = s_coords[-1]

    s_query = np.linspace(0.0, length, ntheta, endpoint=False)
    R_wall = np.interp(s_query, s_coords, coords_open[:, 0])
    Z_wall = np.interp(s_query, s_coords, coords_open[:, 1])

    return R_wall, Z_wall


@pytest.mark.network
def test_extended_wall_identical_inside_vmec(tmp_path: Path) -> None:
    """Verify extended and wall chartmaps give identical R,Z inside VMEC region.

    This is a key integration test that confirms the coordinate systems are
    consistent where they should be (inside VMEC region), and only differ
    where expected (beyond LCFS in the extrapolation region).
    """
    from netCDF4 import Dataset

    from libneo.chartmap import (
        write_chartmap_from_vmec_extended,
        write_chartmap_from_vmec_to_wall,
    )
    from libneo.vmec import VMECGeometry

    wout_url = (
        "https://princetonuniversity.github.io/STELLOPT/examples/wout_ncsx_c09r00_fixed.nc"
    )
    wout = tmp_path / "wout_ncsx.nc"
    _download(wout_url, wout)

    geom = VMECGeometry.from_file(str(wout))

    with Dataset(wout, "r") as ds:
        nfp = int(np.array(ds.variables["nfp"][...]))

    nrho = 33
    ntheta = 65
    nzeta = 5
    rho_lcfs = 0.75
    wall_offset = 0.15

    ext_chartmap = tmp_path / "chartmap_extended.nc"
    write_chartmap_from_vmec_extended(
        wout,
        ext_chartmap,
        nrho=nrho,
        ntheta=ntheta,
        nzeta=nzeta,
        nrho_vmec=int(rho_lcfs * nrho),
        boundary_offset=wall_offset,
    )

    period = 2.0 * np.pi / nfp
    wall_zeta = np.linspace(0.0, period, nzeta, endpoint=False, dtype=float)
    wall_rz = []
    for zeta in wall_zeta:
        R_wall, Z_wall = _create_synthetic_wall(geom, float(zeta), wall_offset)
        wall_rz.append((R_wall, Z_wall))

    wall_chartmap = tmp_path / "chartmap_wall.nc"
    write_chartmap_from_vmec_to_wall(
        wout,
        wall_chartmap,
        wall_rz,
        wall_zeta,
        nrho=nrho,
        ntheta=ntheta,
        rho_lcfs=rho_lcfs,
        num_field_periods=nfp,
    )

    with Dataset(ext_chartmap, "r") as ds:
        rho_ext = np.array(ds.variables["rho"][:], dtype=float)
        x_ext = np.array(ds.variables["x"][:], dtype=float)
        y_ext = np.array(ds.variables["y"][:], dtype=float)
        z_ext = np.array(ds.variables["z"][:], dtype=float)
        units_ext = str(getattr(ds.variables["x"], "units", "cm")).strip().lower()

    with Dataset(wall_chartmap, "r") as ds:
        rho_wall = np.array(ds.variables["rho"][:], dtype=float)
        x_wall = np.array(ds.variables["x"][:], dtype=float)
        y_wall = np.array(ds.variables["y"][:], dtype=float)
        z_wall = np.array(ds.variables["z"][:], dtype=float)
        units_wall = str(getattr(ds.variables["x"], "units", "cm")).strip().lower()

    scale_ext = 0.01 if units_ext == "cm" else 1.0
    scale_wall = 0.01 if units_wall == "cm" else 1.0

    R_ext = np.sqrt((x_ext * scale_ext) ** 2 + (y_ext * scale_ext) ** 2)
    Z_ext = z_ext * scale_ext
    R_wall = np.sqrt((x_wall * scale_wall) ** 2 + (y_wall * scale_wall) ** 2)
    Z_wall = z_wall * scale_wall

    ir_lcfs_ext = int(np.searchsorted(rho_ext, rho_lcfs, side="right") - 1)
    ir_lcfs_wall = int(np.searchsorted(rho_wall, rho_lcfs, side="right") - 1)

    for iz in range(R_ext.shape[0]):
        for ith in range(R_ext.shape[1]):
            for ir in range(min(ir_lcfs_ext, ir_lcfs_wall) + 1):
                R_e = R_ext[iz, ith, ir]
                Z_e = Z_ext[iz, ith, ir]
                R_w = R_wall[iz, ith, ir]
                Z_w = Z_wall[iz, ith, ir]

                diff_R = abs(R_e - R_w)
                diff_Z = abs(Z_e - Z_w)

                assert diff_R < 1.0e-6, (
                    f"R mismatch inside VMEC at iz={iz}, ith={ith}, ir={ir}: "
                    f"ext={R_e:.6f}, wall={R_w:.6f}, diff={diff_R:.2e} m"
                )
                assert diff_Z < 1.0e-6, (
                    f"Z mismatch inside VMEC at iz={iz}, ith={ith}, ir={ir}: "
                    f"ext={Z_e:.6f}, wall={Z_w:.6f}, diff={diff_Z:.2e} m"
                )


@pytest.mark.network
def test_extended_wall_differ_beyond_lcfs(tmp_path: Path) -> None:
    """Verify extended and wall chartmaps differ as expected beyond LCFS.

    This confirms that the extrapolation methods (cubic vs quintic Hermite)
    produce different coordinates beyond the LCFS boundary, which is the
    expected and correct behavior.
    """
    from netCDF4 import Dataset

    from libneo.chartmap import (
        write_chartmap_from_vmec_extended,
        write_chartmap_from_vmec_to_wall,
    )
    from libneo.vmec import VMECGeometry

    wout_url = (
        "https://princetonuniversity.github.io/STELLOPT/examples/wout_ncsx_c09r00_fixed.nc"
    )
    wout = tmp_path / "wout_ncsx.nc"
    _download(wout_url, wout)

    geom = VMECGeometry.from_file(str(wout))

    with Dataset(wout, "r") as ds:
        nfp = int(np.array(ds.variables["nfp"][...]))

    nrho = 65
    ntheta = 65
    nzeta = 5
    rho_lcfs = 0.75
    wall_offset = 0.15

    ext_chartmap = tmp_path / "chartmap_extended.nc"
    write_chartmap_from_vmec_extended(
        wout,
        ext_chartmap,
        nrho=nrho,
        ntheta=ntheta,
        nzeta=nzeta,
        nrho_vmec=int(rho_lcfs * nrho),
        boundary_offset=wall_offset,
    )

    period = 2.0 * np.pi / nfp
    wall_zeta = np.linspace(0.0, period, nzeta, endpoint=False, dtype=float)
    wall_rz = []
    for zeta in wall_zeta:
        R_wall, Z_wall = _create_synthetic_wall(geom, float(zeta), wall_offset)
        wall_rz.append((R_wall, Z_wall))

    wall_chartmap = tmp_path / "chartmap_wall.nc"
    write_chartmap_from_vmec_to_wall(
        wout,
        wall_chartmap,
        wall_rz,
        wall_zeta,
        nrho=nrho,
        ntheta=ntheta,
        rho_lcfs=rho_lcfs,
        num_field_periods=nfp,
    )

    with Dataset(ext_chartmap, "r") as ds:
        x_ext = np.array(ds.variables["x"][:], dtype=float)
        y_ext = np.array(ds.variables["y"][:], dtype=float)
        z_ext = np.array(ds.variables["z"][:], dtype=float)
        units_ext = str(getattr(ds.variables["x"], "units", "cm")).strip().lower()

    with Dataset(wall_chartmap, "r") as ds:
        x_wall = np.array(ds.variables["x"][:], dtype=float)
        y_wall = np.array(ds.variables["y"][:], dtype=float)
        z_wall = np.array(ds.variables["z"][:], dtype=float)
        units_wall = str(getattr(ds.variables["x"], "units", "cm")).strip().lower()

    scale_ext = 0.01 if units_ext == "cm" else 1.0
    scale_wall = 0.01 if units_wall == "cm" else 1.0

    R_ext = np.sqrt((x_ext * scale_ext) ** 2 + (y_ext * scale_ext) ** 2)
    Z_ext = z_ext * scale_ext
    R_wall_arr = np.sqrt((x_wall * scale_wall) ** 2 + (y_wall * scale_wall) ** 2)
    Z_wall_arr = z_wall * scale_wall

    R_ext_outer = R_ext[:, :, -1]
    Z_ext_outer = Z_ext[:, :, -1]
    R_wall_outer = R_wall_arr[:, :, -1]
    Z_wall_outer = Z_wall_arr[:, :, -1]

    max_diff_R = np.max(np.abs(R_ext_outer - R_wall_outer))
    max_diff_Z = np.max(np.abs(Z_ext_outer - Z_wall_outer))

    assert max_diff_R > 0.001 or max_diff_Z > 0.001, (
        f"Extended and wall should differ at rho=1: max_diff_R={max_diff_R:.4f}m, "
        f"max_diff_Z={max_diff_Z:.4f}m"
    )


@pytest.mark.network
def test_extended_wall_visual_comparison(tmp_path: Path) -> None:
    """Generate visual comparison plot of extended vs wall chartmaps."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from netCDF4 import Dataset

    from libneo.chartmap import (
        write_chartmap_from_vmec_extended,
        write_chartmap_from_vmec_to_wall,
    )
    from libneo.vmec import VMECGeometry

    wout_url = (
        "https://princetonuniversity.github.io/STELLOPT/examples/wout_ncsx_c09r00_fixed.nc"
    )
    wout = tmp_path / "wout_ncsx.nc"
    _download(wout_url, wout)

    geom = VMECGeometry.from_file(str(wout))

    with Dataset(wout, "r") as ds:
        nfp = int(np.array(ds.variables["nfp"][...]))

    nrho = 33
    ntheta = 65
    nzeta = 5
    rho_lcfs = 0.75
    wall_offset = 0.15

    ext_chartmap = tmp_path / "chartmap_extended.nc"
    write_chartmap_from_vmec_extended(
        wout,
        ext_chartmap,
        nrho=nrho,
        ntheta=ntheta,
        nzeta=nzeta,
        nrho_vmec=int(rho_lcfs * nrho),
        boundary_offset=wall_offset,
    )

    period = 2.0 * np.pi / nfp
    wall_zeta = np.linspace(0.0, period, nzeta, endpoint=False, dtype=float)
    wall_rz = []
    for zeta in wall_zeta:
        R_wall_curve, Z_wall_curve = _create_synthetic_wall(
            geom, float(zeta), wall_offset
        )
        wall_rz.append((R_wall_curve, Z_wall_curve))

    wall_chartmap = tmp_path / "chartmap_wall.nc"
    write_chartmap_from_vmec_to_wall(
        wout,
        wall_chartmap,
        wall_rz,
        wall_zeta,
        nrho=nrho,
        ntheta=ntheta,
        rho_lcfs=rho_lcfs,
        num_field_periods=nfp,
    )

    def load_chartmap(path):
        with Dataset(path, "r") as ds:
            rho = np.array(ds.variables["rho"][:], dtype=float)
            x = np.array(ds.variables["x"][:], dtype=float)
            y = np.array(ds.variables["y"][:], dtype=float)
            z = np.array(ds.variables["z"][:], dtype=float)
            units = str(getattr(ds.variables["x"], "units", "cm")).strip().lower()
        scale = 0.01 if units == "cm" else 1.0
        R = np.sqrt((x * scale) ** 2 + (y * scale) ** 2)
        Z = z * scale
        return rho, R, Z

    rho_ext, R_ext, Z_ext = load_chartmap(ext_chartmap)
    _, R_wall, Z_wall = load_chartmap(wall_chartmap)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    iz = 0
    for ir in range(0, nrho, 2):
        axes[0].plot(R_ext[iz, :, ir], Z_ext[iz, :, ir], "b-", lw=0.5, alpha=0.7)
        axes[1].plot(R_wall[iz, :, ir], Z_wall[iz, :, ir], "r-", lw=0.5, alpha=0.7)

    axes[0].set_title("Extended (cubic Hermite)")
    axes[0].set_xlabel("R [m]")
    axes[0].set_ylabel("Z [m]")
    axes[0].set_aspect("equal")

    axes[1].set_title("Wall (quintic Hermite)")
    axes[1].set_xlabel("R [m]")
    axes[1].set_ylabel("Z [m]")
    axes[1].set_aspect("equal")

    diff = np.sqrt((R_ext - R_wall) ** 2 + (Z_ext - Z_wall) ** 2) * 1000
    im = axes[2].imshow(
        diff[iz, :, :].T,
        aspect="auto",
        origin="lower",
        extent=[0, 2 * np.pi, float(rho_ext[0]), float(rho_ext[-1])],
        cmap="hot",
    )
    axes[2].axhline(rho_lcfs, color="g", linestyle="--", label=f"rho_lcfs={rho_lcfs}")
    axes[2].set_title("Difference (mm)")
    axes[2].set_xlabel("theta [rad]")
    axes[2].set_ylabel("rho")
    axes[2].legend()
    plt.colorbar(im, ax=axes[2], label="|ΔR,Z| [mm]")

    fig.tight_layout()
    out_path = _figure_path("chartmap_extended_wall_comparison.png")
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
