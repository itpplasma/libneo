"""Tests for VMEC chartmap with cubic Hermite extension beyond LCFS."""
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


def _plot_vmec_extended_grid(
    chartmap_path: Path,
    wout_path: Path,
    out_path: Path,
    ir_lcfs: int,
    izeta: int = 0,
) -> None:
    """Plot R,Z grid showing VMEC interior (blue) and Hermite extension (orange)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from netCDF4 import Dataset

    from libneo.vmec import VMECGeometry

    with Dataset(chartmap_path, "r") as ds:
        rho_grid = np.array(ds.variables["rho"][:], dtype=float)
        theta_grid = np.array(ds.variables["theta"][:], dtype=float)
        zeta_grid = np.array(ds.variables["zeta"][:], dtype=float)
        x = np.array(ds.variables["x"][:], dtype=float)
        y = np.array(ds.variables["y"][:], dtype=float)
        z = np.array(ds.variables["z"][:], dtype=float)
        units = str(getattr(ds.variables["x"], "units", "cm")).strip().lower()
        scale = 0.01 if units == "cm" else 1.0

    R = np.sqrt((x * scale) ** 2 + (y * scale) ** 2)
    Z = z * scale

    geom = VMECGeometry.from_file(str(wout_path))

    fig, ax = plt.subplots(figsize=(8.0, 8.0))

    zeta = float(zeta_grid[izeta])
    nrho = len(rho_grid)
    ntheta = len(theta_grid)

    for ir in range(nrho):
        color = "C0" if ir <= ir_lcfs else "C1"
        lw = 1.5 if ir == ir_lcfs else 0.6
        r_ring = R[izeta, :, ir]
        z_ring = Z[izeta, :, ir]
        ax.plot(
            np.r_[r_ring, r_ring[0]],
            np.r_[z_ring, z_ring[0]],
            color=color,
            linewidth=lw,
            alpha=0.9 if ir == ir_lcfs else 0.7,
        )

    for ith in range(0, ntheta, max(1, ntheta // 24)):
        ax.plot(
            R[izeta, ith, : ir_lcfs + 1],
            Z[izeta, ith, : ir_lcfs + 1],
            color="C0",
            linewidth=0.4,
            alpha=0.5,
        )
        ax.plot(
            R[izeta, ith, ir_lcfs:],
            Z[izeta, ith, ir_lcfs:],
            color="C1",
            linewidth=0.4,
            alpha=0.5,
        )

    theta_dense = np.linspace(0.0, 2.0 * np.pi, 256, endpoint=False)
    R_lcfs, Z_lcfs, _ = geom.coords_s(1.0, theta_dense, zeta, use_asym=True)
    ax.plot(R_lcfs, Z_lcfs, "k--", linewidth=1.0, label="VMEC LCFS")

    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("R [m]")
    ax.set_ylabel("Z [m]")
    ax.set_title(f"VMEC + cubic Hermite extension (zeta={zeta:.3f})")

    legend_elements = [
        Line2D([0], [0], color="C0", linewidth=2, label="VMEC region"),
        Line2D([0], [0], color="C1", linewidth=2, label="Hermite extension"),
        Line2D([0], [0], color="k", linewidth=1, linestyle="--", label="VMEC LCFS"),
    ]
    ax.legend(handles=legend_elements, loc="upper right")

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def _plot_derivative_continuity(
    chartmap_path: Path,
    out_path: Path,
    ir_lcfs: int,
    izeta: int = 0,
) -> None:
    """Plot dR/drho and dZ/drho to verify C1 continuity at LCFS."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from netCDF4 import Dataset

    with Dataset(chartmap_path, "r") as ds:
        rho_grid = np.array(ds.variables["rho"][:], dtype=float)
        x = np.array(ds.variables["x"][:], dtype=float)
        y = np.array(ds.variables["y"][:], dtype=float)
        z = np.array(ds.variables["z"][:], dtype=float)
        units = str(getattr(ds.variables["x"], "units", "cm")).strip().lower()
        scale = 0.01 if units == "cm" else 1.0

    R = np.sqrt((x * scale) ** 2 + (y * scale) ** 2)
    Z = z * scale

    ith_sample = R.shape[1] // 4

    R_line = R[izeta, ith_sample, :]
    Z_line = Z[izeta, ith_sample, :]

    dR_drho = np.gradient(R_line, rho_grid)
    dZ_drho = np.gradient(Z_line, rho_grid)

    fig, axes = plt.subplots(2, 2, figsize=(10, 8))

    axes[0, 0].plot(rho_grid, R_line, "b-", linewidth=1.5)
    axes[0, 0].axvline(rho_grid[ir_lcfs], color="k", linestyle="--", alpha=0.5)
    axes[0, 0].set_xlabel("rho")
    axes[0, 0].set_ylabel("R [m]")
    axes[0, 0].set_title("R vs rho")

    axes[0, 1].plot(rho_grid, Z_line, "b-", linewidth=1.5)
    axes[0, 1].axvline(rho_grid[ir_lcfs], color="k", linestyle="--", alpha=0.5)
    axes[0, 1].set_xlabel("rho")
    axes[0, 1].set_ylabel("Z [m]")
    axes[0, 1].set_title("Z vs rho")

    axes[1, 0].plot(rho_grid, dR_drho, "r-", linewidth=1.5)
    axes[1, 0].axvline(rho_grid[ir_lcfs], color="k", linestyle="--", alpha=0.5)
    axes[1, 0].set_xlabel("rho")
    axes[1, 0].set_ylabel("dR/drho [m]")
    axes[1, 0].set_title("dR/drho (should be continuous at LCFS)")

    axes[1, 1].plot(rho_grid, dZ_drho, "r-", linewidth=1.5)
    axes[1, 1].axvline(rho_grid[ir_lcfs], color="k", linestyle="--", alpha=0.5)
    axes[1, 1].set_xlabel("rho")
    axes[1, 1].set_ylabel("dZ/drho [m]")
    axes[1, 1].set_title("dZ/drho (should be continuous at LCFS)")

    fig.suptitle(f"Derivative continuity check (theta index={ith_sample})")
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


@pytest.mark.network
def test_chartmap_vmec_extended_real_wout(tmp_path: Path) -> None:
    """Test VMEC extended mode with real NCSX data and visual verification."""
    from netCDF4 import Dataset

    from libneo.chartmap import write_chartmap_from_vmec_extended
    from libneo.vmec import VMECGeometry

    wout_url = (
        "https://princetonuniversity.github.io/STELLOPT/examples/wout_ncsx_c09r00_fixed.nc"
    )
    wout = tmp_path / "wout_ncsx.nc"
    _download(wout_url, wout)

    nrho = 33
    nrho_vmec = 25
    ntheta = 65
    nzeta = 9
    boundary_offset = 0.15

    chartmap = tmp_path / "chartmap_vmec_extended.nc"
    write_chartmap_from_vmec_extended(
        wout,
        chartmap,
        nrho=nrho,
        ntheta=ntheta,
        nzeta=nzeta,
        nrho_vmec=nrho_vmec,
        boundary_offset=boundary_offset,
    )

    geom = VMECGeometry.from_file(str(wout))

    with Dataset(chartmap, "r") as ds:
        rho_grid = np.array(ds.variables["rho"][:], dtype=float)
        theta_grid = np.array(ds.variables["theta"][:], dtype=float)
        zeta_grid = np.array(ds.variables["zeta"][:], dtype=float)
        x = np.array(ds.variables["x"][:], dtype=float)
        y = np.array(ds.variables["y"][:], dtype=float)
        z = np.array(ds.variables["z"][:], dtype=float)
        rho_conv = str(ds.getncattr("rho_convention"))
        units = str(getattr(ds.variables["x"], "units", "cm")).strip().lower()
        scale = 0.01 if units == "cm" else 1.0

    assert rho_conv == "vmec_extended"

    R = np.sqrt((x * scale) ** 2 + (y * scale) ** 2)
    Z = z * scale

    ir_lcfs = nrho_vmec - 1
    rho_lcfs = rho_grid[ir_lcfs]

    for iz in [0, len(zeta_grid) // 2]:
        zeta = float(zeta_grid[iz])
        theta = theta_grid

        for ir in range(ir_lcfs + 1):
            rho_val = rho_grid[ir]
            s_val = (float(rho_val) / rho_lcfs) ** 2
            R_vmec, Z_vmec, _ = geom.coords_s(s_val, theta, zeta, use_asym=True)
            assert np.allclose(R_vmec, R[iz, :, ir], rtol=0.0, atol=1.0e-6), (
                f"VMEC region mismatch at iz={iz}, ir={ir}"
            )
            assert np.allclose(Z_vmec, Z[iz, :, ir], rtol=0.0, atol=1.0e-6), (
                f"VMEC region mismatch at iz={iz}, ir={ir}"
            )

    R_lcfs_stored = R[0, :, ir_lcfs]
    R_outer_stored = R[0, :, -1]
    mean_extension = float(np.mean(np.abs(R_outer_stored - R_lcfs_stored)))
    assert mean_extension > 0.05, f"Extension too small: {mean_extension:.3f} m"
    assert mean_extension < 0.5, f"Extension too large: {mean_extension:.3f} m"

    _plot_vmec_extended_grid(
        chartmap,
        wout,
        _figure_path("chartmap_vmec_extended_ncsx.png"),
        ir_lcfs=ir_lcfs,
        izeta=0,
    )
    _plot_vmec_extended_grid(
        chartmap,
        wout,
        _figure_path("chartmap_vmec_extended_ncsx_phi2.png"),
        ir_lcfs=ir_lcfs,
        izeta=len(zeta_grid) // 2,
    )
    _plot_derivative_continuity(
        chartmap,
        _figure_path("chartmap_vmec_extended_derivatives.png"),
        ir_lcfs=ir_lcfs,
        izeta=0,
    )


@pytest.mark.network
def test_chartmap_vmec_extended_c1_continuity(tmp_path: Path) -> None:
    """Test that derivatives are continuous at LCFS (C1 continuity)."""
    from netCDF4 import Dataset

    from libneo.chartmap import write_chartmap_from_vmec_extended

    wout_url = (
        "https://princetonuniversity.github.io/STELLOPT/examples/wout_ncsx_c09r00_fixed.nc"
    )
    wout = tmp_path / "wout_ncsx.nc"
    _download(wout_url, wout)

    nrho = 65
    nrho_vmec = 49
    ntheta = 65
    nzeta = 5

    chartmap = tmp_path / "chartmap_continuity.nc"
    write_chartmap_from_vmec_extended(
        wout,
        chartmap,
        nrho=nrho,
        ntheta=ntheta,
        nzeta=nzeta,
        nrho_vmec=nrho_vmec,
        boundary_offset=0.12,
    )

    with Dataset(chartmap, "r") as ds:
        rho_grid = np.array(ds.variables["rho"][:], dtype=float)
        x = np.array(ds.variables["x"][:], dtype=float)
        y = np.array(ds.variables["y"][:], dtype=float)
        z = np.array(ds.variables["z"][:], dtype=float)
        units = str(getattr(ds.variables["x"], "units", "cm")).strip().lower()
        scale = 0.01 if units == "cm" else 1.0

    R = np.sqrt((x * scale) ** 2 + (y * scale) ** 2)
    Z = z * scale

    ir_lcfs = nrho_vmec - 1

    for iz in range(Z.shape[0]):
        for ith in range(0, Z.shape[1], 8):
            R_line = R[iz, ith, :]
            Z_line = Z[iz, ith, :]

            dR_drho = np.gradient(R_line, rho_grid)
            dZ_drho = np.gradient(Z_line, rho_grid)

            dR_jump = abs(dR_drho[ir_lcfs + 1] - dR_drho[ir_lcfs])
            dZ_jump = abs(dZ_drho[ir_lcfs + 1] - dZ_drho[ir_lcfs])

            dR_scale = max(abs(dR_drho[ir_lcfs]), abs(dR_drho[ir_lcfs + 1]), 0.01)
            dZ_scale = max(abs(dZ_drho[ir_lcfs]), abs(dZ_drho[ir_lcfs + 1]), 0.01)

            assert dR_jump / dR_scale < 2.0, (
                f"dR/drho discontinuity at iz={iz}, ith={ith}: jump={dR_jump:.4f}, "
                f"scale={dR_scale:.4f}, ratio={dR_jump/dR_scale:.2f}"
            )
            assert dZ_jump / dZ_scale < 2.0, (
                f"dZ/drho discontinuity at iz={iz}, ith={ith}: jump={dZ_jump:.4f}, "
                f"scale={dZ_scale:.4f}, ratio={dZ_jump/dZ_scale:.2f}"
            )
