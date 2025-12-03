#!/usr/bin/env python3
"""
Generate Babin (map2disc-based) coordinate overlays for visual inspection.

- GeoFlux vs Babin in the RZ plane from an EQDSK file
- VMEC vs Babin in the RZ plane from a sample VMEC equilibrium
"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path
import ssl
import urllib.request

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from map2disc import map as m2d  # noqa: E402

from libneo.vmec import VMECGeometry  # noqa: E402

GEQDSK_URL = "https://crppwww.epfl.ch/~sauter/benchmark/EQDSK_I"
SSL_CONTEXT = ssl._create_unverified_context()

STELLOPT_WOUT_URL = (
    "https://princetonuniversity.github.io/STELLOPT/examples/wout_ncsx_c09r00_fixed.nc"
)


def ensure_geqdsk(dest_dir: Path) -> Path:
    """Download or reuse a reference GEQDSK file."""
    target = dest_dir / "EQDSK_I.geqdsk"
    if target.exists() and target.stat().st_size > 0:
        return target
    dest_dir.mkdir(parents=True, exist_ok=True)
    tmp_path = target.with_suffix(target.suffix + ".tmp")
    try:
        with urllib.request.urlopen(GEQDSK_URL, context=SSL_CONTEXT) as response, \
                tmp_path.open("wb") as handle:
            handle.write(response.read())
        tmp_path.replace(target)
    except Exception:
        tmp_path.unlink(missing_ok=True)
        raise
    return target


def ensure_vmec_wout(dest_dir: Path) -> Path:
    """Download or reuse a reference VMEC wout file."""
    target = dest_dir / "wout_ncsx_c09r00_fixed.nc"
    if target.exists() and target.stat().st_size > 0:
        return target
    dest_dir.mkdir(parents=True, exist_ok=True)
    tmp_path = target.with_suffix(target.suffix + ".tmp")
    try:
        with urllib.request.urlopen(STELLOPT_WOUT_URL, context=SSL_CONTEXT) as response, \
                tmp_path.open("wb") as handle:
            handle.write(response.read())
        tmp_path.replace(target)
    except Exception:
        tmp_path.unlink(missing_ok=True)
        raise
    return target


def load_geoflux_dump(path: Path) -> dict[str, list[np.ndarray]]:
    """Load GeoFlux coordinate dump (cm) into surfaces/theta groups in meters."""
    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    type_id = data[:, 0].astype(int)
    curve_id = data[:, 1].astype(int)
    r_cm = data[:, 2]
    z_cm = data[:, 3]
    groups: dict[str, list[np.ndarray]] = {"surfaces": [], "thetas": []}

    for idx in np.unique(curve_id[type_id == 0]):
        mask = (type_id == 0) & (curve_id == idx)
        surf = np.stack((r_cm[mask], z_cm[mask]), axis=0) * 1.0e-2
        groups["surfaces"].append(surf)
    for idx in np.unique(curve_id[type_id == 1]):
        mask = (type_id == 1) & (curve_id == idx)
        line = np.stack((r_cm[mask], z_cm[mask]), axis=0) * 1.0e-2
        groups["thetas"].append(line)
    return groups


def make_curve_from_discrete_boundary(boundary: np.ndarray) -> callable:
    """Return a smooth periodic curve(t) from discrete boundary points.

    boundary has shape (2, N) with (R, Z) in meters, sampled uniformly in angle.
    We fit simple Fourier series in theta to obtain a differentiable parametrization
    suitable for map2disc.
    """
    Rb, Zb = boundary
    npts = Rb.size
    if npts < 8:
        raise ValueError("need at least 8 boundary points for Fourier fit")
    theta_nodes = np.linspace(0.0, 2.0 * np.pi, npts, endpoint=False)

    mmax = min(npts // 2, 32)
    m_vals = np.arange(0, mmax + 1, dtype=int)

    aR = np.zeros(mmax + 1)
    bR = np.zeros(mmax + 1)
    aZ = np.zeros(mmax + 1)
    bZ = np.zeros(mmax + 1)

    n_eff = float(npts)
    for j, th in enumerate(theta_nodes):
        cos_mt = np.cos(m_vals * th)
        sin_mt = np.sin(m_vals * th)
        aR = aR + Rb[j] * cos_mt
        bR = bR + Rb[j] * sin_mt
        aZ = aZ + Zb[j] * cos_mt
        bZ = bZ + Zb[j] * sin_mt

    # m = 0 mode
    aR[0] = aR[0] / n_eff
    aZ[0] = aZ[0] / n_eff
    bR[0] = 0.0
    bZ[0] = 0.0
    # higher modes
    scale = 2.0 / n_eff
    aR[1:] = scale * aR[1:]
    bR[1:] = scale * bR[1:]
    aZ[1:] = scale * aZ[1:]
    bZ[1:] = scale * bZ[1:]

    def curve(t):
        t_arr = np.asarray(t, dtype=float)
        t_mod = np.mod(t_arr, 2.0 * np.pi)
        cos_mt = np.cos(np.outer(m_vals, t_mod))
        sin_mt = np.sin(np.outer(m_vals, t_mod))
        R = aR[0] + np.sum(aR[1:, None] * cos_mt[1:, :] + bR[1:, None] * sin_mt[1:, :], axis=0)
        Z = aZ[0] + np.sum(aZ[1:, None] * cos_mt[1:, :] + bZ[1:, None] * sin_mt[1:, :], axis=0)
        return np.stack((R, Z), axis=0)

    return curve


def build_babin_mapping(boundary: np.ndarray, M: int = 10) -> m2d.BoundaryConformingMapping:
    """Construct a Babin (map2disc) mapping for a given RZ boundary."""
    curve = make_curve_from_discrete_boundary(boundary)
    bcm = m2d.BoundaryConformingMapping(curve=curve, M=M, Nt=256, Ng=(128, 128))
    bcm.solve_domain2disk()
    bcm.solve_disk2domain()
    return bcm


def sample_babin_grid(
    bcm: m2d.BoundaryConformingMapping,
    n_surface: int = 5,
    n_surface_samples: int = 180,
    n_theta_lines: int = 6,
    n_radial_samples: int = 120,
) -> dict[str, list[np.ndarray]]:
    """Sample Babin coordinate lines from a map2disc mapping."""
    if n_surface < 2:
        n_surface = 2
    rho_inner = np.linspace(0.1, 0.9, n_surface - 1)
    rho_vals = np.concatenate([rho_inner, [1.0]])
    theta = np.linspace(0.0, 2.0 * np.pi, n_surface_samples + 1)
    surfaces: list[np.ndarray] = []
    thetas: list[np.ndarray] = []

    # Precompute full tensor grid via polar evaluation
    xy_full = bcm.eval_rt_1d(rho_vals, theta)
    R_full = xy_full[0]
    Z_full = xy_full[1]

    for i, _rho in enumerate(rho_vals):
        surfaces.append(np.vstack((R_full[i, :], Z_full[i, :])))

    theta_lines = np.linspace(0.0, 2.0 * np.pi, n_theta_lines, endpoint=False)
    radial = np.linspace(0.0, 1.0, n_radial_samples)
    for ang in theta_lines:
        t_arr = np.full_like(radial, ang)
        xy = bcm.eval_rt(radial, t_arr)
        thetas.append(np.vstack((xy[0], xy[1])))

    return {"surfaces": surfaces, "thetas": thetas}


def plot_overlay(base, other, labels, title, outfile: Path) -> None:
    fig, ax = plt.subplots(figsize=(5.0, 5.0))
    for curve in base["surfaces"]:
        ax.plot(curve[0], curve[1], color="tab:blue", lw=0.9, label=labels[0])
    for curve in base["thetas"]:
        ax.plot(curve[0], curve[1], color="tab:blue", lw=0.7, ls="--")
    for curve in other["surfaces"]:
        ax.plot(curve[0], curve[1], color="tab:red", lw=0.9, label=labels[1])
    for curve in other["thetas"]:
        ax.plot(curve[0], curve[1], color="tab:red", lw=0.7, ls="-.")
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("R [m]")
    ax.set_ylabel("Z [m]")
    ax.set_title(title)
    handles, legend_labels = ax.get_legend_handles_labels()
    if legend_labels:
        ax.legend(loc="best", fontsize=8)
    ax.grid(True, alpha=0.3, linestyle="--")
    fig.tight_layout()
    outfile.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outfile, dpi=180)
    plt.close(fig)


def make_geoflux_overlay(build_dir: Path, out_root: Path) -> Path:
    """GeoFlux vs Babin overlay using EQDSK_I."""
    out_dir = out_root / "geoflux"
    out_dir.mkdir(parents=True, exist_ok=True)
    geqdsk = ensure_geqdsk(out_dir)
    geoflux_dump = out_dir / "geoflux_coords.dat"

    geo_exec = build_dir / "test" / "geoflux_coord_dump.x"
    subprocess.run([str(geo_exec), str(geqdsk), str(geoflux_dump)], check=True)

    geo = load_geoflux_dump(geoflux_dump)
    # Use outermost surface as boundary for Babin mapping
    boundary = geo["surfaces"][-1]
    bcm = build_babin_mapping(boundary)
    babin = sample_babin_grid(bcm)

    outfile = out_dir / "geoflux_vs_babin.png"
    plot_overlay(
        geo,
        babin,
        ("GeoFlux", "Babin"),
        "GeoFlux vs Babin coordinates (EQDSK_I)",
        outfile,
    )
    return outfile


def make_vmec_overlay(out_root: Path) -> Path:
    """VMEC vs Babin overlay using an NCSX reference equilibrium."""
    out_dir = out_root / "vmec"
    out_dir.mkdir(parents=True, exist_ok=True)
    wout_path = ensure_vmec_wout(out_dir)

    geom = VMECGeometry.from_file(str(wout_path))
    ns = geom.rmnc.shape[1]
    s_indices = np.linspace(0, ns - 1, 5, dtype=int)
    s_boundary = ns - 1
    theta = np.linspace(0.0, 2.0 * np.pi, 256, endpoint=False)
    zeta0 = 0.0

    # VMEC flux surfaces
    vmec_surfaces: list[np.ndarray] = []
    for s_idx in s_indices:
        R, Z, _ = geom.coords(s_idx, theta, zeta0, use_asym=True)
        vmec_surfaces.append(np.vstack((R, Z)))

    # VMEC "radial" lines at fixed poloidal angles
    vmec_thetas: list[np.ndarray] = []
    theta_lines = np.linspace(0.0, 2.0 * np.pi, 6, endpoint=False)
    s_line = np.linspace(0, ns - 1, 120)
    for ang in theta_lines:
        R_line = np.empty_like(s_line)
        Z_line = np.empty_like(s_line)
        for i, s_val in enumerate(s_line):
            s_idx = int(round(s_val))
            s_idx = max(0, min(ns - 1, s_idx))
            R_pt, Z_pt, _ = geom.coords(s_idx, np.array([ang]), zeta0, use_asym=True)
            R_line[i] = R_pt[0]
            Z_line[i] = Z_pt[0]
        vmec_thetas.append(np.vstack((R_line, Z_line)))

    vmec = {"surfaces": vmec_surfaces, "thetas": vmec_thetas}

    # Babin mapping based on the outer VMEC surface
    Rb, Zb, _ = geom.coords(s_boundary, theta, zeta0, use_asym=True)
    boundary = np.vstack((Rb, Zb))
    bcm = build_babin_mapping(boundary)
    babin = sample_babin_grid(bcm)

    outfile = out_dir / "vmec_vs_babin.png"
    plot_overlay(
        vmec,
        babin,
        ("VMEC", "Babin"),
        "VMEC vs Babin coordinates (NCSX)",
        outfile,
    )
    return outfile


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", default="build", help="Path to CMake build directory")
    parser.add_argument(
        "--output-dir", default="build/test/babin", help="Root directory for Babin plots"
    )
    args = parser.parse_args()

    build_dir = Path(args.build_dir).resolve()
    out_root = Path(args.output_dir).resolve()
    out_root.mkdir(parents=True, exist_ok=True)

    geoflux_png = make_geoflux_overlay(build_dir, out_root)
    vmec_png = make_vmec_overlay(out_root)
    print(geoflux_png)
    print(vmec_png)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
