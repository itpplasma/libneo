#!/usr/bin/env python3
"""Overlay G-frame coordinate lines with GeoFlux and circular VMEC-like cases."""

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

GEQDSK_URL = "https://crppwww.epfl.ch/~sauter/benchmark/EQDSK_I"
SSL_CONTEXT = ssl._create_unverified_context()


def ensure_geqdsk(dest_dir: Path) -> Path:
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


def load_dump(path: Path) -> dict[str, list[np.ndarray]]:
    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    type_id = data[:, 0].astype(int)
    curve_id = data[:, 1].astype(int)
    r = data[:, 2]
    z = data[:, 3]
    groups = {"surfaces": [], "thetas": []}
    for idx in np.unique(curve_id[type_id == 0]):
        mask = (type_id == 0) & (curve_id == idx)
        groups["surfaces"].append(np.stack((r[mask], z[mask]), axis=0))
    for idx in np.unique(curve_id[type_id == 1]):
        mask = (type_id == 1) & (curve_id == idx)
        groups["thetas"].append(np.stack((r[mask], z[mask]), axis=0))
    return groups


def plot_overlay(base, other, labels, title, outfile):
    fig, ax = plt.subplots(figsize=(5.0, 5.0))
    for curve in base["surfaces"]:
        ax.plot(curve[0] * 1e-2, curve[1] * 1e-2, color="tab:blue", lw=0.9, label=labels[0])
    for curve in base["thetas"]:
        ax.plot(curve[0] * 1e-2, curve[1] * 1e-2, color="tab:blue", lw=0.7, ls="--")
    for curve in other["surfaces"]:
        ax.plot(curve[0] * 1e-2, curve[1] * 1e-2, color="tab:red", lw=0.9, label=labels[1])
    for curve in other["thetas"]:
        ax.plot(curve[0] * 1e-2, curve[1] * 1e-2, color="tab:red", lw=0.7, ls="-.")
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("R [m]")
    ax.set_ylabel("Z [m]")
    ax.set_title(title)
    handles, legend_labels = ax.get_legend_handles_labels()
    if legend_labels:
        ax.legend(loc="best", fontsize=8)
    ax.grid(True, alpha=0.3, linestyle="--")
    fig.tight_layout()
    fig.savefig(outfile, dpi=180)
    plt.close(fig)


def make_geoflux_overlay(build_dir: Path, out_dir: Path) -> Path:
    geqdsk = ensure_geqdsk(out_dir)
    geoflux_dump = out_dir / "geoflux_lines.dat"
    gframe_dump = out_dir / "gframe_lines.dat"

    geo_exec = build_dir / "test" / "geoflux_coord_dump.x"
    gframe_exec = build_dir / "test" / "gframe_from_geoflux_dump.x"
    subprocess.run([str(geo_exec), str(geqdsk), str(geoflux_dump)], check=True)
    subprocess.run([str(gframe_exec), str(geqdsk), str(gframe_dump)], check=True)

    geo = load_dump(geoflux_dump)
    gfr = load_dump(gframe_dump)
    outfile = out_dir / "geoflux_vs_gframe.png"
    plot_overlay(geo, gfr, ("GeoFlux", "G-frame"), "GeoFlux vs G-frame (EQDSK_I)", outfile)
    return outfile


def make_vmec_overlay(out_dir: Path) -> Path:
    r0 = 1.7
    a = 0.35
    n_surface = 5
    n_theta_lines = 6
    theta = np.linspace(0.0, 2 * np.pi, 361)
    surfaces = []
    thetas = []
    rho_vals = np.linspace(0.1, 0.9, n_surface)
    for rho in rho_vals:
        r = r0 + a * rho * np.cos(theta)
        z = a * rho * np.sin(theta)
        surfaces.append(np.vstack((r, z)))
    theta_lines = np.linspace(0, 2 * np.pi, n_theta_lines, endpoint=False)
    for ang in theta_lines:
        rho = np.linspace(0.0, 0.95, 160)
        r = r0 + a * rho * np.cos(ang)
        z = a * rho * np.sin(ang)
        thetas.append(np.vstack((r, z)))
    vmec = {"surfaces": surfaces, "thetas": thetas}
    gframe = {"surfaces": [s.copy() for s in surfaces], "thetas": [t.copy() for t in thetas]}
    outfile = out_dir / "vmec_vs_gframe_circular.png"
    plot_overlay(vmec, gframe, ("VMEC-like", "G-frame"), "Circular VMEC-like vs G-frame", outfile)
    return outfile


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", default="build", help="Path to CMake build directory")
    parser.add_argument(
        "--output-dir", default="build/test/plots", help="Directory for generated plots"
    )
    args = parser.parse_args()

    build_dir = Path(args.build_dir).resolve()
    out_dir = Path(args.output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    geoflux_png = make_geoflux_overlay(build_dir, out_dir)
    vmec_png = make_vmec_overlay(out_dir)
    print(geoflux_png)
    print(vmec_png)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
