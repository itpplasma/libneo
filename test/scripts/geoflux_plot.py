#!/usr/bin/env python3
"""Plot geoflux coordinate lines from the sampled Fortran dump."""

from __future__ import annotations

import argparse
import subprocess
import tempfile
from pathlib import Path
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from run_geoflux_test import ensure_geqdsk


def load_coords(path: Path) -> dict[str, list[np.ndarray]]:
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
        groups["surfaces"].append(np.stack((r_cm[mask], z_cm[mask]), axis=0))
    for idx in np.unique(curve_id[type_id == 1]):
        mask = (type_id == 1) & (curve_id == idx)
        groups["thetas"].append(np.stack((r_cm[mask], z_cm[mask]), axis=0))
    return groups


def make_plot(groups: dict[str, list[np.ndarray]], output: Path) -> None:
    fig, ax = plt.subplots(figsize=(4.5, 4.5))
    for surf in groups["surfaces"]:
        ax.plot(surf[0] * 1.0e-2, surf[1] * 1.0e-2, color="tab:blue", linewidth=0.8)
    for line in groups["thetas"]:
        ax.plot(
            line[0] * 1.0e-2,
            line[1] * 1.0e-2,
            color="tab:orange",
            linestyle="--",
            linewidth=0.6,
        )
    ax.set_xlabel("R [m]")
    ax.set_ylabel("Z [m]")
    ax.set_aspect("equal", adjustable="box")
    ax.set_title("Geoflux coordinate lines")
    ax.grid(True, alpha=0.2)
    fig.tight_layout()
    fig.savefig(output, dpi=150)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="geoflux_plot")
    parser.add_argument("--exe", required=True, help="Path to geoflux dump executable")
    parser.add_argument("--data-dir", required=True, help="Directory for cached GEQDSK files")
    parser.add_argument("--output-dir", required=True, help="Directory for generated plots")
    parser.add_argument("--basename", default="geoflux", help="Base name for output plot")
    args = parser.parse_args(argv)

    exe_path = Path(args.exe).resolve()
    if not exe_path.is_file():
        parser.error(f"executable not found: {exe_path}")

    data_dir = Path(args.data_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    geqdsk_path = ensure_geqdsk(data_dir)

    with tempfile.TemporaryDirectory() as tmpdir:
        dump_path = Path(tmpdir) / "geoflux_coords.dat"
        result = subprocess.run(
            [str(exe_path), str(geqdsk_path), str(dump_path)],
            check=False,
        )
        if result.returncode != 0:
            raise subprocess.CalledProcessError(result.returncode, result.args, result.returncode)
        groups = load_coords(dump_path)

    plot_path = output_dir / f"{args.basename}_coords.png"
    make_plot(groups, plot_path)
    if not plot_path.exists() or plot_path.stat().st_size == 0:
        raise RuntimeError(f"plot was not created: {plot_path}")
    print(plot_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
