#!/usr/bin/env python3
"""Generate quick-look plots from a public GEQDSK sample."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402  (after Agg selection)

from run_geoflux_test import ensure_geqdsk, require_testing_enabled  # reuse downloader


def parse_geqdsk(path: Path):
    with path.open("r", encoding="utf-8") as handle:
        header = handle.readline()
        parts = header.split()
        if len(parts) < 3:
            raise ValueError("invalid GEQDSK header")
        nw = int(parts[-2])
        nh = int(parts[-1])

        def read_values(count: int, width: int, converter=float) -> list:
            values: list = []
            while len(values) < count:
                line = handle.readline()
                if not line:
                    raise ValueError("unexpected EOF in GEQDSK file")
                for start in range(0, len(line), width):
                    chunk = line[start:start + width]
                    if not chunk.strip():
                        continue
                    values.append(converter(chunk))
                    if len(values) == count:
                        break
            return values

        rdim, zdim, _rcentr, rleft, zmid = read_values(5, 16)
        _rmaxis, _zmaxis, simag, sibry, _bcentr = read_values(5, 16)
        current_line = read_values(5, 16)
        _current = current_line[0]
        # skip remaining placeholder entries
        _ = read_values(5, 16)

        fpol = np.array(read_values(nw, 16), dtype=float)
        pres = np.array(read_values(nw, 16), dtype=float)
        ffprim = np.array(read_values(nw, 16), dtype=float)
        pprime = np.array(read_values(nw, 16), dtype=float)
        psirz = np.array(read_values(nw * nh, 16), dtype=float)
        qpsi = np.array(read_values(nw, 16), dtype=float)

        nbbbs, limitr = read_values(2, 5, int)
        # consume contour definitions (pairs of values)
        _ = read_values(nbbbs * 2, 16)
        _ = read_values(limitr * 2, 16)

    psi_axis = min(simag, sibry)
    psi_sep = max(simag, sibry)
    psinorm = (psirz - psi_axis) / (psi_sep - psi_axis)
    psinorm = psinorm.reshape((nw, nh), order="F")

    rz_r = rleft + np.linspace(0.0, rdim, nw)
    rz_z = zmid - 0.5 * zdim + np.linspace(0.0, zdim, nh)

    return {
        "R": rz_r,
        "Z": rz_z,
        "psinorm": psinorm,
        "fpol": fpol,
        "pres": pres,
        "ffprim": ffprim,
        "pprime": pprime,
        "qpsi": qpsi,
    }


def make_plots(data, output_dir: Path, basename: str = "geqdsk") -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []

    r, z = data["R"], data["Z"]
    psinorm = data["psinorm"]
    R_grid, Z_grid = np.meshgrid(r, z, indexing="ij")

    fig, ax = plt.subplots(figsize=(4.0, 4.0))
    levels = np.linspace(0.0, 1.0, 11)
    cs = ax.contour(R_grid * 1e-2, Z_grid * 1e-2, psinorm, levels=levels)
    ax.clabel(cs, fmt="%.2f", fontsize=6)
    ax.set_xlabel("R [m]")
    ax.set_ylabel("Z [m]")
    ax.set_title("Normalized poloidal flux")
    ax.set_aspect("equal", adjustable="box")
    contour_path = output_dir / f"{basename}_contours.png"
    fig.tight_layout()
    fig.savefig(contour_path, dpi=150)
    plt.close(fig)
    paths.append(contour_path)

    profiles_fig, profiles_ax = plt.subplots(figsize=(4.0, 3.0))
    psi = np.linspace(0.0, 1.0, data["fpol"].size)
    profiles_ax.plot(psi, data["fpol"] / data["fpol"][0], label="fpol/fpol(0)")
    profiles_ax.plot(psi, data["pres"] / max(abs(data["pres"])), label="pres")
    profiles_ax.legend(fontsize=7)
    profiles_ax.set_xlabel("Normalized psi")
    profiles_ax.set_ylabel("Scaled value")
    profiles_ax.set_title("GEQDSK profiles (scaled)")
    profiles_ax.grid(True, alpha=0.3)
    profiles_path = output_dir / f"{basename}_profiles.png"
    profiles_fig.tight_layout()
    profiles_fig.savefig(profiles_path, dpi=150)
    plt.close(profiles_fig)
    paths.append(profiles_path)

    return paths


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="geqdsk_plot")
    parser.add_argument("--data-dir", required=True, help="directory for cached GEQDSK files")
    parser.add_argument("--output-dir", required=True, help="directory for generated plots")
    parser.add_argument("--basename", default="geqdsk", help="base name for output plots")
    args = parser.parse_args(argv)

    require_testing_enabled("geqdsk plotting test")

    data_dir = Path(args.data_dir).resolve()
    geqdsk_path = ensure_geqdsk(data_dir)
    data = parse_geqdsk(geqdsk_path)

    output_dir = Path(args.output_dir).resolve()
    paths = make_plots(data, output_dir, basename=args.basename)

    for path in paths:
        if not path.exists() or path.stat().st_size == 0:
            raise RuntimeError(f"missing plot: {path}")
        print(path)

    return 0


if __name__ == "__main__":
    sys.exit(main())
