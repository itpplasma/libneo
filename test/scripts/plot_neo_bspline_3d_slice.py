#!/usr/bin/env python3
import argparse
import pathlib

import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot 3D neo_bspline LSQ fit error on a mid-plane slice.",
    )
    parser.add_argument(
        "--data",
        type=pathlib.Path,
        default=pathlib.Path("bspline_3d_lsq_grid.dat"),
        help="Path to 3D grid data file written by test_neo_bspline_3d.x",
    )
    parser.add_argument(
        "--output",
        type=pathlib.Path,
        default=pathlib.Path("bspline_3d_lsq_slice.png"),
        help="Output PNG file.",
    )
    args = parser.parse_args()

    data = np.loadtxt(args.data)
    x = data[:, 0]
    y = data[:, 1]
    z = data[:, 2]
    f_true = data[:, 3]
    f_fit = data[:, 4]
    err = f_fit - f_true

    z_levels = np.unique(z)
    z_mid = z_levels[len(z_levels) // 2]
    mask = np.isclose(z, z_mid)

    fig, ax = plt.subplots(figsize=(6, 4))
    sc = ax.scatter(x[mask], y[mask], c=err[mask], cmap="coolwarm", s=8, edgecolor="none")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(f"neo_bspline 3D LSQ mid-plane error (z = {z_mid:.3f})")
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("error")
    fig.tight_layout()
    fig.savefig(args.output, dpi=150)


if __name__ == "__main__":
    main()

