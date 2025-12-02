#!/usr/bin/env python3
import argparse
import pathlib

import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot 2D neo_bspline LSQ fit error on full grid."
    )
    parser.add_argument(
        "--data",
        type=pathlib.Path,
        default=pathlib.Path("bspline_2d_lsq_grid.dat"),
        help="Path to grid data file written by test_neo_bspline_2d.x",
    )
    parser.add_argument(
        "--output",
        type=pathlib.Path,
        default=pathlib.Path("bspline_2d_lsq_grid.png"),
        help="Output PNG file.",
    )
    args = parser.parse_args()

    data = np.loadtxt(args.data)
    x = data[:, 0]
    y = data[:, 1]
    f_true = data[:, 2]
    f_fit = data[:, 3]
    err = f_fit - f_true

    fig, ax = plt.subplots(figsize=(6, 4))
    sc = ax.scatter(x, y, c=err, cmap="coolwarm", s=8, edgecolor="none")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("neo_bspline 2D LSQ full-grid error (fit - analytic)")
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("error")
    fig.tight_layout()
    fig.savefig(args.output, dpi=150)


if __name__ == "__main__":
    main()

