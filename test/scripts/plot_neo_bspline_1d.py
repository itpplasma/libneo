#!/usr/bin/env python3
import argparse
import pathlib

import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot 1D neo_bspline LSQ fit vs analytic reference."
    )
    parser.add_argument(
        "--data",
        type=pathlib.Path,
        default=pathlib.Path("bspline_1d_lsq.dat"),
        help="Path to data file written by test_neo_bspline_1d.x",
    )
    parser.add_argument(
        "--output",
        type=pathlib.Path,
        default=pathlib.Path("bspline_1d_lsq.png"),
        help="Output PNG file for the plot.",
    )
    args = parser.parse_args()

    data = np.loadtxt(args.data)
    x = data[:, 0]
    f_true = data[:, 1]
    f_fit = data[:, 2]

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(x, f_true, label="analytic", color="k", linewidth=1.5)
    ax.plot(x, f_fit, label="neo_bspline fit", color="C0", linestyle="--")
    ax.set_xlabel("x")
    ax.set_ylabel("f(x)")
    ax.set_title("neo_bspline 1D LSQ fit")
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.output, dpi=150)


if __name__ == "__main__":
    main()

