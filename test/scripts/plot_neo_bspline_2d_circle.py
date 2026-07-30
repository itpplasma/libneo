#!/usr/bin/env python3
import argparse
import pathlib

import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot 2D neo_bspline LSQ fit vs analytic along a circle."
    )
    parser.add_argument(
        "--data",
        type=pathlib.Path,
        default=pathlib.Path("bspline_2d_lsq_circle.dat"),
        help="Path to circle data file written by test_neo_bspline_2d.x",
    )
    parser.add_argument(
        "--output",
        type=pathlib.Path,
        default=pathlib.Path("bspline_2d_lsq_circle.png"),
        help="Output PNG file.",
    )
    args = parser.parse_args()

    data = np.loadtxt(args.data)
    theta = data[:, 0]
    f_true = data[:, 3]
    f_fit = data[:, 4]

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(theta, f_true, label="analytic", color="k", linewidth=1.5)
    ax.plot(theta, f_fit, label="neo_bspline fit", color="C1", linestyle="--")
    ax.set_xlabel("theta")
    ax.set_ylabel("f(theta)")
    ax.set_title("neo_bspline 2D LSQ fit along circle")
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.output, dpi=150)


if __name__ == "__main__":
    main()

