#!/usr/bin/env python3
import argparse
import pathlib

import matplotlib.pyplot as plt
import numpy as np


def load_bench(path: pathlib.Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    data = np.loadtxt(path)
    n = data[:, 0]
    t_create_interp = data[:, 1]
    t_eval_interp = data[:, 2]
    t_create_bs = data[:, 3]
    t_eval_bs = data[:, 4]
    return n, t_create_interp, t_eval_interp, t_create_bs, t_eval_bs


def plot_dim(path: pathlib.Path, title: str, output: pathlib.Path) -> None:
    n, t_ci, t_ei, t_cb, t_eb = load_bench(path)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.loglog(n, t_ci, "o-", label="interp create")
    ax.loglog(n, t_ei, "o--", label="interp eval")
    ax.loglog(n, t_cb, "s-", label="neo_bspline create")
    ax.loglog(n, t_eb, "s--", label="neo_bspline eval")

    ax.set_xlabel("number of data points")
    ax.set_ylabel("time [s]")
    ax.set_title(title)
    ax.legend()
    ax.grid(True, which="both", ls=":")

    fig.tight_layout()
    fig.savefig(output, dpi=150)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot runtime benchmarks: interpolate vs neo_bspline.",
    )
    parser.add_argument(
        "--data-dir",
        type=pathlib.Path,
        default=pathlib.Path("."),
        help="Directory with bench_neo_bspline_*.dat files.",
    )
    parser.add_argument(
        "--output-dir",
        type=pathlib.Path,
        default=pathlib.Path("."),
        help="Directory for output PNG files.",
    )
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    plot_dim(
        args.data_dir / "bench_neo_bspline_1d.dat",
        "1D: interpolate vs neo_bspline",
        args.output_dir / "bench_neo_bspline_1d_runtime.png",
    )
    plot_dim(
        args.data_dir / "bench_neo_bspline_2d.dat",
        "2D: interpolate vs neo_bspline",
        args.output_dir / "bench_neo_bspline_2d_runtime.png",
    )
    plot_dim(
        args.data_dir / "bench_neo_bspline_3d.dat",
        "3D: interpolate vs neo_bspline",
        args.output_dir / "bench_neo_bspline_3d_runtime.png",
    )


if __name__ == "__main__":
    main()

