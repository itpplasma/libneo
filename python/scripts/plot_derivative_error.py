#!/usr/bin/env python3
"""Plot derivative error series produced by test_coil_tools_vector_potential_derivs.

The script expects text files with two columns: coordinate and relative error.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import os

import matplotlib

if os.environ.get("MPLBACKEND") is None:  # pragma: no cover
    matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


def _plot_series(data_path: Path, output_path: Path, xlabel: str, ylabel: str, title: str, dpi: int) -> None:
    data = np.loadtxt(data_path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.size == 0:
        raise ValueError(f"No data found in {data_path}")

    x = data[:, 0]
    y = data[:, 1]

    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(7.0, 4.5))
    ax.plot(x, y)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(output_path, dpi=dpi)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--r-input", type=Path, required=True, help="Data file for dAphi/dR error")
    parser.add_argument("--r-output", type=Path, required=True, help="PNG output for dAphi/dR error")
    parser.add_argument("--z-input", type=Path, required=True, help="Data file for dAphi/dZ error")
    parser.add_argument("--z-output", type=Path, required=True, help="PNG output for dAphi/dZ error")
    parser.add_argument("--xlabel-r", default="R [cm]")
    parser.add_argument("--xlabel-z", default="Z [cm]")
    parser.add_argument("--xlabel-x", default="R [cm]")
    parser.add_argument("--xlabel-y", default="R [cm]")
    parser.add_argument("--ylabel", default="Relative error")
    parser.add_argument("--title-r", default="Relative error |Δ(dAφ/dR)|")
    parser.add_argument("--title-z", default="Relative error |Δ(dAφ/dZ)|")
    parser.add_argument("--title-x", default="Relative error |Δ(dAφ/dx)|")
    parser.add_argument("--title-y", default="Relative error |Δ(dAφ/dy)|")
    parser.add_argument("--x-input", type=Path, help="Data file for dAphi/dx error")
    parser.add_argument("--x-output", type=Path, help="PNG output for dAphi/dx error")
    parser.add_argument("--y-input", type=Path, help="Data file for dAphi/dy error")
    parser.add_argument("--y-output", type=Path, help="PNG output for dAphi/dy error")
    parser.add_argument("--dpi", type=int, default=150)
    args = parser.parse_args()

    _plot_series(args.r_input, args.r_output, args.xlabel_r, args.ylabel, args.title_r, args.dpi)
    _plot_series(args.z_input, args.z_output, args.xlabel_z, args.ylabel, args.title_z, args.dpi)
    if args.x_input and args.x_output:
        _plot_series(args.x_input, args.x_output, args.xlabel_x, args.ylabel, args.title_x, args.dpi)
    if args.y_input and args.y_output:
        _plot_series(args.y_input, args.y_output, args.xlabel_y, args.ylabel, args.title_y, args.dpi)


if __name__ == "__main__":  # pragma: no cover
    main()
