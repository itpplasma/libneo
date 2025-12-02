#!/usr/bin/env python
"""Plot masks and errors used in interpolate masked tests."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def load_mask(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    data = np.loadtxt(path)
    if data.ndim != 2 or data.shape[1] != 3:
        raise RuntimeError(f"Expected 3 columns in {path}, got shape {data.shape}")
    x = data[:, 0]
    y = data[:, 1]
    m = data[:, 2]
    x_unique = np.unique(x)
    y_unique = np.unique(y)
    nx = x_unique.size
    ny = y_unique.size
    if nx * ny != x.size:
        raise RuntimeError(f"Cannot reshape data from {path} to a regular grid")
    x_grid = x_unique
    y_grid = y_unique
    mask_grid = m.reshape(ny, nx)
    return x_grid, y_grid, mask_grid


def plot_mask(x, y, mask, title: str, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(4.0, 4.0))
    xx, yy = np.meshgrid(x, y, indexing="xy")
    im = ax.pcolormesh(xx, yy, mask, cmap="Greys", shading="nearest")
    ax.set_aspect("equal")
    ax.set_xlabel("x1")
    ax.set_ylabel("x2")
    ax.set_title(title)
    fig.colorbar(im, ax=ax, label="mask")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def load_eval(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    data = np.loadtxt(path)
    if data.ndim != 2 or data.shape[1] != 6:
        raise RuntimeError(f"Expected 6 columns in {path}, got shape {data.shape}")
    x = data[:, 0]
    y = data[:, 1]
    expected = data[:, 2]
    actual = data[:, 3]
    error = data[:, 4]
    valid = data[:, 5] > 0.5
    x_unique = np.unique(x)
    y_unique = np.unique(y)
    nx = x_unique.size
    ny = y_unique.size
    if nx * ny != x.size:
        raise RuntimeError(f"Cannot reshape data from {path} to a regular grid")
    x_grid = x_unique
    y_grid = y_unique
    error_grid = error.reshape(ny, nx)
    valid_grid = valid.reshape(ny, nx)
    error_plot = np.where(valid_grid, error_grid, np.nan)
    return x_grid, y_grid, error_plot, valid_grid


def plot_error(x, y, error, title: str, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(4.0, 4.0))
    xx, yy = np.meshgrid(x, y, indexing="xy")
    vmax = np.nanmax(np.abs(error))
    if not np.isfinite(vmax) or vmax == 0.0:
        vmax = 1.0
    im = ax.pcolormesh(
        xx,
        yy,
        error,
        cmap="coolwarm",
        shading="nearest",
        vmin=-vmax,
        vmax=vmax,
    )
    ax.set_aspect("equal")
    ax.set_xlabel("x1")
    ax.set_ylabel("x2")
    ax.set_title(title)
    fig.colorbar(im, ax=ax, label="f_interp - f_exact")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Plot interpolate masked-ring and masked-sphere masks and errors."
        )
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path("build") / "test",
        help="Directory with interpolate_mask_*.dat files.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("build") / "test",
        help="Directory where PNGs will be written.",
    )
    args = parser.parse_args()

    data_dir = args.data_dir
    out_dir = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    ring_dat = data_dir / "interpolate_mask_2d_ring.dat"
    ring_eval = data_dir / "interpolate_mask_2d_ring_eval.dat"
    shell_dat = data_dir / "interpolate_mask_3d_shell_slice.dat"
    shell_eval = data_dir / "interpolate_mask_3d_shell_slice_eval.dat"

    if ring_dat.exists():
        x_ring, y_ring, m_ring = load_mask(ring_dat)
        plot_mask(
            x_ring,
            y_ring,
            m_ring,
            "interpolate 2D masked ring",
            out_dir / "interpolate_mask_2d_ring.png",
        )
        if ring_eval.exists():
            x_re, y_re, err_ring, _ = load_eval(ring_eval)
            plot_error(
                x_re,
                y_re,
                err_ring,
                "interpolate 2D masked ring error",
                out_dir / "interpolate_mask_2d_ring_error.png",
            )

    if shell_dat.exists():
        x_shell, y_shell, m_shell = load_mask(shell_dat)
        plot_mask(
            x_shell,
            y_shell,
            m_shell,
            "interpolate 3D masked sphere slice",
            out_dir / "interpolate_mask_3d_shell_slice.png",
        )
        if shell_eval.exists():
            x_se, y_se, err_shell, _ = load_eval(shell_eval)
            plot_error(
                x_se,
                y_se,
                err_shell,
                "interpolate 3D masked sphere slice error",
                out_dir / "interpolate_mask_3d_shell_slice_error.png",
            )


if __name__ == "__main__":
    main()
