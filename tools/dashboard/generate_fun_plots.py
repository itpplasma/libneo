#!/usr/bin/env python3
"""Generate playful PNG plots for the CI dashboard artifacts.

The plots are deterministic so that successive runs produce comparable
artifacts while still being visually interesting. All outputs are written
under ``build/test/dashboard`` so the existing artifact collection pattern
can pick them up without further changes.
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


OUTPUT_ROOT = Path("build/test/dashboard")


def _ensure_output_dir() -> None:
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)


def _save_figure(fig: plt.Figure, name: str) -> None:
    path = OUTPUT_ROOT / name
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plasma_spiral() -> None:
    """Create a synthetic magnetic surface style spiral plot."""

    theta = np.linspace(0.0, 6.0 * np.pi, 2048)
    radius = np.linspace(0.1, 1.0, theta.size)
    perturb = 0.05 * np.sin(12.0 * theta)

    x = (radius + perturb) * np.cos(theta)
    y = (radius + perturb) * np.sin(theta)

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.plot(x, y, color="#c751e5", linewidth=2.0)
    ax.set_title("Synthetic Flux Surface", fontsize=14)
    ax.set_xlabel("X [a.u.]")
    ax.set_ylabel("Y [a.u.]")
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.25)

    _save_figure(fig, "flux_surface.png")


def standing_wave_heatmap() -> None:
    """Render a 2D standing-wave pattern as a heatmap."""

    xs = np.linspace(-math.pi, math.pi, 256)
    ys = np.linspace(-math.pi, math.pi, 256)
    xx, yy = np.meshgrid(xs, ys)
    field = np.sin(3.0 * xx) * np.cos(4.0 * yy)

    fig, ax = plt.subplots(figsize=(6, 5))
    mesh = ax.pcolormesh(xx, yy, field, cmap="magma")
    fig.colorbar(mesh, ax=ax, label="Amplitude")
    ax.set_title("Standing Wave Interference")
    ax.set_xlabel("x [rad]")
    ax.set_ylabel("y [rad]")

    _save_figure(fig, "standing_wave.png")


def poincare_scatter() -> None:
    """Produce a pseudo Poincaré section using a simple twist map."""

    n_points = 6000
    theta = np.linspace(0.0, 2.0 * np.pi, n_points, endpoint=False)
    radius = 0.25 + 0.55 * np.sin(5.0 * theta)

    twist = 0.35
    phi = np.mod(theta + twist * radius, 2.0 * np.pi)

    x = radius * np.cos(phi)
    y = radius * np.sin(phi)

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(x, y, s=4.0, c=theta, cmap="viridis", alpha=0.8, linewidths=0)
    ax.set_title("Twist-Map Poincaré Section")
    ax.set_xlabel("X [a.u.]")
    ax.set_ylabel("Y [a.u.]")
    ax.set_aspect("equal")
    ax.grid(False)

    _save_figure(fig, "poincare_section.png")


def main() -> None:
    _ensure_output_dir()
    plasma_spiral()
    standing_wave_heatmap()
    poincare_scatter()


if __name__ == "__main__":
    main()
