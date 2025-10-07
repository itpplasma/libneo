#!/usr/bin/env python3
"""Validate axisymmetric Biot-Savart outputs against analytic loop solution."""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path
from typing import Tuple

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[2]
PYTHON_DIR = REPO_ROOT / "python"
BUILD_DIR = REPO_ROOT / "build"
SCRIPTS_DIR = PYTHON_DIR / "scripts"

if PYTHON_DIR.exists() and str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))
if BUILD_DIR.exists() and str(BUILD_DIR) not in sys.path:
    sys.path.insert(0, str(BUILD_DIR))
if SCRIPTS_DIR.exists() and str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from plot_biotsavart_fourier import (
    _load_mode_from_bnvac,
    _load_mode_from_anvac,
    _tensordot_currents,
    _evaluate_modes_at_point,
    FOURIER_REFERENCE_CURRENT,
)
from _magfie import compute_field_from_gpec_file

MU0 = 4.0e-7 * math.pi


def read_currents(path: Path) -> np.ndarray:
    data = np.loadtxt(path, dtype=float)
    return np.atleast_1d(data)


def read_coil_radius(path: Path) -> Tuple[float, float]:
    with path.open("r") as fh:
        header = fh.readline().split()
        if len(header) < 4:
            raise ValueError(f"Malformed coil file header in {path}")
        # Read first point
        first = fh.readline().split()
        if len(first) != 3:
            raise ValueError(f"Malformed coil coordinates in {path}")
        x, y, z = map(float, first)
    radius_m = math.hypot(x, y)
    z_offset_m = z
    return radius_m, z_offset_m


def compute_axis_field_fourier(grid, BnR, Bnphi, BnZ, z_cm: np.ndarray) -> np.ndarray:
    results = np.zeros_like(z_cm, dtype=float)
    positive_radii = grid.R[grid.R > 0.0]
    fallback_r = float(positive_radii[0]) if positive_radii.size > 0 else 1.0e-6
    for idx, z in enumerate(z_cm):
        Bx, By, Bz = _evaluate_modes_at_point(BnR, Bnphi, BnZ, 0, grid.R, grid.Z, 0.0, 0.0, z)
        if not np.isfinite(Bz):
            Bz_samples = []
            for sign in (-1.0, 1.0):
                Bx_eps, By_eps, Bz_eps = _evaluate_modes_at_point(
                    BnR, Bnphi, BnZ, 0, grid.R, grid.Z, sign * fallback_r, 0.0, z
                )
                if np.isfinite(Bz_eps):
                    Bz_samples.append(Bz_eps)
            if Bz_samples:
                Bz = float(np.mean(Bz_samples))
        results[idx] = float(np.real(Bz)) if np.isfinite(Bz) else float("nan")
    return results


def compute_axis_field_direct(coil_file: Path, coil_currents: np.ndarray, z_cm: np.ndarray) -> np.ndarray:
    x = np.zeros_like(z_cm)
    y = np.zeros_like(z_cm)
    _, _, Bz = compute_field_from_gpec_file(str(coil_file), coil_currents, x, y, z_cm)
    return Bz


def analytic_loop_field(radius_m: float, currents: np.ndarray, z_cm: np.ndarray) -> np.ndarray:
    total_current = float(np.sum(currents))
    z_m = z_cm / 100.0
    denom = (radius_m ** 2 + z_m ** 2) ** 1.5
    B_tesla = MU0 * total_current * radius_m ** 2 / (2.0 * denom)
    return B_tesla * 1.0e4  # Tesla -> Gauss


def main() -> int:
    parser = argparse.ArgumentParser(description="Check axisymmetric coil response")
    parser.add_argument("--reference", required=True, type=Path, help="HDF5 Fourier output file")
    parser.add_argument("--vector", required=True, type=Path, help="NetCDF vector potential output")
    parser.add_argument("--coil", required=True, type=Path, help="Input GPEC coil file")
    parser.add_argument("--currents", required=True, type=Path, help="Coil current file")
    parser.add_argument("--summary", required=True, type=Path, help="Summary text output")
    parser.add_argument("--plot", type=Path, help="Optional PNG output showing axis comparison")
    parser.add_argument("--axis-range", type=float, default=60.0, help="Axis half-length in cm")
    parser.add_argument("--samples", type=int, default=181, help="Number of samples along axis")
    parser.add_argument("--abs-tol", type=float, default=5.0e-4, help="Absolute tolerance (Gauss)")
    parser.add_argument("--rel-tol", type=float, default=5.0e-2, help="Relative tolerance (fraction)")
    args = parser.parse_args()

    radius_m, z_offset_m = read_coil_radius(args.coil)
    currents = read_currents(args.currents)

    if abs(z_offset_m) > 1e-9:
        print("WARNING: coil not centered on z=0; analytic comparison assumes zero offset", file=sys.stderr)

    mode_fourier = _load_mode_from_bnvac(args.reference, ntor=0)
    mode_vector, _, _, _ = _load_mode_from_anvac(args.vector, ntor=0)

    weights = currents * FOURIER_REFERENCE_CURRENT
    BnR_ref = _tensordot_currents(weights, mode_fourier.BnR)
    Bnphi_ref = _tensordot_currents(weights, mode_fourier.Bnphi)
    BnZ_ref = _tensordot_currents(weights, mode_fourier.BnZ)

    BnR_vec = _tensordot_currents(weights, mode_vector.BnR)
    Bnphi_vec = _tensordot_currents(weights, mode_vector.Bnphi)
    BnZ_vec = _tensordot_currents(weights, mode_vector.BnZ)

    z_cm = np.linspace(-args.axis_range, args.axis_range, args.samples, dtype=float)

    B_ref_axis = compute_axis_field_fourier(mode_fourier.grid, BnR_ref, Bnphi_ref, BnZ_ref, z_cm)
    B_vec_axis = compute_axis_field_fourier(mode_vector.grid, BnR_vec, Bnphi_vec, BnZ_vec, z_cm)
    B_direct_axis = compute_axis_field_direct(args.coil, currents, z_cm)
    B_analytic_axis = analytic_loop_field(radius_m, currents, z_cm)

    abs_diff_ref = np.abs(B_ref_axis - B_analytic_axis)
    rel_diff_ref = abs_diff_ref / np.maximum(np.abs(B_analytic_axis), 1e-12)

    abs_diff_vec = np.abs(B_vec_axis - B_analytic_axis)
    rel_diff_vec = abs_diff_vec / np.maximum(np.abs(B_analytic_axis), 1e-12)

    abs_diff_direct = np.abs(B_direct_axis - B_analytic_axis)
    rel_diff_direct = abs_diff_direct / np.maximum(np.abs(B_analytic_axis), 1e-12)

    mask_ref = np.isfinite(abs_diff_ref)
    mask_vec = np.isfinite(abs_diff_vec)
    mask_direct = np.isfinite(abs_diff_direct)

    max_abs_ref = float(np.max(abs_diff_ref[mask_ref])) if np.any(mask_ref) else float("nan")
    max_rel_ref = float(np.max(rel_diff_ref[mask_ref])) if np.any(mask_ref) else float("nan")
    max_abs_vec = float(np.max(abs_diff_vec[mask_vec])) if np.any(mask_vec) else float("nan")
    max_rel_vec = float(np.max(rel_diff_vec[mask_vec])) if np.any(mask_vec) else float("nan")
    max_abs_direct = float(np.max(abs_diff_direct[mask_direct])) if np.any(mask_direct) else float("nan")
    max_rel_direct = float(np.max(rel_diff_direct[mask_direct])) if np.any(mask_direct) else float("nan")

    summary = (
        "Axisymmetric loop validation\n"
        f"  Radius (m):            {radius_m:.6f}\n"
        f"  Samples:               {args.samples}\n"
        f"  Max |B_ref - analytic|:   {max_abs_ref:.6e} G\n"
        f"  Max rel diff ref:        {max_rel_ref:.6e}\n"
        f"  Max |B_vec - analytic|:   {max_abs_vec:.6e} G\n"
        f"  Max rel diff vec:        {max_rel_vec:.6e}\n"
        f"  Max |B_dir - analytic|:   {max_abs_direct:.6e} G\n"
        f"  Max rel diff direct:     {max_rel_direct:.6e}\n"
    )

    args.summary.write_text(summary)

    if args.plot is not None:
      fig, (ax_field, ax_err) = plt.subplots(2, 1, sharex=True, figsize=(7, 6))

      colors = {
          'analytic': '#000000',
          'fourier': '#1f77b4',
          'vector': '#ff7f0e',
          'direct': '#2ca02c',
      }

      ax_field.plot(z_cm, B_analytic_axis, label='Analytic', linewidth=2, color=colors['analytic'])
      ax_field.plot(z_cm, B_ref_axis, label='Fourier', linestyle='-', color=colors['fourier'])
      ax_field.plot(z_cm, B_vec_axis, label='Vector', linestyle='--', color=colors['vector'])
      finite_direct = np.isfinite(B_direct_axis)
      ax_field.plot(z_cm[finite_direct], B_direct_axis[finite_direct], label='Direct', linestyle=':', color=colors['direct'])
      ax_field.set_ylabel('B_parallel [G]')
      ax_field.set_title('Axis field comparison')
      ax_field.legend(loc='best')
      ax_field.grid(True, alpha=0.3)

      eps = 1e-20
      ax_err.semilogy(z_cm, np.maximum(abs_diff_ref, eps), label='|Fourier - analytic|', linestyle='-', color=colors['fourier'])
      ax_err.semilogy(z_cm, np.maximum(abs_diff_vec, eps), label='|Vector - analytic|', linestyle='--', color=colors['vector'])
      if np.any(finite_direct):
        ax_err.semilogy(z_cm[finite_direct], np.maximum(abs_diff_direct[finite_direct], eps), label='|Direct - analytic|', linestyle=':', color=colors['direct'])
      ax_err.set_ylabel('|ΔB| [G]')
      ax_err.set_xlabel('z [cm]')
      ax_err.set_title('Absolute error along axis')
      ax_err.grid(True, which='both', alpha=0.3)
      ax_err.legend(loc='best')

      fig.tight_layout()
      fig.savefig(args.plot, dpi=150)
      plt.close(fig)

    if not np.isfinite(max_abs_ref) or not np.isfinite(max_rel_ref):
        print(summary, file=sys.stderr)
        return 1

    if max_abs_ref > args.abs_tol or max_rel_ref > args.rel_tol:
        print(summary, file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
