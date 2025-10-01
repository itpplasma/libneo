#!/usr/bin/env python3
"""Validate axisymmetric Biot-Savart outputs against analytic loop solution."""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path
from typing import Tuple

import numpy as np

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

try:
    from compare_superposition import (
        load_bnvac_modes,
        load_anvac_modes,
        superpose_modes,
        evaluate_modes_at_point,
        FOURIER_REFERENCE_CURRENT,
        compute_field_from_gpec_file,
    )
except ImportError as exc:  # pragma: no cover - hard failure on missing deps
    print(f"ERROR: unable to import compare_superposition helpers: {exc}", file=sys.stderr)
    sys.exit(1)

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


def compute_axis_field_fourier(ref_grid, ntor_vals, BnR, Bnphi, BnZ, z_cm: np.ndarray) -> np.ndarray:
    results = np.zeros_like(z_cm, dtype=float)
    for idx, z in enumerate(z_cm):
        x = 0.0
        y = 0.0
        Bx, By, Bz = evaluate_modes_at_point(BnR, Bnphi, BnZ, ntor_vals, ref_grid.R, ref_grid.Z, x, y, z)
        results[idx] = float(np.real(Bz))
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
    parser.add_argument("--axis-range", type=float, default=60.0, help="Axis half-length in cm")
    parser.add_argument("--samples", type=int, default=181, help="Number of samples along axis")
    parser.add_argument("--abs-tol", type=float, default=5.0e-4, help="Absolute tolerance (Gauss)")
    parser.add_argument("--rel-tol", type=float, default=5.0e-2, help="Relative tolerance (fraction)")
    args = parser.parse_args()

    radius_m, z_offset_m = read_coil_radius(args.coil)
    currents = read_currents(args.currents)

    if abs(z_offset_m) > 1e-9:
        print("WARNING: coil not centered on z=0; analytic comparison assumes zero offset", file=sys.stderr)

    ref_grid, ntor_vals_ref, BnR_modes_ref, Bnphi_modes_ref, BnZ_modes_ref = load_bnvac_modes(str(args.reference))
    vec_grid, ntor_vals_vec, BnR_modes_vec, Bnphi_modes_vec, BnZ_modes_vec = load_anvac_modes(str(args.vector))

    if not np.array_equal(ntor_vals_ref, ntor_vals_vec):
        raise RuntimeError("ntor spectra differ between reference and vector datasets")

    BnR_ref = superpose_modes(BnR_modes_ref, currents, FOURIER_REFERENCE_CURRENT)
    Bnphi_ref = superpose_modes(Bnphi_modes_ref, currents, FOURIER_REFERENCE_CURRENT)
    BnZ_ref = superpose_modes(BnZ_modes_ref, currents, FOURIER_REFERENCE_CURRENT)

    BnR_vec = superpose_modes(BnR_modes_vec, currents, FOURIER_REFERENCE_CURRENT)
    Bnphi_vec = superpose_modes(Bnphi_modes_vec, currents, FOURIER_REFERENCE_CURRENT)
    BnZ_vec = superpose_modes(BnZ_modes_vec, currents, FOURIER_REFERENCE_CURRENT)

    z_cm = np.linspace(-args.axis_range, args.axis_range, args.samples, dtype=float)

    B_ref_axis = compute_axis_field_fourier(ref_grid, ntor_vals_ref, BnR_ref, Bnphi_ref, BnZ_ref, z_cm)
    B_vec_axis = compute_axis_field_fourier(vec_grid, ntor_vals_vec, BnR_vec, Bnphi_vec, BnZ_vec, z_cm)
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

    if not np.isfinite(max_abs_ref) or not np.isfinite(max_rel_ref):
        print(summary, file=sys.stderr)
        return 1

    if max_abs_ref > args.abs_tol or max_rel_ref > args.rel_tol:
        print(summary, file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
