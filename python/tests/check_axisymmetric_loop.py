#!/usr/bin/env python3
"""Validate axisymmetric Biot-Savart outputs against an analytic loop."""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path
from typing import Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[2]
PYTHON_DIR = REPO_ROOT / "python"

if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))

from libneo.biotsavart_fourier import (
    field_divfree,
    gauge_Anvac,
    read_Anvac_fourier_all,
    read_Bnvac_fourier_all,
    reconstruct_field_from_modes,
    spline_gauged_Anvac,
)
from _magfie import compute_field_from_gpec_file

MU0 = 4.0e-7 * math.pi
FOURIER_REFERENCE_CURRENT = 0.1


def read_currents(path: Path) -> np.ndarray:
    data = np.loadtxt(path, dtype=float)
    return np.atleast_1d(data)


def read_coil_radius(path: Path) -> Tuple[float, float]:
    with path.open("r") as fh:
        header = fh.readline().split()
        if len(header) < 4:
            raise ValueError(f"Malformed coil file header in {path}")
        first = fh.readline().split()
        if len(first) != 3:
            raise ValueError(f"Malformed coil coordinates in {path}")
        x, y, z = map(float, first)
    radius_m = math.hypot(x, y)
    z_offset_m = z
    return radius_m, z_offset_m


def compute_axis_field_direct(coil_file: Path, coil_currents: np.ndarray, z_cm: np.ndarray) -> np.ndarray:
    x = np.zeros_like(z_cm)
    y = np.zeros_like(z_cm)
    _, _, Bz = compute_field_from_gpec_file(str(coil_file), coil_currents, x, y, z_cm)
    return Bz


def analytic_loop_field(radius_m: float, currents: np.ndarray, z_cm: np.ndarray) -> np.ndarray:
    total_current = float(np.sum(currents))
    z_m = z_cm / 100.0
    denom = (radius_m**2 + z_m**2) ** 1.5
    B_tesla = MU0 * total_current * radius_m**2 / (2.0 * denom)
    return B_tesla * 1.0e4


def load_fourier_sum(path: Path, currents: np.ndarray):
    grid_ref, modes, BnR, Bnphi, BnZ = read_Bnvac_fourier_all(str(path))
    weights = currents * FOURIER_REFERENCE_CURRENT
    BnR_sum = np.tensordot(weights, BnR, axes=(0, 1))
    Bnphi_sum = np.tensordot(weights, Bnphi, axes=(0, 1))
    BnZ_sum = np.tensordot(weights, BnZ, axes=(0, 1))
    return grid_ref, modes, BnR_sum, Bnphi_sum, BnZ_sum


def load_anvac_sum(path: Path, modes_ref: np.ndarray, grid_ref, currents: np.ndarray):
    grid_raw, modes, AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ = read_Anvac_fourier_all(str(path))
    if not np.array_equal(modes, modes_ref):
        raise RuntimeError("Mode numbers differ between Fourier and Anvac data")

    weights = currents * FOURIER_REFERENCE_CURRENT
    nmodes = modes.size
    BnR_sum = np.empty((nmodes, grid_ref.nR, grid_ref.nZ), dtype=complex)
    Bnphi_sum = np.empty_like(BnR_sum)
    BnZ_sum = np.empty_like(BnR_sum)

    for idx, ntor in enumerate(modes):
        gauged_AnR, gauged_AnZ = gauge_Anvac(
            grid_raw,
            AnR[idx],
            Anphi[idx],
            AnZ[idx],
            dAnphi_dR[idx],
            dAnphi_dZ[idx],
            ntor=ntor,
        )
        if ntor != 0:
            gauged_AnR = gauged_AnR - 0.5j * Anphi[idx]
        spline = spline_gauged_Anvac(
            grid_raw,
            gauged_AnR,
            gauged_AnZ,
            ntor=ntor,
            Anphi=Anphi[idx],
        )
        BnR_mode, Bnphi_mode, BnZ_mode = field_divfree(spline, grid_ref.R, grid_ref.Z, ntor=ntor)
        if BnR_mode.ndim == 2:
            BnR_mode = BnR_mode[None, :, :]
            Bnphi_mode = Bnphi_mode[None, :, :]
            BnZ_mode = BnZ_mode[None, :, :]
        BnR_sum[idx] = np.tensordot(weights, BnR_mode, axes=(0, 0))
        Bnphi_sum[idx] = np.tensordot(weights, Bnphi_mode, axes=(0, 0))
        BnZ_sum[idx] = np.tensordot(weights, BnZ_mode, axes=(0, 0))

    return BnR_sum, Bnphi_sum, BnZ_sum


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

    grid_ref, modes, BnR_ref, Bnphi_ref, BnZ_ref = load_fourier_sum(args.reference, currents)
    BnR_vec, Bnphi_vec, BnZ_vec = load_anvac_sum(args.vector, modes, grid_ref, currents)

    z_cm = np.linspace(-args.axis_range, args.axis_range, args.samples, dtype=float)
    zero = np.zeros_like(z_cm)

    Bx_ref, By_ref, Bz_ref = reconstruct_field_from_modes(
        BnR_ref, Bnphi_ref, BnZ_ref, modes, grid_ref.R, grid_ref.Z, zero, zero, z_cm
    )
    Bx_vec, By_vec, Bz_vec = reconstruct_field_from_modes(
        BnR_vec, Bnphi_vec, BnZ_vec, modes, grid_ref.R, grid_ref.Z, zero, zero, z_cm
    )

    B_ref_axis = Bz_ref
    B_vec_axis = Bz_vec
    B_direct_axis = compute_axis_field_direct(args.coil, currents, z_cm)
    B_analytic_axis = analytic_loop_field(radius_m, currents, z_cm)

    abs_diff_ref = np.abs(B_ref_axis - B_analytic_axis)
    rel_diff_ref = abs_diff_ref / np.maximum(np.abs(B_analytic_axis), 1e-12)

    abs_diff_vec = np.abs(B_vec_axis - B_analytic_axis)
    rel_diff_vec = abs_diff_vec / np.maximum(np.abs(B_analytic_axis), 1e-12)

    abs_diff_direct = np.abs(B_direct_axis - B_analytic_axis)
    rel_diff_direct = abs_diff_direct / np.maximum(np.abs(B_analytic_axis), 1e-12)

    max_abs_ref = float(np.max(abs_diff_ref))
    max_rel_ref = float(np.max(rel_diff_ref))
    max_abs_vec = float(np.max(abs_diff_vec))
    max_rel_vec = float(np.max(rel_diff_vec))
    max_abs_direct = float(np.max(abs_diff_direct))
    max_rel_direct = float(np.max(rel_diff_direct))

    summary = (
        "Axisymmetric loop validation\n"
        f"  Radius (m):            {radius_m:.6f}\n"
        f"  Samples:               {args.samples}\n"
        f"  Max abs error (Fourier): {max_abs_ref:.6e} G\n"
        f"  Max rel error (Fourier): {max_rel_ref:.6%}\n"
        f"  Max abs error (Anvac):   {max_abs_vec:.6e} G\n"
        f"  Max rel error (Anvac):   {max_rel_vec:.6%}\n"
        f"  Max abs error (Direct):  {max_abs_direct:.6e} G\n"
        f"  Max rel error (Direct):  {max_rel_direct:.6%}\n"
    )

    args.summary.write_text(summary)

    if args.plot is not None:
        fig, ax = plt.subplots(figsize=(7.5, 4.5), layout="constrained")
        ax.plot(z_cm, B_analytic_axis, label="Analytic", linestyle="--")
        ax.plot(z_cm, B_ref_axis, label="Fourier")
        ax.plot(z_cm, B_vec_axis, label="Anvac")
        ax.plot(z_cm, B_direct_axis, label="Direct (segments)")
        ax.set_xlabel("z [cm]")
        ax.set_ylabel("B_z [G]")
        ax.legend(loc="best")
        fig.savefig(args.plot, dpi=150, bbox_inches="tight")
        plt.close(fig)

    tol_abs = args.abs_tol
    tol_rel = args.rel_tol
    if max_abs_ref > tol_abs or max_rel_ref > tol_rel:
        raise RuntimeError("Fourier axisymmetric field exceeds tolerance")
    if max_abs_vec > tol_abs or max_rel_vec > tol_rel:
        raise RuntimeError("Anvac axisymmetric field exceeds tolerance")
    if max_abs_direct > tol_abs or max_rel_direct > tol_rel:
        raise RuntimeError("Direct axisymmetric field exceeds tolerance")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
