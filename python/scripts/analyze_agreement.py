#!/usr/bin/env python3
"""Detailed analysis of agreement between GPEC Fourier and coil_tools implementations."""

import sys
import numpy as np
from pathlib import Path

try:
    from libneo.biotsavart_fourier import (
        read_Bnvac_fourier,
        read_Anvac_fourier,
        gauged_Anvac_from_Bnvac,
        gauge_Anvac,
        spline_gauged_Anvac,
        field_divfree,
    )
except ImportError:
    print("Error: libneo module not found")
    sys.exit(1)

# Read files
ref_file = Path("test/magfie/test_data/aug_reference.h5")
test_file = Path("test/magfie/test_data/aug_test.nc")

print("=" * 70)
print("AGREEMENT ANALYSIS: GPEC Fourier vs coil_tools vector_potential")
print("=" * 70)

ref_grid, BnR_ref, Bnphi_ref, BnZ_ref = read_Bnvac_fourier(str(ref_file), ntor=2)
test_grid, AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ = read_Anvac_fourier(str(test_file), ntor=2)

print(f"\nReference (GPEC Fourier) grid: {ref_grid.nR} × {ref_grid.nZ}")
print(f"Test (coil_tools) grid:        {test_grid.nR} × {test_grid.nZ}")
print(f"Number of coils: {BnR_ref.shape[0]}")

# Create gauged vector potentials from both
ref_spl = spline_gauged_Anvac(
    ref_grid,
    *gauged_Anvac_from_Bnvac(ref_grid, BnR_ref, BnZ_ref, ntor=2),
    ntor=2
)

test_spl = spline_gauged_Anvac(
    test_grid,
    *gauge_Anvac(test_grid, AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ, ntor=2),
    ntor=2
)

# Evaluate on same grid
R = np.linspace(ref_grid.R_min, ref_grid.R_max, 2 * ref_grid.nR - 1)
Z = np.linspace(ref_grid.Z_min, ref_grid.Z_max, 2 * ref_grid.nZ - 1)

print(f"\nEvaluation grid: {len(R)} × {len(Z)}")

BnR_ref_eval, Bnphi_ref_eval, BnZ_ref_eval = field_divfree(ref_spl, R, Z, ntor=2)
BnR_test_eval, Bnphi_test_eval, BnZ_test_eval = field_divfree(test_spl, R, Z, ntor=2)

# Compute field magnitudes
B_ref_mag = np.sqrt((BnR_ref_eval * np.conj(BnR_ref_eval) +
                     Bnphi_ref_eval * np.conj(Bnphi_ref_eval) +
                     BnZ_ref_eval * np.conj(BnZ_ref_eval)).real)

B_test_mag = np.sqrt((BnR_test_eval * np.conj(BnR_test_eval) +
                      Bnphi_test_eval * np.conj(Bnphi_test_eval) +
                      BnZ_test_eval * np.conj(BnZ_test_eval)).real)

# Compute relative errors
relative_error = np.abs(B_test_mag - B_ref_mag) / (B_ref_mag + 1e-15)

print("\n" + "=" * 70)
print("FIELD MAGNITUDE COMPARISON")
print("=" * 70)

ncoil = BnR_ref_eval.shape[0]
for kcoil in range(ncoil):
    ref_mag = B_ref_mag[kcoil]
    test_mag = B_test_mag[kcoil]
    rel_err = relative_error[kcoil]

    # Statistics
    mean_ref = np.mean(ref_mag)
    mean_test = np.mean(test_mag)
    median_err = np.median(rel_err)
    mean_err = np.mean(rel_err)
    max_err = np.max(rel_err)
    p95_err = np.percentile(rel_err, 95)

    print(f"\nCoil {kcoil + 1}:")
    print(f"  Mean |B| (GPEC):      {mean_ref:.6e}")
    print(f"  Mean |B| (coil_tools): {mean_test:.6e}")
    print(f"  Relative error (median): {median_err*100:.4f}%")
    print(f"  Relative error (mean):   {mean_err*100:.4f}%")
    print(f"  Relative error (95th):   {p95_err*100:.4f}%")
    print(f"  Relative error (max):    {max_err*100:.4f}%")

# Overall statistics
print("\n" + "=" * 70)
print("OVERALL STATISTICS (all coils)")
print("=" * 70)
print(f"Median relative error: {np.median(relative_error)*100:.4f}%")
print(f"Mean relative error:   {np.mean(relative_error)*100:.4f}%")
print(f"95th percentile error: {np.percentile(relative_error, 95)*100:.4f}%")
print(f"Maximum error:         {np.max(relative_error)*100:.4f}%")

# Check component-wise
print("\n" + "=" * 70)
print("COMPONENT-WISE COMPARISON (Sample coil)")
print("=" * 70)
kcoil = 0
for comp_name, ref_comp, test_comp in [
    ("B_R", BnR_ref_eval[kcoil], BnR_test_eval[kcoil]),
    ("B_phi", Bnphi_ref_eval[kcoil], Bnphi_test_eval[kcoil]),
    ("B_Z", BnZ_ref_eval[kcoil], BnZ_test_eval[kcoil]),
]:
    ref_abs = np.abs(ref_comp)
    test_abs = np.abs(test_comp)
    rel_err = np.abs(test_abs - ref_abs) / (ref_abs + 1e-15)

    print(f"\n{comp_name}:")
    print(f"  Mean |{comp_name}| (GPEC):      {np.mean(ref_abs):.6e}")
    print(f"  Mean |{comp_name}| (coil_tools): {np.mean(test_abs):.6e}")
    print(f"  Median relative error: {np.median(rel_err)*100:.4f}%")
    print(f"  Max relative error:    {np.max(rel_err)*100:.4f}%")
