#!/usr/bin/env python3
"""Inspect the single-coil reference after summing every Fourier mode."""

import sys
sys.path.insert(0, 'python')

import numpy as np

from libneo.biotsavart_fourier import read_Bnvac_fourier_all, reconstruct_field_from_modes

REFERENCE_H5 = 'build/test/magfie/single_coil/single_reference.h5'
CURRENT_FILE = 'test/magfie/test_data/single_coil_currents.txt'
PREFactor = 0.1

grid, mode_numbers, BnR, Bnphi, BnZ = read_Bnvac_fourier_all(REFERENCE_H5)
currents = np.loadtxt(CURRENT_FILE)
weights = currents * PREFactor

BnR_sum = np.tensordot(weights, BnR, axes=(0, 1))
Bnphi_sum = np.tensordot(weights, Bnphi, axes=(0, 1))
BnZ_sum = np.tensordot(weights, BnZ, axes=(0, 1))

print("=" * 70)
print("Reference Bnvac (all modes summed with weights)")
print("=" * 70)
print(f"Modes stored: {mode_numbers}")
print(f"Grid R=[{grid.R[0]:.1f}, {grid.R[-1]:.1f}] cm, Z=[{grid.Z[0]:.1f}, {grid.Z[-1]:.1f}] cm")

R_center = 197.2682697
Y_center = 72.00853957
Z_center = 78.0
Bx, By, Bz = reconstruct_field_from_modes(
    BnR_sum,
    Bnphi_sum,
    BnZ_sum,
    mode_numbers,
    grid.R,
    grid.Z,
    R_center,
    Y_center,
    Z_center,
)

print(f"\nTotal field at coil centre (weighted sum over modes):")
print(f"  Bx = {Bx:.6e} G")
print(f"  By = {By:.6e} G")
print(f"  Bz = {Bz:.6e} G")
print(f"  |B| = {np.sqrt(Bx**2 + By**2 + Bz**2):.6e} G")

if abs(By) > 1e-4:
    print("\nNon-zero By confirms toroidal component persists after mode summation (discretised coil).")
else:
    print("\nBy ≈ 0: check discretisation assumptions.")
