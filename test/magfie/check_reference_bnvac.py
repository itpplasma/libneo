#!/usr/bin/env python3
"""Check what the reference Bnvac contains for ntor=0"""

import sys
sys.path.insert(0, 'python')

import numpy as np
from libneo.biotsavart_fourier import read_Bnvac_fourier

# Load reference Bnvac for ntor=0
grid, BnR, Bnphi, BnZ = read_Bnvac_fourier('build/test/magfie/single_coil/single_reference.h5', ntor=0)

print("="*70)
print("Reference Bnvac (ntor=0) from H5 file")
print("="*70)
print(f"\nGrid: R=[{grid.R[0]:.1f}, {grid.R[-1]:.1f}] cm, Z=[{grid.Z[0]:.1f}, {grid.Z[-1]:.1f}] cm")
print(f"      nR={grid.nR}, nZ={grid.nZ}, ncoil={BnR.shape[0]}")

print(f"\nBnR statistics:")
print(f"  Real: min={BnR.real.min():.6e}, max={BnR.real.max():.6e} G")
print(f"  Imag: min={BnR.imag.min():.6e}, max={BnR.imag.max():.6e} G")

print(f"\nBnphi statistics:")
print(f"  Real: min={Bnphi.real.min():.6e}, max={Bnphi.real.max():.6e} G")
print(f"  Imag: min={Bnphi.imag.min():.6e}, max={Bnphi.imag.max():.6e} G")

print(f"\nBnZ statistics:")
print(f"  Real: min={BnZ.real.min():.6e}, max={BnZ.real.max():.6e} G")
print(f"  Imag: min={BnZ.imag.min():.6e}, max={BnZ.imag.max():.6e} G")

# Check at coil center
R_center = np.hypot(197.2682697, 72.00853957)
Z_center = 78.0

iR = np.argmin(np.abs(grid.R - R_center))
iZ = np.argmin(np.abs(grid.Z - Z_center))

print(f"\n" + "="*70)
print(f"At coil center (nearest grid point)")
print("="*70)
print(f"Grid point: R={grid.R[iR]:.2f} cm, Z={grid.Z[iZ]:.2f} cm")
print(f"\nReference Bnvac (ntor=0):")
print(f"  BnR[{iR},{iZ}]   = {BnR[0, iR, iZ].real:.6e} G")
print(f"  Bnphi[{iR},{iZ}] = {Bnphi[0, iR, iZ].real:.6e} G")
print(f"  BnZ[{iR},{iZ}]   = {BnZ[0, iR, iZ].real:.6e} G")
print(f"  |B| = {np.sqrt(BnR[0,iR,iZ].real**2 + Bnphi[0,iR,iZ].real**2 + BnZ[0,iR,iZ].real**2):.6e} G")

print(f"\n" + "="*70)
print("ANALYSIS")
print("="*70)
print("If Bnphi is zero or very small → reference represents true n=0 Fourier mode of B")
print("If Bnphi is large → reference has same issue as our computation")
print("\nExpected for true axisymmetric n=0 mode:")
print("  - Bnphi should be exactly zero everywhere")
print("  - BnR and BnZ should be smooth and axisymmetric")
