#!/usr/bin/env python3
"""Debug field components at coil center to find sign error"""

import sys
sys.path.insert(0, 'python')

import numpy as np
from scipy.interpolate import RectBivariateSpline
import netCDF4 as nc

from libneo.biotsavart_fourier import grid_t, field_divfree

# Load data
ds = nc.Dataset('build/test/magfie/single_coil/single_test.nc', 'r')

grid = grid_t()
grid.R = np.array(ds['R'][:])
grid.Z = np.array(ds['Z'][:])

ntor = 0
AnR_real = ds['AnR_real'][0, :, :, ntor].T
Anphi_real = ds['Anphi_real'][0, :, :, ntor].T
AnZ_real = ds['AnZ_real'][0, :, :, ntor].T

print("Coil center location: (197.27, 72.01, 78.0) cm")
print("In cylindrical: R = sqrt(197.27² + 72.01²) = 210.00 cm, Z = 78.0 cm")

# Find nearest grid point to coil center
R_center = np.hypot(197.2682697, 72.00853957)
Z_center = 78.0

iR = np.argmin(np.abs(grid.R - R_center))
iZ = np.argmin(np.abs(grid.Z - Z_center))

print(f"\nNearest grid point: R = {grid.R[iR]:.2f} cm, Z = {grid.Z[iZ]:.2f} cm")
print(f"Indices: iR={iR}, iZ={iZ}")

# Vector potential at this point
print(f"\nVector potential at grid point:")
print(f"  AnR   = {AnR_real[iR, iZ]:.6e}")
print(f"  Anphi = {Anphi_real[iR, iZ]:.6e}")
print(f"  AnZ   = {AnZ_real[iR, iZ]:.6e}")

# Now compute B using manual derivatives
# For ntor=0: BR = -∂Aphi/∂Z, BZ = ∂Aphi/∂R + Aphi/R

# Manual finite difference
if iZ > 0 and iZ < len(grid.Z) - 1:
    dZ = grid.Z[iZ+1] - grid.Z[iZ-1]
    dAphi_dZ_manual = (Anphi_real[iR, iZ+1] - Anphi_real[iR, iZ-1]) / dZ
else:
    dAphi_dZ_manual = np.nan

if iR > 0 and iR < len(grid.R) - 1:
    dR = grid.R[iR+1] - grid.R[iR-1]
    dAphi_dR_manual = (Anphi_real[iR+1, iZ] - Anphi_real[iR-1, iZ]) / dR
else:
    dAphi_dR_manual = np.nan

print(f"\nManual finite difference derivatives:")
print(f"  ∂Aphi/∂R = {dAphi_dR_manual:.6e}")
print(f"  ∂Aphi/∂Z = {dAphi_dZ_manual:.6e}")

# Compute B components
BR_manual = -dAphi_dZ_manual
BZ_manual = dAphi_dR_manual + Anphi_real[iR, iZ] / grid.R[iR]

print(f"\nMagnetic field (manual calculation):")
print(f"  BR   = -∂Aphi/∂Z = {BR_manual:.6e} G")
print(f"  BZ   = ∂Aphi/∂R + Aphi/R = {dAphi_dR_manual:.6e} + {Anphi_real[iR, iZ]/grid.R[iR]:.6e} = {BZ_manual:.6e} G")

# Now use field_divfree and compare
AnR = np.zeros((1, len(grid.R), len(grid.Z)), dtype=complex)
Anphi = np.zeros((1, len(grid.R), len(grid.Z)), dtype=complex)
AnZ = np.zeros((1, len(grid.R), len(grid.Z)), dtype=complex)

AnR[0, :, :].real = AnR_real
Anphi[0, :, :].real = Anphi_real
AnZ[0, :, :].real = AnZ_real

spl = {
    'AnR_Re': [RectBivariateSpline(grid.R, grid.Z, AnR_real, kx=5, ky=5)],
    'AnR_Im': [RectBivariateSpline(grid.R, grid.Z, AnR[0].imag, kx=5, ky=5)],
    'AnZ_Re': [RectBivariateSpline(grid.R, grid.Z, AnZ_real, kx=5, ky=5)],
    'AnZ_Im': [RectBivariateSpline(grid.R, grid.Z, AnZ[0].imag, kx=5, ky=5)],
    'Anphi_Re': [RectBivariateSpline(grid.R, grid.Z, Anphi_real, kx=5, ky=5)],
    'Anphi_Im': [RectBivariateSpline(grid.R, grid.Z, Anphi[0].imag, kx=5, ky=5)],
}

BnR_field, Bnphi_field, BnZ_field = field_divfree(spl, grid.R, grid.Z, ntor=0)

print(f"\nMagnetic field from field_divfree:")
print(f"  BnR[{iR},{iZ}]   = {BnR_field[iR, iZ].real:.6e} G")
print(f"  Bnphi[{iR},{iZ}] = {Bnphi_field[iR, iZ].real:.6e} G")
print(f"  BnZ[{iR},{iZ}]   = {BnZ_field[iR, iZ].real:.6e} G")

print(f"\nComparison:")
print(f"  BR manual vs field_divfree: {BR_manual:.6e} vs {BnR_field[iR, iZ].real:.6e}")
print(f"  BZ manual vs field_divfree: {BZ_manual:.6e} vs {BnZ_field[iR, iZ].real:.6e}")

# Check stored derivatives from Fortran
dAphi_dR_stored = ds['dAnphi_dR_real'][0, :, :, ntor].T
dAphi_dZ_stored = ds['dAnphi_dZ_real'][0, :, :, ntor].T

print(f"\nStored derivatives from Fortran:")
print(f"  dAnphi/dR[{iR},{iZ}] = {dAphi_dR_stored[iR, iZ]:.6e}")
print(f"  dAnphi/dZ[{iR},{iZ}] = {dAphi_dZ_stored[iR, iZ]:.6e}")

print(f"\nUsing stored derivatives:")
BR_stored = -dAphi_dZ_stored[iR, iZ]
BZ_stored = dAphi_dR_stored[iR, iZ] + Anphi_real[iR, iZ] / grid.R[iR]
print(f"  BR = {BR_stored:.6e} G")
print(f"  BZ = {BZ_stored:.6e} G")

ds.close()

print(f"\n" + "="*70)
print("ANALYSIS")
print("="*70)
print("For a circular loop at R~210cm, Z~78cm with 1A current:")
print("- Field should point mostly in +Z direction (along loop normal)")
print("- Magnitude should be ~0.018 G at center")
print("\nIf we're getting negative values or wrong signs, check:")
print("1. Sign of curl formulas in field_divfree")
print("2. Sign conventions in Fortran derivatives")
print("3. Current direction in coil file")
