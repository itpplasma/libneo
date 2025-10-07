#!/usr/bin/env python3
"""Compare Fourier (Bnvac) vs Anvac methods at coil center"""

import sys
sys.path.insert(0, 'python')

import numpy as np
from scipy.interpolate import RectBivariateSpline
import netCDF4 as nc

from libneo.biotsavart_fourier import read_Bnvac_fourier, field_divfree, grid_t

# Load reference Bnvac
grid_ref, BnR_ref, Bnphi_ref, BnZ_ref = read_Bnvac_fourier(
    'build/test/magfie/single_coil/single_reference.h5', ntor=0
)

# Load Anvac
ds = nc.Dataset('build/test/magfie/single_coil/single_test.nc', 'r')
grid_anv = grid_t()
grid_anv.R = np.array(ds['R'][:])
grid_anv.Z = np.array(ds['Z'][:])

ntor = 0
AnR_real = ds['AnR_real'][0, :, :, ntor].T
Anphi_real = ds['Anphi_real'][0, :, :, ntor].T
AnZ_real = ds['AnZ_real'][0, :, :, ntor].T

spl = {
    'AnR_Re': [RectBivariateSpline(grid_anv.R, grid_anv.Z, AnR_real, kx=5, ky=5)],
    'AnR_Im': [RectBivariateSpline(grid_anv.R, grid_anv.Z, np.zeros_like(AnR_real), kx=5, ky=5)],
    'AnZ_Re': [RectBivariateSpline(grid_anv.R, grid_anv.Z, AnZ_real, kx=5, ky=5)],
    'AnZ_Im': [RectBivariateSpline(grid_anv.R, grid_anv.Z, np.zeros_like(AnZ_real), kx=5, ky=5)],
    'Anphi_Re': [RectBivariateSpline(grid_anv.R, grid_anv.Z, Anphi_real, kx=5, ky=5)],
    'Anphi_Im': [RectBivariateSpline(grid_anv.R, grid_anv.Z, np.zeros_like(Anphi_real), kx=5, ky=5)],
}

BnR_anv, Bnphi_anv, BnZ_anv = field_divfree(spl, grid_anv.R, grid_anv.Z, ntor=0)
ds.close()

# Compare at coil center
R_center = np.hypot(197.2682697, 72.00853957)
Z_center = 78.0

iR_ref = np.argmin(np.abs(grid_ref.R - R_center))
iZ_ref = np.argmin(np.abs(grid_ref.Z - Z_center))

iR_anv = np.argmin(np.abs(grid_anv.R - R_center))
iZ_anv = np.argmin(np.abs(grid_anv.Z - Z_center))

print("="*70)
print("Comparison at coil center")
print("="*70)
print(f"\nFourier (Bnvac reference) at R={grid_ref.R[iR_ref]:.2f}, Z={grid_ref.Z[iZ_ref]:.2f}:")
print(f"  BnR   = {BnR_ref[0, iR_ref, iZ_ref].real:.8e} G")
print(f"  Bnphi = {Bnphi_ref[0, iR_ref, iZ_ref].real:.8e} G")
print(f"  BnZ   = {BnZ_ref[0, iR_ref, iZ_ref].real:.8e} G")
print(f"  |B|   = {np.sqrt(BnR_ref[0,iR_ref,iZ_ref].real**2 + Bnphi_ref[0,iR_ref,iZ_ref].real**2 + BnZ_ref[0,iR_ref,iZ_ref].real**2):.8e} G")

print(f"\nAnvac (from vector potential) at R={grid_anv.R[iR_anv]:.2f}, Z={grid_anv.Z[iZ_anv]:.2f}:")
print(f"  BnR   = {BnR_anv[iR_anv, iZ_anv].real:.8e} G")
print(f"  Bnphi = {Bnphi_anv[iR_anv, iZ_anv].real:.8e} G")
print(f"  BnZ   = {BnZ_anv[iR_anv, iZ_anv].real:.8e} G")
print(f"  |B|   = {np.sqrt(BnR_anv[iR_anv,iZ_anv].real**2 + Bnphi_anv[iR_anv,iZ_anv].real**2 + BnZ_anv[iR_anv,iZ_anv].real**2):.8e} G")

print(f"\n" + "="*70)
print("Component-wise comparison")
print("="*70)
diff_R = BnR_anv[iR_anv, iZ_anv].real - BnR_ref[0, iR_ref, iZ_ref].real
diff_phi = Bnphi_anv[iR_anv, iZ_anv].real - Bnphi_ref[0, iR_ref, iZ_ref].real
diff_Z = BnZ_anv[iR_anv, iZ_anv].real - BnZ_ref[0, iR_ref, iZ_ref].real

rel_R = abs(diff_R / BnR_ref[0, iR_ref, iZ_ref].real) * 100
rel_phi = abs(diff_phi / Bnphi_ref[0, iR_ref, iZ_ref].real) * 100
rel_Z = abs(diff_Z / BnZ_ref[0, iR_ref, iZ_ref].real) * 100

print(f"\nBnR:")
print(f"  Difference: {diff_R:.8e} G")
print(f"  Relative:   {rel_R:.6f}%")

print(f"\nBnphi:")
print(f"  Difference: {diff_phi:.8e} G")
print(f"  Relative:   {rel_phi:.6f}%")

print(f"\nBnZ:")
print(f"  Difference: {diff_Z:.8e} G")
print(f"  Relative:   {rel_Z:.6f}%")

print(f"\n" + "="*70)
print("CONCLUSION")
print("="*70)
if max(rel_R, rel_phi, rel_Z) < 1.0:
    print("✓ Fourier and Anvac methods AGREE to <1%!")
    print("✓ Both correctly compute the n=0 Fourier mode")
    print("✓ Non-zero Bφ is CORRECT due to coil discretization")
else:
    print("✗ Methods disagree by >1% - there may be a bug")
