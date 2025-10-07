#!/usr/bin/env python3
"""Debug B field computed from Aphi"""

import numpy as np
import netCDF4 as nc
import sys
sys.path.insert(0, 'python')

from libneo.biotsavart_fourier import field_divfree, grid_t

# Load the test file
ds = nc.Dataset('build/test/magfie/single_coil/single_test.nc', 'r')

# Get grid
grid = grid_t()
grid.R = np.array(ds['R'][:])
grid.Z = np.array(ds['Z'][:])
grid.nR = len(grid.R)
grid.nZ = len(grid.Z)

print(f"Grid: R=[{grid.R[0]:.1f}, {grid.R[-1]:.1f}] cm, Z=[{grid.Z[0]:.1f}, {grid.Z[-1]:.1f}] cm")
print(f"      nR={grid.nR}, nZ={grid.nZ}")
print()

# Load Anvac for ntor=0
ntor = 0
ncoil = 1

AnR = np.empty((ncoil, grid.nR, grid.nZ), dtype=complex)
Anphi = np.empty((ncoil, grid.nR, grid.nZ), dtype=complex)
AnZ = np.empty((ncoil, grid.nR, grid.nZ), dtype=complex)

AnR[0, :, :].real = ds['AnR_real'][0, :, :, ntor].T
AnR[0, :, :].imag = ds['AnR_imag'][0, :, :, ntor].T
Anphi[0, :, :].real = ds['Anphi_real'][0, :, :, ntor].T
Anphi[0, :, :].imag = ds['Anphi_imag'][0, :, :, ntor].T
AnZ[0, :, :].real = ds['AnZ_real'][0, :, :, ntor].T
AnZ[0, :, :].imag = ds['AnZ_imag'][0, :, :, ntor].T

ds.close()

print("Vector potential statistics (ntor=0, single coil):")
print(f"  AnR:   {AnR.real.min():.6e} to {AnR.real.max():.6e}")
print(f"  Anphi: {Anphi.real.min():.6e} to {Anphi.real.max():.6e}")
print(f"  AnZ:   {AnZ.real.min():.6e} to {AnZ.real.max():.6e}")
print()

# (Don't need stored derivatives for this test)

# For ntor=0, gauging doesn't change anything
gauged_AnR = AnR
gauged_AnZ = AnZ

# Create splines including Aphi
from scipy.interpolate import RectBivariateSpline

spl = {
    'AnR_Re': [RectBivariateSpline(grid.R, grid.Z, AnR[0].real, kx=5, ky=5)],
    'AnR_Im': [RectBivariateSpline(grid.R, grid.Z, AnR[0].imag, kx=5, ky=5)],
    'AnZ_Re': [RectBivariateSpline(grid.R, grid.Z, AnZ[0].real, kx=5, ky=5)],
    'AnZ_Im': [RectBivariateSpline(grid.R, grid.Z, AnZ[0].imag, kx=5, ky=5)],
    'Anphi_Re': [RectBivariateSpline(grid.R, grid.Z, Anphi[0].real, kx=5, ky=5)],
    'Anphi_Im': [RectBivariateSpline(grid.R, grid.Z, Anphi[0].imag, kx=5, ky=5)],
}

# Evaluate B at grid points
BnR, Bnphi, BnZ = field_divfree(spl, grid.R, grid.Z, ntor=0)

print("Magnetic field statistics (ntor=0, single coil):")
print(f"  BnR:   {BnR.real.min():.6e} to {BnR.real.max():.6e} Gauss")
print(f"  Bnphi: {Bnphi.real.min():.6e} to {Bnphi.real.max():.6e} Gauss")
print(f"  BnZ:   {BnZ.real.min():.6e} to {BnZ.real.max():.6e} Gauss")
print()

# Check at coil center
coil_center = np.array([197.2682697, 72.00853957, 78.0])  # cm
R_center = np.hypot(coil_center[0], coil_center[1])
Z_center = coil_center[2]

print(f"Evaluating at coil center: R={R_center:.2f} cm, Z={Z_center:.2f} cm")

# Find nearest grid point
iR = np.argmin(np.abs(grid.R - R_center))
iZ = np.argmin(np.abs(grid.Z - Z_center))

print(f"  Nearest grid point: R={grid.R[iR]:.2f} cm, Z={grid.Z[iZ]:.2f} cm")
print(f"  BnR   = {BnR[iR, iZ].real:.6e} Gauss")
print(f"  Bnphi = {Bnphi[iR, iZ].real:.6e} Gauss")
print(f"  BnZ   = {BnZ[iR, iZ].real:.6e} Gauss")
print(f"  |B|   = {np.sqrt(BnR[iR, iZ].real**2 + Bnphi[iR, iZ].real**2 + BnZ[iR, iZ].real**2):.6e} Gauss")
