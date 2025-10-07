#!/usr/bin/env python3
"""Debug axis field evaluation"""

import numpy as np
import netCDF4 as nc

# Load the test file
ds = nc.Dataset('build/test/magfie/single_coil/single_test.nc', 'r')

print("NetCDF variables:")
for var in ds.variables:
    print(f"  {var}: {ds[var].shape}")
print()

# Get grid (in cm)
R = ds['R'][:]
Z = ds['Z'][:]
print(f"R grid: {R.min():.2f} to {R.max():.2f} cm, {len(R)} points")
print(f"Z grid: {Z.min():.2f} to {Z.max():.2f} cm, {len(Z)} points")
print()

# Get fields for ntor=0
# Variables are (coil, Z, R, ntor) - note the axis order!
AnR_real = ds['AnR_real'][:]
AnR_imag = ds['AnR_imag'][:]
AnZ_real = ds['AnZ_real'][:]
AnZ_imag = ds['AnZ_imag'][:]
print(f"AnR_real shape: {AnR_real.shape}")
print()

# Check ntor=0 mode (last axis is ntor)
ntor_idx = 0
AnR_n0 = AnR_real[0, :, :, ntor_idx] + 1j * AnR_imag[0, :, :, ntor_idx]  # (nZ, nR)
AnZ_n0 = AnZ_real[0, :, :, ntor_idx] + 1j * AnZ_imag[0, :, :, ntor_idx]

print(f"AnR(ntor=0) statistics:")
print(f"  Real: min={np.real(AnR_n0).min():.6e}, max={np.real(AnR_n0).max():.6e}")
print(f"  Imag: min={np.imag(AnR_n0).min():.6e}, max={np.imag(AnR_n0).max():.6e}")
print()

print(f"AnZ(ntor=0) statistics:")
print(f"  Real: min={np.real(AnZ_n0).min():.6e}, max={np.real(AnZ_n0).max():.6e}")
print(f"  Imag: min={np.imag(AnZ_n0).min():.6e}, max={np.imag(AnZ_n0).max():.6e}")
print()

# Check if we have Aphi
if 'Anphi_real' in ds.variables:
    Anphi_real = ds['Anphi_real'][:]
    Anphi_imag = ds['Anphi_imag'][:]
    Anphi_n0 = Anphi_real[0, :, :, ntor_idx] + 1j * Anphi_imag[0, :, :, ntor_idx]
    print(f"Anphi(ntor=0) statistics:")
    print(f"  Real: min={np.real(Anphi_n0).min():.6e}, max={np.real(Anphi_n0).max():.6e}")
    print(f"  Imag: min={np.imag(Anphi_n0).min():.6e}, max={np.imag(Anphi_n0).max():.6e}")
else:
    print("WARNING: Anphi not found in NetCDF file!")

ds.close()
