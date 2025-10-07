#!/usr/bin/env python3
"""Debug Bφ = ∂AR/∂Z - ∂AZ/∂R computation"""

import sys
sys.path.insert(0, 'python')

import numpy as np
import netCDF4 as nc

# Load data
ds = nc.Dataset('build/test/magfie/single_coil/single_test.nc', 'r')

grid_R = np.array(ds['R'][:])
grid_Z = np.array(ds['Z'][:])

ntor = 0
AnR_real = ds['AnR_real'][0, :, :, ntor].T  # (nR, nZ)
AnZ_real = ds['AnZ_real'][0, :, :, ntor].T

# Coil center
R_center = np.hypot(197.2682697, 72.00853957)
Z_center = 78.0

iR = np.argmin(np.abs(grid_R - R_center))
iZ = np.argmin(np.abs(grid_Z - Z_center))

print("="*70)
print("Computing Bφ = ∂AR/∂Z - ∂AZ/∂R at coil center")
print("="*70)
print(f"\nGrid point: R={grid_R[iR]:.2f} cm, Z={grid_Z[iZ]:.2f} cm (indices {iR},{iZ})")
print(f"\nVector potential components:")
print(f"  AR = {AnR_real[iR, iZ]:.6e}")
print(f"  AZ = {AnZ_real[iR, iZ]:.6e}")

# Compute derivatives using centered finite differences
if iZ > 0 and iZ < len(grid_Z) - 1:
    dZ = grid_Z[iZ+1] - grid_Z[iZ-1]
    dAR_dZ = (AnR_real[iR, iZ+1] - AnR_real[iR, iZ-1]) / dZ
else:
    dAR_dZ = np.nan

if iR > 0 and iR < len(grid_R) - 1:
    dR = grid_R[iR+1] - grid_R[iR-1]
    dAZ_dR = (AnZ_real[iR+1, iZ] - AnZ_real[iR-1, iZ]) / dR
else:
    dAZ_dR = np.nan

print(f"\nDerivatives:")
print(f"  ∂AR/∂Z = {dAR_dZ:.6e}")
print(f"  ∂AZ/∂R = {dAZ_dR:.6e}")

Bphi = dAR_dZ - dAZ_dR
print(f"\nBφ = ∂AR/∂Z - ∂AZ/∂R = {dAR_dZ:.6e} - {dAZ_dR:.6e} = {Bphi:.6e} G")

print(f"\n" + "="*70)
print("Analysis")
print("="*70)
print(f"\nFor a perfectly axisymmetric coil:")
print(f"  - AR and AZ should vary smoothly with R and Z")
print(f"  - But Bφ = ∂AR/∂Z - ∂AZ/∂R should be ZERO everywhere")
print(f"  - (Because ∇×A has only BR and BZ components for axisymmetric fields)")
print(f"\nMeasured Bφ = {Bphi:.6e} G")
print(f"This is NOT zero! This indicates:")
print(f"  1. The coil is not perfectly axisymmetric (has 3D structure)")
print(f"  2. OR there are numerical errors in AR, AZ from Fortran")
print(f"  3. OR the derivatives are not being computed correctly")

# Check magnitude relative to other components
from scipy.interpolate import RectBivariateSpline
Anphi_real = ds['Anphi_real'][0, :, :, ntor].T

spl_Anphi = RectBivariateSpline(grid_R, grid_Z, Anphi_real, kx=5, ky=5)
dAnphi_dR = spl_Anphi(grid_R[iR], grid_Z[iZ], dx=1).item()
dAnphi_dZ = spl_Anphi(grid_R[iR], grid_Z[iZ], dy=1).item()

BR = -dAnphi_dZ
BZ = dAnphi_dR + Anphi_real[iR, iZ] / grid_R[iR]

print(f"\nOther field components at this point:")
print(f"  BR = {BR:.6e} G")
print(f"  Bφ = {Bphi:.6e} G")
print(f"  BZ = {BZ:.6e} G")
print(f"\nRatios:")
print(f"  |Bφ|/|BR| = {abs(Bphi)/abs(BR):.2f}")
print(f"  |Bφ|/|BZ| = {abs(Bphi)/abs(BZ):.2f}")
print(f"\nBφ is {abs(Bphi)/abs(BR):.1f}x larger than BR!")
print(f"This is NOT expected for an axisymmetric coil.")

ds.close()

print(f"\n" + "="*70)
print("HYPOTHESIS")
print("="*70)
print("The coil geometry in single_coil.dat is a CIRCULAR LOOP,")
print("but it's discretized into 129 segments. This discretization")
print("creates small non-axisymmetric perturbations that produce")
print("a non-zero Bφ component.")
print("\nTo verify: compute the axisymmetric n=0 Fourier mode from")
print("the 3D field, which should remove the Bφ component.")
