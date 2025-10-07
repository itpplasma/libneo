#!/usr/bin/env python3
"""Debug script to trace where the 10x error occurs in axis plotting"""

import sys
sys.path.insert(0, 'python')

import numpy as np
from numpy import asarray, linspace, outer, pi, hypot, arctan2, cos, sin
from numpy.linalg import norm
import netCDF4 as nc
from scipy.interpolate import RectBivariateSpline

from libneo.biotsavart_fourier import grid_t, field_divfree

# ============================================================================
# Step 1: Verify analytic formula computation
# ============================================================================
print("="*70)
print("STEP 1: Verify analytic formula computation")
print("="*70)

origin = asarray([197.2682697, 72.00853957, 78.0])  # cm
normal = asarray([0.35045589, -0.45058614, 0.82106808])
direction = normal / norm(normal)
coil_radius = 35.0  # cm
axis_range = 60.0  # cm
axis_samples = 181
current_ampere = 1.0  # A

s_vals = linspace(-axis_range, axis_range, axis_samples)
mu0 = 4e-7 * pi
radius_m = coil_radius / 100.0
s_vals_m = s_vals / 100.0

B_analytic = mu0 * current_ampere * radius_m**2 / (2.0 * (radius_m**2 + s_vals_m**2)**1.5)
B_analytic *= 1.0e4  # Tesla -> Gauss

center_idx = axis_samples // 2
print(f"\nAxis parameters:")
print(f"  Origin: {origin} cm")
print(f"  Direction: {direction}")
print(f"  Coil radius: {coil_radius} cm")
print(f"  Current: {current_ampere} A")
print(f"\nAnalytic formula results:")
print(f"  B_analytic at center (s=0): {B_analytic[center_idx]:.8e} G")
print(f"  B_analytic max: {B_analytic.max():.8e} G")
print(f"  Expected at center: 0.01795196 G")

# ============================================================================
# Step 2: Load and evaluate computed fields
# ============================================================================
print("\n" + "="*70)
print("STEP 2: Load NetCDF data and evaluate computed fields")
print("="*70)

ds = nc.Dataset('build/test/magfie/single_coil/single_test.nc', 'r')

grid = grid_t()
grid.R = np.array(ds['R'][:])
grid.Z = np.array(ds['Z'][:])
grid.nR = len(grid.R)
grid.nZ = len(grid.Z)

print(f"\nGrid: R=[{grid.R[0]:.1f}, {grid.R[-1]:.1f}] cm, Z=[{grid.Z[0]:.1f}, {grid.Z[-1]:.1f}] cm")

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

# Create splines
spl = {
    'AnR_Re': [RectBivariateSpline(grid.R, grid.Z, AnR[0].real, kx=5, ky=5)],
    'AnR_Im': [RectBivariateSpline(grid.R, grid.Z, AnR[0].imag, kx=5, ky=5)],
    'AnZ_Re': [RectBivariateSpline(grid.R, grid.Z, AnZ[0].real, kx=5, ky=5)],
    'AnZ_Im': [RectBivariateSpline(grid.R, grid.Z, AnZ[0].imag, kx=5, ky=5)],
    'Anphi_Re': [RectBivariateSpline(grid.R, grid.Z, Anphi[0].real, kx=5, ky=5)],
    'Anphi_Im': [RectBivariateSpline(grid.R, grid.Z, Anphi[0].imag, kx=5, ky=5)],
}

# Evaluate B at grid points
BnR_grid, Bnphi_grid, BnZ_grid = field_divfree(spl, grid.R, grid.Z, ntor=0)

print(f"\nMagnetic field on grid:")
print(f"  BnR:   {BnR_grid.real.min():.6e} to {BnR_grid.real.max():.6e} G")
print(f"  Bnphi: {Bnphi_grid.real.min():.6e} to {Bnphi_grid.real.max():.6e} G")
print(f"  BnZ:   {BnZ_grid.real.min():.6e} to {BnZ_grid.real.max():.6e} G")

# ============================================================================
# Step 3: Evaluate computed fields along axis
# ============================================================================
print("\n" + "="*70)
print("STEP 3: Evaluate computed fields along axis points")
print("="*70)

points = origin[np.newaxis, :] + outer(s_vals, direction)
x_eval = points[:, 0]
y_eval = points[:, 1]
z_eval = points[:, 2]

B_computed = np.zeros(axis_samples)

# Manual interpolation and transformation
for idx in range(axis_samples):
    x, y, z = x_eval[idx], y_eval[idx], z_eval[idx]
    R_val = hypot(x, y)
    phi_val = arctan2(y, x)

    # Bilinear interpolation
    BR = spl['AnR_Re'][0](R_val, z).item() + 1j * spl['AnR_Im'][0](R_val, z).item()
    Bphi = spl['Anphi_Re'][0](R_val, z).item() + 1j * spl['Anphi_Im'][0](R_val, z).item()
    BZ = spl['AnZ_Re'][0](R_val, z).item() + 1j * spl['AnZ_Im'][0](R_val, z).item()

    # For ntor=0, take real part
    BR_total = BR.real
    Bphi_total = Bphi.real
    BZ_total = BZ.real

    # Transform to Cartesian
    cosphi = cos(phi_val)
    sinphi = sin(phi_val)
    Bx = BR_total * cosphi - Bphi_total * sinphi
    By = BR_total * sinphi + Bphi_total * cosphi

    # Project onto axis direction
    B_computed[idx] = Bx * direction[0] + By * direction[1] + BZ_total * direction[2]

print(f"\nComputed field along axis:")
print(f"  B_computed at center: {B_computed[center_idx]:.8e} G")
print(f"  B_computed max: {B_computed.max():.8e} G")
print(f"  B_computed min: {B_computed.min():.8e} G")

# ============================================================================
# Step 4: Compare analytic vs computed
# ============================================================================
print("\n" + "="*70)
print("STEP 4: Compare analytic vs computed")
print("="*70)

ratio_at_center = B_computed[center_idx] / B_analytic[center_idx]
ratio_max = B_computed.max() / B_analytic.max()

print(f"\nAt center (s=0):")
print(f"  Analytic: {B_analytic[center_idx]:.8e} G")
print(f"  Computed: {B_computed[center_idx]:.8e} G")
print(f"  Ratio (computed/analytic): {ratio_at_center:.6f}")
print(f"\nMaximum values:")
print(f"  Analytic max: {B_analytic.max():.8e} G")
print(f"  Computed max: {B_computed.max():.8e} G")
print(f"  Ratio: {ratio_max:.6f}")

# ============================================================================
# Step 5: Check what the plot would show
# ============================================================================
print("\n" + "="*70)
print("STEP 5: Expected plot values")
print("="*70)

print(f"\nIf plotted correctly:")
print(f"  Analytic curve peak should be at: ~{B_analytic.max():.5f} G")
print(f"  Computed curve peak should be at: ~{B_computed.max():.5f} G")
print(f"\nActual plot shows:")
print(f"  Analytic curve peak at: ~0.00175 G (10x too small!)")
print(f"  Computed curve peak at: ~0.0009 G")
print(f"\nConclusions:")
print(f"  1. Analytic formula computes correct value: {B_analytic.max():.8e} G")
print(f"  2. But plot shows 10x smaller: ~{B_analytic.max()/10:.8e} G")
print(f"  3. Computed field is {ratio_at_center:.2f}x of analytic (should be ~1.0)")
print(f"  4. There are TWO separate bugs:")
print(f"     a) Analytic curve plotted with 10x error")
print(f"     b) Computed field is ~{1/ratio_at_center:.1f}x too small")
