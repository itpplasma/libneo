#!/usr/bin/env python3
"""Debug cylindrical to Cartesian transformation at coil center"""

import sys
sys.path.insert(0, 'python')

import numpy as np
from numpy import hypot, arctan2, cos, sin

# Coil center in Cartesian coordinates (cm)
x_center = 197.2682697
y_center = 72.00853957
z_center = 78.0

# Convert to cylindrical
R_center = hypot(x_center, y_center)
phi_center = arctan2(y_center, x_center)
Z_center = z_center

print("="*70)
print("Coil center coordinates")
print("="*70)
print(f"Cartesian: x={x_center:.2f}, y={y_center:.2f}, z={z_center:.2f} cm")
print(f"Cylindrical: R={R_center:.2f}, φ={phi_center:.6f} rad = {np.degrees(phi_center):.2f}°, Z={Z_center:.2f} cm")

# Field in cylindrical at this point (from previous debug)
BR_cyl = 1.053591e-03  # G
Bphi_cyl = -9.458670e-03  # G  (NOTE: This is unusually large and negative!)
BZ_cyl = 5.030582e-03  # G

print(f"\n" + "="*70)
print("Magnetic field in cylindrical coordinates")
print("="*70)
print(f"BR   = {BR_cyl:.6e} G")
print(f"Bφ   = {Bphi_cyl:.6e} G  ← NOTE: Large and negative!")
print(f"BZ   = {BZ_cyl:.6e} G")
print(f"|B|  = {np.sqrt(BR_cyl**2 + Bphi_cyl**2 + BZ_cyl**2):.6e} G")

# Transform to Cartesian
# B_x = B_R * cos(φ) - B_φ * sin(φ)
# B_y = B_R * sin(φ) + B_φ * cos(φ)
# B_z = B_Z

cosphi = cos(phi_center)
sinphi = sin(phi_center)

Bx_cart = BR_cyl * cosphi - Bphi_cyl * sinphi
By_cart = BR_cyl * sinphi + Bphi_cyl * cosphi
Bz_cart = BZ_cyl

print(f"\n" + "="*70)
print("Transformation to Cartesian")
print("="*70)
print(f"cos(φ) = {cosphi:.6f}")
print(f"sin(φ) = {sinphi:.6f}")
print(f"\nBx = BR*cos(φ) - Bφ*sin(φ)")
print(f"   = {BR_cyl:.6e}*{cosphi:.6f} - ({Bphi_cyl:.6e})*{sinphi:.6f}")
print(f"   = {BR_cyl*cosphi:.6e} - ({Bphi_cyl*sinphi:.6e})")
print(f"   = {Bx_cart:.6e} G")
print(f"\nBy = BR*sin(φ) + Bφ*cos(φ)")
print(f"   = {BR_cyl:.6e}*{sinphi:.6f} + ({Bphi_cyl:.6e})*{cosphi:.6f}")
print(f"   = {BR_cyl*sinphi:.6e} + ({Bphi_cyl*cosphi:.6e})")
print(f"   = {By_cart:.6e} G")
print(f"\nBz = BZ = {Bz_cart:.6e} G")
print(f"\n|B_cart| = {np.sqrt(Bx_cart**2 + By_cart**2 + Bz_cart**2):.6e} G")

# Project onto axis direction
normal = np.array([0.35045589, -0.45058614, 0.82106808])
direction = normal / np.linalg.norm(normal)

B_cart = np.array([Bx_cart, By_cart, Bz_cart])
B_parallel = np.dot(B_cart, direction)

print(f"\n" + "="*70)
print("Project onto axis direction")
print("="*70)
print(f"Axis direction: n̂ = ({direction[0]:.6f}, {direction[1]:.6f}, {direction[2]:.6f})")
print(f"B·n̂ = Bx*nx + By*ny + Bz*nz")
print(f"    = {Bx_cart:.6e}*{direction[0]:.6f} + {By_cart:.6e}*{direction[1]:.6f} + {Bz_cart:.6e}*{direction[2]:.6f}")
print(f"    = {Bx_cart*direction[0]:.6e} + {By_cart*direction[1]:.6e} + {Bz_cart*direction[2]:.6e}")
print(f"    = {B_parallel:.6e} G")

print(f"\n" + "="*70)
print("ANALYSIS")
print("="*70)
print(f"\nExpected field magnitude at center: ~0.018 G")
print(f"Computed B_parallel: {B_parallel:.6e} G")
print(f"Ratio: {B_parallel/0.018:.3f}")
print(f"\nPROBLEM: Bφ component is HUGE ({Bphi_cyl:.6e} G) compared to BR and BZ!")
print(f"For an axisymmetric coil, Bφ should be ZERO or very small.")
print(f"This suggests:")
print(f"  1. Bφ formula in field_divfree is wrong")
print(f"  2. OR dAnphi/dR, dAnphi/dZ formulas have errors")
print(f"  3. OR Fortran is computing wrong derivatives")
print(f"\nRecall: For ntor=0, Bφ = ∂AR/∂Z - ∂AZ/∂R")
print(f"Let's check if this is being computed correctly...")
