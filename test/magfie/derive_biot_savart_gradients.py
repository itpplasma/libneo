#!/usr/bin/env python3
"""
Derive analytical gradients for Biot-Savart vector potential using SymPy.

The vector potential from a straight line segment is:
    A = (dl/|dl|) * log((1+e)/(1-e))

where:
    dl = segment vector (from coil point i to coil point f)
    r = evaluation position
    r_i = vector from evaluation point to segment start
    r_f = vector from evaluation point to segment end
    e = eccentricity = |dl|/(|r_i| + |r_f|)

We need ∂A/∂r (derivative w.r.t. evaluation position, NOT coil position!)
"""

import sympy as sp
from sympy import symbols, sqrt, log, simplify, diff, Matrix, latex

# Define symbolic variables
# Evaluation position
rx, ry, rz = symbols('r_x r_y r_z', real=True)
r = Matrix([rx, ry, rz])

# Coil segment endpoints
ix, iy, iz = symbols('i_x i_y i_z', real=True)  # Start point of segment
fx, fy, fz = symbols('f_x f_y f_z', real=True)  # End point of segment
coil_i = Matrix([ix, iy, iz])
coil_f = Matrix([fx, fy, fz])

# Segment vector (coil coordinate, not evaluation coordinate)
dl = coil_f - coil_i
dl_mag = sqrt(dl.dot(dl))

# Vectors from evaluation point to segment endpoints
# NOTE: These point FROM evaluation TO coil, so r_i = coil_i - r
r_i = coil_i - r
r_f = coil_f - r

# Distances
dist_i = sqrt(r_i.dot(r_i))
dist_f = sqrt(r_f.dot(r_f))

# Eccentricity
e = dl_mag / (dist_i + dist_f)

# Log term
log_term = log((1 + e) / (1 - e))

# Vector potential contribution (for one Cartesian component)
# A_x component contribution from this segment
A_x = (dl[0] / dl_mag) * log_term

print("=" * 80)
print("SYMBOLIC FORMULAS")
print("=" * 80)
print("\nSegment vector (dl):")
print(f"  dl = {dl.T}")
print(f"  |dl| = {dl_mag}")

print("\nVectors from evaluation point to coil:")
print(f"  r_i = coil_i - r = {r_i.T}")
print(f"  r_f = coil_f - r = {r_f.T}")

print("\nEccentricity:")
print(f"  e = |dl| / (|r_i| + |r_f|)")

print("\nVector potential A_x component:")
print(f"  A_x = (dl_x / |dl|) * log((1+e)/(1-e))")

print("\n" + "=" * 80)
print("COMPUTING GRADIENT: ∂A_x/∂r")
print("=" * 80)

# Compute gradient of A_x with respect to evaluation position
grad_Ax_rx = diff(A_x, rx)
grad_Ax_ry = diff(A_x, ry)
grad_Ax_rz = diff(A_x, rz)

print("\n∂A_x/∂r_x:")
print(grad_Ax_rx)

print("\n∂A_x/∂r_y:")
print(grad_Ax_ry)

print("\n∂A_x/∂r_z:")
print(grad_Ax_rz)

print("\n" + "=" * 80)
print("SIMPLIFYING...")
print("=" * 80)

# Try to simplify
grad_Ax_rx_simp = simplify(grad_Ax_rx)
grad_Ax_ry_simp = simplify(grad_Ax_ry)
grad_Ax_rz_simp = simplify(grad_Ax_rz)

print("\nSimplified ∂A_x/∂r_x:")
print(grad_Ax_rx_simp)

print("\nSimplified ∂A_x/∂r_y:")
print(grad_Ax_ry_simp)

print("\nSimplified ∂A_x/∂r_z:")
print(grad_Ax_rz_simp)

print("\n" + "=" * 80)
print("CHECKING CURRENT IMPLEMENTATION")
print("=" * 80)

print("\nCurrent code computes:")
print("  common_gradient_term = (r_i/|r_i| + r_f/|r_f|) / (|r_i|*|r_f| + r_i·r_f)")
print("  grad_AX = dl_x * common_gradient_term")
print("\nThis gives grad_AX as a 3-vector [∂A_x/∂x, ∂A_x/∂y, ∂A_x/∂z]")

# Current implementation formula
common_term = (r_i / dist_i + r_f / dist_f) / (dist_i * dist_f + r_i.dot(r_f))

grad_Ax_current = dl[0] * common_term

print("\nCurrent formula for grad_A_x:")
for i, comp in enumerate(['x', 'y', 'z']):
    print(f"  ∂A_x/∂{comp} = dl_x * common_term[{i}]")
    print(f"           = {grad_Ax_current[i]}")

print("\n" + "=" * 80)
print("COMPARING SYMPY vs CURRENT IMPLEMENTATION")
print("=" * 80)

# Compare
diff_x = simplify(grad_Ax_rx_simp - grad_Ax_current[0])
diff_y = simplify(grad_Ax_ry_simp - grad_Ax_current[1])
diff_z = simplify(grad_Ax_rz_simp - grad_Ax_current[2])

print("\nDifference in ∂A_x/∂r_x:")
print(f"  SymPy - Current = {diff_x}")

print("\nDifference in ∂A_x/∂r_y:")
print(f"  SymPy - Current = {diff_y}")

print("\nDifference in ∂A_x/∂r_z:")
print(f"  SymPy - Current = {diff_z}")

if diff_x == 0 and diff_y == 0 and diff_z == 0:
    print("\n✓ FORMULAS MATCH!")
else:
    print("\n✗ FORMULAS DO NOT MATCH!")
    print("\nCorrect gradient formula should be:")
    print(f"  ∂A_x/∂r_x = {grad_Ax_rx_simp}")
    print(f"  ∂A_x/∂r_y = {grad_Ax_ry_simp}")
    print(f"  ∂A_x/∂r_z = {grad_Ax_rz_simp}")

print("\n" + "=" * 80)
print("DERIVATIVE OF LOG TERM")
print("=" * 80)

# Focus on the derivative of the log term
print("\nLet's examine ∂/∂r_x [log((1+e)/(1-e))]")

de_drx = diff(e, rx)
print(f"\n∂e/∂r_x = {simplify(de_drx)}")

dlog_de = diff(log((1+e)/(1-e)), e)
print(f"\n∂/∂e [log((1+e)/(1-e))] = {simplify(dlog_de)}")

dlog_drx = diff(log_term, rx)
print(f"\n∂/∂r_x [log((1+e)/(1-e))] = {simplify(dlog_drx)}")
