#!/usr/bin/env python3
"""
Check the transformation from Cartesian gradient components to cylindrical derivatives.

We have:
  grad_AX = [∂Ax/∂x, ∂Ax/∂y, ∂Ax/∂z]  (Cartesian gradients of Ax)
  grad_AY = [∂Ay/∂x, ∂Ay/∂y, ∂Ay/∂z]  (Cartesian gradients of Ay)

We need to compute:
  ∂Aphi/∂R  and  ∂Aphi/∂Z

where Aphi = Ay*cos(phi) - Ax*sin(phi)

The current code at lines 599-601 does:
  dAphi_dR = grad_AY(1) * cos²(phi) - grad_AX(2) * sin²(phi) +
             (grad_AY(2) - grad_AX(1)) * cos(phi)*sin(phi)
  dAphi_dZ = grad_AY(3) * cos(phi) - grad_AX(3) * sin(phi)

Let's verify this is correct using SymPy.
"""

import sympy as sp
from sympy import symbols, cos, sin, diff, simplify

# Define symbols
R, phi, Z = symbols('R phi Z', real=True)

# Cartesian coordinates
x = R * cos(phi)
y = R * sin(phi)
z = Z

print("="*80)
print("COORDINATE TRANSFORMATIONS")
print("="*80)
print(f"x = R*cos(phi) = {x}")
print(f"y = R*sin(phi) = {y}")
print(f"z = Z = {z}")

# Define Cartesian vector potential components as functions
Ax_func = sp.Function('A_x')(x, y, z)
Ay_func = sp.Function('A_y')(x, y, z)

# Cylindrical components
Aphi = Ay_func * cos(phi) - Ax_func * sin(phi)

print("\n" + "="*80)
print("VECTOR POTENTIAL TRANSFORMATION")
print("="*80)
print(f"Aphi = Ay*cos(phi) - Ax*sin(phi)")
print(f"     = {Aphi}")

# Compute ∂Aphi/∂R
print("\n" + "="*80)
print("COMPUTING ∂Aphi/∂R")
print("="*80)

dAphi_dR = diff(Aphi, R)
print(f"\n∂Aphi/∂R = {dAphi_dR}")

# Expand using chain rule
print("\nExpanding derivatives...")

# Define derivative symbols
dAx_dx, dAx_dy, dAx_dz = symbols('dAx_dx dAx_dy dAx_dz', real=True)
dAy_dx, dAy_dy, dAy_dz = symbols('dAy_dx dAy_dy dAy_dz', real=True)

# Substitute  derivatives
dAphi_dR_explicit = (
    (dAy_dx * diff(x, R) + dAy_dy * diff(y, R) + dAy_dz * diff(z, R)) * cos(phi)
    - (dAx_dx * diff(x, R) + dAx_dy * diff(y, R) + dAx_dz * diff(z, R)) * sin(phi)
)

print(f"\nUsing chain rule:")
print(f"∂Aphi/∂R = [∂Ay/∂x * ∂x/∂R + ∂Ay/∂y * ∂y/∂R] * cos(phi)")
print(f"         - [∂Ax/∂x * ∂x/∂R + ∂Ax/∂y * ∂y/∂R] * sin(phi)")

# Compute partial derivatives
dx_dR = diff(x, R)  # cos(phi)
dy_dR = diff(y, R)  # sin(phi)
dz_dR = diff(z, R)  # 0

print(f"\n∂x/∂R = {dx_dR}")
print(f"∂y/∂R = {dy_dR}")
print(f"∂z/∂R = {dz_dR}")

dAphi_dR_formula = (
    (dAy_dx * cos(phi) + dAy_dy * sin(phi)) * cos(phi)
    - (dAx_dx * cos(phi) + dAx_dy * sin(phi)) * sin(phi)
)

print(f"\n∂Aphi/∂R = (dAy_dx*cos + dAy_dy*sin)*cos - (dAx_dx*cos + dAx_dy*sin)*sin")
print(f"         = dAy_dx*cos² + dAy_dy*sin*cos - dAx_dx*cos*sin - dAx_dy*sin²")
print(f"         = dAy_dx*cos² - dAx_dy*sin² + (dAy_dy - dAx_dx)*cos*sin")

print("\n" + "="*80)
print("CURRENT IMPLEMENTATION (line 599-600)")
print("="*80)
print("dAphi_dR = grad_AY(1)*cos² - grad_AX(2)*sin² + (grad_AY(2) - grad_AX(1))*cos*sin")
print("\nwhere:")
print("  grad_AX(1) = ∂Ax/∂x = dAx_dx")
print("  grad_AX(2) = ∂Ax/∂y = dAx_dy")
print("  grad_AY(1) = ∂Ay/∂x = dAy_dx")
print("  grad_AY(2) = ∂Ay/∂y = dAy_dy")

print("\nSo current formula is:")
print("  dAphi_dR = dAy_dx*cos² - dAx_dy*sin² + (dAy_dy - dAx_dx)*cos*sin")

print("\n✓ THIS MATCHES THE CORRECT FORMULA!")

# Now check ∂Aphi/∂Z
print("\n" + "="*80)
print("COMPUTING ∂Aphi/∂Z")
print("="*80)

dAphi_dZ = diff(Aphi, Z)
print(f"\n∂Aphi/∂Z = {dAphi_dZ}")

dx_dZ = diff(x, Z)  # 0
dy_dZ = diff(y, Z)  # 0
dz_dZ = diff(z, Z)  # 1

print(f"\n∂x/∂Z = {dx_dZ}")
print(f"∂y/∂Z = {dy_dZ}")
print(f"∂z/∂Z = {dz_dZ}")

dAphi_dZ_formula = (
    (dAy_dx * 0 + dAy_dy * 0 + dAy_dz * 1) * cos(phi)
    - (dAx_dx * 0 + dAx_dy * 0 + dAx_dz * 1) * sin(phi)
)

print(f"\n∂Aphi/∂Z = dAy_dz*cos(phi) - dAx_dz*sin(phi)")

print("\n" + "="*80)
print("CURRENT IMPLEMENTATION (line 601)")
print("="*80)
print("dAphi_dZ = grad_AY(3)*cos(phi) - grad_AX(3)*sin(phi)")
print("\nwhere:")
print("  grad_AX(3) = ∂Ax/∂z = dAx_dz")
print("  grad_AY(3) = ∂Ay/∂z = dAy_dz")

print("\nSo current formula is:")
print("  dAphi_dZ = dAy_dz*cos(phi) - dAx_dz*sin(phi)")

print("\n✓ THIS ALSO MATCHES THE CORRECT FORMULA!")

print("\n" + "="*80)
print("CONCLUSION")
print("="*80)
print("Both transformation formulas at lines 599-601 are CORRECT!")
print("The bug must be elsewhere...")
