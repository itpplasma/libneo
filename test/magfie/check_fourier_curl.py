#!/usr/bin/env python3
"""Check curl formulas for Fourier modes with physical cylindrical components"""

import sympy as sp
from sympy import symbols, cos, sin, exp, I, diff, simplify, re

# Cylindrical coordinates and mode number
R, phi, Z, n = symbols('R phi Z n', real=True)
n = symbols('n', integer=True, nonzero=True)

# Fourier mode amplitudes (complex)
AR_re, AR_im = symbols('AR_re AR_im', real=True)
Aphi_re, Aphi_im = symbols('Aphi_re Aphi_im', real=True)
AZ_re, AZ_im = symbols('AZ_re AZ_im', real=True)

AR_tilde = AR_re + I * AR_im
Aphi_tilde = Aphi_re + I * Aphi_im
AZ_tilde = AZ_re + I * AZ_im

# Physical vector potential components for mode n
# A(R,φ,Z) = 2*Re[Ã(R,Z) * e^(inφ)]
AR = 2 * re(AR_tilde * exp(I * n * phi))
Aphi = 2 * re(Aphi_tilde * exp(I * n * phi))
AZ = 2 * re(AZ_tilde * exp(I * n * phi))

print("Physical vector potential components:")
print(f"AR = {AR}")
print(f"Aphi = {Aphi}")
print(f"AZ = {AZ}")
print()

# Curl in cylindrical coordinates for physical components:
# B_R = (1/R)∂A_Z/∂φ - ∂A_φ/∂Z
# B_φ = ∂A_R/∂Z - ∂A_Z/∂R
# B_Z = (1/R)∂(R*A_φ)/∂R - (1/R)∂A_R/∂φ

BR = diff(AZ, phi) / R - diff(Aphi, Z)
Bphi = diff(AR, Z) - diff(AZ, R)
BZ = diff(R * Aphi, R) / R - diff(AR, phi) / R

print("Computing curl B = ∇×A...")
print()

# Simplify
BR_simplified = simplify(BR)
Bphi_simplified = simplify(Bphi)
BZ_simplified = simplify(BZ)

print("B_R =", BR_simplified)
print()
print("B_φ =", Bphi_simplified)
print()
print("B_Z =", BZ_simplified)
print()

# Extract coefficients of cos(n*phi) and sin(n*phi) for each component
# B = 2*Re[B̃ * e^(inφ)] = 2*(B̃_re*cos(nφ) - B̃_im*sin(nφ))

print("\n" + "="*70)
print("Checking if for ntor≠0, when Aphi_tilde = 0:")
print("="*70)

# Substitute Aphi_tilde = 0 (which is what we store for ntor≠0)
BR_no_aphi = simplify(BR.subs([(Aphi_re, 0), (Aphi_im, 0)]))
BZ_no_aphi = simplify(BZ.subs([(Aphi_re, 0), (Aphi_im, 0)]))

print(f"\nB_R (with Aphi=0) = {BR_no_aphi}")
print(f"B_Z (with Aphi=0) = {BZ_no_aphi}")

# Check if this matches the Fourier formula:
# B_R = 2*Re[i*n/R * AZ_tilde * e^(inφ)]
# B_Z = 2*Re[-i*n/R * AR_tilde * e^(inφ)]

expected_BR = 2 * re(I * n / R * AZ_tilde * exp(I * n * phi))
expected_BZ = 2 * re(-I * n / R * AR_tilde * exp(I * n * phi))

print(f"\nExpected B_R (Fourier) = {simplify(expected_BR)}")
print(f"Expected B_Z (Fourier) = {simplify(expected_BZ)}")

diff_BR = simplify(BR_no_aphi - expected_BR)
diff_BZ = simplify(BZ_no_aphi - expected_BZ)

print(f"\nDifference B_R: {diff_BR}")
print(f"Difference B_Z: {diff_BZ}")

if diff_BR == 0 and diff_BZ == 0:
    print("\n✓ Fourier formulas are CORRECT when Aphi=0 for ntor≠0")
else:
    print("\n✗ Fourier formulas are WRONG!")
