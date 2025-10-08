#!/usr/bin/env python3
"""Symbolically verify the curl(A) formulas used for Fourier harmonics."""

import sympy as sp

def main() -> int:
    R, Z, phi, n = sp.symbols('R Z phi n', real=True)
    AR_tilde = sp.Function('AR_tilde')(R, Z)
    Aphi_tilde = sp.Function('Aphi_tilde')(R, Z)
    AZ_tilde = sp.Function('AZ_tilde')(R, Z)

    exp_factor = sp.exp(sp.I * n * phi)

    A_R = AR_tilde * exp_factor
    A_phi = Aphi_tilde * exp_factor
    A_Z = AZ_tilde * exp_factor

    BR = sp.simplify(sp.diff(A_Z, phi) / R - sp.diff(A_phi, Z))
    Bphi = sp.simplify(sp.diff(A_R, Z) - sp.diff(A_Z, R))
    BZ = sp.simplify(sp.diff(R * A_phi, R) / R - sp.diff(A_R, phi) / R)

    expected_BR = (sp.I * n / R) * AZ_tilde * exp_factor - sp.diff(Aphi_tilde, Z) * exp_factor
    expected_Bphi = (sp.diff(AR_tilde, Z) - sp.diff(AZ_tilde, R)) * exp_factor
    expected_BZ = (
        sp.diff(Aphi_tilde, R) * exp_factor + Aphi_tilde * exp_factor / R
        - sp.I * n * AR_tilde * exp_factor / R
    )

    diff_BR = sp.simplify(BR - expected_BR)
    diff_Bphi = sp.simplify(Bphi - expected_Bphi)
    diff_BZ = sp.simplify(BZ - expected_BZ)

    if any(expr != 0 for expr in (diff_BR, diff_Bphi, diff_BZ)):
        raise SystemExit(
            f"Curl identities failed\n"
            f"  B_R difference: {diff_BR}\n"
            f"  B_phi difference: {diff_Bphi}\n"
            f"  B_Z difference: {diff_BZ}\n"
        )

    axis_BR = sp.simplify(sp.diff(Aphi_tilde, Z).subs(n, 0))
    axis_Bphi = sp.simplify((sp.diff(AR_tilde, Z) - sp.diff(AZ_tilde, R)).subs(n, 0))
    axis_BZ = sp.simplify((sp.diff(Aphi_tilde, R) + Aphi_tilde / R).subs(n, 0))

    assert axis_BR == sp.diff(Aphi_tilde, Z)
    assert axis_Bphi == sp.diff(AR_tilde, Z) - sp.diff(AZ_tilde, R)
    assert axis_BZ == sp.diff(Aphi_tilde, R) + Aphi_tilde / R

    return 0


if __name__ == '__main__':
    raise SystemExit(main())
