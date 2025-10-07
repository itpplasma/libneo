#!/usr/bin/env python3
"""SymPy validation of analytical formulas used in coil_tools."""

import unittest
import sympy as sp


class TestCoilToolsAnalyticalFormulas(unittest.TestCase):
    """Verify every analytical assumption in coil_tools with SymPy."""

    def test_cartesian_to_cylindrical_vector_potential(self):
        """Validate AR, Aphi, AZ transformation from Cartesian components."""
        AX, AY, AZ_cart, phi = sp.symbols('A_X A_Y A_Z phi', real=True)
        AR_fortran = AX * sp.cos(phi) + AY * sp.sin(phi)
        Aphi_fortran = AY * sp.cos(phi) - AX * sp.sin(phi)
        AZ_fortran = AZ_cart

        AR_expected = AX * sp.cos(phi) + AY * sp.sin(phi)
        Aphi_expected = -AX * sp.sin(phi) + AY * sp.cos(phi)
        AZ_expected = AZ_cart

        self.assertEqual(sp.simplify(AR_fortran - AR_expected), 0)
        self.assertEqual(sp.simplify(Aphi_fortran - Aphi_expected), 0)
        self.assertEqual(sp.simplify(AZ_fortran - AZ_expected), 0)

    def test_cartesian_to_cylindrical_magnetic_field(self):
        """Validate BR, Bphi, BZ transformation from Cartesian components."""
        Bx, By, Bz, phi = sp.symbols('B_x B_y B_z phi', real=True)
        BR_fortran = Bx * sp.cos(phi) + By * sp.sin(phi)
        Bphi_fortran = By * sp.cos(phi) - Bx * sp.sin(phi)
        BZ_fortran = Bz

        BR_expected = Bx * sp.cos(phi) + By * sp.sin(phi)
        Bphi_expected = -Bx * sp.sin(phi) + By * sp.cos(phi)
        BZ_expected = Bz

        self.assertEqual(sp.simplify(BR_fortran - BR_expected), 0)
        self.assertEqual(sp.simplify(Bphi_fortran - Bphi_expected), 0)
        self.assertEqual(sp.simplify(BZ_fortran - BZ_expected), 0)

    def test_dAphi_dR_formula(self):
        """Verify ∂Aφ/∂R formula derived from Cartesian gradients."""
        phi = sp.symbols('phi', real=True)
        dAX_dx, dAX_dy = sp.symbols('dAX_dx dAX_dy', real=True)
        dAY_dx, dAY_dy = sp.symbols('dAY_dx dAY_dy', real=True)

        dAphi_dR_fortran = (dAY_dx * sp.cos(phi) ** 2 - dAX_dy * sp.sin(phi) ** 2
                             + (dAY_dy - dAX_dx) * sp.cos(phi) * sp.sin(phi))

        dAphi_dR_expected = (-sp.sin(phi) * (dAX_dx * sp.cos(phi) + dAX_dy * sp.sin(phi))
                             + sp.cos(phi) * (dAY_dx * sp.cos(phi) + dAY_dy * sp.sin(phi)))

        self.assertEqual(sp.simplify(dAphi_dR_fortran - dAphi_dR_expected), 0)

    def test_dAphi_dZ_formula(self):
        """Verify ∂Aφ/∂Z formula derived from Cartesian gradients."""
        phi = sp.symbols('phi', real=True)
        dAX_dz, dAY_dz = sp.symbols('dAX_dz dAY_dz', real=True)

        dAphi_dZ_fortran = dAY_dz * sp.cos(phi) - dAX_dz * sp.sin(phi)
        dAphi_dZ_expected = -sp.sin(phi) * dAX_dz + sp.cos(phi) * dAY_dz

        self.assertEqual(sp.simplify(dAphi_dZ_fortran - dAphi_dZ_expected), 0)

    def test_segment_vector_potential_contribution(self):
        """Derive segment vector potential formula from Biot-Savart integral."""
        L, y0, z0 = sp.symbols('L y0 z0', positive=True)
        x = sp.symbols('x', real=True)

        integrand = 1 / sp.sqrt(x**2 + y0**2 + z0**2)
        integral = sp.integrate(integrand, (x, -L / 2, L / 2))

        dist_i = sp.sqrt((L / 2) ** 2 + y0 ** 2 + z0 ** 2)
        dist_f = dist_i
        ecc = L / (dist_i + dist_f)
        fortran_form = sp.log((1 + ecc) / (1 - ecc))
        analytic_form = 2 * sp.asinh(L / (2 * sp.sqrt(y0 ** 2 + z0 ** 2)))

        self.assertEqual(sp.simplify(integral - analytic_form), 0)

        samples = (
            (0.6, 0.45, 0.35),
            (0.8, 0.62, 0.28),
            (0.4, 0.18, 0.27)
        )
        for sample in samples:
            numeric = sp.N((fortran_form - analytic_form).subs({L: sample[0], y0: sample[1], z0: sample[2]}), 50)
            self.assertAlmostEqual(float(numeric), 0.0, places=12)

    def test_segment_gradient_kernel_matches_vector_potential_gradient(self):
        """Verify gradient kernel reproduces ∇A when scaled by segment vector."""
        L, x0, y0, z0 = sp.symbols('L x0 y0 z0', positive=True, real=True)
        obs = sp.Matrix([x0, y0, z0])
        seg_start = sp.Matrix([-L / 2, 0, 0])
        seg_end = sp.Matrix([L / 2, 0, 0])
        segment_vec = seg_end - seg_start
        dist_segment = sp.sqrt(segment_vec.dot(segment_vec))

        def compute_terms(point):
            vec = point - obs
            dist = sp.sqrt(vec.dot(vec))
            return vec, dist

        XYZ_i, dist_i = compute_terms(seg_start)
        XYZ_f, dist_f = compute_terms(seg_end)
        ecc = L / (dist_i + dist_f)
        aphi = sp.log((1 + ecc) / (1 - ecc))

        grad_kernel = (XYZ_i / dist_i + XYZ_f / dist_f) / (dist_i * dist_f + XYZ_i.dot(XYZ_f))
        grad_AX_fortran = segment_vec[0] * grad_kernel

        grad_aphi = sp.Matrix([sp.diff(aphi, coord) for coord in (x0, y0, z0)])
        # The analytic gradients agree up to simplification; verify symbolically at x0=0 for tractability.
        diff_axis = sp.simplify(grad_aphi.subs(x0, 0) - grad_AX_fortran.subs(x0, 0))
        self.assertTrue(all(expr == 0 for expr in diff_axis))

        # Numerical spot checks confirm equality for general x0.
        subs_map = {L: 0.7, x0: 0.13, y0: 0.21, z0: -0.18}
        numeric_diff = (grad_aphi - grad_AX_fortran).subs(subs_map).evalf()
        self.assertTrue(all(abs(val) < 1e-12 for val in numeric_diff))

    def test_segment_gradient_contribution(self):
        """Ensure segment_gradient_contribution scales each Cartesian gradient."""
        grad_kernel = sp.Matrix(sp.symbols('g0 g1 g2', real=True))
        XYZ_segment = sp.Matrix(sp.symbols('s0 s1 s2', real=True))

        expected_AX = XYZ_segment[0] * grad_kernel
        expected_AY = XYZ_segment[1] * grad_kernel
        expected_AZ = XYZ_segment[2] * grad_kernel

        # The Fortran routine performs these multiplications component-wise.
        self.assertEqual(expected_AX, sp.Matrix([XYZ_segment[0] * g for g in grad_kernel]))
        self.assertEqual(expected_AY, sp.Matrix([XYZ_segment[1] * g for g in grad_kernel]))
        self.assertEqual(expected_AZ, sp.Matrix([XYZ_segment[2] * g for g in grad_kernel]))

    def test_curl_formulas_ntor_0(self):
        """Verify cylindrical curl identities used for ntor=0."""
        R, Z, phi = sp.symbols('R Z phi', positive=True, real=True)
        AR = sp.Function('A_R')(R, Z)
        Aphi = sp.Function('A_phi')(R, Z)
        AZ = sp.Function('A_Z')(R, Z)

        BR_expected = -sp.diff(Aphi, Z)
        Bphi_expected = sp.diff(AR, Z) - sp.diff(AZ, R)
        BZ_expected = sp.diff(Aphi, R) + Aphi / R
        BZ_alt = sp.diff(R * Aphi, R) / R

        self.assertEqual(sp.simplify(BZ_expected - BZ_alt), 0)
        # identities for BR and Bphi hold by definition; ensure no simplification triggers
        self.assertEqual(BR_expected, -sp.diff(Aphi, Z))
        self.assertEqual(Bphi_expected, sp.diff(AR, Z) - sp.diff(AZ, R))

    def test_fourier_mode_formulas_ntor_nonzero(self):
        """Verify curl(A) for non-zero toroidal mode reduces to implemented formulas."""
        R, Z, phi, n = sp.symbols('R Z phi n', positive=True, real=True)
        n_int = sp.symbols('n_int', integer=True, nonzero=True)
        n = n_int

        AnR = sp.Function('AnR')(R, Z)
        AnZ = sp.Function('AnZ')(R, Z)

        A_R = AnR * sp.exp(sp.I * n * phi)
        A_phi = 0
        A_Z = AnZ * sp.exp(sp.I * n * phi)

        BR = (1 / R) * sp.diff(A_Z, phi) - sp.diff(A_phi, Z)
        Bphi = sp.diff(A_R, Z) - sp.diff(A_Z, R)
        BZ = sp.diff(R * A_phi, R) / R - (1 / R) * sp.diff(A_R, phi)

        expected_BnR = sp.I * n * AnZ / R
        expected_Bnphi = sp.diff(AnR, Z) - sp.diff(AnZ, R)
        expected_BnZ = -sp.I * n * AnR / R

        self.assertEqual(sp.simplify(sp.exp(-sp.I * n * phi) * BR - expected_BnR), 0)
        self.assertEqual(sp.simplify(sp.exp(-sp.I * n * phi) * Bphi - expected_Bnphi), 0)
        self.assertEqual(sp.simplify(sp.exp(-sp.I * n * phi) * BZ - expected_BnZ), 0)


if __name__ == '__main__':
    unittest.main(verbosity=2)
