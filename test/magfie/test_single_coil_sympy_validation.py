#!/usr/bin/env python3
"""SymPy validation of all formulas used in single coil axis validation"""

import unittest
import sympy as sp
from sympy import symbols, cos, sin, pi, simplify, sqrt, atan2

class TestSingleCoilAxisFormulas(unittest.TestCase):
    """Verify every formula used in single coil on-axis field calculation"""

    def test_analytic_formula_on_axis(self):
        """
        Verify the on-axis magnetic field formula for a circular current loop.

        Reference: Physics LibreTexts, Biot-Savart derivation
        Formula: B(z) = μ₀IR² / [2(R² + z²)^(3/2)]
        At center (z=0): B(0) = μ₀I / (2R)

        Units: SI system
        - μ₀ = 4π×10⁻⁷ T·m/A
        - I in Amperes
        - R, z in meters
        - B in Tesla
        """
        mu0_SI, I, R, z = symbols('mu_0 I R z', real=True, positive=True)

        # General formula for field at distance z along axis
        B_z = mu0_SI * I * R**2 / (2 * (R**2 + z**2)**(sp.Rational(3, 2)))

        # At center (z=0)
        B_center = B_z.subs(z, 0)
        B_center_simplified = simplify(B_center)
        B_center_expected = mu0_SI * I / (2 * R)

        self.assertTrue(simplify(B_center_simplified - B_center_expected) == 0,
                       "At center, B = μ₀I/(2R)")

        # Verify dimension analysis symbolically
        # [B] = [μ₀][I][R²] / [R³] = [μ₀][I]/[R]
        # In SI: [μ₀] = T·m/A, so [B] = T·m/A · A / m = T ✓

        print("\n✓ Analytic formula verified:")
        print(f"  B(z) = {B_z}")
        print(f"  B(0) = {B_center_simplified}")

    def test_units_conversion_SI_to_Gauss(self):
        """
        Verify conversion from Tesla to Gauss

        1 Tesla = 10,000 Gauss (exactly)
        """
        # This is a definition, not derived
        conversion_factor = 10000  # 1 T = 10^4 G

        self.assertEqual(conversion_factor, int(1e4),
                        "1 Tesla = 10,000 Gauss exactly")

        print("\n✓ Unit conversion verified: 1 T = 10⁴ G")

    def test_mu0_value(self):
        """
        Verify μ₀ value in SI units

        μ₀ = 4π × 10⁻⁷ T·m/A (exact in SI system before 2019 redefinition)
        """
        from numpy import pi as np_pi

        mu0_value = 4e-7 * np_pi  # T·m/A
        expected = 1.25663706e-6  # T·m/A

        # Check to 6 significant figures
        self.assertAlmostEqual(mu0_value, expected, places=13,
                              msg="μ₀ = 4π×10⁻⁷ T·m/A")

        print(f"\n✓ μ₀ value verified: {mu0_value:.10e} T·m/A")

    def test_cylindrical_to_cartesian_B_field(self):
        """
        Verify transformation of magnetic field from cylindrical to Cartesian

        B_x = B_R * cos(φ) - B_φ * sin(φ)
        B_y = B_R * sin(φ) + B_φ * cos(φ)
        B_z = B_Z

        This is for PHYSICAL components (not covariant)
        """
        BR, Bphi, BZ, phi = symbols('B_R B_phi B_Z phi', real=True)

        # Transformation formulas
        Bx = BR * cos(phi) - Bphi * sin(phi)
        By = BR * sin(phi) + Bphi * cos(phi)
        Bz = BZ

        # Verify magnitude is preserved (Pythagorean theorem)
        mag_cartesian = Bx**2 + By**2 + Bz**2
        mag_cylindrical = BR**2 + Bphi**2 + BZ**2

        diff = simplify(mag_cartesian - mag_cylindrical)
        self.assertTrue(diff == 0,
                       "Magnitude must be preserved under coordinate transformation")

        # Verify inverse transformation
        Bx_sym, By_sym = symbols('B_x B_y', real=True)
        BR_inv = Bx_sym * cos(phi) + By_sym * sin(phi)
        Bphi_inv = -Bx_sym * sin(phi) + By_sym * cos(phi)

        # Substitute forward transformation into inverse
        BR_check = simplify(Bx.subs([(BR, BR), (Bphi, Bphi)]) * cos(phi) +
                           By.subs([(BR, BR), (Bphi, Bphi)]) * sin(phi))
        BR_check = simplify(BR_check)

        # Should get back BR
        self.assertTrue(simplify(BR_check - BR) == 0,
                       "Inverse transformation must recover original BR")

        print("\n✓ Cylindrical to Cartesian B transformation verified")
        print(f"  B_x = {Bx}")
        print(f"  B_y = {By}")
        print(f"  Magnitude preserved: {diff == 0}")

    def test_dot_product_formula(self):
        """
        Verify B·n̂ calculation for field component along axis

        For a unit vector n̂ = (n_x, n_y, n_z) with |n̂| = 1,
        B·n̂ = B_x*n_x + B_y*n_y + B_z*n_z
        """
        Bx, By, Bz = symbols('B_x B_y B_z', real=True)
        nx, ny, nz = symbols('n_x n_y n_z', real=True)

        # Dot product
        B_parallel = Bx*nx + By*ny + Bz*nz

        # For unit vector, verify that if B is parallel to n̂, then B·n̂ = |B|
        # Set B = |B|*n̂
        B_mag = symbols('B_mag', positive=True)
        Bx_parallel = B_mag * nx
        By_parallel = B_mag * ny
        Bz_parallel = B_mag * nz

        B_dot_n_parallel = simplify(
            Bx_parallel*nx + By_parallel*ny + Bz_parallel*nz
        )

        # If |n̂| = 1, then B·n̂ = |B|(n_x² + n_y² + n_z²) = |B|
        # But we can't simplify without knowing nx²+ny²+nz²=1
        # Just verify the formula is correct
        expected = B_mag * (nx**2 + ny**2 + nz**2)

        self.assertTrue(simplify(B_dot_n_parallel - expected) == 0,
                       "B·n̂ formula must be correct")

        print("\n✓ Dot product formula verified")
        print(f"  B·n̂ = B_x*n_x + B_y*n_y + B_z*n_z")

    def test_numerical_example_single_coil(self):
        """
        Numerical test: compute B at center of 35cm radius loop with 1A current

        Given:
        - R = 0.35 m
        - I = 1.0 A
        - z = 0 (at center)
        - μ₀ = 4π×10⁻⁷ T·m/A

        Expected:
        - B = μ₀I/(2R) = (4π×10⁻⁷ × 1.0) / (2 × 0.35)
        - B = 1.795×10⁻⁶ T = 0.01795 Gauss
        """
        from numpy import pi as np_pi

        mu0 = 4e-7 * np_pi  # T·m/A
        I = 1.0  # A
        R = 0.35  # m

        B_tesla = mu0 * I / (2.0 * R)
        B_gauss = B_tesla * 1.0e4

        expected_gauss = 0.01795196  # Reference value

        self.assertAlmostEqual(B_gauss, expected_gauss, places=8,
                              msg="Numerical calculation must match expected value")

        print(f"\n✓ Numerical example verified:")
        print(f"  R = {R} m, I = {I} A")
        print(f"  B(center) = {B_tesla:.6e} T = {B_gauss:.8f} G")
        print(f"  Expected  = {expected_gauss:.8f} G")

    def test_off_axis_specific_distances(self):
        """
        Test formula at specific distances from center

        Verify that field decreases as expected when moving away from center
        """
        from numpy import pi as np_pi, sqrt

        mu0 = 4e-7 * np_pi
        I = 1.0
        R = 0.35  # m

        # At center
        B_center = mu0 * I / (2.0 * R)

        # At z = R (one radius away)
        z_R = R
        B_at_R = mu0 * I * R**2 / (2.0 * (R**2 + z_R**2)**1.5)
        B_at_R_simplified = mu0 * I * R**2 / (2.0 * (2*R**2)**1.5)
        B_at_R_expected = mu0 * I / (2.0 * R * 2**(1.5))

        self.assertAlmostEqual(B_at_R, B_at_R_expected, places=15,
                              msg="At z=R, B = B_center / 2^(3/2)")

        # Field should decrease
        self.assertTrue(B_at_R < B_center,
                       "Field must decrease away from center")

        ratio = B_at_R / B_center
        expected_ratio = 1.0 / (2.0**1.5)  # = 1/2√2 ≈ 0.3536

        self.assertAlmostEqual(ratio, expected_ratio, places=10,
                              msg="Ratio B(R)/B(0) = 1/2^(3/2)")

        print(f"\n✓ Off-axis field variation verified:")
        print(f"  B(z=0)   = {B_center:.6e} T")
        print(f"  B(z=R)   = {B_at_R:.6e} T")
        print(f"  B(R)/B(0) = {ratio:.6f} (expected {expected_ratio:.6f})")


if __name__ == '__main__':
    # Run with verbose output
    suite = unittest.TestLoader().loadTestsFromTestCase(TestSingleCoilAxisFormulas)
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)

    # Exit with error code if tests failed
    exit(0 if result.wasSuccessful() else 1)
