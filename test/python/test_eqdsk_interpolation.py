"""
Tests for EQDSK psi interpolation and flux coordinate conversion.
"""

import pytest
import numpy as np

from libneo.eqdsk import eqdsk_file


TEST_EQDSK = "test/resources/input_efit_file.dat"


class TestPsiInterpolation:

    @pytest.fixture
    def eq(self):
        return eqdsk_file(TEST_EQDSK)

    def test_psi_at_magnetic_axis_equals_psi_axis(self, eq):
        R_axis = eq.Rpsi0
        Z_axis = eq.Zpsi0
        psi = eq.psi_at_rz(R_axis, Z_axis)
        i_r = np.argmin(np.abs(eq.R - R_axis))
        i_z = np.argmin(np.abs(eq.Z - Z_axis))
        psi_grid_at_axis = eq.PsiVs[i_z, i_r]
        assert np.isclose(psi, psi_grid_at_axis, rtol=1e-2)

    def test_psi_at_grid_point_matches_grid_value(self, eq):
        i_r = eq.nrgr // 2
        i_z = eq.nzgr // 2
        R_grid = eq.R[i_r]
        Z_grid = eq.Z[i_z]
        psi_interp = eq.psi_at_rz(R_grid, Z_grid)
        psi_grid = eq.PsiVs[i_z, i_r]
        assert np.isclose(psi_interp, psi_grid, rtol=1e-4)

    def test_psi_at_rz_vectorized(self, eq):
        R_arr = np.array([eq.R[10], eq.R[20], eq.R[30]])
        Z_arr = np.array([eq.Z[50], eq.Z[60], eq.Z[70]])
        psi = eq.psi_at_rz(R_arr, Z_arr)
        assert psi.shape == (3,)
        for i in range(3):
            psi_scalar = eq.psi_at_rz(R_arr[i], Z_arr[i])
            assert np.isclose(psi[i], psi_scalar)

    def test_psi_at_rz_2d_grid(self, eq):
        R_arr = eq.R[10:13]
        Z_arr = eq.Z[50:53]
        psi = eq.psi_at_rz(R_arr, Z_arr, grid=True)
        assert psi.shape == (3, 3)


class TestNormalizedFluxCoordinate:

    @pytest.fixture
    def eq(self):
        return eqdsk_file(TEST_EQDSK)

    def test_spol_at_axis_is_zero(self, eq):
        R_axis = eq.Rpsi0
        Z_axis = eq.Zpsi0
        s_pol = eq.spol_at_rz(R_axis, Z_axis)
        assert np.isclose(s_pol, 0.0, atol=1e-2)

    def test_spol_at_lcfs_is_one(self, eq):
        R_lcfs = eq.Lcfs[0, 0]
        Z_lcfs = eq.Lcfs[1, 0]
        s_pol = eq.spol_at_rz(R_lcfs, Z_lcfs)
        assert np.isclose(s_pol, 1.0, atol=0.05)

    def test_spol_at_rz_vectorized(self, eq):
        R_arr = np.array([eq.R[30], eq.R[35], eq.R[40]])
        Z_arr = np.array([eq.Z[60], eq.Z[65], eq.Z[64]])
        s_pol = eq.spol_at_rz(R_arr, Z_arr)
        assert s_pol.shape == (3,)
        assert np.all(s_pol >= 0.0)

    def test_spol_monotonic_from_axis_outward(self, eq):
        R_axis = eq.Rpsi0
        Z_axis = eq.Zpsi0
        R_line = np.linspace(R_axis, eq.R[-5], 20)
        Z_line = np.full_like(R_line, Z_axis)
        s_pol = eq.spol_at_rz(R_line, Z_line)
        assert np.all(np.diff(s_pol) >= -1e-6)


class TestGeometricAngle:

    @pytest.fixture
    def eq(self):
        return eqdsk_file(TEST_EQDSK)

    def test_theta_at_outboard_midplane_is_zero(self, eq):
        R_out = eq.Rpsi0 + 0.2
        Z_mid = eq.Zpsi0
        theta = eq.theta_geometric_at_rz(R_out, Z_mid)
        assert np.isclose(theta, 0.0, atol=0.05)

    def test_theta_at_top_is_pi_half(self, eq):
        R_axis = eq.Rpsi0
        Z_top = eq.Zpsi0 + 0.5
        theta = eq.theta_geometric_at_rz(R_axis, Z_top)
        assert np.isclose(theta, np.pi / 2, atol=0.05)

    def test_theta_at_bottom_is_minus_pi_half(self, eq):
        R_axis = eq.Rpsi0
        Z_bot = eq.Zpsi0 - 0.5
        theta = eq.theta_geometric_at_rz(R_axis, Z_bot)
        assert np.isclose(theta, -np.pi / 2, atol=0.05)

    def test_theta_vectorized(self, eq):
        R_arr = np.array([eq.Rpsi0 + 0.2, eq.Rpsi0, eq.Rpsi0 - 0.1])
        Z_arr = np.array([eq.Zpsi0, eq.Zpsi0 + 0.3, eq.Zpsi0])
        theta = eq.theta_geometric_at_rz(R_arr, Z_arr)
        assert theta.shape == (3,)


class TestFluxCoordinateConversion:

    @pytest.fixture
    def eq(self):
        return eqdsk_file(TEST_EQDSK)

    def test_rz_to_flux_returns_spol_and_theta(self, eq):
        R = eq.Rpsi0 + 0.2
        Z = eq.Zpsi0
        s_pol, theta = eq.rz_to_flux_coords(R, Z)
        assert s_pol > 0.0
        assert np.isclose(theta, 0.0, atol=0.05)

    def test_rz_to_flux_vectorized(self, eq):
        R_arr = np.array([eq.Rpsi0 + 0.1, eq.Rpsi0 + 0.2])
        Z_arr = np.array([eq.Zpsi0, eq.Zpsi0 + 0.1])
        s_pol, theta = eq.rz_to_flux_coords(R_arr, Z_arr)
        assert s_pol.shape == (2,)
        assert theta.shape == (2,)
