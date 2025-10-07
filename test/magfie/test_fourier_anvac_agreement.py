#!/usr/bin/env python3
"""
Persistent test: Verify Fourier and Anvac methods agree for ntor=0

This test validates that computing B from curl(A) using the n=0 Fourier mode
of the vector potential gives the same result as the n=0 Fourier mode of B
computed directly by the Fortran code.

Both methods correctly handle the discretization of the coil geometry.
"""

import unittest
import sys
import os
from pathlib import Path

# Add python module to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / 'python'))

import numpy as np
from scipy.interpolate import RectBivariateSpline
import netCDF4 as nc

from libneo.biotsavart_fourier import read_Bnvac_fourier, field_divfree, grid_t


class TestFourierAnvacAgreement(unittest.TestCase):
    """Test that Fourier and Anvac methods agree for single coil ntor=0"""

    @classmethod
    def setUpClass(cls):
        """Load test data files"""
        cls.test_dir = Path(__file__).parent.parent.parent / 'build' / 'test' / 'magfie' / 'single_coil'

        cls.h5_file = cls.test_dir / 'single_reference.h5'
        cls.nc_file = cls.test_dir / 'single_test.nc'

        if not cls.h5_file.exists():
            raise unittest.SkipTest(f"Reference file not found: {cls.h5_file}")
        if not cls.nc_file.exists():
            raise unittest.SkipTest(f"Test file not found: {cls.nc_file}")

        # Load Fourier (Bnvac) reference
        cls.grid_ref, cls.BnR_ref, cls.Bnphi_ref, cls.BnZ_ref = read_Bnvac_fourier(
            str(cls.h5_file), ntor=0
        )

        # Load Anvac and compute B
        ds = nc.Dataset(str(cls.nc_file), 'r')
        cls.grid_anv = grid_t()
        cls.grid_anv.R = np.array(ds['R'][:])
        cls.grid_anv.Z = np.array(ds['Z'][:])

        ntor = 0
        AnR_real = ds['AnR_real'][0, :, :, ntor].T
        Anphi_real = ds['Anphi_real'][0, :, :, ntor].T
        AnZ_real = ds['AnZ_real'][0, :, :, ntor].T

        spl = {
            'AnR_Re': [RectBivariateSpline(cls.grid_anv.R, cls.grid_anv.Z, AnR_real, kx=5, ky=5)],
            'AnR_Im': [RectBivariateSpline(cls.grid_anv.R, cls.grid_anv.Z, np.zeros_like(AnR_real), kx=5, ky=5)],
            'AnZ_Re': [RectBivariateSpline(cls.grid_anv.R, cls.grid_anv.Z, AnZ_real, kx=5, ky=5)],
            'AnZ_Im': [RectBivariateSpline(cls.grid_anv.R, cls.grid_anv.Z, np.zeros_like(AnZ_real), kx=5, ky=5)],
            'Anphi_Re': [RectBivariateSpline(cls.grid_anv.R, cls.grid_anv.Z, Anphi_real, kx=5, ky=5)],
            'Anphi_Im': [RectBivariateSpline(cls.grid_anv.R, cls.grid_anv.Z, np.zeros_like(Anphi_real), kx=5, ky=5)],
        }

        cls.BnR_anv, cls.Bnphi_anv, cls.BnZ_anv = field_divfree(spl, cls.grid_anv.R, cls.grid_anv.Z, ntor=0)
        ds.close()

        # Find grid point nearest to coil center
        R_center = np.hypot(197.2682697, 72.00853957)
        Z_center = 78.0
        cls.iR_ref = np.argmin(np.abs(cls.grid_ref.R - R_center))
        cls.iZ_ref = np.argmin(np.abs(cls.grid_ref.Z - Z_center))
        cls.iR_anv = np.argmin(np.abs(cls.grid_anv.R - R_center))
        cls.iZ_anv = np.argmin(np.abs(cls.grid_anv.Z - Z_center))

    def test_BnR_agreement(self):
        """Test that BnR agrees between Fourier and Anvac methods"""
        BnR_fourier = self.BnR_ref[0, self.iR_ref, self.iZ_ref].real
        BnR_anvac = self.BnR_anv[self.iR_anv, self.iZ_anv].real

        rel_diff = abs((BnR_anvac - BnR_fourier) / BnR_fourier) * 100

        self.assertLess(rel_diff, 1.0,
                       f"BnR differs by {rel_diff:.3f}% (>{1.0}%)")

        print(f"\n  BnR: Fourier={BnR_fourier:.6e}, Anvac={BnR_anvac:.6e}, diff={rel_diff:.3f}%")

    def test_Bnphi_agreement(self):
        """Test that Bnphi agrees between Fourier and Anvac methods"""
        Bnphi_fourier = self.Bnphi_ref[0, self.iR_ref, self.iZ_ref].real
        Bnphi_anvac = self.Bnphi_anv[self.iR_anv, self.iZ_anv].real

        rel_diff = abs((Bnphi_anvac - Bnphi_fourier) / Bnphi_fourier) * 100

        self.assertLess(rel_diff, 1.0,
                       f"Bnphi differs by {rel_diff:.3f}% (>1.0%)")

        print(f"  Bnphi: Fourier={Bnphi_fourier:.6e}, Anvac={Bnphi_anvac:.6e}, diff={rel_diff:.3f}%")

    def test_BnZ_agreement(self):
        """Test that BnZ agrees between Fourier and Anvac methods"""
        BnZ_fourier = self.BnZ_ref[0, self.iR_ref, self.iZ_ref].real
        BnZ_anvac = self.BnZ_anv[self.iR_anv, self.iZ_anv].real

        rel_diff = abs((BnZ_anvac - BnZ_fourier) / BnZ_fourier) * 100

        self.assertLess(rel_diff, 1.0,
                       f"BnZ differs by {rel_diff:.3f}% (>1.0%)")

        print(f"  BnZ: Fourier={BnZ_fourier:.6e}, Anvac={BnZ_anvac:.6e}, diff={rel_diff:.3f}%")

    def test_Bnphi_nonzero(self):
        """Test that Bnphi is non-zero (expected for discretized coil)"""
        Bnphi_fourier = self.Bnphi_ref[0, self.iR_ref, self.iZ_ref].real
        Bnphi_anvac = self.Bnphi_anv[self.iR_anv, self.iZ_anv].real

        # Both should be significantly non-zero (>1e-4 G)
        self.assertGreater(abs(Bnphi_fourier), 1e-4,
                          "Bnphi from Fourier should be non-zero for discretized coil")
        self.assertGreater(abs(Bnphi_anvac), 1e-4,
                          "Bnphi from Anvac should be non-zero for discretized coil")

        print(f"  ✓ Bnphi is non-zero as expected for 129-segment discretized coil")
        print(f"    (Bφ={Bnphi_fourier:.6e} G from discretization effects)")

    def test_magnitude_agreement(self):
        """Test that |B| agrees between methods"""
        BnR_f = self.BnR_ref[0, self.iR_ref, self.iZ_ref].real
        Bnphi_f = self.Bnphi_ref[0, self.iR_ref, self.iZ_ref].real
        BnZ_f = self.BnZ_ref[0, self.iR_ref, self.iZ_ref].real
        mag_fourier = np.sqrt(BnR_f**2 + Bnphi_f**2 + BnZ_f**2)

        BnR_a = self.BnR_anv[self.iR_anv, self.iZ_anv].real
        Bnphi_a = self.Bnphi_anv[self.iR_anv, self.iZ_anv].real
        BnZ_a = self.BnZ_anv[self.iR_anv, self.iZ_anv].real
        mag_anvac = np.sqrt(BnR_a**2 + Bnphi_a**2 + BnZ_a**2)

        rel_diff = abs((mag_anvac - mag_fourier) / mag_fourier) * 100

        self.assertLess(rel_diff, 1.0,
                       f"|B| differs by {rel_diff:.3f}% (>1.0%)")

        print(f"  |B|: Fourier={mag_fourier:.6e}, Anvac={mag_anvac:.6e}, diff={rel_diff:.3f}%")


if __name__ == '__main__':
    # Run with verbose output
    suite = unittest.TestLoader().loadTestsFromTestCase(TestFourierAnvacAgreement)
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)

    if result.wasSuccessful():
        print("\n" + "="*70)
        print("✅ ALL TESTS PASSED")
        print("="*70)
        print("Fourier and Anvac methods agree to <1% for ntor=0")
        print("Non-zero Bφ is correct for discretized 129-segment coil")
        print("="*70)

    exit(0 if result.wasSuccessful() else 1)
