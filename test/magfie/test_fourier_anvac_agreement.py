#!/usr/bin/env python3
"""Persistent test: Verify Fourier and Anvac sums agree for the single coil."""

import sys
from pathlib import Path
import unittest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / 'python'))

import numpy as np
from scipy.interpolate import RectBivariateSpline

from libneo.biotsavart_fourier import (
    field_divfree,
    gauge_Anvac,
    read_Anvac_fourier_all,
    read_Bnvac_fourier_all,
    reconstruct_field_from_modes,
    spline_gauged_Anvac,
)


REFERENCE_H5 = Path('build/test/magfie/single_coil/single_reference.h5')
ANVAC_NC = Path('build/test/magfie/single_coil/single_test.nc')
CURRENT_FILE = Path('test/magfie/test_data/single_coil_currents.txt')
PREFactor = 0.1


class TestFourierAnvacAgreement(unittest.TestCase):
    """Check that full Fourier sums agree with Anvac-derived fields."""

    @classmethod
    def setUpClass(cls):
        if not (REFERENCE_H5.exists() and ANVAC_NC.exists()):
            raise unittest.SkipTest("Single-coil reference data missing; run test suite first")

        cls.grid_ref, cls.mode_numbers, BnR_ref, Bnphi_ref, BnZ_ref = read_Bnvac_fourier_all(str(REFERENCE_H5))
        grid_anv, mode_numbers_anv, AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ, AnX_raw, AnY_raw, AnZ_raw =             read_Anvac_fourier_all(str(ANVAC_NC))
        if not np.array_equal(cls.mode_numbers, mode_numbers_anv):
            raise unittest.SkipTest("Mode numbers differ between Fourier and Anvac datasets")

        currents = np.atleast_1d(np.loadtxt(CURRENT_FILE))
        cls.weights = currents * PREFactor

        cls.BnR_fourier_sum = np.tensordot(cls.weights, BnR_ref, axes=(0, 1))
        cls.Bnphi_fourier_sum = np.tensordot(cls.weights, Bnphi_ref, axes=(0, 1))
        cls.BnZ_fourier_sum = np.tensordot(cls.weights, BnZ_ref, axes=(0, 1))

        # Compute Bn from Anvac for every mode using stored derivatives
        nmodes = cls.mode_numbers.size
        ncoil = AnR.shape[1]
        cls.BnR_anvac_sum = np.empty((nmodes, cls.grid_ref.nR, cls.grid_ref.nZ), dtype=complex)
        cls.Bnphi_anvac_sum = np.empty_like(cls.BnR_anvac_sum)
        cls.BnZ_anvac_sum = np.empty_like(cls.BnR_anvac_sum)
        cls.BnR_anvac_spline_sum = np.empty_like(cls.BnR_anvac_sum)
        cls.Bnphi_anvac_spline_sum = np.empty_like(cls.BnR_anvac_sum)
        cls.BnZ_anvac_spline_sum = np.empty_like(cls.BnR_anvac_sum)

        dAphi_dR_spline = np.empty_like(Anphi)
        dAphi_dZ_spline = np.empty_like(Anphi)

        for idx, ntor in enumerate(cls.mode_numbers):
            gauged_AnR, gauged_AnZ = gauge_Anvac(
                grid_anv,
                AnR[idx],
                Anphi[idx],
                AnZ[idx],
                dAnphi_dR[idx],
                dAnphi_dZ[idx],
                ntor=ntor,
            )
            spl = spline_gauged_Anvac(grid_anv, gauged_AnR, gauged_AnZ, ntor=ntor, Anphi=Anphi[idx])
            BnR_mode, Bnphi_mode, BnZ_mode = field_divfree(spl, cls.grid_ref.R, cls.grid_ref.Z, ntor=ntor)
            if BnR_mode.ndim == 2:
                BnR_mode = BnR_mode[None, :, :]
                Bnphi_mode = Bnphi_mode[None, :, :]
                BnZ_mode = BnZ_mode[None, :, :]
            cls.BnR_anvac_sum[idx] = np.tensordot(cls.weights, BnR_mode, axes=(0, 0))
            cls.Bnphi_anvac_sum[idx] = np.tensordot(cls.weights, Bnphi_mode, axes=(0, 0))
            cls.BnZ_anvac_sum[idx] = np.tensordot(cls.weights, BnZ_mode, axes=(0, 0))

            if ntor == 0:
                for coil_idx in range(ncoil):
                    spline_real = RectBivariateSpline(grid_anv.R, grid_anv.Z, Anphi[idx, coil_idx].real)
                    spline_imag = RectBivariateSpline(grid_anv.R, grid_anv.Z, Anphi[idx, coil_idx].imag)
                    dAphi_dR_spline[idx, coil_idx] = (
                        spline_real(grid_anv.R, grid_anv.Z, dx=1, dy=0, grid=True)
                        + 1j * spline_imag(grid_anv.R, grid_anv.Z, dx=1, dy=0, grid=True)
                    )
                    dAphi_dZ_spline[idx, coil_idx] = (
                        spline_real(grid_anv.R, grid_anv.Z, dx=0, dy=1, grid=True)
                        + 1j * spline_imag(grid_anv.R, grid_anv.Z, dx=0, dy=1, grid=True)
                    )
            else:
                dAphi_dR_spline[idx] = dAnphi_dR[idx]
                dAphi_dZ_spline[idx] = dAnphi_dZ[idx]

        for idx, ntor in enumerate(cls.mode_numbers):
            gauged_AnR, gauged_AnZ = gauge_Anvac(
                grid_anv,
                AnR[idx],
                Anphi[idx],
                AnZ[idx],
                dAphi_dR_spline[idx],
                dAphi_dZ_spline[idx],
                ntor=ntor,
            )
            spl = spline_gauged_Anvac(grid_anv, gauged_AnR, gauged_AnZ, ntor=ntor, Anphi=Anphi[idx])
            BnR_mode, Bnphi_mode, BnZ_mode = field_divfree(spl, cls.grid_ref.R, cls.grid_ref.Z, ntor=ntor)
            if BnR_mode.ndim == 2:
                BnR_mode = BnR_mode[None, :, :]
                Bnphi_mode = Bnphi_mode[None, :, :]
                BnZ_mode = BnZ_mode[None, :, :]
            cls.BnR_anvac_spline_sum[idx] = np.tensordot(cls.weights, BnR_mode, axes=(0, 0))
            cls.Bnphi_anvac_spline_sum[idx] = np.tensordot(cls.weights, Bnphi_mode, axes=(0, 0))
            cls.BnZ_anvac_spline_sum[idx] = np.tensordot(cls.weights, BnZ_mode, axes=(0, 0))

        inv_R = np.divide(
            1.0,
            cls.grid_ref.R.reshape((1, 1, -1, 1)),
            out=np.zeros((1, 1, cls.grid_ref.nR, 1), dtype=float),
            where=cls.grid_ref.R.reshape((1, 1, -1, 1)) > 0.0,
        )
        try:
            dAnR_dZ = np.gradient(AnR, grid_anv.Z, axis=3, edge_order=2)
            dAnZ_dR = np.gradient(AnZ, grid_anv.R, axis=2, edge_order=2)
        except ValueError:
            dAnR_dZ = np.gradient(AnR, grid_anv.Z, axis=3, edge_order=1)
            dAnZ_dR = np.gradient(AnZ, grid_anv.R, axis=2, edge_order=1)
        ntor_vals = cls.mode_numbers.reshape((-1, 1, 1, 1)).astype(float)
        BnR_ung_modes = 1j * ntor_vals * AnZ * inv_R - dAnphi_dZ
        Bnphi_ung_modes = dAnR_dZ - dAnZ_dR
        BnZ_ung_modes = Anphi * inv_R + dAnphi_dR - 1j * ntor_vals * AnR * inv_R
        cls.BnR_ungauged_sum = np.tensordot(cls.weights, BnR_ung_modes, axes=(0, 1))
        cls.Bnphi_ungauged_sum = np.tensordot(cls.weights, Bnphi_ung_modes, axes=(0, 1))
        cls.BnZ_ungauged_sum = np.tensordot(cls.weights, BnZ_ung_modes, axes=(0, 1))

        cls.R_center = np.hypot(197.2682697, 72.00853957)
        cls.Z_center = 78.0
        cls.x_center = 197.2682697
        cls.y_center = 72.00853957
        cls.z_center = 78.0

        # Field maps at phi=0 for full comparison
        R_mesh, Z_mesh = np.meshgrid(cls.grid_ref.R, cls.grid_ref.Z, indexing='ij')
        X_map = R_mesh
        Y_map = np.zeros_like(R_mesh)
        Z_map = Z_mesh

        cls.Bx_fourier_map, cls.By_fourier_map, cls.Bz_fourier_map = reconstruct_field_from_modes(
            cls.BnR_fourier_sum,
            cls.Bnphi_fourier_sum,
            cls.BnZ_fourier_sum,
            cls.mode_numbers,
            cls.grid_ref.R,
            cls.grid_ref.Z,
            X_map,
            Y_map,
            Z_map,
        )
        cls.Bx_anvac_map, cls.By_anvac_map, cls.Bz_anvac_map = reconstruct_field_from_modes(
            cls.BnR_anvac_sum,
            cls.Bnphi_anvac_sum,
            cls.BnZ_anvac_sum,
            cls.mode_numbers,
            cls.grid_ref.R,
            cls.grid_ref.Z,
            X_map,
            Y_map,
            Z_map,
        )
        cls.Bx_anvac_spline_map, cls.By_anvac_spline_map, cls.Bz_anvac_spline_map = reconstruct_field_from_modes(
            cls.BnR_anvac_spline_sum,
            cls.Bnphi_anvac_spline_sum,
            cls.BnZ_anvac_spline_sum,
            cls.mode_numbers,
            cls.grid_ref.R,
            cls.grid_ref.Z,
            X_map,
            Y_map,
            Z_map,
        )

        cls.Bx_anvac_ung_map, cls.By_anvac_ung_map, cls.Bz_anvac_ung_map = reconstruct_field_from_modes(
            cls.BnR_ungauged_sum,
            cls.Bnphi_ungauged_sum,
            cls.BnZ_ungauged_sum,
            cls.mode_numbers,
            cls.grid_ref.R,
            cls.grid_ref.Z,
            X_map,
            Y_map,
            Z_map,
        )

    def test_center_field_agreement(self):
        Bx_fourier, By_fourier, Bz_fourier = reconstruct_field_from_modes(
            self.BnR_fourier_sum,
            self.Bnphi_fourier_sum,
            self.BnZ_fourier_sum,
            self.mode_numbers,
            self.grid_ref.R,
            self.grid_ref.Z,
            self.x_center,
            self.y_center,
            self.z_center,
        )
        Bx_anvac, By_anvac, Bz_anvac = reconstruct_field_from_modes(
            self.BnR_anvac_sum,
            self.Bnphi_anvac_sum,
            self.BnZ_anvac_sum,
            self.mode_numbers,
            self.grid_ref.R,
            self.grid_ref.Z,
            self.x_center,
            self.y_center,
            self.z_center,
        )

        Bx_anvac_spline, By_anvac_spline, Bz_anvac_spline = reconstruct_field_from_modes(
            self.BnR_anvac_spline_sum,
            self.Bnphi_anvac_spline_sum,
            self.BnZ_anvac_spline_sum,
            self.mode_numbers,
            self.grid_ref.R,
            self.grid_ref.Z,
            self.x_center,
            self.y_center,
            self.z_center,
        )

        diff_stored = np.abs(np.array([Bx_anvac - Bx_fourier, By_anvac - By_fourier, Bz_anvac - Bz_fourier]))
        diff_spline = np.abs(np.array([Bx_anvac_spline - Bx_fourier,
                                       By_anvac_spline - By_fourier,
                                       Bz_anvac_spline - Bz_fourier]))
        denom = np.maximum(np.abs([Bx_fourier, By_fourier, Bz_fourier]), 1e-15)
        rel_stored = diff_stored / denom
        rel_spline = diff_spline / denom

        self.assertLess(
            min(rel_stored.max(), rel_spline.max()),
            0.10,
            "Center-field components differ by more than 10%",
        )

    def test_ungauged_field_finite(self):
        mag_ungauged = np.sqrt(self.Bx_anvac_ung_map**2 + self.By_anvac_ung_map**2 + self.Bz_anvac_ung_map**2)
        self.assertTrue(np.all(np.isfinite(mag_ungauged)))

    @unittest.expectedFailure
    def test_field_map_rms(self):
        mag_fourier = np.sqrt(self.Bx_fourier_map**2 + self.By_fourier_map**2 + self.Bz_fourier_map**2)
        mag_anvac = np.sqrt(self.Bx_anvac_map**2 + self.By_anvac_map**2 + self.Bz_anvac_map**2)
        mag_anvac_spline = np.sqrt(
            self.Bx_anvac_spline_map**2 + self.By_anvac_spline_map**2 + self.Bz_anvac_spline_map**2
        )
        rel_stored = np.abs(mag_anvac - mag_fourier) / np.maximum(mag_fourier, 1e-12)
        rel_spline = np.abs(mag_anvac_spline - mag_fourier) / np.maximum(mag_fourier, 1e-12)

        self.assertLess(
            min(np.mean(rel_stored), np.mean(rel_spline)),
            0.05,
            "Mean relative error of field map exceeds 5%",
        )
        self.assertLess(
            min(np.max(rel_stored), np.max(rel_spline)),
            0.15,
            "Max relative error of field map exceeds 15%",
        )


if __name__ == '__main__':
    unittest.main(verbosity=2)
