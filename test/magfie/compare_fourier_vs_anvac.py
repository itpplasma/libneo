#!/usr/bin/env python3
"""Compare coil_tools Fourier vs Anvac totals at coil center using all modes."""

import sys
sys.path.insert(0, 'python')

import numpy as np

from libneo.biotsavart_fourier import (
    field_divfree,
    gauge_Anvac,
    read_Anvac_fourier_all,
    read_Bnvac_fourier_all,
    reconstruct_field_from_modes,
    spline_gauged_Anvac,
)


REFERENCE_H5 = 'build/test/magfie/single_coil/single_reference.h5'
ANVAC_NC = 'build/test/magfie/single_coil/single_test.nc'
CURRENT_FILE = 'test/magfie/test_data/single_coil_currents.txt'
PREFactor = 0.1


def _compute_anvac_modes():
    grid_raw, mode_numbers, AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ, AnX_raw, AnY_raw, AnZ_raw =         read_Anvac_fourier_all(ANVAC_NC)
    nmodes = mode_numbers.size
    ncoil = AnR.shape[1]
    BnR_modes = np.empty_like(AnR)
    Bnphi_modes = np.empty_like(AnR)
    BnZ_modes = np.empty_like(AnR)

    for idx, ntor in enumerate(mode_numbers):
        gauged_AnR, gauged_AnZ = gauge_Anvac(
            grid_raw,
            AnR[idx],
            Anphi[idx],
            AnZ[idx],
            dAnphi_dR[idx],
            dAnphi_dZ[idx],
            ntor=ntor,
        )
        spl = spline_gauged_Anvac(grid_raw, gauged_AnR, gauged_AnZ, ntor=ntor, Anphi=Anphi[idx])
        BnR_tmp, Bnphi_tmp, BnZ_tmp = field_divfree(spl, grid_raw.R, grid_raw.Z, ntor=ntor)
        BnR_modes[idx] = BnR_tmp
        Bnphi_modes[idx] = Bnphi_tmp
        BnZ_modes[idx] = BnZ_tmp

    inv_R = np.divide(
        1.0,
        grid_raw.R.reshape((1, 1, -1, 1)),
        out=np.zeros((1, 1, grid_raw.nR, 1), dtype=float),
        where=grid_raw.R.reshape((1, 1, -1, 1)) > 0.0,
    )
    try:
        dAnR_dZ = np.gradient(AnR, grid_raw.Z, axis=3, edge_order=2)
        dAnZ_dR = np.gradient(AnZ, grid_raw.R, axis=2, edge_order=2)
    except ValueError:
        dAnR_dZ = np.gradient(AnR, grid_raw.Z, axis=3, edge_order=1)
        dAnZ_dR = np.gradient(AnZ, grid_raw.R, axis=2, edge_order=1)
    ntor_vals = mode_numbers.reshape((-1, 1, 1, 1)).astype(float)
    BnR_ung = 1j * ntor_vals * AnZ * inv_R - dAnphi_dZ
    Bnphi_ung = dAnR_dZ - dAnZ_dR
    BnZ_ung = Anphi * inv_R + dAnphi_dR - 1j * ntor_vals * AnR * inv_R

    return grid_raw, mode_numbers, BnR_modes, Bnphi_modes, BnZ_modes, BnR_ung, Bnphi_ung, BnZ_ung


def main() -> None:
    grid_ref_raw, mode_numbers, BnR_ref, Bnphi_ref, BnZ_ref = read_Bnvac_fourier_all(REFERENCE_H5)
    (
        grid_anv_raw,
        mode_numbers_anv,
        BnR_anv_modes,
        Bnphi_anv_modes,
        BnZ_anv_modes,
        BnR_ung_modes,
        Bnphi_ung_modes,
        BnZ_ung_modes,
    ) = _compute_anvac_modes()

    if not np.array_equal(mode_numbers, mode_numbers_anv):
        raise SystemExit("Mode sets differ between Fourier and Anvac data")

    currents = np.atleast_1d(np.loadtxt(CURRENT_FILE))
    weights = currents * PREFactor

    BnR_ref_sum = np.tensordot(weights, BnR_ref, axes=(0, 1))
    Bnphi_ref_sum = np.tensordot(weights, Bnphi_ref, axes=(0, 1))
    BnZ_ref_sum = np.tensordot(weights, BnZ_ref, axes=(0, 1))

    BnR_anv_sum = np.tensordot(weights, BnR_anv_modes, axes=(0, 1))
    Bnphi_anv_sum = np.tensordot(weights, Bnphi_anv_modes, axes=(0, 1))
    BnZ_anv_sum = np.tensordot(weights, BnZ_anv_modes, axes=(0, 1))
    BnR_ung_sum = np.tensordot(weights, BnR_ung_modes, axes=(0, 1))
    Bnphi_ung_sum = np.tensordot(weights, Bnphi_ung_modes, axes=(0, 1))
    BnZ_ung_sum = np.tensordot(weights, BnZ_ung_modes, axes=(0, 1))

    # Evaluate total field at coil centre (in cm)
    x_c, y_c, z_c = 197.2682697, 72.00853957, 78.0
    Bx_ref, By_ref, Bz_ref = reconstruct_field_from_modes(
        BnR_ref_sum,
        Bnphi_ref_sum,
        BnZ_ref_sum,
        mode_numbers,
        grid_ref_raw.R,
        grid_ref_raw.Z,
        x_c,
        y_c,
        z_c,
    )
    Bx_anv, By_anv, Bz_anv = reconstruct_field_from_modes(
        BnR_anv_sum,
        Bnphi_anv_sum,
        BnZ_anv_sum,
        mode_numbers,
        grid_anv_raw.R,
        grid_anv_raw.Z,
        x_c,
        y_c,
        z_c,
    )
    Bx_ung, By_ung, Bz_ung = reconstruct_field_from_modes(
        BnR_ung_sum,
        Bnphi_ung_sum,
        BnZ_ung_sum,
        mode_numbers,
        grid_anv_raw.R,
        grid_anv_raw.Z,
        x_c,
        y_c,
        z_c,
    )

    B_mag_ref = np.sqrt(Bx_ref**2 + By_ref**2 + Bz_ref**2)
    B_mag_anv = np.sqrt(Bx_anv**2 + By_anv**2 + Bz_anv**2)
    B_mag_ung = np.sqrt(Bx_ung**2 + By_ung**2 + Bz_ung**2)

    print("=" * 70)
    print("Total field comparison at coil centre (all modes)")
    print("=" * 70)
    print(f"Fourier total B (G):   ({Bx_ref:.6e}, {By_ref:.6e}, {Bz_ref:.6e}) -> |B|={B_mag_ref:.6e}")
    print(f"Anvac total B (G):     ({Bx_anv:.6e}, {By_anv:.6e}, {Bz_anv:.6e}) -> |B|={B_mag_anv:.6e}")
    print(f"Anvac ungauged B (G): ({Bx_ung:.6e}, {By_ung:.6e}, {Bz_ung:.6e}) -> |B|={B_mag_ung:.6e}")

    diff_vec = np.array([Bx_anv - Bx_ref, By_anv - By_ref, Bz_anv - Bz_ref])
    rel_err = np.abs(diff_vec) / np.maximum(np.abs([Bx_ref, By_ref, Bz_ref]), 1e-15)
    diff_ung = np.array([Bx_ung - Bx_ref, By_ung - By_ref, Bz_ung - Bz_ref])
    rel_err_ung = np.abs(diff_ung) / np.maximum(np.abs([Bx_ref, By_ref, Bz_ref]), 1e-15)
    print("\nComponent-wise relative errors (stored):")
    print(f"  Bx: {rel_err[0]*100:.6f}%")
    print(f"  By: {rel_err[1]*100:.6f}%")
    print(f"  Bz: {rel_err[2]*100:.6f}%")
    print("\nComponent-wise relative errors (ungauged):")
    print(f"  Bx: {rel_err_ung[0]*100:.6f}%")
    print(f"  By: {rel_err_ung[1]*100:.6f}%")
    print(f"  Bz: {rel_err_ung[2]*100:.6f}%")
    print(f"\nMagnitude relative error (stored): {(abs(B_mag_anv - B_mag_ref)/max(B_mag_ref,1e-15))*100:.6f}%")
    print(f"Magnitude relative error (ungauged): {(abs(B_mag_ung - B_mag_ref)/max(B_mag_ref,1e-15))*100:.6f}%")

    if np.all(rel_err < 1e-2):
        print("\n✓ Fourier and gauged Anvac agree to better than 1% for the full field sum.")
    else:
        print("\n✗ Detected >1% discrepancy in full-field comparison (gauged).")

    if np.all(rel_err_ung < 1e-2):
        print("✓ Ungauged reconstruction matches Fourier to better than 1%.")
    else:
        print("✗ Ungauged reconstruction deviates by more than 1%.")


if __name__ == '__main__':
    main()
