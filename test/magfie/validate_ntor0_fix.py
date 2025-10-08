#!/usr/bin/env python3
"""Validate ntor=0 fix and produce comparison plots."""

import numpy as np
from numpy import pi, cos, sin
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import sys

def biot_savart_direct(coil_xyz, R, phi, Z, current=1.0):
    """Compute field at (R,phi,Z) using direct Biot-Savart summation."""
    mu0 = 4e-7 * pi
    x_obs = R * cos(phi)
    y_obs = R * sin(phi)
    z_obs = Z
    r_obs = np.array([x_obs, y_obs, z_obs])

    B_cart = np.zeros(3)
    nseg = coil_xyz.shape[1]
    for i in range(nseg):
        i_next = (i + 1) % nseg
        r_i = coil_xyz[:, i]
        r_f = coil_xyz[:, i_next]
        dl = r_f - r_i
        r_vec_i = r_obs - r_i
        r_vec_f = r_obs - r_f
        dist_i = np.linalg.norm(r_vec_i)
        dist_f = np.linalg.norm(r_vec_f)

        if dist_i > 1e-10 and dist_f > 1e-10:
            cross_prod = np.cross(dl, r_vec_i)
            denom = dist_i * dist_f * (dist_i * dist_f + np.dot(r_vec_i, r_vec_f))
            B_cart += (mu0 * current / (4*pi)) * cross_prod * (dist_i + dist_f) / denom

    B_SI = B_cart
    B_gauss = B_SI * 1e4  # Tesla to Gauss

    # Convert to cylindrical
    Bx, By, Bz = B_gauss
    BR = Bx * cos(phi) + By * sin(phi)
    Bphi = -Bx * sin(phi) + By * cos(phi)
    BZ = Bz

    return BR, Bphi, BZ

def main():
    from libneo.biotsavart_fourier import (
        read_Anvac_fourier, gauge_Anvac, spline_gauged_Anvac, field_divfree
    )

    # Load coil geometry
    coil_xyz = np.loadtxt('tilted_coil.dat', skiprows=1).T

    # Load ntor=0 Fourier modes
    grid, AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ = read_Anvac_fourier(
        'tilted_coil_Anvac.nc', ntor=0
    )

    # Apply gauge transformation
    gauged_AnR, gauged_AnZ = gauge_Anvac(
        grid, AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ, ntor=0
    )

    # Spline including Anphi
    spl = spline_gauged_Anvac(grid, gauged_AnR, gauged_AnZ, ntor=0, Anphi=Anphi)

    # Compute field from Fourier modes
    BnR_fourier, Bnphi_fourier, BnZ_fourier = field_divfree(
        spl, grid.R, grid.Z, ntor=0
    )

    # Test points along Z-axis at R = 2.1 m
    R_test = 2.1
    phi_test = 0.0
    Z_vals = np.linspace(-1.0, 2.0, 30)

    BR_direct_list = []
    Bphi_direct_list = []
    BZ_direct_list = []
    BR_fourier_list = []
    Bphi_fourier_list = []
    BZ_fourier_list = []

    for Z_test in Z_vals:
        # Direct Biot-Savart
        BR_d, Bphi_d, BZ_d = biot_savart_direct(coil_xyz, R_test, phi_test, Z_test)
        BR_direct_list.append(BR_d)
        Bphi_direct_list.append(Bphi_d)
        BZ_direct_list.append(BZ_d)

        # Fourier reconstruction
        idx_R = np.argmin(np.abs(grid.R - R_test))
        idx_Z = np.argmin(np.abs(grid.Z - Z_test))
        BR_fourier_list.append(BnR_fourier[0, idx_R, idx_Z].real)
        Bphi_fourier_list.append(Bnphi_fourier[0, idx_R, idx_Z].real)
        BZ_fourier_list.append(BnZ_fourier[0, idx_R, idx_Z].real)

    # Create comparison plot
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    axes[0].plot(Z_vals, BR_direct_list, 'ko-', label='Direct Biot-Savart', linewidth=2, markersize=4)
    axes[0].plot(Z_vals, BR_fourier_list, 'r^--', label='Fourier (ntor=0)', linewidth=2, markersize=6)
    axes[0].set_xlabel('Z (m)', fontsize=12)
    axes[0].set_ylabel('$B_R$ (G)', fontsize=12)
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    axes[0].set_title('Radial Component')

    axes[1].plot(Z_vals, Bphi_direct_list, 'ko-', label='Direct Biot-Savart', linewidth=2, markersize=4)
    axes[1].plot(Z_vals, Bphi_fourier_list, 'r^--', label='Fourier (ntor=0)', linewidth=2, markersize=6)
    axes[1].set_xlabel('Z (m)', fontsize=12)
    axes[1].set_ylabel('$B_\\phi$ (G)', fontsize=12)
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    axes[1].set_title('Toroidal Component')

    axes[2].plot(Z_vals, BZ_direct_list, 'ko-', label='Direct Biot-Savart', linewidth=2, markersize=4)
    axes[2].plot(Z_vals, BZ_fourier_list, 'r^--', label='Fourier (ntor=0)', linewidth=2, markersize=6)
    axes[2].set_xlabel('Z (m)', fontsize=12)
    axes[2].set_ylabel('$B_Z$ (G)', fontsize=12)
    axes[2].legend()
    axes[2].grid(True, alpha=0.3)
    axes[2].set_title('Vertical Component')

    plt.tight_layout()
    plt.savefig('ntor0_validation.png', dpi=150)
    print("Saved plot to ntor0_validation.png")

    # Compute relative errors
    BR_direct_arr = np.array(BR_direct_list)
    Bphi_direct_arr = np.array(Bphi_direct_list)
    BZ_direct_arr = np.array(BZ_direct_list)
    BR_fourier_arr = np.array(BR_fourier_list)
    Bphi_fourier_arr = np.array(Bphi_fourier_list)
    BZ_fourier_arr = np.array(BZ_fourier_list)

    # Avoid division by very small numbers
    mask = np.abs(BR_direct_arr) > 1e-5
    if mask.sum() > 0:
        err_BR = np.mean(np.abs(BR_fourier_arr[mask] - BR_direct_arr[mask]) / np.abs(BR_direct_arr[mask])) * 100
    else:
        err_BR = 0

    mask = np.abs(Bphi_direct_arr) > 1e-5
    err_Bphi = np.mean(np.abs(Bphi_fourier_arr[mask] - Bphi_direct_arr[mask]) / np.abs(Bphi_direct_arr[mask])) * 100

    mask = np.abs(BZ_direct_arr) > 1e-5
    if mask.sum() > 0:
        err_BZ = np.mean(np.abs(BZ_fourier_arr[mask] - BZ_direct_arr[mask]) / np.abs(BZ_direct_arr[mask])) * 100
    else:
        err_BZ = 0

    print(f"\nRelative errors (mean where |B_direct| > 1e-5 G):")
    print(f"  BR:   {err_BR:.2f}%")
    print(f"  Bphi: {err_Bphi:.2f}%")
    print(f"  BZ:   {err_BZ:.2f}%")

    # Test passes if all errors < 2%
    tolerance = 2.0
    if err_BR < tolerance and err_Bphi < tolerance and err_BZ < tolerance:
        print(f"\n✓ TEST PASSED: ntor=0 reconstruction agrees with direct Biot-Savart (< {tolerance}%)")
        return 0
    else:
        print(f"\n✗ TEST FAILED: ntor=0 reconstruction has large errors (> {tolerance}%)")
        return 1

if __name__ == '__main__':
    sys.exit(main())
