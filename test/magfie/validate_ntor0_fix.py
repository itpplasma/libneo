#!/usr/bin/env python3
"""Validate ntor=0 fix against direct Biot-Savart along tilted coil axis."""

import numpy as np
from numpy import pi, cos, sin, sqrt
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import sys

def circular_coil_analytical(s_vals, R_center, Z_center, tilt_rad, coil_radius, current=1.0):
    """Analytical field on axis of tilted circular coil in CGS units.

    For a circular coil of radius a with current I, the field on axis at distance z from center is:
    B = (2*pi*I*a^2) / (a^2 + z^2)^(3/2)  (in CGS: I in abamperes, distances in cm)

    For tilted coil: transform observation points to coil frame, compute field along coil axis,
    then transform back to cylindrical coordinates.
    """
    BR_list = []
    Bphi_list = []
    BZ_list = []

    for s in s_vals:
        # Position along tilted coil axis in lab frame
        R_obs = R_center + s * cos(tilt_rad)
        Z_obs = Z_center + s * sin(tilt_rad)

        # Distance along coil axis (in coil frame, this is the z coordinate)
        z_coil = s

        # Analytical field magnitude along coil axis (pointing along axis direction)
        a = coil_radius
        B_mag = (2.0 * pi * current * a**2) / (a**2 + z_coil**2)**(1.5)

        # Field direction is along the coil axis (tilted in R-Z plane)
        # Coil axis direction: (cos(tilt), 0, sin(tilt)) in cylindrical coords
        # This gives B_R and B_Z components, B_phi = 0 by symmetry
        BR_list.append(-B_mag * cos(tilt_rad))  # Negative because field points opposite to axis direction on one side
        Bphi_list.append(0.0)
        BZ_list.append(-B_mag * sin(tilt_rad))

    return np.array(BR_list), np.array(Bphi_list), np.array(BZ_list)

def main():
    from libneo.biotsavart_fourier import (
        read_Anvac_fourier, gauge_Anvac, spline_gauged_Anvac, field_divfree
    )

    # Coil parameters (from generate_tilted_coil.py)
    R_center = 2.0
    Z_center = 0.8
    coil_radius = 2.0
    tilt_deg = 30.0
    tilt_rad = np.deg2rad(tilt_deg)
    current = 1.0  # abampere

    # Tilted coil axis: passes through center at angle to Z-axis
    # Points along axis: (R, Z) = (R_center + s*cos(tilt), Z_center + s*sin(tilt))
    s_vals = np.linspace(-2.0, 2.0, 30)  # Distance along tilted axis
    R_vals = R_center + s_vals * cos(tilt_rad)
    Z_vals = Z_center + s_vals * sin(tilt_rad)
    phi_test = 0.0  # At phi=0 plane

    # Load nmax from NetCDF to know how many modes to sum
    import netCDF4
    with netCDF4.Dataset('tilted_coil_Anvac.nc') as nc:
        nphi = int(nc.variables['nphi'][:])
        ntor_array = nc.variables['ntor'][:]
        nmax = int(np.max(ntor_array))

    # Sum ALL Fourier modes from ntor=0 to ntor=nmax for full field reconstruction
    BR_fourier_full = None
    Bphi_fourier_full = None
    BZ_fourier_full = None

    for ntor in range(0, nmax + 1):
        grid, AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ = read_Anvac_fourier(
            'tilted_coil_Anvac.nc', ntor=ntor
        )

        # Apply gauge transformation
        gauged_AnR, gauged_AnZ = gauge_Anvac(
            grid, AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ, ntor=ntor
        )

        # Spline including Anphi
        spl = spline_gauged_Anvac(grid, gauged_AnR, gauged_AnZ, ntor=ntor, Anphi=Anphi)

        # Compute field from this Fourier mode
        BnR, Bnphi, BnZ = field_divfree(spl, grid.R, grid.Z, ntor=ntor)

        # Sum contributions (for ntor>0, need to evaluate cos(ntor*phi) and sin(ntor*phi) at phi=0)
        # At phi=0: cos(ntor*0) = 1, sin(ntor*0) = 0, so only real part contributes
        if BR_fourier_full is None:
            BR_fourier_full = BnR.real
            Bphi_fourier_full = Bnphi.real
            BZ_fourier_full = BnZ.real
        else:
            BR_fourier_full += BnR.real
            Bphi_fourier_full += Bnphi.real
            BZ_fourier_full += BnZ.real

    print(f"Summed Fourier modes ntor=0 to ntor={nmax} for full field reconstruction")

    # Compute analytical field on coil axis
    BR_analytical, Bphi_analytical, BZ_analytical = circular_coil_analytical(
        s_vals, R_center, Z_center, tilt_rad, coil_radius, current=current
    )

    # Extract Fourier reconstruction along tilted axis from the full summed field
    BR_fourier_list = []
    Bphi_fourier_list = []
    BZ_fourier_list = []

    for R, Z in zip(R_vals, Z_vals):
        idx_R = np.argmin(np.abs(grid.R - R))
        idx_Z = np.argmin(np.abs(grid.Z - Z))
        # Handle squeezed arrays (single coil case)
        if BR_fourier_full.ndim == 2:
            BR_f = BR_fourier_full[idx_R, idx_Z]
            Bphi_f = Bphi_fourier_full[idx_R, idx_Z]
            BZ_f = BZ_fourier_full[idx_R, idx_Z]
        else:
            BR_f = BR_fourier_full[0, idx_R, idx_Z]
            Bphi_f = Bphi_fourier_full[0, idx_R, idx_Z]
            BZ_f = BZ_fourier_full[0, idx_R, idx_Z]

        BR_fourier_list.append(BR_f)
        Bphi_fourier_list.append(Bphi_f)
        BZ_fourier_list.append(BZ_f)

    BR_fourier_arr = np.array(BR_fourier_list)
    Bphi_fourier_arr = np.array(Bphi_fourier_list)
    BZ_fourier_arr = np.array(BZ_fourier_list)

    # Create comparison plot
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    axes[0].plot(s_vals, BR_analytical, 'ko-', label='Analytical (on-axis)', linewidth=2, markersize=4)
    axes[0].plot(s_vals, BR_fourier_arr, 'r^--', label='Fourier Sum (all modes)', linewidth=2, markersize=6)
    axes[0].set_xlabel('Distance along coil axis (m)', fontsize=12)
    axes[0].set_ylabel('$B_R$ (G)', fontsize=12)
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    axes[0].set_title('Radial Component')

    axes[1].plot(s_vals, Bphi_analytical, 'ko-', label='Analytical (on-axis)', linewidth=2, markersize=4)
    axes[1].plot(s_vals, Bphi_fourier_arr, 'r^--', label='Fourier Sum (all modes)', linewidth=2, markersize=6)
    axes[1].set_xlabel('Distance along coil axis (m)', fontsize=12)
    axes[1].set_ylabel('$B_\\phi$ (G)', fontsize=12)
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    axes[1].set_title('Toroidal Component')

    axes[2].plot(s_vals, BZ_analytical, 'ko-', label='Analytical (on-axis)', linewidth=2, markersize=4)
    axes[2].plot(s_vals, BZ_fourier_arr, 'r^--', label='Fourier Sum (all modes)', linewidth=2, markersize=6)
    axes[2].set_xlabel('Distance along coil axis (m)', fontsize=12)
    axes[2].set_ylabel('$B_Z$ (G)', fontsize=12)
    axes[2].legend()
    axes[2].grid(True, alpha=0.3)
    axes[2].set_title('Vertical Component')

    plt.tight_layout()
    plt.savefig('ntor0_validation.png', dpi=150)
    print("\nSaved plot to ntor0_validation.png")

    # Compute relative errors
    mask = np.abs(BR_analytical) > 1e-5
    if mask.sum() > 0:
        err_BR = np.mean(np.abs(BR_fourier_arr[mask] - BR_analytical[mask]) / np.abs(BR_analytical[mask])) * 100
    else:
        err_BR = 0

    mask = np.abs(Bphi_analytical) > 1e-5
    if mask.sum() > 0:
        err_Bphi = np.mean(np.abs(Bphi_fourier_arr[mask] - Bphi_analytical[mask]) / np.abs(Bphi_analytical[mask])) * 100
    else:
        err_Bphi = 0

    mask = np.abs(BZ_analytical) > 1e-5
    if mask.sum() > 0:
        err_BZ = np.mean(np.abs(BZ_fourier_arr[mask] - BZ_analytical[mask]) / np.abs(BZ_analytical[mask])) * 100
    else:
        err_BZ = 0

    print(f"\nRelative errors along coil axis (mean where |B_analytical| > 1e-5 G):")
    print(f"  BR:   {err_BR:.2f}%")
    print(f"  Bphi: {err_Bphi:.2f}%")
    print(f"  BZ:   {err_BZ:.2f}%")

    # Test criterion: Full Fourier sum should match analytical solution within 5%
    max_err = max(err_BR, err_Bphi, err_BZ)

    if max_err < 5.0:
        print(f"\n✓ TEST PASSED: Full Fourier reconstruction matches analytical on-axis field (max error {max_err:.2f}%)")
        return 0
    else:
        print(f"\n✗ TEST FAILED: Fourier reconstruction does not match analytical field (max error {max_err:.2f}%)")
        return 1

if __name__ == '__main__':
    sys.exit(main())
