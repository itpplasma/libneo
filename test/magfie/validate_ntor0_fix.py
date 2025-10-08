#!/usr/bin/env python3
"""Validate ntor=0 fix against direct Biot-Savart along tilted coil axis."""

import numpy as np
from numpy import pi, cos, sin, sqrt
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import sys

from libneo.biotsavart_fourier import (
    read_Anvac_fourier, gauge_Anvac, spline_gauged_Anvac
)

def coil_orientation_vectors(phi0, tilt_theta_deg, tilt_psi_deg):
    """Construct in-plane basis and axis vectors for the tilted coil."""
    theta = np.deg2rad(tilt_theta_deg)
    psi = np.deg2rad(tilt_psi_deg)

    e_R = np.array([cos(phi0), sin(phi0), 0.0])
    e_phi = np.array([-sin(phi0), cos(phi0), 0.0])
    e_Z = np.array([0.0, 0.0, 1.0])

    n_vec = sin(theta) * cos(psi) * e_R + sin(theta) * sin(psi) * e_phi + cos(theta) * e_Z
    n_vec = n_vec / np.linalg.norm(n_vec)

    candidate = np.cross(n_vec, e_Z)
    if np.linalg.norm(candidate) < 1e-12:
        candidate = np.cross(n_vec, e_R)
    u_vec = candidate / np.linalg.norm(candidate)
    v_vec = np.cross(n_vec, u_vec)
    return u_vec, v_vec, n_vec


def circular_coil_analytical(s_vals, center_xyz, axis_hat, coil_radius, current=1.0):
    """Analytical magnetic field along the coil axis for a circular coil."""
    s_vals = np.asarray(s_vals, dtype=float)
    axis_hat = axis_hat / np.linalg.norm(axis_hat)
    axis_points = center_xyz[:, None] + axis_hat[:, None] * s_vals[None, :]
    phi_points = np.mod(np.arctan2(axis_points[1, :], axis_points[0, :]), 2.0 * pi)

    B_mag = (2.0 * pi * current * coil_radius**2) / (coil_radius**2 + s_vals**2) ** 1.5
    B_cart = -axis_hat[:, None] * B_mag[None, :]

    cosphi = np.cos(phi_points)
    sinphi = np.sin(phi_points)
    BR = B_cart[0, :] * cosphi + B_cart[1, :] * sinphi
    Bphi = -B_cart[0, :] * sinphi + B_cart[1, :] * cosphi
    BZ = B_cart[2, :]
    return BR, Bphi, BZ, axis_points, phi_points

def main():
    # Coil parameters (from generate_tilted_coil.py)
    R_center = 2.0
    phi_center = 0.35
    Z_center = 0.8
    coil_radius = 2.0
    tilt_theta = 30.0
    tilt_psi = 35.0
    current = 1.0  # abampere

    s_vals = np.linspace(-2.0, 2.0, 30)  # Distance along tilted axis
    center_xyz = np.array([R_center * cos(phi_center), R_center * sin(phi_center), Z_center])
    _, _, axis_hat = coil_orientation_vectors(phi_center, tilt_theta, tilt_psi)
    BR_analytical_full, Bphi_analytical_full, BZ_analytical_full, axis_points, phi_vals = (
        circular_coil_analytical(s_vals, center_xyz, axis_hat, coil_radius, current=current)
    )
    X_axis, Y_axis, Z_axis = axis_points
    R_vals = np.sqrt(X_axis**2 + Y_axis**2)

    # Load nmax from NetCDF to know how many modes to sum
    import netCDF4
    with netCDF4.Dataset('tilted_coil_Anvac.nc') as nc:
        nphi = int(nc.variables['nphi'][:])
        ntor_array = nc.variables['ntor'][:]
        nmax = int(np.max(ntor_array))

    # Read ntor=0 data once to obtain grid and cache vectors
    grid0, AnR0, Anphi0, AnZ0, dAnphi_dR0, dAnphi_dZ0 = read_Anvac_fourier(
        'tilted_coil_Anvac.nc', ntor=0
    )

    valid_mask = (
        (R_vals >= grid0.R_min) & (R_vals <= grid0.R_max) &
        (Z_axis >= grid0.Z_min) & (Z_axis <= grid0.Z_max)
    )
    if not np.all(valid_mask):
        clipped = np.count_nonzero(~valid_mask)
        print(f"Warning: {clipped} evaluation points lie outside spline domain and will be skipped.")

    R_eval = R_vals[valid_mask]
    Z_eval = Z_axis[valid_mask]
    s_eval = s_vals[valid_mask]
    phi_eval = phi_vals[valid_mask]

    if R_eval.size == 0:
        print("No evaluation points fall inside the spline grid. Aborting.")
        return 1

    BR_analytical = BR_analytical_full[valid_mask]
    Bphi_analytical = Bphi_analytical_full[valid_mask]
    BZ_analytical = BZ_analytical_full[valid_mask]

    # Sum ALL Fourier modes from ntor=0 to ntor=nmax for full field reconstruction
    npts = R_eval.size
    BR_fourier_gauged = np.zeros(npts)
    Bphi_fourier_gauged = np.zeros(npts)
    BZ_fourier_gauged = np.zeros(npts)
    BR_fourier_ungauged = np.zeros(npts)
    Bphi_fourier_ungauged = np.zeros(npts)
    BZ_fourier_ungauged = np.zeros(npts)

    def evaluate_mode_at_points(spl, R_points, Z_points, ntor):
        """Return summed coil contribution of Fourier mode at specific points."""
        from numpy import sum as np_sum

        R_points = np.asarray(R_points, dtype=float)
        Z_points = np.asarray(Z_points, dtype=float)
        npts_local = R_points.size
        ncoil = len(spl['AnR_Re'])
        BnR = np.empty((ncoil, npts_local), dtype=complex)
        Bnphi = np.empty_like(BnR)
        BnZ = np.empty_like(BnR)
        has_Anphi = 'Anphi_Re' in spl

        for kcoil in range(ncoil):
            AnR_re = spl['AnR_Re'][kcoil].ev(R_points, Z_points)
            AnR_im = spl['AnR_Im'][kcoil].ev(R_points, Z_points)
            AnZ_re = spl['AnZ_Re'][kcoil].ev(R_points, Z_points)
            AnZ_im = spl['AnZ_Im'][kcoil].ev(R_points, Z_points)
            dAnR_dZ_re = spl['AnR_Re'][kcoil].ev(R_points, Z_points, dy=1)
            dAnR_dZ_im = spl['AnR_Im'][kcoil].ev(R_points, Z_points, dy=1)
            dAnZ_dR_re = spl['AnZ_Re'][kcoil].ev(R_points, Z_points, dx=1)
            dAnZ_dR_im = spl['AnZ_Im'][kcoil].ev(R_points, Z_points, dx=1)

            AnR = AnR_re + 1j * AnR_im
            AnZ = AnZ_re + 1j * AnZ_im
            dAnR_dZ = dAnR_dZ_re + 1j * dAnR_dZ_im
            dAnZ_dR = dAnZ_dR_re + 1j * dAnZ_dR_im

            if has_Anphi:
                Anphi_re = spl['Anphi_Re'][kcoil].ev(R_points, Z_points)
                Anphi_im = spl['Anphi_Im'][kcoil].ev(R_points, Z_points)
                dAnphi_dR_re = spl['Anphi_Re'][kcoil].ev(R_points, Z_points, dx=1)
                dAnphi_dR_im = spl['Anphi_Im'][kcoil].ev(R_points, Z_points, dx=1)
                dAnphi_dZ_re = spl['Anphi_Re'][kcoil].ev(R_points, Z_points, dy=1)
                dAnphi_dZ_im = spl['Anphi_Im'][kcoil].ev(R_points, Z_points, dy=1)

                Anphi = Anphi_re + 1j * Anphi_im
                dAnphi_dR = dAnphi_dR_re + 1j * dAnphi_dR_im
                dAnphi_dZ = dAnphi_dZ_re + 1j * dAnphi_dZ_im

                BnR[kcoil, :] = 1j * ntor * AnZ / R_points - dAnphi_dZ
                Bnphi[kcoil, :] = dAnR_dZ - dAnZ_dR
                BnZ[kcoil, :] = dAnphi_dR + Anphi / R_points - 1j * ntor * AnR / R_points
            else:
                if ntor == 0:
                    raise ValueError("For ntor=0, Aphi splines must be provided")
                BnR[kcoil, :] = 1j * ntor * AnZ / R_points
                Bnphi[kcoil, :] = dAnR_dZ - dAnZ_dR
                BnZ[kcoil, :] = -1j * ntor * AnR / R_points

        return (
            np_sum(BnR, axis=0),
            np_sum(Bnphi, axis=0),
            np_sum(BnZ, axis=0),
        )

    data_cache = {
        0: (grid0, AnR0, Anphi0, AnZ0, dAnphi_dR0, dAnphi_dZ0),
    }

    for ntor in range(0, nmax + 1):
        if ntor in data_cache:
            grid, AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ = data_cache[ntor]
        else:
            grid, AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ = read_Anvac_fourier(
                'tilted_coil_Anvac.nc', ntor=ntor
            )

        # Ungauged (direct Biot-Savart gauge)
        spl_ungauged = spline_gauged_Anvac(grid, AnR, AnZ, ntor=ntor, Anphi=Anphi)
        BnR_ungauged, Bnphi_ungauged, BnZ_ungauged = evaluate_mode_at_points(
            spl_ungauged, R_eval, Z_eval, ntor
        )

        # Apply gauge transformation
        gauged_AnR, gauged_AnZ = gauge_Anvac(
            grid, AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ, ntor=ntor
        )

        # Spline for gauged fields (keep ntor=0 ungauged by reusing original Anphi)
        spl_gauged = spline_gauged_Anvac(
            grid,
            gauged_AnR,
            gauged_AnZ,
            ntor=ntor,
            Anphi=Anphi if ntor == 0 else None,
        )

        # Compute field from this Fourier mode
        BnR_gauged, Bnphi_gauged, BnZ_gauged = evaluate_mode_at_points(
            spl_gauged, R_eval, Z_eval, ntor
        )

        phase = np.exp(1j * ntor * phi_eval)
        BR_fourier_ungauged += np.real(BnR_ungauged * phase)
        Bphi_fourier_ungauged += np.real(Bnphi_ungauged * phase)
        BZ_fourier_ungauged += np.real(BnZ_ungauged * phase)

        if ntor == 0:
            BR_fourier_gauged += np.real(BnR_ungauged * phase)
            Bphi_fourier_gauged += np.real(Bnphi_ungauged * phase)
            BZ_fourier_gauged += np.real(BnZ_ungauged * phase)
        else:
            BR_fourier_gauged += np.real(BnR_gauged * phase)
            Bphi_fourier_gauged += np.real(Bnphi_gauged * phase)
            BZ_fourier_gauged += np.real(BnZ_gauged * phase)

    print(f"Summed Fourier modes ntor=0 to ntor={nmax} for full field reconstruction")

    # Compute analytical field on coil axis
    # Create comparison plot
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    axes[0].plot(s_eval, BR_analytical, 'ko-', label='Analytical (on-axis)', linewidth=2, markersize=4)
    axes[0].plot(s_eval, BR_fourier_ungauged, 'r^--', label='Fourier Sum (ungauged)', linewidth=2, markersize=6)
    axes[0].plot(s_eval, BR_fourier_gauged, 'bs-.', label='Fourier Sum (gauged n>0)', linewidth=2, markersize=5)
    axes[0].set_xlabel('Distance along coil axis (m)', fontsize=12)
    axes[0].set_ylabel('$B_R$ (G)', fontsize=12)
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    axes[0].set_title('Radial Component')

    axes[1].plot(s_eval, Bphi_analytical, 'ko-', label='Analytical (on-axis)', linewidth=2, markersize=4)
    axes[1].plot(s_eval, Bphi_fourier_ungauged, 'r^--', label='Fourier Sum (ungauged)', linewidth=2, markersize=6)
    axes[1].plot(s_eval, Bphi_fourier_gauged, 'bs-.', label='Fourier Sum (gauged n>0)', linewidth=2, markersize=5)
    axes[1].set_xlabel('Distance along coil axis (m)', fontsize=12)
    axes[1].set_ylabel('$B_\\phi$ (G)', fontsize=12)
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    axes[1].set_title('Toroidal Component')

    axes[2].plot(s_eval, BZ_analytical, 'ko-', label='Analytical (on-axis)', linewidth=2, markersize=4)
    axes[2].plot(s_eval, BZ_fourier_ungauged, 'r^--', label='Fourier Sum (ungauged)', linewidth=2, markersize=6)
    axes[2].plot(s_eval, BZ_fourier_gauged, 'bs-.', label='Fourier Sum (gauged n>0)', linewidth=2, markersize=5)
    axes[2].set_xlabel('Distance along coil axis (m)', fontsize=12)
    axes[2].set_ylabel('$B_Z$ (G)', fontsize=12)
    axes[2].legend()
    axes[2].grid(True, alpha=0.3)
    axes[2].set_title('Vertical Component')

    plt.tight_layout()
    plt.savefig('ntor0_validation.png', dpi=150)
    print("\nSaved plot to ntor0_validation.png")

    # Compute relative errors
    def relative_errors(reconstructed):
        mask = np.abs(BR_analytical) > 1e-5
        err_BR = (
            np.mean(np.abs(reconstructed[0][mask] - BR_analytical[mask]) / np.abs(BR_analytical[mask])) * 100
            if mask.sum() > 0 else 0
        )

        mask = np.abs(Bphi_analytical) > 1e-5
        err_Bphi = (
            np.mean(np.abs(reconstructed[1][mask] - Bphi_analytical[mask]) / np.abs(Bphi_analytical[mask])) * 100
            if mask.sum() > 0 else 0
        )

        mask = np.abs(BZ_analytical) > 1e-5
        err_BZ = (
            np.mean(np.abs(reconstructed[2][mask] - BZ_analytical[mask]) / np.abs(BZ_analytical[mask])) * 100
            if mask.sum() > 0 else 0
        )
        return err_BR, err_Bphi, err_BZ

    err_gauged = relative_errors(
        (BR_fourier_gauged, Bphi_fourier_gauged, BZ_fourier_gauged)
    )
    err_ungauged = relative_errors(
        (BR_fourier_ungauged, Bphi_fourier_ungauged, BZ_fourier_ungauged)
    )

    print(f"\nRelative errors along coil axis (mean where |B_analytical| > 1e-5 G):")
    print(
        "  Gauged   -> BR: {:.2f}%, Bphi: {:.2f}%, BZ: {:.2f}%".format(*err_gauged)
    )
    print(
        "  Ungauged -> BR: {:.2f}%, Bphi: {:.2f}%, BZ: {:.2f}%".format(*err_ungauged)
    )

    max_err_gauged = max(err_gauged)
    max_err_ungauged = max(err_ungauged)
    max_err = max(max_err_gauged, max_err_ungauged)

    if max_err < 5.0:
        print(
            "\n✓ TEST PASSED: Fourier reconstructions (gauged and ungauged) match analytical on-axis field "
            f"(worst-case error {max_err:.2f}%)"
        )
        return 0
    else:
        print(
            "\n✗ TEST FAILED: Fourier reconstructions deviate from analytical field "
            f"(worst-case error {max_err:.2f}%)"
        )
        return 1

if __name__ == '__main__':
    sys.exit(main())
