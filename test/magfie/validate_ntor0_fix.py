#!/usr/bin/env python3
"""Validate ntor=0 fix: ensure formulas don't give zero field."""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import sys

def main():
    from libneo.biotsavart_fourier import (
        read_Anvac_fourier, gauge_Anvac, spline_gauged_Anvac, field_divfree
    )

    # Load ntor=0 Fourier modes from vector potential
    grid, AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ = read_Anvac_fourier(
        'tilted_coil_Anvac.nc', ntor=0
    )

    # Apply gauge transformation
    gauged_AnR, gauged_AnZ = gauge_Anvac(
        grid, AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ, ntor=0
    )

    # Spline including Anphi
    spl = spline_gauged_Anvac(grid, gauged_AnR, gauged_AnZ, ntor=0, Anphi=Anphi)

    # Compute field from Fourier modes using the new ntor=0 formulas
    BnR_new, Bnphi_new, BnZ_new = field_divfree(
        spl, grid.R, grid.Z, ntor=0
    )

    # Test that the fix produces non-zero fields
    # The bug was that the original formulas gave exactly zero for B_R and B_Z when ntor=0
    print(f"\nField statistics:")
    print(f"  BR:   max={np.abs(BnR_new).max():.3e} G, mean={np.abs(BnR_new).mean():.3e} G")
    print(f"  Bphi: max={np.abs(Bnphi_new).max():.3e} G, mean={np.abs(Bnphi_new).mean():.3e} G")
    print(f"  BZ:   max={np.abs(BnZ_new).max():.3e} G, mean={np.abs(BnZ_new).mean():.3e} G")

    # Test points along vertical line at R = 2.1 m
    R_test = 2.1
    Z_vals = np.linspace(-1.0, 2.0, 30)

    BR_new_list = []
    Bphi_new_list = []
    BZ_new_list = []

    for Z_test in Z_vals:
        # Fourier reconstruction at test point
        idx_R = np.argmin(np.abs(grid.R - R_test))
        idx_Z = np.argmin(np.abs(grid.Z - Z_test))
        # Handle squeezed arrays (single coil case)
        if BnR_new.ndim == 2:
            BR_f = BnR_new[idx_R, idx_Z].real
            Bphi_f = Bnphi_new[idx_R, idx_Z].real
            BZ_f = BnZ_new[idx_R, idx_Z].real
        else:
            BR_f = BnR_new[0, idx_R, idx_Z].real
            Bphi_f = Bnphi_new[0, idx_R, idx_Z].real
            BZ_f = BnZ_new[0, idx_R, idx_Z].real

        BR_new_list.append(BR_f)
        Bphi_new_list.append(Bphi_f)
        BZ_new_list.append(BZ_f)

    BR_new_arr = np.array(BR_new_list)
    Bphi_new_arr = np.array(Bphi_new_list)
    BZ_new_arr = np.array(BZ_new_list)

    # Create visualization plot
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    axes[0].plot(Z_vals, BR_new_arr, 'r^-', label='ntor=0 (with fix)', linewidth=2, markersize=6)
    axes[0].axhline(0, color='k', linestyle='--', alpha=0.3)
    axes[0].set_xlabel('Z (m)', fontsize=12)
    axes[0].set_ylabel('$B_R$ (G)', fontsize=12)
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    axes[0].set_title('Radial Component')

    axes[1].plot(Z_vals, Bphi_new_arr, 'r^-', label='ntor=0 (with fix)', linewidth=2, markersize=6)
    axes[1].axhline(0, color='k', linestyle='--', alpha=0.3)
    axes[1].set_xlabel('Z (m)', fontsize=12)
    axes[1].set_ylabel('$B_\\phi$ (G)', fontsize=12)
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    axes[1].set_title('Toroidal Component')

    axes[2].plot(Z_vals, BZ_new_arr, 'r^-', label='ntor=0 (with fix)', linewidth=2, markersize=6)
    axes[2].axhline(0, color='k', linestyle='--', alpha=0.3)
    axes[2].set_xlabel('Z (m)', fontsize=12)
    axes[2].set_ylabel('$B_Z$ (G)', fontsize=12)
    axes[2].legend()
    axes[2].grid(True, alpha=0.3)
    axes[2].set_title('Vertical Component')

    plt.tight_layout()
    plt.savefig('ntor0_validation.png', dpi=150)
    print("\nSaved plot to ntor0_validation.png")

    # Test criterion: fields must be non-trivial (not all zero)
    # The bug caused B_R and B_Z to be exactly zero when using the wrong formulas
    BR_nonzero = np.abs(BR_new_arr).max() > 1e-3
    BZ_nonzero = np.abs(BZ_new_arr).max() > 1e-3

    print(f"\nValidation checks:")
    print(f"  B_R is non-zero: {BR_nonzero} (max = {np.abs(BR_new_arr).max():.3e} G)")
    print(f"  B_Z is non-zero: {BZ_nonzero} (max = {np.abs(BZ_new_arr).max():.3e} G)")
    print(f"  B_φ magnitude: {np.abs(Bphi_new_arr).max():.3e} G (expected ~0 for axisymmetric coil)")

    if BR_nonzero and BZ_nonzero:
        print(f"\n✓ TEST PASSED: ntor=0 formulas produce non-zero B_R and B_Z components")
        return 0
    else:
        print(f"\n✗ TEST FAILED: ntor=0 fields are still zero (bug not fixed)")
        return 1

if __name__ == '__main__':
    sys.exit(main())
