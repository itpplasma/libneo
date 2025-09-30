#!/usr/bin/env python3
"""
Compare total real-space magnetic field from GPEC Fourier vs coil_tools.

This script:
1. Sums over ALL toroidal modes (n=0 to nmax)
2. Applies coil currents
3. Reconstructs real-space field at phi=0
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize

try:
    from libneo.biotsavart_fourier import (
        read_Bnvac_fourier,
        read_Anvac_fourier,
        gauged_Anvac_from_Bnvac,
        gauge_Anvac,
        spline_gauged_Anvac,
        field_divfree,
    )
except ImportError:
    print("Error: libneo module not found. Please install the Python package.")
    print("Run: pip install -e . from the repository root")
    sys.exit(1)


def read_currents(currents_file):
    """Read coil currents from file."""
    with open(currents_file, 'r') as f:
        line = f.readline().strip()
        currents = np.array([float(x) for x in line.split()])
    return currents


def compute_field_magnitude(BnR, Bnphi, BnZ):
    """Compute log10 of field magnitude squared."""
    return np.log10((BnR * np.conj(BnR) +
                     Bnphi * np.conj(Bnphi) +
                     BnZ * np.conj(BnZ)).real + 1e-20)


def main():
    parser = argparse.ArgumentParser(
        description="Compare total real-space field from GPEC Fourier vs coil_tools"
    )
    parser.add_argument("reference", type=Path, help="HDF5 file from GPEC Fourier")
    parser.add_argument("test", type=Path, help="NetCDF file from coil_tools")
    parser.add_argument("currents", type=Path, help="Coil currents file")
    parser.add_argument("-o", "--output", type=Path, default="realspace_comparison.png",
                       help="Output plot filename (default: realspace_comparison.png)")
    parser.add_argument("--nmax", type=int, default=8,
                       help="Maximum toroidal mode number (default: 8)")
    parser.add_argument("--phi", type=float, default=0.0,
                       help="Toroidal angle in radians for reconstruction (default: 0)")

    args = parser.parse_args()

    # Validate input files
    for f in [args.reference, args.test, args.currents]:
        if not f.exists():
            print(f"Error: File not found: {f}")
            sys.exit(1)

    # Read currents
    print(f"Reading currents from: {args.currents}")
    currents = read_currents(args.currents)
    print(f"  Number of currents: {len(currents)}")
    print(f"  Current range: [{np.min(currents):.1f}, {np.max(currents):.1f}] A")

    print(f"\nReconstructing real-space field at phi = {args.phi:.3f} rad")
    print(f"Using toroidal modes n = 0 to {args.nmax}")

    # Read first mode to get grid info
    ref_grid_0, _, _, _ = read_Bnvac_fourier(str(args.reference), ntor=0)
    test_grid_0, _, _, _, _, _ = read_Anvac_fourier(str(args.test), ntor=0)

    # Create evaluation grid
    nR_interp = 2 * ref_grid_0.nR - 1
    nZ_interp = 2 * ref_grid_0.nZ - 1
    R = np.linspace(ref_grid_0.R_min, ref_grid_0.R_max, nR_interp)
    Z = np.linspace(ref_grid_0.Z_min, ref_grid_0.Z_max, nZ_interp)

    print(f"Evaluation grid: {len(R)} × {len(Z)}")

    # Initialize accumulator for real-space field
    B_ref_realspace = np.zeros((3, len(R), len(Z)), dtype=complex)
    B_test_realspace = np.zeros((3, len(R), len(Z)), dtype=complex)

    # Loop over all toroidal modes (skip n=0 as it requires different gauge)
    # Note: n=0 (axisymmetric) uses A_phi, not A_R/A_Z gauge
    for n in range(1, args.nmax + 1):
        print(f"\nProcessing mode n={n}...")

        # Read reference (GPEC Fourier)
        ref_grid, BnR, Bnphi, BnZ = read_Bnvac_fourier(str(args.reference), ntor=n)

        # Read test (coil_tools)
        test_grid, AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ = read_Anvac_fourier(str(args.test), ntor=n)

        # Create splines
        ref_spl = spline_gauged_Anvac(
            ref_grid,
            *gauged_Anvac_from_Bnvac(ref_grid, BnR, BnZ, ntor=n),
            ntor=n
        )

        test_spl = spline_gauged_Anvac(
            test_grid,
            *gauge_Anvac(test_grid, AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ, ntor=n),
            ntor=n
        )

        # Evaluate on grid
        BnR_ref, Bnphi_ref, BnZ_ref = field_divfree(ref_spl, R, Z, ntor=n)
        BnR_test, Bnphi_test, BnZ_test = field_divfree(test_spl, R, Z, ntor=n)

        # Sum over coils with currents
        BnR_ref_n = np.sum(currents[:, np.newaxis, np.newaxis] * BnR_ref, axis=0)
        Bnphi_ref_n = np.sum(currents[:, np.newaxis, np.newaxis] * Bnphi_ref, axis=0)
        BnZ_ref_n = np.sum(currents[:, np.newaxis, np.newaxis] * BnZ_ref, axis=0)

        BnR_test_n = np.sum(currents[:, np.newaxis, np.newaxis] * BnR_test, axis=0)
        Bnphi_test_n = np.sum(currents[:, np.newaxis, np.newaxis] * Bnphi_test, axis=0)
        BnZ_test_n = np.sum(currents[:, np.newaxis, np.newaxis] * BnZ_test, axis=0)

        # Fourier reconstruction: B(phi) = sum_n [ B_n * exp(i*n*phi) ]
        phase_factor = np.exp(1j * n * args.phi)

        B_ref_realspace[0] += BnR_ref_n * phase_factor
        B_ref_realspace[1] += Bnphi_ref_n * phase_factor
        B_ref_realspace[2] += BnZ_ref_n * phase_factor

        B_test_realspace[0] += BnR_test_n * phase_factor
        B_test_realspace[1] += Bnphi_test_n * phase_factor
        B_test_realspace[2] += BnZ_test_n * phase_factor

        print(f"  Mode n={n} contribution: |B_ref|_max = {np.max(np.abs([BnR_ref_n, Bnphi_ref_n, BnZ_ref_n])):.3e}")

    # Take real part (imaginary should be small due to symmetry)
    B_ref_real = B_ref_realspace.real
    B_test_real = B_test_realspace.real

    print(f"\nMax imaginary part (should be small):")
    print(f"  Reference: {np.max(np.abs(B_ref_realspace.imag)):.3e}")
    print(f"  Test:      {np.max(np.abs(B_test_realspace.imag)):.3e}")

    # Compute field magnitudes
    B_ref_mag = np.sqrt(np.sum(B_ref_real**2, axis=0))
    B_test_mag = np.sqrt(np.sum(B_test_real**2, axis=0))

    log_B2_ref = np.log10(B_ref_mag**2 + 1e-20)
    log_B2_test = np.log10(B_test_mag**2 + 1e-20)

    # Compute difference
    diff = log_B2_test - log_B2_ref

    # Statistics
    print("\nReal-space field comparison statistics:")
    relative_error = np.abs(B_test_mag - B_ref_mag) / (B_ref_mag + 1e-15)
    print(f"  Mean |B| (GPEC):      {np.mean(B_ref_mag):.6e}")
    print(f"  Mean |B| (coil_tools): {np.mean(B_test_mag):.6e}")
    print(f"  Median relative error: {np.median(relative_error)*100:.4f}%")
    print(f"  Mean relative error:   {np.mean(relative_error)*100:.4f}%")
    print(f"  95th percentile error: {np.percentile(relative_error, 95)*100:.4f}%")
    print(f"  Max relative error:    {np.max(relative_error)*100:.4f}%")
    print(f"  Log difference range: [{np.min(diff):.3f}, {np.max(diff):.3f}]")

    # Create comparison plot
    print(f"\nCreating comparison plot...")

    fig = plt.figure(figsize=(15, 5), layout='constrained')
    axs = fig.subplots(1, 3)

    # Shared colorbar range for absolute fields
    vmin_abs = min(np.min(log_B2_ref), np.min(log_B2_test))
    vmax_abs = max(np.max(log_B2_ref), np.max(log_B2_test))
    norm_abs = Normalize(vmin=vmin_abs, vmax=vmax_abs)

    extent = [R[0], R[-1], Z[0], Z[-1]]

    # Plot GPEC reference
    im0 = axs[0].imshow(log_B2_ref.T, origin='lower', cmap='magma',
                        extent=extent, norm=norm_abs, aspect='auto', interpolation='bilinear')
    axs[0].set_title(f'GPEC Fourier\n(n=1 to {args.nmax}, φ={args.phi:.2f} rad)', fontsize=11, fontweight='bold')
    axs[0].set_xlabel('R [cm]')
    axs[0].set_ylabel('Z [cm]')

    # Plot coil_tools test
    im1 = axs[1].imshow(log_B2_test.T, origin='lower', cmap='magma',
                        extent=extent, norm=norm_abs, aspect='auto', interpolation='bilinear')
    axs[1].set_title(f'coil_tools vector_potential\n(n=1 to {args.nmax}, φ={args.phi:.2f} rad)', fontsize=11, fontweight='bold')
    axs[1].set_xlabel('R [cm]')
    axs[1].set_ylabel('Z [cm]')

    # Plot difference
    diff_max = max(abs(np.min(diff)), abs(np.max(diff)))
    im2 = axs[2].imshow(diff.T, origin='lower', cmap='RdBu_r',
                        extent=extent, vmin=-diff_max, vmax=diff_max, aspect='auto', interpolation='bilinear')
    axs[2].set_title('Difference (Test - Reference)', fontsize=11, fontweight='bold')
    axs[2].set_xlabel('R [cm]')
    axs[2].set_ylabel('Z [cm]')

    # Colorbars
    cbar0 = fig.colorbar(im0, ax=axs[:2], location='bottom', fraction=0.05, pad=0.08)
    cbar0.set_label(r'$\log_{10} |B|^{2}$ (real space)', fontsize=11)

    cbar2 = fig.colorbar(im2, ax=axs[2], location='bottom', fraction=0.05, pad=0.08)
    cbar2.set_label(r'$\Delta \log_{10} |B|^{2}$', fontsize=11)

    # Add text with mode info
    fig.suptitle(f'Real-space field reconstruction (n=0 to {args.nmax})', fontsize=13, fontweight='bold')

    print(f"Saving plot to {args.output}")
    fig.savefig(args.output, dpi=150, bbox_inches='tight')
    print(f"Done! Plot saved to {args.output}")


if __name__ == "__main__":
    main()
