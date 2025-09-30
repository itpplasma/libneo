#!/usr/bin/env python3
"""
Compare total superposed magnetic field from GPEC Fourier vs coil_tools.

This script computes the full superposition of all coils weighted by their
currents, which is the physically relevant comparison.
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
        description="Compare total superposed field from GPEC Fourier vs coil_tools"
    )
    parser.add_argument("reference", type=Path, help="HDF5 file from GPEC Fourier")
    parser.add_argument("test", type=Path, help="NetCDF file from coil_tools")
    parser.add_argument("currents", type=Path, help="Coil currents file")
    parser.add_argument("-o", "--output", type=Path, default="superposition_comparison.png",
                       help="Output plot filename (default: superposition_comparison.png)")
    parser.add_argument("--ntor", type=int, default=2,
                       help="Toroidal mode number (default: 2)")

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

    # Read field data
    print(f"\nReading reference (GPEC Fourier): {args.reference}")
    ref_grid, BnR, Bnphi, BnZ = read_Bnvac_fourier(str(args.reference), ntor=args.ntor)

    print(f"Reading test (coil_tools vector_potential): {args.test}")
    test_grid, AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ = read_Anvac_fourier(
        str(args.test), ntor=args.ntor
    )

    ncoil = len(currents)
    if BnR.shape[0] != ncoil:
        print(f"Error: Number of coils ({BnR.shape[0]}) != number of currents ({ncoil})")
        sys.exit(1)

    # Use raw Fourier data directly (no spline interpolation)
    print("\nUsing raw Fourier data on native grid...")

    # Reference has B directly
    BnR_ref = BnR
    # Compute Bnphi from div B = 0: ∂(R B_R)/∂R + ∂(R B_Z)/∂Z + i n B_phi = 0
    Bnphi_ref = np.zeros_like(BnR)
    BnZ_ref = BnZ

    # Test: compute B from curl of A using stored derivatives
    # B_R = (1/R) ∂A_Z/∂φ - ∂A_φ/∂Z = (i n / R) A_Z - ∂A_φ/∂Z
    # B_φ = ∂A_R/∂Z - ∂A_Z/∂R
    # B_Z = (1/R)[∂(R A_φ)/∂R - ∂A_R/∂φ] = (1/R)[A_φ + R ∂A_φ/∂R - i n A_R]
    R_grid = np.linspace(test_grid.R_min, test_grid.R_max, test_grid.nR)
    Z_grid = np.linspace(test_grid.Z_min, test_grid.Z_max, test_grid.nZ)
    R_mesh = R_grid[np.newaxis, :, np.newaxis]  # (1, nR, 1) for shape (ncoil, nR, nZ)

    BnR_test = (1j * args.ntor / R_mesh) * AnZ - dAnphi_dZ
    # For B_phi, need ∂A_R/∂Z and ∂A_Z/∂R - compute with finite differences
    dAnR_dZ = np.zeros_like(AnR)
    dAnZ_dR = np.zeros_like(AnZ)
    dZ = (test_grid.Z_max - test_grid.Z_min) / (test_grid.nZ - 1)
    dR = (test_grid.R_max - test_grid.R_min) / (test_grid.nR - 1)
    dAnR_dZ[:, :, 1:-1] = (AnR[:, :, 2:] - AnR[:, :, :-2]) / (2 * dZ)
    dAnZ_dR[:, 1:-1, :] = (AnZ[:, 2:, :] - AnZ[:, :-2, :]) / (2 * dR)
    Bnphi_test = dAnR_dZ - dAnZ_dR
    BnZ_test = (Anphi + R_mesh * dAnphi_dR - 1j * args.ntor * AnR) / R_mesh

    # Compute superposition weighted by currents
    print("\nComputing superposition with coil currents...")

    # Sum over coils weighted by currents
    BnR_ref_total = np.sum(currents[:, np.newaxis, np.newaxis] * BnR_ref, axis=0)
    Bnphi_ref_total = np.sum(currents[:, np.newaxis, np.newaxis] * Bnphi_ref, axis=0)
    BnZ_ref_total = np.sum(currents[:, np.newaxis, np.newaxis] * BnZ_ref, axis=0)

    BnR_test_total = np.sum(currents[:, np.newaxis, np.newaxis] * BnR_test, axis=0)
    Bnphi_test_total = np.sum(currents[:, np.newaxis, np.newaxis] * Bnphi_test, axis=0)
    BnZ_test_total = np.sum(currents[:, np.newaxis, np.newaxis] * BnZ_test, axis=0)

    # Compute field magnitudes
    B_ref_mag = np.sqrt((BnR_ref_total * np.conj(BnR_ref_total) +
                         Bnphi_ref_total * np.conj(Bnphi_ref_total) +
                         BnZ_ref_total * np.conj(BnZ_ref_total)).real)
    B_test_mag = np.sqrt((BnR_test_total * np.conj(BnR_test_total) +
                          Bnphi_test_total * np.conj(Bnphi_test_total) +
                          BnZ_test_total * np.conj(BnZ_test_total)).real)

    log_B_ref = np.log10(B_ref_mag + 1e-20)
    log_B_test = np.log10(B_test_mag + 1e-20)

    # Compute difference
    diff = log_B_test - log_B_ref

    # Statistics
    print("\nSuperposition comparison statistics:")

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

    fig = plt.figure(figsize=(15, 4), layout='constrained')
    axs = fig.subplots(1, 3)

    # Shared colorbar range for absolute fields
    vmin_abs = min(np.min(log_B_ref), np.min(log_B_test))
    vmax_abs = max(np.max(log_B_ref), np.max(log_B_test))
    norm_abs = Normalize(vmin=vmin_abs, vmax=vmax_abs)

    extent = [R_grid[0], R_grid[-1], Z_grid[0], Z_grid[-1]]

    # Plot GPEC reference
    # field is (nR, nZ), need .T so R is horizontal and Z is vertical in imshow
    im0 = axs[0].imshow(log_B_ref.T, origin='lower', cmap='magma',
                        extent=extent, norm=norm_abs, aspect='auto', interpolation='bilinear')
    axs[0].set_title('GPEC Fourier (Reference)', fontsize=12, fontweight='bold')
    axs[0].set_xlabel('R [cm]')
    axs[0].set_ylabel('Z [cm]')

    # Plot coil_tools test
    im1 = axs[1].imshow(log_B_test.T, origin='lower', cmap='magma',
                        extent=extent, norm=norm_abs, aspect='auto', interpolation='bilinear')
    axs[1].set_title('coil_tools vector_potential (Test)', fontsize=12, fontweight='bold')
    axs[1].set_xlabel('R [cm]')
    axs[1].set_ylabel('Z [cm]')

    # Plot difference
    diff_max = max(abs(np.min(diff)), abs(np.max(diff)))
    im2 = axs[2].imshow(diff.T, origin='lower', cmap='RdBu_r',
                        extent=extent, vmin=-diff_max, vmax=diff_max, aspect='auto', interpolation='bilinear')
    axs[2].set_title('Difference (Test - Reference)', fontsize=12, fontweight='bold')
    axs[2].set_xlabel('R [cm]')
    axs[2].set_ylabel('Z [cm]')

    # Colorbars
    cbar0 = fig.colorbar(im0, ax=axs[:2], location='bottom', fraction=0.05, pad=0.08)
    cbar0.set_label(r'$\log_{10} |\vec{B}_{n=2}|$', fontsize=11)

    cbar2 = fig.colorbar(im2, ax=axs[2], location='bottom', fraction=0.05, pad=0.08)
    cbar2.set_label(r'$\Delta \log_{10} |B|$', fontsize=11)

    print(f"Saving plot to {args.output}")
    fig.savefig(args.output, dpi=150, bbox_inches='tight')
    print(f"Done! Plot saved to {args.output}")


if __name__ == "__main__":
    main()
