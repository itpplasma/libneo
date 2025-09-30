#!/usr/bin/env python3
"""
Compare three Biot-Savart implementations:
1. GPEC Fourier (direct B_n computation)
2. coil_tools vector_potential (A_n then curl)
3. Direct real-space Biot-Savart (on-the-fly evaluation)

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
    from _magfie import compute_field_from_gpec_file
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


def compute_direct_biotsavart_field(coil_files, coil_currents, R_grid, Z_grid, ntor=2, n_phi=64):
    """Compute n=ntor Fourier mode from direct Biot-Savart by evaluating at multiple phi."""
    nR, nZ = len(R_grid), len(Z_grid)

    # Create phi grid for Fourier decomposition
    phi_grid = np.linspace(0, 2*np.pi, n_phi, endpoint=False)

    # Storage for fields at all phi values
    BR_all = np.zeros((nR, nZ, n_phi))
    Bphi_all = np.zeros((nR, nZ, n_phi))
    BZ_all = np.zeros((nR, nZ, n_phi))

    # Split currents between coil files
    n_coils_per_file = len(coil_currents) // len(coil_files)

    # Evaluate field at each toroidal angle
    for i_phi, phi in enumerate(phi_grid):
        R_mesh, Z_mesh = np.meshgrid(R_grid, Z_grid, indexing='ij')
        x_eval = R_mesh.flatten() * np.cos(phi)
        y_eval = R_mesh.flatten() * np.sin(phi)
        z_eval = Z_mesh.flatten()

        # Evaluate field for each coil file (bu, bl) and sum
        for i, coil_file in enumerate(coil_files):
            currents_for_file = coil_currents[i*n_coils_per_file:(i+1)*n_coils_per_file]

            Bx, By, Bz = compute_field_from_gpec_file(
                str(coil_file), currents_for_file, x_eval, y_eval, z_eval
            )

            # Convert Cartesian to cylindrical
            Bx = Bx.reshape(nR, nZ)
            By = By.reshape(nR, nZ)
            Bz = Bz.reshape(nR, nZ)

            BR_all[:, :, i_phi] += Bx * np.cos(phi) + By * np.sin(phi)
            Bphi_all[:, :, i_phi] += -Bx * np.sin(phi) + By * np.cos(phi)
            BZ_all[:, :, i_phi] += Bz

    # Extract n=ntor Fourier mode using FFT
    # For each (R,Z) point, compute FFT over phi
    # FFT convention: B(phi) = sum_n B_n * exp(i*n*phi)
    # So B_n = (1/N) * sum_phi B(phi) * exp(-i*n*phi)
    BnR = np.fft.fft(BR_all, axis=2)[:, :, ntor] / n_phi
    Bnphi = np.fft.fft(Bphi_all, axis=2)[:, :, ntor] / n_phi
    BnZ = np.fft.fft(BZ_all, axis=2)[:, :, ntor] / n_phi

    return BnR, Bnphi, BnZ


def compute_field_magnitude(BnR, Bnphi, BnZ):
    """Compute log10 of field magnitude squared."""
    return np.log10((BnR * np.conj(BnR) +
                     Bnphi * np.conj(Bnphi) +
                     BnZ * np.conj(BnZ)).real + 1e-20)


def main():
    parser = argparse.ArgumentParser(
        description="Compare three Biot-Savart implementations"
    )
    parser.add_argument("reference", type=Path, help="HDF5 file from GPEC Fourier")
    parser.add_argument("test", type=Path, help="NetCDF file from coil_tools")
    parser.add_argument("currents", type=Path, help="Coil currents file")
    parser.add_argument("coil_upper", type=Path, help="Upper coil file (e.g., aug_bu.dat)")
    parser.add_argument("coil_lower", type=Path, help="Lower coil file (e.g., aug_bl.dat)")
    parser.add_argument("-o", "--output", type=Path, default="superposition_comparison.png",
                       help="Output plot filename (default: superposition_comparison.png)")
    parser.add_argument("--ntor", type=int, default=2,
                       help="Toroidal mode number (default: 2)")

    args = parser.parse_args()

    # Validate input files
    for f in [args.reference, args.test, args.currents, args.coil_upper, args.coil_lower]:
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

    # Sum over coils weighted by currents (Fourier mode n=2)
    BnR_ref_total = np.sum(currents[:, np.newaxis, np.newaxis] * BnR_ref, axis=0)
    Bnphi_ref_total = np.sum(currents[:, np.newaxis, np.newaxis] * Bnphi_ref, axis=0)
    BnZ_ref_total = np.sum(currents[:, np.newaxis, np.newaxis] * BnZ_ref, axis=0)

    BnR_test_total = np.sum(currents[:, np.newaxis, np.newaxis] * BnR_test, axis=0)
    Bnphi_test_total = np.sum(currents[:, np.newaxis, np.newaxis] * Bnphi_test, axis=0)
    BnZ_test_total = np.sum(currents[:, np.newaxis, np.newaxis] * BnZ_test, axis=0)

    # Compute field magnitudes for Fourier methods
    B_ref_mag = np.sqrt((BnR_ref_total * np.conj(BnR_ref_total) +
                         Bnphi_ref_total * np.conj(Bnphi_ref_total) +
                         BnZ_ref_total * np.conj(BnZ_ref_total)).real)
    B_test_mag = np.sqrt((BnR_test_total * np.conj(BnR_test_total) +
                          Bnphi_test_total * np.conj(Bnphi_test_total) +
                          BnZ_test_total * np.conj(BnZ_test_total)).real)

    # Compute direct Biot-Savart field n=2 mode via Fourier decomposition
    print("\nComputing direct Biot-Savart field (on-the-fly with Fourier decomposition)...")
    BnR_direct, Bnphi_direct, BnZ_direct = compute_direct_biotsavart_field(
        [args.coil_upper, args.coil_lower], currents, R_grid, Z_grid, ntor=args.ntor
    )
    B_direct_mag = np.sqrt((BnR_direct * np.conj(BnR_direct) +
                            Bnphi_direct * np.conj(Bnphi_direct) +
                            BnZ_direct * np.conj(BnZ_direct)).real)

    log_B_ref = np.log10(B_ref_mag + 1e-20)
    log_B_test = np.log10(B_test_mag + 1e-20)
    log_B_direct = np.log10(B_direct_mag + 1e-20)

    # Compute pairwise differences
    diff_1_2 = log_B_test - log_B_ref
    diff_1_3 = log_B_direct - log_B_ref
    diff_2_3 = log_B_direct - log_B_test

    # Statistics
    print("\nThree-way comparison statistics (Fourier mode n=2):")
    print(f"  Mean |B| (GPEC Fourier):      {np.mean(B_ref_mag):.6e}  (max: {np.max(B_ref_mag):.2f}, min: {np.min(B_ref_mag):.2f})")
    print(f"  Mean |B| (vector_potential):  {np.mean(B_test_mag):.6e}  (max: {np.max(B_test_mag):.2f}, min: {np.min(B_test_mag):.2f})")
    print(f"  Mean |B| (direct Biot-Savart): {np.mean(B_direct_mag):.6e}  (max: {np.max(B_direct_mag):.2f}, min: {np.min(B_direct_mag):.2f})")

    rel_err_1_2 = np.abs(B_test_mag - B_ref_mag) / (B_ref_mag + 1e-15)
    rel_err_1_3 = np.abs(B_direct_mag - B_ref_mag) / (B_ref_mag + 1e-15)

    print(f"\nGPEC vs vector_potential:")
    print(f"  Median relative error: {np.median(rel_err_1_2)*100:.4f}%")
    print(f"  Mean relative error:   {np.mean(rel_err_1_2)*100:.4f}%")

    print(f"\nGPEC vs direct Biot-Savart:")
    print(f"  Median relative error: {np.median(rel_err_1_3)*100:.4f}%")
    print(f"  Mean relative error:   {np.mean(rel_err_1_3)*100:.4f}%")

    # Create 2x2 comparison plot
    print(f"\nCreating comparison plot...")

    fig = plt.figure(figsize=(12, 10), layout='constrained')
    axs = fig.subplots(2, 2)

    # Shared colorbar range for absolute fields
    vmin_abs = min(np.min(log_B_ref), np.min(log_B_test), np.min(log_B_direct))
    vmax_abs = max(np.max(log_B_ref), np.max(log_B_test), np.max(log_B_direct))
    norm_abs = Normalize(vmin=vmin_abs, vmax=vmax_abs)

    extent = [R_grid[0], R_grid[-1], Z_grid[0], Z_grid[-1]]

    # Plot 1: GPEC Fourier
    im0 = axs[0, 0].imshow(log_B_ref.T, origin='lower', cmap='magma',
                           extent=extent, norm=norm_abs, aspect='auto', interpolation='bilinear')
    axs[0, 0].set_title('1. GPEC Fourier', fontsize=11, fontweight='bold')
    axs[0, 0].set_xlabel('R [cm]')
    axs[0, 0].set_ylabel('Z [cm]')

    # Plot 2: vector_potential
    im1 = axs[0, 1].imshow(log_B_test.T, origin='lower', cmap='magma',
                           extent=extent, norm=norm_abs, aspect='auto', interpolation='bilinear')
    axs[0, 1].set_title('2. vector_potential', fontsize=11, fontweight='bold')
    axs[0, 1].set_xlabel('R [cm]')
    axs[0, 1].set_ylabel('Z [cm]')

    # Plot 3: Direct Biot-Savart
    im2 = axs[1, 0].imshow(log_B_direct.T, origin='lower', cmap='magma',
                           extent=extent, norm=norm_abs, aspect='auto', interpolation='bilinear')
    axs[1, 0].set_title('3. Direct Biot-Savart', fontsize=11, fontweight='bold')
    axs[1, 0].set_xlabel('R [cm]')
    axs[1, 0].set_ylabel('Z [cm]')

    # Plot 4: Max difference
    diff_max = max(abs(np.min(diff_1_2)), abs(np.max(diff_1_2)),
                   abs(np.min(diff_1_3)), abs(np.max(diff_1_3)))
    im3 = axs[1, 1].imshow(diff_1_2.T, origin='lower', cmap='RdBu_r',
                           extent=extent, vmin=-diff_max, vmax=diff_max, aspect='auto', interpolation='bilinear')
    axs[1, 1].set_title('Difference: 2 - 1', fontsize=11, fontweight='bold')
    axs[1, 1].set_xlabel('R [cm]')
    axs[1, 1].set_ylabel('Z [cm]')

    # Colorbars
    cbar_abs = fig.colorbar(im0, ax=axs[:, :].ravel().tolist()[:3], location='bottom', fraction=0.05, pad=0.08)
    cbar_abs.set_label(r'$\log_{10} |\vec{B}_{n=2}|$', fontsize=10)

    cbar_diff = fig.colorbar(im3, ax=axs[1, 1], location='bottom', fraction=0.05, pad=0.08)
    cbar_diff.set_label(r'$\Delta \log_{10} |B|$', fontsize=10)

    print(f"Saving plot to {args.output}")
    fig.savefig(args.output, dpi=150, bbox_inches='tight')
    print(f"Done! Plot saved to {args.output}")


if __name__ == "__main__":
    main()
