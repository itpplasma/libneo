#!/usr/bin/env python3
"""
Compare GPEC Fourier Biot-Savart vs coil_tools vector_potential implementation.

This script reads the reference HDF5 file (from GPEC Fourier) and the test
NetCDF file (from coil_tools vector_potential), computes the magnetic field
from both, and plots a comparison.
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


def compute_field_magnitude(BnR, Bnphi, BnZ):
    """Compute log10 of field magnitude squared."""
    return np.log10((BnR * np.conj(BnR) +
                     Bnphi * np.conj(Bnphi) +
                     BnZ * np.conj(BnZ)).real)


def main():
    parser = argparse.ArgumentParser(
        description="Compare GPEC Fourier vs coil_tools vector_potential Biot-Savart"
    )
    parser.add_argument("reference", type=Path, help="HDF5 file from GPEC Fourier")
    parser.add_argument("test", type=Path, help="NetCDF file from coil_tools vector_potential")
    parser.add_argument("-o", "--output", type=Path, default="comparison.png",
                       help="Output plot filename (default: comparison.png)")
    parser.add_argument("--ntor", type=int, default=2,
                       help="Toroidal mode number (default: 2)")
    parser.add_argument("--resolution-factor", type=float, default=2.0,
                       help="Resolution increase factor for interpolation (default: 2.0)")

    args = parser.parse_args()

    # Validate input files
    if not args.reference.exists():
        print(f"Error: Reference file not found: {args.reference}")
        sys.exit(1)
    if not args.test.exists():
        print(f"Error: Test file not found: {args.test}")
        sys.exit(1)

    print(f"Reading reference (GPEC Fourier): {args.reference}")
    ref_grid, BnR, Bnphi, BnZ = read_Bnvac_fourier(str(args.reference), ntor=args.ntor)

    print(f"Reading test (coil_tools vector_potential): {args.test}")
    test_grid, AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ = read_Anvac_fourier(
        str(args.test), ntor=args.ntor
    )

    # Create splines for both methods
    print("Creating splines for reference...")
    ref_spl = spline_gauged_Anvac(
        ref_grid,
        *gauged_Anvac_from_Bnvac(ref_grid, BnR, BnZ, ntor=args.ntor),
        ntor=args.ntor
    )

    print("Creating splines for test...")
    test_spl = spline_gauged_Anvac(
        test_grid,
        *gauge_Anvac(test_grid, AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ, ntor=args.ntor),
        ntor=args.ntor
    )

    # Verify grids match
    if not np.allclose(ref_grid.R, test_grid.R):
        print("Warning: R grids do not match exactly")
    if not np.allclose(ref_grid.Z, test_grid.Z):
        print("Warning: Z grids do not match exactly")

    # Create higher resolution evaluation grid
    nR_interp = int((2 * ref_grid.nR - 1) * args.resolution_factor / 2)
    nZ_interp = int((2 * ref_grid.nZ - 1) * args.resolution_factor / 2)
    R = np.linspace(ref_grid.R_min, ref_grid.R_max, nR_interp)
    Z = np.linspace(ref_grid.Z_min, ref_grid.Z_max, nZ_interp)

    print(f"Evaluating fields on {len(R)}x{len(Z)} grid...")

    # Compute fields for both methods
    ncoil = len(ref_spl['AnR_Re'])
    log_Bn2 = np.empty((2, ncoil, R.size, Z.size))

    for k, (name, spl) in enumerate([("reference", ref_spl), ("test", test_spl)]):
        print(f"  Computing {name} fields...")
        BnR_eval, Bnphi_eval, BnZ_eval = field_divfree(spl, R, Z, ntor=args.ntor)
        log_Bn2[k] = compute_field_magnitude(BnR_eval, Bnphi_eval, BnZ_eval)

    # Create comparison plot with 4x4 layout (2 columns per coil, multiple rows)
    print(f"Creating comparison plot for {ncoil} coils...")

    ncols = 4  # 2 columns per coil (reference + test)
    nrows = int(np.ceil(ncoil / 2))  # 2 coils per row

    fig = plt.figure(figsize=(12, 3 * nrows), layout='constrained')
    axs = fig.subplots(nrows, ncols)

    # Ensure axs is 2D even if nrows=1
    if nrows == 1:
        axs = axs.reshape(1, -1)

    # Use consistent colorbar range
    vmin = np.amin(log_Bn2)
    vmax = np.amax(log_Bn2)
    norm = Normalize(vmin=vmin, vmax=vmax)

    for kcoil in range(ncoil):
        row = kcoil // 2
        col_offset = (kcoil % 2) * 2

        # Plot reference
        im_ref = axs[row, col_offset].imshow(
            log_Bn2[0, kcoil].T,
            origin='lower',
            cmap='magma',
            extent=[R[0], R[-1], Z[0], Z[-1]],
            norm=norm,
            aspect='auto'
        )
        axs[row, col_offset].set_title(f"Coil {kcoil + 1} (GPEC)", fontsize=10)
        axs[row, col_offset].set_xlabel("R [cm]")
        axs[row, col_offset].set_ylabel("Z [cm]")

        # Plot test
        im_test = axs[row, col_offset + 1].imshow(
            log_Bn2[1, kcoil].T,
            origin='lower',
            cmap='magma',
            extent=[R[0], R[-1], Z[0], Z[-1]],
            norm=norm,
            aspect='auto'
        )
        axs[row, col_offset + 1].set_title(f"Coil {kcoil + 1} (coil_tools)", fontsize=10)
        axs[row, col_offset + 1].set_xlabel("R [cm]")
        axs[row, col_offset + 1].set_ylabel("Z [cm]")

    # Hide unused subplots if ncoil is odd
    if ncoil % 2 == 1:
        axs[-1, -2].axis('off')
        axs[-1, -1].axis('off')

    # Add colorbar
    cbar = fig.colorbar(im_ref, ax=axs, location='bottom', fraction=0.02, pad=0.04)
    cbar.set_label(r'$\log_{10} |\vec{B}_{n}|^{2}$', fontsize=12)

    print(f"Saving plot to {args.output}")
    fig.savefig(args.output, dpi=150, bbox_inches='tight')
    print(f"Done! Plot saved to {args.output}")

    # Compute and report statistics
    print("\nComparison statistics:")
    for kcoil in range(ncoil):
        ref_field = log_Bn2[0, kcoil]
        test_field = log_Bn2[1, kcoil]
        diff = ref_field - test_field
        rms_diff = np.sqrt(np.mean(diff**2))
        max_diff = np.max(np.abs(diff))
        print(f"  Coil {kcoil + 1}: RMS diff = {rms_diff:.6f}, max diff = {max_diff:.6f}")


if __name__ == "__main__":
    main()
