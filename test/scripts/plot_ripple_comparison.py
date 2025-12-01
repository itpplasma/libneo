#!/usr/bin/env python3
"""
Plot 2D comparison of Bphi with and without TF ripple perturbation.

This script reads CSV data from the Fortran test and creates pcolormesh plots
showing the toroidal field component on a flux surface in (phi, Z) coordinates.
"""
import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def load_bphi_data(csv_path):
    """Load Bphi data from CSV file."""
    data = np.loadtxt(csv_path, delimiter=',', skiprows=1)
    phi = data[:, 0]
    Z = data[:, 1]
    Bphi = data[:, 2]
    return phi, Z, Bphi


def reshape_to_grid(phi, Z, Bphi):
    """Reshape 1D arrays to 2D grid for pcolormesh."""
    phi_unique = np.unique(phi)
    Z_unique = np.unique(Z)
    nphi = len(phi_unique)
    nZ = len(Z_unique)

    phi_grid = phi.reshape(nZ, nphi)
    Z_grid = Z.reshape(nZ, nphi)
    Bphi_grid = Bphi.reshape(nZ, nphi)

    return phi_grid, Z_grid, Bphi_grid


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', type=Path, default=Path.cwd(),
                        help='Directory containing CSV files')
    parser.add_argument('--output-dir', type=Path, default=None,
                        help='Directory for output plots (default: data-dir)')
    args = parser.parse_args()

    data_dir = args.data_dir.resolve()
    output_dir = (args.output_dir or data_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    csv_no_ripple = data_dir / 'bphi_no_ripple.csv'
    csv_ripple = data_dir / 'bphi_with_ripple.csv'

    if not csv_no_ripple.exists():
        print(f"Error: {csv_no_ripple} not found", file=sys.stderr)
        return 1

    if not csv_ripple.exists():
        print(f"Error: {csv_ripple} not found", file=sys.stderr)
        return 1

    print(f"Loading data from {data_dir}")
    phi_no, Z_no, Bphi_no = load_bphi_data(csv_no_ripple)
    phi_rip, Z_rip, Bphi_rip = load_bphi_data(csv_ripple)

    phi_grid_no, Z_grid_no, Bphi_grid_no = reshape_to_grid(phi_no, Z_no, Bphi_no)
    phi_grid_rip, Z_grid_rip, Bphi_grid_rip = reshape_to_grid(phi_rip, Z_rip, Bphi_rip)

    phi_deg_no = np.rad2deg(phi_grid_no)
    phi_deg_rip = np.rad2deg(phi_grid_rip)

    vmin = min(Bphi_grid_no.min(), Bphi_grid_rip.min())
    vmax = max(Bphi_grid_no.max(), Bphi_grid_rip.max())

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    im0 = axes[0].pcolormesh(phi_deg_no, Z_grid_no, Bphi_grid_no,
                             shading='auto', vmin=vmin, vmax=vmax, cmap='viridis')
    axes[0].set_xlabel('Toroidal angle φ [deg]', fontsize=12)
    axes[0].set_ylabel('Z [m]', fontsize=12)
    axes[0].set_title('B_φ without ripple', fontsize=14)
    axes[0].set_aspect('auto')
    plt.colorbar(im0, ax=axes[0], label='B_φ [T]')

    im1 = axes[1].pcolormesh(phi_deg_rip, Z_grid_rip, Bphi_grid_rip,
                             shading='auto', vmin=vmin, vmax=vmax, cmap='viridis')
    axes[1].set_xlabel('Toroidal angle φ [deg]', fontsize=12)
    axes[1].set_ylabel('Z [m]', fontsize=12)
    axes[1].set_title('B_φ with 9-coil ripple (δ₀=0.10)', fontsize=14)
    axes[1].set_aspect('auto')
    plt.colorbar(im1, ax=axes[1], label='B_φ [T]')

    delta_Bphi = Bphi_grid_rip - Bphi_grid_no
    delta_max = np.abs(delta_Bphi).max()

    im2 = axes[2].pcolormesh(phi_deg_rip, Z_grid_rip, delta_Bphi,
                             shading='auto', vmin=-delta_max, vmax=delta_max,
                             cmap='RdBu_r')
    axes[2].set_xlabel('Toroidal angle φ [deg]', fontsize=12)
    axes[2].set_ylabel('Z [m]', fontsize=12)
    axes[2].set_title('Ripple perturbation ΔB_φ', fontsize=14)
    axes[2].set_aspect('auto')
    plt.colorbar(im2, ax=axes[2], label='ΔB_φ [T]')

    fig.tight_layout()

    output_png = output_dir / 'bphi_ripple_comparison.png'
    fig.savefig(output_png, dpi=150, bbox_inches='tight')
    print(f"Saved plot to {output_png}")

    ripple_fraction = (delta_Bphi.max() - delta_Bphi.min()) / Bphi_grid_no.mean()
    print(f"\nRipple statistics:")
    print(f"  Mean B_φ (no ripple): {Bphi_grid_no.mean():.4f} T")
    print(f"  ΔB_φ range: [{delta_Bphi.min():.6f}, {delta_Bphi.max():.6f}] T")
    print(f"  Peak-to-peak ripple: {ripple_fraction*100:.2f}%")
    print(f"  Expected for 9 coils: ~{360/9:.1f}° periodicity")

    plt.close(fig)

    return 0


if __name__ == '__main__':
    sys.exit(main())
