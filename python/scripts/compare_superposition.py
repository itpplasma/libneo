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


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[2]
PYTHON_DIR = REPO_ROOT / "python"
BUILD_DIR = REPO_ROOT / "build"
VENV_SITE_PACKAGES = None

def find_site_packages(base: Path):
    lib_dir = base / "lib"
    if not lib_dir.exists():
        return None
    for sub in lib_dir.iterdir():
        if not sub.is_dir():
            continue
        candidate = sub / "site-packages"
        if candidate.exists():
            return candidate
    return None

for possible in (REPO_ROOT / ".venv", REPO_ROOT.parent / ".venv"):
    if possible.exists():
        site_packages = find_site_packages(possible)
        if site_packages is not None:
            VENV_SITE_PACKAGES = site_packages
            break

for candidate in (PYTHON_DIR, BUILD_DIR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

if VENV_SITE_PACKAGES is not None and str(VENV_SITE_PACKAGES) not in sys.path:
    sys.path.insert(0, str(VENV_SITE_PACKAGES))

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
except ImportError as exc:
    print("Error: required Python modules for libneo are missing:", exc)
    print("Ensure the libneo package and runtime dependencies (h5py, netCDF4) are installed.")
    sys.exit(1)


def read_currents(currents_file):
    """Read coil currents from file."""
    with open(currents_file, 'r') as f:
        line = f.readline().strip()
        currents = np.array([float(x) for x in line.split()])
    return currents


def read_gpec_header(filename: Path):
    """Return (ncoil, nseg, nwind) from the header of a GPEC coil file."""
    with open(filename, 'r') as f:
        header = f.readline().split()
    if len(header) < 4:
        raise ValueError(f"Malformed GPEC header in {filename}")
    return int(header[0]), int(header[2]), float(header[3])


def compute_direct_biotsavart_field(coil_files, coil_currents, coils_per_file,
                                    R_grid, Z_grid, ntor=2, n_phi=128):
    """Compute n=ntor Fourier mode from direct Biot-Savart by evaluating at multiple phi."""
    nR, nZ = len(R_grid), len(Z_grid)

    # Create phi grid for Fourier decomposition
    phi_grid = np.linspace(0, 2*np.pi, n_phi, endpoint=False)

    # Storage for fields at all phi values
    BR_all = np.zeros((nR, nZ, n_phi))
    Bphi_all = np.zeros((nR, nZ, n_phi))
    BZ_all = np.zeros((nR, nZ, n_phi))

    # Evaluate field at each toroidal angle
    for i_phi, phi in enumerate(phi_grid):
        R_mesh, Z_mesh = np.meshgrid(R_grid, Z_grid, indexing='ij')
        x_eval = R_mesh.flatten() * np.cos(phi)
        y_eval = R_mesh.flatten() * np.sin(phi)
        z_eval = Z_mesh.flatten()

        # Evaluate field for each coil file (bu, bl) and sum
        offset = 0
        for i, coil_file in enumerate(coil_files):
            n_coils = coils_per_file[i]
            currents_for_file = coil_currents[offset:offset + n_coils]
            offset += n_coils

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


def bilinear_interpolate(R_grid, Z_grid, field, R_val, Z_val):
    """Simple bilinear interpolation on an (nR, nZ) grid."""
    R_min, R_max = R_grid[0], R_grid[-1]
    Z_min, Z_max = Z_grid[0], Z_grid[-1]

    if (R_val < R_min) or (R_val > R_max) or (Z_val < Z_min) or (Z_val > Z_max):
        return np.nan

    i_hi = np.searchsorted(R_grid, R_val)
    if i_hi == 0:
        i0, i1, t = 0, 0, 0.0
    elif i_hi >= R_grid.size:
        i0, i1, t = R_grid.size - 1, R_grid.size - 1, 0.0
    else:
        i0, i1 = i_hi - 1, i_hi
        denom = R_grid[i1] - R_grid[i0]
        t = 0.0 if denom == 0 else (R_val - R_grid[i0]) / denom

    j_hi = np.searchsorted(Z_grid, Z_val)
    if j_hi == 0:
        j0, j1, u = 0, 0, 0.0
    elif j_hi >= Z_grid.size:
        j0, j1, u = Z_grid.size - 1, Z_grid.size - 1, 0.0
    else:
        j0, j1 = j_hi - 1, j_hi
        denom = Z_grid[j1] - Z_grid[j0]
        u = 0.0 if denom == 0 else (Z_val - Z_grid[j0]) / denom

    f00 = field[i0, j0]
    f01 = field[i0, j1]
    f10 = field[i1, j0]
    f11 = field[i1, j1]

    return ((1 - t) * (1 - u) * f00 +
            (1 - t) * u * f01 +
            t * (1 - u) * f10 +
            t * u * f11)


def validate_axis_response(args, coil_files, coil_currents, coils_per_file,
                           R_grid, Z_grid,
                           BnR_ref_total, Bnphi_ref_total, BnZ_ref_total,
                           BnR_test_total, Bnphi_test_total, BnZ_test_total):
    if args.axis_origin is None or args.axis_normal is None or args.coil_radius is None:
        return None

    axis_origin_cm = np.array(args.axis_origin, dtype=float)
    axis_normal = np.array(args.axis_normal, dtype=float)
    norm = np.linalg.norm(axis_normal)
    if norm == 0.0:
        raise ValueError("Axis normal vector must be non-zero")
    axis_normal /= norm

    s_vals = np.linspace(-args.axis_range, args.axis_range, args.axis_samples)
    sample_points_cm = axis_origin_cm[np.newaxis, :] + np.outer(s_vals, axis_normal)

    x_eval = sample_points_cm[:, 0]
    y_eval = sample_points_cm[:, 1]
    z_eval = sample_points_cm[:, 2]

    Bx_direct = np.zeros_like(x_eval)
    By_direct = np.zeros_like(y_eval)
    Bz_direct = np.zeros_like(z_eval)

    offset = 0
    for i, coil_file in enumerate(coil_files):
        n_coils = coils_per_file[i]
        currents_slice = coil_currents[offset:offset + n_coils]
        offset += n_coils
        Bx, By, Bz = compute_field_from_gpec_file(
            str(coil_file), currents_slice, x_eval, y_eval, z_eval
        )
        Bx_direct += Bx
        By_direct += By
        Bz_direct += Bz

    B_parallel_direct = Bx_direct * axis_normal[0] + By_direct * axis_normal[1] + Bz_direct * axis_normal[2]

    # Fourier-based evaluations (bilinear interpolation on cylindrical grid)
    B_parallel_gpec = np.zeros_like(x_eval)
    B_parallel_vec = np.zeros_like(x_eval)

    for idx, (x_val, y_val, z_val) in enumerate(zip(x_eval, y_eval, z_eval)):
        R_val = np.hypot(x_val, y_val)
        phi_val = np.arctan2(y_val, x_val)

        BR_ref = bilinear_interpolate(R_grid, Z_grid, BnR_ref_total, R_val, z_val)
        Bphi_ref = bilinear_interpolate(R_grid, Z_grid, Bnphi_ref_total, R_val, z_val)
        BZ_ref = bilinear_interpolate(R_grid, Z_grid, BnZ_ref_total, R_val, z_val)

        BR_vec = bilinear_interpolate(R_grid, Z_grid, BnR_test_total, R_val, z_val)
        Bphi_vec = bilinear_interpolate(R_grid, Z_grid, Bnphi_test_total, R_val, z_val)
        BZ_vec = bilinear_interpolate(R_grid, Z_grid, BnZ_test_total, R_val, z_val)

        cosphi = np.cos(phi_val)
        sinphi = np.sin(phi_val)

        if not (np.isnan(BR_ref) or np.isnan(Bphi_ref) or np.isnan(BZ_ref)):
            Bx_ref = (BR_ref * cosphi - Bphi_ref * sinphi).real
            By_ref = (BR_ref * sinphi + Bphi_ref * cosphi).real
            Bz_ref = BZ_ref.real
            B_parallel_gpec[idx] = (
                Bx_ref * axis_normal[0] + By_ref * axis_normal[1] + Bz_ref * axis_normal[2]
            )
        else:
            B_parallel_gpec[idx] = np.nan

        if not (np.isnan(BR_vec) or np.isnan(Bphi_vec) or np.isnan(BZ_vec)):
            Bx_vec = (BR_vec * cosphi - Bphi_vec * sinphi).real
            By_vec = (BR_vec * sinphi + Bphi_vec * cosphi).real
            Bz_vec = BZ_vec.real
            B_parallel_vec[idx] = (
                Bx_vec * axis_normal[0] + By_vec * axis_normal[1] + Bz_vec * axis_normal[2]
            )
        else:
            B_parallel_vec[idx] = np.nan

    # Analytical reference (SI -> Gauss)
    mu0 = 4e-7 * np.pi
    radius_m = args.coil_radius / 100.0
    s_vals_m = s_vals / 100.0
    current_ampere = np.sum(coil_currents)
    B_analytic = mu0 * current_ampere * radius_m**2 / (2.0 * (radius_m**2 + s_vals_m**2)**1.5)
    B_analytic *= 1.0e4  # Tesla -> Gauss

    abs_err_direct = np.abs(B_parallel_direct - B_analytic)
    rel_err_direct = abs_err_direct / np.maximum(np.abs(B_analytic), 1e-20)

    abs_err_gpec = np.abs(B_parallel_gpec - B_analytic)
    rel_err_gpec = abs_err_gpec / np.maximum(np.abs(B_analytic), 1e-20)

    abs_err_vec = np.abs(B_parallel_vec - B_analytic)
    rel_err_vec = abs_err_vec / np.maximum(np.abs(B_analytic), 1e-20)

    return {
        's_vals_cm': s_vals,
        'analytic': B_analytic,
        'direct': B_parallel_direct,
        'gpec': B_parallel_gpec,
        'vector': B_parallel_vec,
        'abs_err_direct': abs_err_direct,
        'rel_err_direct': rel_err_direct,
        'abs_err_gpec': abs_err_gpec,
        'rel_err_gpec': rel_err_gpec,
        'abs_err_vector': abs_err_vec,
        'rel_err_vector': rel_err_vec,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Compare three Biot-Savart implementations"
    )
    parser.add_argument("reference", type=Path, help="HDF5 file from GPEC Fourier")
    parser.add_argument("test", type=Path, help="NetCDF file from coil_tools")
    parser.add_argument("currents", type=Path, help="Coil currents file")
    parser.add_argument("coil_files", type=Path, nargs='+',
                        help="One or more coil geometry files in GPEC format")
    parser.add_argument("-o", "--output", type=Path, default="superposition_comparison.png",
                        help="Output plot filename (default: superposition_comparison.png)")
    parser.add_argument("--ntor", type=int, default=2,
                        help="Toroidal mode number (default: 2)")
    parser.add_argument("--axis-origin", type=float, nargs=3,
                        help="Axis origin (cm) for analytic validation, e.g., x y z")
    parser.add_argument("--axis-normal", type=float, nargs=3,
                        help="Axis unit normal vector (dimensionless)")
    parser.add_argument("--coil-radius", type=float,
                        help="Coil radius in cm for analytic validation")
    parser.add_argument("--axis-range", type=float, default=100.0,
                        help="Half-length of axis line segment for validation (cm, default 100)")
    parser.add_argument("--axis-samples", type=int, default=101,
                        help="Number of samples along axis (default 101)")

    args = parser.parse_args()

    # Validate input files
    required = [args.reference, args.test, args.currents]
    required.extend(args.coil_files)
    for f in required:
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
    coils_per_file = [read_gpec_header(path)[0] for path in args.coil_files]
    if sum(coils_per_file) != ncoil:
        print("Error: Number of currents does not match total coils defined in geometry files.")
        sys.exit(1)
    if BnR.shape[0] != ncoil:
        print(f"Error: Number of coils ({BnR.shape[0]}) in reference file != currents ({ncoil})")
        sys.exit(1)

    # Use raw Fourier data directly (no spline interpolation)
    print("\nUsing raw Fourier data on native grid...")

    # Reference has B directly
    BnR_ref = BnR
    Bnphi_ref = Bnphi
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
        args.coil_files, currents, coils_per_file, R_grid, Z_grid, ntor=args.ntor
    )
    B_direct_mag = np.sqrt((BnR_direct * np.conj(BnR_direct) +
                            Bnphi_direct * np.conj(Bnphi_direct) +
                            BnZ_direct * np.conj(BnZ_direct)).real)

    axis_result = validate_axis_response(
        args,
        args.coil_files,
        currents,
        coils_per_file,
        R_grid,
        Z_grid,
        BnR_ref_total,
        Bnphi_ref_total,
        BnZ_ref_total,
        BnR_test_total,
        Bnphi_test_total,
        BnZ_test_total,
    )
    if axis_result is not None:
        print("\nAxis validation against analytical solution:")

        def _nanmax(arr):
            finite = np.isfinite(arr)
            if not np.any(finite):
                return float('nan')
            return np.nanmax(arr[finite])

        print(f"  Max |B_GPEC - B_analytic|: {_nanmax(axis_result['abs_err_gpec']):.3e} G")
        print(f"  Max rel error GPEC: {_nanmax(axis_result['rel_err_gpec']) * 100.0:.4f}%")
        print(f"  Max |B_vector - B_analytic|: {_nanmax(axis_result['abs_err_vector']):.3e} G")
        print(f"  Max rel error vector: {_nanmax(axis_result['rel_err_vector']) * 100.0:.4f}%")
        print(f"  Max |B_direct - B_analytic|: {_nanmax(axis_result['abs_err_direct']):.3e} G")
        print(f"  Max rel error direct: {_nanmax(axis_result['rel_err_direct']) * 100.0:.4f}%")

        axis_fig = plt.figure(figsize=(6, 4), layout='constrained')
        ax = axis_fig.add_subplot(111)
        ax.plot(axis_result['s_vals_cm'], axis_result['analytic'], label='Analytical', linestyle='--', linewidth=2)
        ax.plot(axis_result['s_vals_cm'], axis_result['gpec'], label='GPEC Fourier', linewidth=1.4)
        ax.plot(axis_result['s_vals_cm'], axis_result['vector'], label='vector_potential', linewidth=1.4)
        ax.plot(axis_result['s_vals_cm'], axis_result['direct'], label='Direct Biot-Savart', linewidth=1.4)
        ax.set_xlabel('Axis coordinate s [cm]')
        ax.set_ylabel('B_parallel [Gauss]')
        ax.set_title('Coil-axis validation')
        ax.legend(loc='best')
        axis_output = args.output.with_name(args.output.stem + '_axis.png')
        axis_fig.savefig(axis_output, dpi=150, bbox_inches='tight')
        print(f"  Axis comparison plot saved to {axis_output}")

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
