#!/usr/bin/env python3
"""Validate ntor=0 reconstruction against analytic and direct Biot-Savart fields."""

from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # Non-interactive backend
import matplotlib.pyplot as plt
import netCDF4
import numpy as np
from scipy.interpolate import RectBivariateSpline, RegularGridInterpolator

SCRIPT_DIR = Path(__file__).resolve().parent
BUILD_DIR = Path.cwd()
COIL_FILE = BUILD_DIR / "tilted_coil.dat"
FOURIER_FILE = BUILD_DIR / "tilted_coil_Anvac.nc"
BNVAC_FILE = BUILD_DIR / "tilted_coil_Bnvac.h5"
BVAC_FILE = BUILD_DIR / "tilted_coil_Bvac.dat"
AXIS_SOLVER = BUILD_DIR / "compute_axis_biot_savart.x"
CLIGHT = 2.99792458e10  # cm / s

sys.path.append(str(SCRIPT_DIR))
from generate_tilted_coil import _plane_basis  # noqa: E402


def circular_coil_analytical(
    s_vals: np.ndarray,
    center_xyz: np.ndarray,
    axis_hat: np.ndarray,
    coil_radius: float,
    current: float,
):
    """Analytical circular coil field along its axis (CGS units)."""
    axis_hat = axis_hat / np.linalg.norm(axis_hat)
    pts = center_xyz[:, None] + axis_hat[:, None] * s_vals[None, :]
    R = np.sqrt(pts[0, :] ** 2 + pts[1, :] ** 2)
    phi = np.mod(np.arctan2(pts[1, :], pts[0, :]), 2.0 * np.pi)

    scale = 100.0  # meters -> centimeters
    a_cm = coil_radius * scale
    s_cm = s_vals * scale
    B_mag = (2.0 * np.pi * current * a_cm**2) / (a_cm**2 + s_cm**2) ** 1.5
    B_cart = axis_hat[:, None] * B_mag[None, :]

    cosphi = np.cos(phi)
    sinphi = np.sin(phi)
    BR = B_cart[0, :] * cosphi + B_cart[1, :] * sinphi
    Bphi = -B_cart[0, :] * sinphi + B_cart[1, :] * cosphi
    BZ = B_cart[2, :]
    return BR, Bphi, BZ, pts, phi


def evaluate_mode_at_points(spl: dict, R_points: np.ndarray, Z_points: np.ndarray, ntor: int):
    """Evaluate Fourier mode at specific (R, Z) points for all coils."""
    from numpy import sum as np_sum

    R_points = np.asarray(R_points, dtype=float)
    Z_points = np.asarray(Z_points, dtype=float)

    npts = R_points.size
    ncoil = len(spl["AnR_Re"])
    BnR = np.empty((ncoil, npts), dtype=complex)
    Bnphi = np.empty_like(BnR)
    BnZ = np.empty_like(BnR)
    has_Anphi = "Anphi_Re" in spl

    for kcoil in range(ncoil):
        AnR = spl["AnR_Re"][kcoil].ev(R_points, Z_points) + 1j * spl["AnR_Im"][kcoil].ev(R_points, Z_points)
        AnZ = spl["AnZ_Re"][kcoil].ev(R_points, Z_points) + 1j * spl["AnZ_Im"][kcoil].ev(R_points, Z_points)
        dAnR_dZ = spl["AnR_Re"][kcoil].ev(R_points, Z_points, dy=1) + 1j * spl["AnR_Im"][kcoil].ev(
            R_points, Z_points, dy=1
        )
        dAnZ_dR = spl["AnZ_Re"][kcoil].ev(R_points, Z_points, dx=1) + 1j * spl["AnZ_Im"][kcoil].ev(
            R_points, Z_points, dx=1
        )

        if has_Anphi:
            Anphi = spl["Anphi_Re"][kcoil].ev(R_points, Z_points) + 1j * spl["Anphi_Im"][kcoil].ev(
                R_points, Z_points
            )
            dAnphi_dR = spl["Anphi_Re"][kcoil].ev(R_points, Z_points, dx=1) + 1j * spl["Anphi_Im"][kcoil].ev(
                R_points, Z_points, dx=1
            )
            dAnphi_dZ = spl["Anphi_Re"][kcoil].ev(R_points, Z_points, dy=1) + 1j * spl["Anphi_Im"][kcoil].ev(
                R_points, Z_points, dy=1
            )

            BnR[kcoil, :] = 1j * ntor * AnZ / R_points - dAnphi_dZ
            Bnphi[kcoil, :] = dAnR_dZ - dAnZ_dR
            BnZ[kcoil, :] = dAnphi_dR + Anphi / R_points - 1j * ntor * AnR / R_points
        else:
            if ntor == 0:
                raise ValueError("For ntor=0, Anphi splines must be provided")
            BnR[kcoil, :] = 1j * ntor * AnZ / R_points
            Bnphi[kcoil, :] = dAnR_dZ - dAnZ_dR
            BnZ[kcoil, :] = -1j * ntor * AnR / R_points

    return np_sum(BnR, axis=0), np_sum(Bnphi, axis=0), np_sum(BnZ, axis=0)


def compute_direct_biot_savart(axis_points: np.ndarray) -> np.ndarray:
    """Call the Fortran biotsavart_field solver for direct B along axis."""
    if not AXIS_SOLVER.exists():
        raise FileNotFoundError(f"Missing executable {AXIS_SOLVER}")
    if not COIL_FILE.exists():
        raise FileNotFoundError(f"Missing coil file {COIL_FILE}")

    scale = 100.0  # meters to centimeters

    coil_data = np.loadtxt(COIL_FILE, skiprows=1)
    if coil_data.ndim == 1:
        coil_data = coil_data[np.newaxis, :]
    scaled_coil = coil_data.copy()
    scaled_coil[:, 0:3] *= scale
    scaled_coil[:, 3] *= CLIGHT  # convert abampere to statampere

    with tempfile.NamedTemporaryFile("w", dir=BUILD_DIR, suffix=".coil", delete=False) as coil_f:
        coil_path = Path(coil_f.name)
        coil_f.write(f"{scaled_coil.shape[0]}\n")
        for row in scaled_coil:
            coil_f.write(f"{row[0]:.18e} {row[1]:.18e} {row[2]:.18e} {row[3]:.18e}\n")

    with tempfile.NamedTemporaryFile("w", dir=BUILD_DIR, suffix=".axis", delete=False) as axis_f:
        axis_path = Path(axis_f.name)
        axis_f.write(f"{axis_points.shape[1]}\n")
        for col in axis_points.T:
            axis_f.write(
                f"{(col[0]*scale):.18e} {(col[1]*scale):.18e} {(col[2]*scale):.18e}\n"
            )

    with tempfile.NamedTemporaryFile("w", dir=BUILD_DIR, suffix=".field", delete=False) as out_f:
        out_path = Path(out_f.name)

    try:
        subprocess.run(
            [str(AXIS_SOLVER), str(coil_path), str(axis_path), str(out_path)],
            cwd=BUILD_DIR,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        data = np.loadtxt(out_path, dtype=float)
    finally:
        coil_path.unlink(missing_ok=True)
        axis_path.unlink(missing_ok=True)
        out_path.unlink(missing_ok=True)

    if data.ndim == 1:
        data = data[np.newaxis, :]
    return data


def read_Bvac_nemov(path: Path):
    with open(path, "r", encoding="utf-8") as f:
        header = f.readline().split()
        if len(header) < 4:
            raise ValueError("Malformed Bvac header")
        nR, nphi_plus1, nZ, _ = map(int, header[:4])
        Rmin, Rmax = map(float, f.readline().split())
        phimin, phimax = map(float, f.readline().split())
        Zmin, Zmax = map(float, f.readline().split())
        data = np.loadtxt(f, dtype=float)
    expected = nR * nZ * nphi_plus1
    if data.shape[0] != expected:
        raise ValueError(f"Bvac file row mismatch: expected {expected}, got {data.shape[0]}")
    raw = data.reshape((nR, nphi_plus1, nZ, 3), order="C")
    components = np.moveaxis(raw, -1, 0)  # (3, nR, nphi+1, nZ)
    components = components[:, :, :-1, :]  # drop duplicate phi
    components = np.transpose(components, (0, 2, 3, 1))  # (3, nphi, nZ, nR)
    R_grid = np.linspace(Rmin, Rmax, nR)
    Z_grid = np.linspace(Zmin, Zmax, nZ)
    phi_grid = np.linspace(phimin, phimax, nphi_plus1)[:-1]
    return R_grid, phi_grid, Z_grid, components


def eval_complex_spline(x_grid, y_grid, values, x_points, y_points):
    spline_re = RectBivariateSpline(x_grid, y_grid, values.real, kx=3, ky=3)
    spline_im = RectBivariateSpline(x_grid, y_grid, values.imag, kx=3, ky=3)
    return spline_re.ev(x_points, y_points) + 1j * spline_im.ev(x_points, y_points)


def main() -> int:
    from libneo.biotsavart_fourier import (
        field_divfree,
        gauge_Anvac,
        read_Anvac_fourier,
        read_Bnvac_fourier,
        spline_gauged_Anvac,
    )

    if not FOURIER_FILE.exists():
        print(f"Missing Fourier data file {FOURIER_FILE}")
        return 1

    # Coil geometry (matches generate_tilted_coil.py)
    R_center = 2.0
    phi_center = 0.35
    Z_center = 0.8
    coil_radius = 2.0
    tilt_theta = 30.0
    tilt_psi = 35.0
    current = 1.0  # abampere

    center_xyz = np.array([R_center * np.cos(phi_center), R_center * np.sin(phi_center), Z_center])
    _, _, axis_hat = _plane_basis(phi_center, tilt_theta, tilt_psi)

    s_vals = np.linspace(-2.0, 2.0, 30)
    BR_analytical_full, Bphi_analytical_full, BZ_analytical_full, axis_points, phi_vals = circular_coil_analytical(
        s_vals, center_xyz, axis_hat, coil_radius, current
    )
    R_vals = np.sqrt(axis_points[0, :] ** 2 + axis_points[1, :] ** 2)
    Z_vals = axis_points[2, :]
    scale = 100.0  # m -> cm
    R_vals_cm = R_vals * scale
    Z_vals_cm = Z_vals * scale

    with netCDF4.Dataset(FOURIER_FILE) as nc:
        ntor_array = nc.variables["ntor"][:]
        nmax = int(np.max(ntor_array))

    grid0, AnR0, Anphi0, AnZ0, dAnphi_dR0, dAnphi_dZ0 = read_Anvac_fourier(str(FOURIER_FILE), ntor=0)

    valid_mask = (
        (R_vals_cm >= grid0.R_min)
        & (R_vals_cm <= grid0.R_max)
        & (Z_vals_cm >= grid0.Z_min)
        & (Z_vals_cm <= grid0.Z_max)
    )
    if not np.any(valid_mask):
        print("No axis points fall inside the spline domain.")
        return 1

    if not np.all(valid_mask):
        skipped = np.count_nonzero(~valid_mask)
        print(f"Warning: {skipped} axis points lie outside the spline grid and will be ignored.")

    R_eval_cm = R_vals_cm[valid_mask]
    Z_eval_cm = Z_vals_cm[valid_mask]
    phi_eval = phi_vals[valid_mask]
    s_eval = s_vals[valid_mask]
    axis_eval = axis_points[:, valid_mask]

    BR_analytical = BR_analytical_full[valid_mask]
    Bphi_analytical = Bphi_analytical_full[valid_mask]
    BZ_analytical = BZ_analytical_full[valid_mask]

    npts = R_eval_cm.size
    BR_fourier_ungauged = np.zeros(npts)
    Bphi_fourier_ungauged = np.zeros(npts)
    BZ_fourier_ungauged = np.zeros(npts)
    BR_fourier_gauged = np.zeros(npts)
    Bphi_fourier_gauged = np.zeros(npts)
    BZ_fourier_gauged = np.zeros(npts)
    Bn2_fourier_ungauged = None
    Bn2_fourier_gauged = None
    BR_bnvac = Bphi_bnvac = BZ_bnvac = None
    BR_bvac = Bphi_bvac = BZ_bvac = None
    has_bnvac = BNVAC_FILE.exists()
    has_bvac = BVAC_FILE.exists()

    cache = {0: (grid0, AnR0, Anphi0, AnZ0, dAnphi_dR0, dAnphi_dZ0)}

    ntor_target = 2
    for ntor in range(0, nmax + 1):
        if ntor in cache:
            grid, AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ = cache[ntor]
        else:
            grid, AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ = read_Anvac_fourier(str(FOURIER_FILE), ntor=ntor)
            cache[ntor] = (grid, AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ)

        spl_ungauged = spline_gauged_Anvac(grid, AnR, AnZ, ntor=ntor, Anphi=Anphi)
        BnR_ungauged, Bnphi_ungauged, BnZ_ungauged = evaluate_mode_at_points(
            spl_ungauged, R_eval_cm, Z_eval_cm, ntor
        )

        gauged_AnR, gauged_AnZ = gauge_Anvac(grid, AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ, ntor=ntor)
        spl_gauged = spline_gauged_Anvac(
            grid,
            gauged_AnR,
            gauged_AnZ,
            ntor=ntor,
            Anphi=Anphi if ntor == 0 else None,
        )
        BnR_gauged, Bnphi_gauged, BnZ_gauged = evaluate_mode_at_points(spl_gauged, R_eval_cm, Z_eval_cm, ntor)

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

        if ntor == ntor_target:
            Bn2_fourier_ungauged = BnR_ungauged, Bnphi_ungauged, BnZ_ungauged
            Bn2_fourier_gauged = BnR_gauged, Bnphi_gauged, BnZ_gauged

    Bn2_bnvac = None
    if has_bnvac:
        BR_bnvac = np.zeros(npts)
        Bphi_bnvac = np.zeros(npts)
        BZ_bnvac = np.zeros(npts)
        for ntor in range(0, nmax + 1):
            grid_b, BnR, Bnphi, BnZ = read_Bnvac_fourier(str(BNVAC_FILE), ntor=ntor)
            total_R = np.zeros(npts, dtype=complex)
            total_phi = np.zeros(npts, dtype=complex)
            total_Z = np.zeros(npts, dtype=complex)
            for kcoil in range(BnR.shape[0]):
                total_R += eval_complex_spline(grid_b.R, grid_b.Z, BnR[kcoil], R_eval_cm, Z_eval_cm)
                total_phi += eval_complex_spline(grid_b.R, grid_b.Z, Bnphi[kcoil], R_eval_cm, Z_eval_cm)
                total_Z += eval_complex_spline(grid_b.R, grid_b.Z, BnZ[kcoil], R_eval_cm, Z_eval_cm)
            phase = np.exp(1j * ntor * phi_eval)
            BR_bnvac += np.real(total_R * phase)
            Bphi_bnvac += np.real(total_phi * phase)
            BZ_bnvac += np.real(total_Z * phase)
            if ntor == ntor_target:
                Bn2_bnvac = total_R, total_phi, total_Z

    Bn2_bvac = None
    if has_bvac:
        R_grid_bvac, phi_grid_bvac, Z_grid_bvac, Bvac_components = read_Bvac_nemov(BVAC_FILE)
        phi_mod_eval = np.mod(phi_eval, 2.0 * np.pi)
        phi_wrap_eval = phi_mod_eval.copy()
        mask_wrap_eval = phi_wrap_eval >= phi_grid_bvac[-1]
        phi_wrap_eval[mask_wrap_eval] -= 2.0 * np.pi
        points_eval = np.column_stack((phi_wrap_eval, Z_eval_cm, R_eval_cm))
        interpolators_eval = [
            RegularGridInterpolator(
                (phi_grid_bvac, Z_grid_bvac, R_grid_bvac),
                Bvac_components[i],
                bounds_error=False,
                fill_value=None,
            )
            for i in range(3)
        ]
        BR_bvac = interpolators_eval[0](points_eval)
        Bphi_bvac = interpolators_eval[1](points_eval)
        BZ_bvac = interpolators_eval[2](points_eval)

        # Compute Fourier harmonic n=2 directly from grid data
        nphi_bvac = phi_grid_bvac.size
        accum_R = np.zeros(npts, dtype=complex)
        accum_phi = np.zeros(npts, dtype=complex)
        accum_Z = np.zeros(npts, dtype=complex)
        for idx, phi_sample in enumerate(phi_grid_bvac):
            values_R = RegularGridInterpolator(
                (Z_grid_bvac, R_grid_bvac), Bvac_components[0, idx, :, :], bounds_error=False, fill_value=None
            )(np.column_stack((Z_eval_cm, R_eval_cm)))
            values_phi = RegularGridInterpolator(
                (Z_grid_bvac, R_grid_bvac), Bvac_components[1, idx, :, :], bounds_error=False, fill_value=None
            )(np.column_stack((Z_eval_cm, R_eval_cm)))
            values_Z = RegularGridInterpolator(
                (Z_grid_bvac, R_grid_bvac), Bvac_components[2, idx, :, :], bounds_error=False, fill_value=None
            )(np.column_stack((Z_eval_cm, R_eval_cm)))
            accum_R += values_R * np.exp(-1j * ntor_target * phi_sample)
            accum_phi += values_phi * np.exp(-1j * ntor_target * phi_sample)
            accum_Z += values_Z * np.exp(-1j * ntor_target * phi_sample)
        Bn2_bvac = accum_R / nphi_bvac, accum_phi / nphi_bvac, accum_Z / nphi_bvac

    direct_cart = compute_direct_biot_savart(axis_eval)
    if direct_cart.shape[0] != npts:
        raise RuntimeError("Direct Biot-Savart solver returned unexpected sample count.")
    cosphi = np.cos(phi_eval)
    sinphi = np.sin(phi_eval)
    BR_direct = direct_cart[:, 0] * cosphi + direct_cart[:, 1] * sinphi
    Bphi_direct = -direct_cart[:, 0] * sinphi + direct_cart[:, 1] * cosphi
    BZ_direct = direct_cart[:, 2]

    # Prepare n=2 datasets
    def real_part(data_tuple):
        return [np.real(comp) if data_tuple is not None else None for comp in data_tuple] if data_tuple else (None, None, None)

    BRn2_fourier_ung, Bphin2_fourier_ung, BZn2_fourier_ung = real_part(Bn2_fourier_ungauged)
    BRn2_fourier_gauged, Bphin2_fourier_gauged, BZn2_fourier_gauged = real_part(Bn2_fourier_gauged)
    BRn2_bnvac, Bphin2_bnvac, BZn2_bnvac = real_part(Bn2_bnvac)
    BRn2_bvac, Bphin2_bvac, BZn2_bvac = real_part(Bn2_bvac)

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    axes[0].plot(s_eval, BR_analytical, "k-", label="Analytical", linewidth=2)
    axes[0].plot(s_eval, BR_fourier_ungauged, "r^--", label="Fourier (ungauged)", linewidth=2, markersize=5)
    axes[0].plot(s_eval, BR_fourier_gauged, "bs-.", label="Fourier (gauged n>0)", linewidth=2, markersize=5)
    if has_bnvac:
        axes[0].plot(s_eval, BR_bnvac, "mD:", label="Fourier (Bnvac)", linewidth=1.8, markersize=4)
    if has_bvac:
        axes[0].plot(s_eval, BR_bvac, "co-", label="Bvac (grid)", linewidth=1.8, markersize=4)
    axes[0].plot(s_eval, BR_direct, "gx-", label="Direct Biot-Savart", linewidth=1.8, markersize=5)
    axes[0].set_xlabel("Axis coordinate s (m)")
    axes[0].set_ylabel("$B_R$ (G)")
    axes[0].grid(True, alpha=0.3)
    axes[0].set_title("Radial Component")

    axes[1].plot(s_eval, Bphi_analytical, "k-", label="Analytical", linewidth=2)
    axes[1].plot(s_eval, Bphi_fourier_ungauged, "r^--", linewidth=2, markersize=5)
    axes[1].plot(s_eval, Bphi_fourier_gauged, "bs-.", linewidth=2, markersize=5)
    if has_bnvac:
        axes[1].plot(s_eval, Bphi_bnvac, "mD:", linewidth=1.8, markersize=4)
    if has_bvac:
        axes[1].plot(s_eval, Bphi_bvac, "co-", linewidth=1.8, markersize=4)
    axes[1].plot(s_eval, Bphi_direct, "gx-", linewidth=1.8, markersize=5)
    axes[1].set_xlabel("Axis coordinate s (m)")
    axes[1].set_ylabel("$B_\\phi$ (G)")
    axes[1].grid(True, alpha=0.3)
    axes[1].set_title("Toroidal Component")

    axes[2].plot(s_eval, BZ_analytical, "k-", label="Analytical", linewidth=2)
    axes[2].plot(s_eval, BZ_fourier_ungauged, "r^--", linewidth=2, markersize=5)
    axes[2].plot(s_eval, BZ_fourier_gauged, "bs-.", linewidth=2, markersize=5)
    if has_bnvac:
        axes[2].plot(s_eval, BZ_bnvac, "mD:", linewidth=1.8, markersize=4)
    if has_bvac:
        axes[2].plot(s_eval, BZ_bvac, "co-", linewidth=1.8, markersize=4)
    axes[2].plot(s_eval, BZ_direct, "gx-", linewidth=1.8, markersize=5)
    axes[2].set_xlabel("Axis coordinate s (m)")
    axes[2].set_ylabel("$B_Z$ (G)")
    axes[2].grid(True, alpha=0.3)
    axes[2].set_title("Vertical Component")

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig("ntor0_validation.png", dpi=150)
    print("\nSaved plot to ntor0_validation.png")

    fig2, axes2 = plt.subplots(1, 3, figsize=(16, 5))
    axes2[0].set_title("Radial Component (n=2)")
    axes2[1].set_title("Toroidal Component (n=2)")
    axes2[2].set_title("Vertical Component (n=2)")

    def plot_n2(axis, label, color, marker, series):
        if series is None:
            return
        axis.plot(s_eval, series, color + marker, label=label, linewidth=2, markersize=5)

    plot_n2(axes2[0], "Fourier (ungauged)", "r", "^--", BRn2_fourier_ung)
    plot_n2(axes2[0], "Fourier (gauged)", "b", "s-.", BRn2_fourier_gauged)
    plot_n2(axes2[0], "Fourier Bnvac", "m", "D:", BRn2_bnvac)
    plot_n2(axes2[0], "Bvac grid", "c", "o-", BRn2_bvac)

    plot_n2(axes2[1], "Fourier (ungauged)", "r", "^--", Bphin2_fourier_ung)
    plot_n2(axes2[1], "Fourier (gauged)", "b", "s-.", Bphin2_fourier_gauged)
    plot_n2(axes2[1], "Fourier Bnvac", "m", "D:", Bphin2_bnvac)
    plot_n2(axes2[1], "Bvac grid", "c", "o-", Bphin2_bvac)

    plot_n2(axes2[2], "Fourier (ungauged)", "r", "^--", BZn2_fourier_ung)
    plot_n2(axes2[2], "Fourier (gauged)", "b", "s-.", BZn2_fourier_gauged)
    plot_n2(axes2[2], "Fourier Bnvac", "m", "D:", BZn2_bnvac)
    plot_n2(axes2[2], "Bvac grid", "c", "o-", BZn2_bvac)

    for ax in axes2:
        ax.set_xlabel("Axis coordinate s (m)")
        ax.set_ylabel("Re[B^{(n=2)}] (G)")
        ax.grid(True, alpha=0.3)

    handles2, labels2 = axes2[0].get_legend_handles_labels()
    fig2.legend(handles2, labels2, loc="upper center", ncol=2)
    fig2.tight_layout(rect=(0, 0, 1, 0.92))
    fig2.savefig("ntor0_mode2.png", dpi=150)
    print("Saved plot to ntor0_mode2.png")

    print("\nField samples at first valid axis point (s = {:.3f} m):".format(s_eval[0]))
    print(
        "  Analytical      : BR = {:+.3e}, Bphi = {:+.3e}, BZ = {:+.3e}".format(
            BR_analytical[0], Bphi_analytical[0], BZ_analytical[0]
        )
    )
    print(
        "  Fourier ungauged: BR = {:+.3e}, Bphi = {:+.3e}, BZ = {:+.3e}".format(
            BR_fourier_ungauged[0], Bphi_fourier_ungauged[0], BZ_fourier_ungauged[0]
        )
    )
    print(
        "  Fourier gauged  : BR = {:+.3e}, Bphi = {:+.3e}, BZ = {:+.3e}".format(
            BR_fourier_gauged[0], Bphi_fourier_gauged[0], BZ_fourier_gauged[0]
        )
    )
    print(
        "  Direct Biot-Savart: BR = {:+.3e}, Bphi = {:+.3e}, BZ = {:+.3e}".format(
            BR_direct[0], Bphi_direct[0], BZ_direct[0]
        )
    )
    if BR_bnvac is not None:
        print(
            "  Fourier Bnvac   : BR = {:+.3e}, Bphi = {:+.3e}, BZ = {:+.3e}".format(
                BR_bnvac[0], Bphi_bnvac[0], BZ_bnvac[0]
            )
        )
    if BR_bvac is not None:
        print(
            "  Bvac (grid)     : BR = {:+.3e}, Bphi = {:+.3e}, BZ = {:+.3e}".format(
                BR_bvac[0], Bphi_bvac[0], BZ_bvac[0]
            )
        )

    def rel_error(ref, val):
        mask = np.abs(ref) > 1e-5
        if np.count_nonzero(mask) == 0:
            return 0.0
        return float(np.mean(np.abs(val[mask] - ref[mask]) / np.abs(ref[mask])) * 100.0)
    error_sets = [
        ("Fourier ungauged", BR_fourier_ungauged, Bphi_fourier_ungauged, BZ_fourier_ungauged),
        ("Fourier gauged", BR_fourier_gauged, Bphi_fourier_gauged, BZ_fourier_gauged),
        ("Direct Biot-Savart", BR_direct, Bphi_direct, BZ_direct),
    ]
    if BR_bnvac is not None:
        error_sets.append(("Fourier Bnvac", BR_bnvac, Bphi_bnvac, BZ_bnvac))
    if BR_bvac is not None:
        error_sets.append(("Bvac grid", BR_bvac, Bphi_bvac, BZ_bvac))

    print("\nRelative errors vs analytical solution (mean where |B| > 1e-5 G):")
    errors = []
    for label, BR_set, Bphi_set, BZ_set in error_sets:
        err = (
            rel_error(BR_analytical, BR_set),
            rel_error(Bphi_analytical, Bphi_set),
            rel_error(BZ_analytical, BZ_set),
        )
        errors.append((label, err))
        print("  {label:17s}-> BR: {0:.2f}%, Bphi: {1:.2f}%, BZ: {2:.2f}%".format(*err, label=label))

    max_err = max(max(err) for _, err in errors)

    if max_err < 5.0:
        print(
            "\n✓ TEST PASSED: All reconstructions (Fourier ungauged/gauged and direct) match analytical field "
            f"(worst-case error {max_err:.2f}%)"
        )
        return 0

    print(
        "\n✗ TEST FAILED: At least one reconstruction deviates more than 5% from the analytical field "
        f"(worst-case error {max_err:.2f}%)"
    )
    if Bn2_bvac is not None:
        def rel_error_complex(ref, val):
            mask = np.abs(ref) > 1e-5
            if np.count_nonzero(mask) == 0:
                return 0.0
            return float(np.mean(np.abs(val[mask] - ref[mask]) / np.abs(ref[mask])) * 100.0)

        print("\nRelative errors for n=2 harmonic (reference = Bvac grid):")
        for label, dataset in (
            ("Fourier ungauged", Bn2_fourier_ungauged),
            ("Fourier gauged", Bn2_fourier_gauged),
            ("Fourier Bnvac", Bn2_bnvac),
        ):
            if dataset is None:
                continue
            err = (
                rel_error_complex(Bn2_bvac[0], dataset[0]),
                rel_error_complex(Bn2_bvac[1], dataset[1]),
                rel_error_complex(Bn2_bvac[2], dataset[2]),
            )
            print("  {label:17s}-> BR: {0:.2f}%, Bphi: {1:.2f}%, BZ: {2:.2f}%".format(*err, label=label))

    return 1


if __name__ == "__main__":
    sys.exit(main())
