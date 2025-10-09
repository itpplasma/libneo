#!/usr/bin/env python3
"""Validate tilted coil vacuum fields against analytic and direct Biot-Savart solutions."""

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
from dataclasses import dataclass
from matplotlib.colors import LogNorm
from scipy.interpolate import RectBivariateSpline, RegularGridInterpolator
from typing import Dict, Iterable, List, Tuple

SCRIPT_DIR = Path(__file__).resolve().parent
BUILD_DIR = Path.cwd()
COIL_FILE = BUILD_DIR / "tilted_coil.dat"
FOURIER_FILE = BUILD_DIR / "tilted_coil_Anvac.nc"
BNVAC_FILE = BUILD_DIR / "tilted_coil_Bnvac.h5"
BVAC_FILE = BUILD_DIR / "tilted_coil_Bvac.dat"
AXIS_SOLVER = BUILD_DIR / "tilted_coil_axis_field.x"
CLIGHT = 2.99792458e10  # cm / s
MODE_MAX = 2
RELATIVE_TOLERANCE_PERCENT = 6.0
XFAIL_LABELS = {"Fourier gauged", "Fourier (gauged n>0)"}

sys.path.append(str(SCRIPT_DIR))
from tilted_coil_geometry import _plane_basis  # noqa: E402

@dataclass
class VectorComponents:
    radial: np.ndarray
    toroidal: np.ndarray
    vertical: np.ndarray

    def magnitude(self) -> np.ndarray:
        return np.sqrt(self.radial**2 + self.toroidal**2 + self.vertical**2)

DEFAULT_STYLE = {"linewidth": 2}
FIELD_STYLES: Dict[str, Dict[str, float]] = {
    "Analytical": {"color": "k", "linestyle": "-"},
    "Fourier (ungauged)": {"color": "r", "linestyle": "--", "marker": "^", "markersize": 5},
    "Fourier (gauged n>0)": {"color": "b", "linestyle": "-.", "marker": "s", "markersize": 5},
    "Fourier (gauged)": {"color": "b", "linestyle": "-.", "marker": "s", "markersize": 5},
    "Fourier Bnvac": {"color": "m", "linestyle": ":", "marker": "D", "markersize": 4},
    "Bvac grid": {"color": "c", "linestyle": "-", "marker": "o", "markersize": 4},
    "Direct Biot-Savart": {"color": "g", "linestyle": "-", "marker": "x", "markersize": 5},
    "An (ungauged)": {"color": "r", "linestyle": "--", "marker": "^", "markersize": 5},
    "An (gauged)": {"color": "b", "linestyle": "-.", "marker": "s", "markersize": 5},
    "An from B (Fourier)": {"color": "g", "linestyle": "-", "marker": "x", "markersize": 5},
    "An from B (Bnvac)": {"color": "m", "linestyle": ":", "marker": "D", "markersize": 4},
    "An from B (Bvac)": {"color": "c", "linestyle": "-", "marker": "o", "markersize": 4},
}

def get_plot_style(label: str) -> Dict[str, float]:
    style = DEFAULT_STYLE.copy()
    style.update(FIELD_STYLES.get(label, {}))
    return style

def mean_relative_error(reference: np.ndarray, values: np.ndarray, threshold: float = 1e-5) -> float:
    mask = np.abs(reference) > threshold
    if np.count_nonzero(mask) == 0:
        return 0.0
    return float(np.mean(np.abs(values[mask] - reference[mask]) / np.abs(reference[mask])) * 100.0)

def plot_component_series(
    filename: str,
    s_eval: np.ndarray,
    component_sets: List[Tuple[str, VectorComponents]],
    ylabels: Tuple[str, str, str],
    title_prefix: str,
) -> None:
    if not component_sets:
        return

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    components = (
        ("radial", ylabels[0], f"{title_prefix} Radial Component"),
        ("toroidal", ylabels[1], f"{title_prefix} Toroidal Component"),
        ("vertical", ylabels[2], f"{title_prefix} Vertical Component"),
    )

    handles = None
    labels = None
    for ax, (attr, ylabel, title) in zip(axes, components):
        for label, vec in component_sets:
            data = getattr(vec, attr)
            if data is None or np.all(np.isnan(data)):
                continue
            ax.plot(s_eval, data, **get_plot_style(label), label=label)
        ax.set_xlabel("Axis coordinate s (m)")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.grid(True, alpha=0.3)
        if handles is None:
            handles, labels = ax.get_legend_handles_labels()

    if handles:
        fig.legend(handles, labels, loc="upper center", ncol=2)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(filename, dpi=150)
    print(f"Saved plot to {filename}")

def plot_mode_series(
    filename: str,
    mode_index: int,
    s_eval: np.ndarray,
    component_sets: List[Tuple[str, VectorComponents]],
    ylabels: Tuple[str, str, str],
    title_prefix: str,
) -> None:
    title = f"{title_prefix} (n={mode_index})"
    plot_component_series(filename, s_eval, component_sets, ylabels, title_prefix=title)

def plot_plane_magnitudes(
    plane_variants: List[Tuple[str, np.ndarray, np.ndarray, np.ndarray]], filename: str
) -> None:
    if not plane_variants:
        return

    fig, axes = plt.subplots(
        1, len(plane_variants), figsize=(5.5 * len(plane_variants), 5), constrained_layout=True
    )
    if len(plane_variants) == 1:
        axes = [axes]

    direct_data = next((pd for pd in plane_variants if pd[0] == "Direct Biot-Savart"), None)
    if direct_data is not None:
        vmax = float(np.max(direct_data[3]))
        positive = direct_data[3][direct_data[3] > 0]
        vmin = float(np.min(positive)) if positive.size else vmax * 1e-12
    else:
        vmax = max(float(np.max(pd[3])) for pd in plane_variants)
        positive = np.concatenate(
            [pd[3][pd[3] > 0] for pd in plane_variants if np.any(pd[3] > 0)]
        )
        vmin = float(np.min(positive)) if positive.size else vmax * 1e-12

    vmin = max(vmin, 1e-12, vmax * 1e-12)
    norm = LogNorm(vmin=vmin, vmax=max(vmax, vmin * 10.0))

    pcm = None
    for ax, (label, R_mesh, Z_mesh, mag) in zip(axes, plane_variants):
        safe_mag = np.where(mag > vmin, mag, vmin)
        pcm = ax.pcolormesh(
            R_mesh,
            Z_mesh,
            safe_mag,
            shading="auto",
            cmap="viridis",
            norm=norm,
        )
        ax.set_title(label)
        ax.set_xlabel("R (m)")
        if ax is axes[0]:
            ax.set_ylabel("Z (m)")
        ax.set_aspect("equal")

    fig.colorbar(pcm, ax=axes, label="|B| (G)")
    fig.savefig(filename, dpi=150)
    print(f"Saved plot to {filename}")

def components_from_tuple(data: Tuple[np.ndarray, ...]) -> VectorComponents:
    def _to_array(values):
        if values is None:
            return None
        return np.asarray(values, dtype=float)

    radial = _to_array(data[0])
    toroidal = _to_array(data[1]) if len(data) > 1 else None
    vertical = _to_array(data[2]) if len(data) > 2 else None
    return VectorComponents(
        radial=radial,
        toroidal=toroidal,
        vertical=vertical,
    )

def components_from_totals(array: np.ndarray) -> VectorComponents:
    return VectorComponents(array[0, :], array[1, :], array[2, :])

def components_from_radial_vertical(radial: np.ndarray, vertical: np.ndarray) -> VectorComponents:
    filler = np.full_like(radial, np.nan)
    return VectorComponents(radial, filler, vertical)

def summarize_component_errors(
    analytical: VectorComponents,
    comparison_sets: List[Tuple[str, VectorComponents]],
    tolerance_percent: float,
    xfail_labels: Iterable[str],
    threshold: float = 1e-5,
) -> Tuple[float, List[Tuple[str, Tuple[float, float, float]]], List[Tuple[str, Tuple[float, float, float]]]]:
    results = []
    for label, data in comparison_sets:
        errs = (
            mean_relative_error(analytical.radial, data.radial, threshold),
            mean_relative_error(analytical.toroidal, data.toroidal, threshold),
            mean_relative_error(analytical.vertical, data.vertical, threshold),
        )
        results.append((label, errs))

    non_xfail = [max(err) for label, err in results if label not in xfail_labels]
    max_err = max(non_xfail) if non_xfail else 0.0
    return max_err, results, [(label, err) for label, err in results if label in xfail_labels]

def report_harmonic_errors(
    header: str,
    mode_indices: Iterable[int],
    reference_modes: Dict[int, Tuple[np.ndarray, np.ndarray, np.ndarray]],
    datasets: List[Tuple[str, Dict[int, Tuple[np.ndarray, np.ndarray, np.ndarray]]]],
    xfail_labels: Iterable[str],
    threshold: float = 1e-5,
) -> None:
    xfail_set = set(xfail_labels)
    for n in mode_indices:
        if n not in reference_modes:
            continue
        ref_components = components_from_tuple(reference_modes[n])
        print(f"\n{header.format(n=n)}")
        for label, data_dict in datasets:
            if n not in data_dict:
                continue
            comp = components_from_tuple(data_dict[n])
            errs = (
                mean_relative_error(ref_components.radial, comp.radial, threshold),
                mean_relative_error(ref_components.toroidal, comp.toroidal, threshold),
                mean_relative_error(ref_components.vertical, comp.vertical, threshold),
            )
            suffix = " [xfail]" if label in xfail_set else ""
            print(
                "  {label:17s}-> BR: {0:.2f}%, Bphi: {1:.2f}%, BZ: {2:.2f}%{suffix}".format(
                    *errs, label=label, suffix=suffix
                )
            )

def report_vector_potential_comparisons(
    mode_indices: Iterable[int],
    reference_modes: Dict[int, Tuple[np.ndarray, np.ndarray, np.ndarray]],
    reconstructions: List[Tuple[str, Dict[int, Tuple[np.ndarray, np.ndarray]]]],
    xfail_labels: Iterable[str],
    threshold: float = 1e-6,
) -> None:
    xfail_set = set(xfail_labels)
    print("\nRelative errors for gauged An components (reference = Fourier gauged)")

    for n in mode_indices:
        ref = reference_modes.get(n)
        if ref is None:
            continue
        ref_vec = components_from_tuple(ref)
        entries = [f"  n={n}"]
        for label, data_dict in reconstructions:
            if n not in data_dict:
                continue
            comp_R, comp_Z = data_dict[n]
            comp_vec = components_from_radial_vertical(np.asarray(comp_R), np.asarray(comp_Z))
            err_R = mean_relative_error(ref_vec.radial, comp_vec.radial, threshold)
            err_Z = mean_relative_error(ref_vec.vertical, comp_vec.vertical, threshold)
            suffix = " [xfail]" if label in xfail_set else ""
            entries.append(f"{label}: AR {err_R:.2f}%, AZ {err_Z:.2f}%{suffix}")
        if len(entries) > 1:
            print("; ".join(entries))

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

def evaluate_vector_potential_spline(spl: dict, R_points: np.ndarray, Z_points: np.ndarray):
    ncoil = len(spl["AnR_Re"])
    total_R = np.zeros(R_points.size, dtype=complex)
    total_Z = np.zeros(R_points.size, dtype=complex)
    total_phi = np.zeros(R_points.size, dtype=complex) if "Anphi_Re" in spl else None

    for kcoil in range(ncoil):
        total_R += spl["AnR_Re"][kcoil].ev(R_points, Z_points) + 1j * spl["AnR_Im"][kcoil].ev(R_points, Z_points)
        total_Z += spl["AnZ_Re"][kcoil].ev(R_points, Z_points) + 1j * spl["AnZ_Im"][kcoil].ev(R_points, Z_points)
        if total_phi is not None:
            total_phi += spl["Anphi_Re"][kcoil].ev(R_points, Z_points) + 1j * spl["Anphi_Im"][kcoil].ev(R_points, Z_points)

    return total_R, total_phi, total_Z

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

    # Coil geometry (matches tilted_coil_geometry.py)
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
    nphi_fourier = getattr(grid0, "nphi", None)

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
    An_total_ung = np.zeros((3, npts))
    An_total_gauged = np.zeros((3, npts))
    fourier_modes_ungauged = {}
    fourier_modes_gauged = {}
    fourier_modes_ungauged_A = {}
    fourier_modes_gauged_A = {}
    BR_bnvac = Bphi_bnvac = BZ_bnvac = None
    BR_bvac = Bphi_bvac = BZ_bvac = None
    has_bnvac = BNVAC_FILE.exists()
    has_bvac = BVAC_FILE.exists()

    cache = {0: (grid0, AnR0, Anphi0, AnZ0, dAnphi_dR0, dAnphi_dZ0)}

    def reconstruction_weight(ntor: int) -> float:
        """Return scaling so real reconstruction includes ±n harmonics."""
        if ntor == 0:
            return 1.0
        if nphi_fourier is not None and nphi_fourier > 0:
            if nphi_fourier % 2 == 0:
                nyquist = nphi_fourier // 2
                if ntor == nyquist and nyquist <= nmax:
                    return 1.0
        return 2.0

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
        AnR_ung, Anphi_ung, AnZ_ung = evaluate_vector_potential_spline(spl_ungauged, R_eval_cm, Z_eval_cm)

        gauged_AnR, gauged_AnZ = gauge_Anvac(grid, AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ, ntor=ntor)
        spl_gauged = spline_gauged_Anvac(
            grid,
            gauged_AnR,
            gauged_AnZ,
            ntor=ntor,
            Anphi=Anphi if ntor == 0 else None,
        )
        BnR_gauged, Bnphi_gauged, BnZ_gauged = evaluate_mode_at_points(spl_gauged, R_eval_cm, Z_eval_cm, ntor)
        AnR_gauged, Anphi_gauged, AnZ_gauged = evaluate_vector_potential_spline(spl_gauged, R_eval_cm, Z_eval_cm)

        phase = np.exp(1j * ntor * phi_eval)
        weight = reconstruction_weight(ntor)
        BR_fourier_ungauged += weight * np.real(BnR_ungauged * phase)
        Bphi_fourier_ungauged += weight * np.real(Bnphi_ungauged * phase)
        BZ_fourier_ungauged += weight * np.real(BnZ_ungauged * phase)
        An_total_ung[0, :] += weight * np.real(AnR_ung * phase)
        if Anphi_ung is not None:
            An_total_ung[1, :] += weight * np.real(Anphi_ung * phase)
        An_total_ung[2, :] += weight * np.real(AnZ_ung * phase)

        if ntor == 0:
            BR_fourier_gauged += weight * np.real(BnR_ungauged * phase)
            Bphi_fourier_gauged += weight * np.real(Bnphi_ungauged * phase)
            BZ_fourier_gauged += weight * np.real(BnZ_ungauged * phase)
            An_total_gauged[0, :] += weight * np.real(AnR_ung * phase)
            if Anphi_ung is not None:
                An_total_gauged[1, :] += weight * np.real(Anphi_ung * phase)
            An_total_gauged[2, :] += weight * np.real(AnZ_ung * phase)
        else:
            BR_fourier_gauged += weight * np.real(BnR_gauged * phase)
            Bphi_fourier_gauged += weight * np.real(Bnphi_gauged * phase)
            BZ_fourier_gauged += weight * np.real(BnZ_gauged * phase)
            An_total_gauged[0, :] += weight * np.real(AnR_gauged * phase)
            if Anphi_gauged is not None:
                An_total_gauged[1, :] += weight * np.real(Anphi_gauged * phase)
            An_total_gauged[2, :] += weight * np.real(AnZ_gauged * phase)

        if ntor <= MODE_MAX:
            fourier_modes_ungauged[ntor] = (BnR_ungauged.copy(), Bnphi_ungauged.copy(), BnZ_ungauged.copy())
            fourier_modes_gauged[ntor] = (BnR_gauged.copy(), Bnphi_gauged.copy(), BnZ_gauged.copy())
            fourier_modes_ungauged_A[ntor] = (AnR_ung.copy(), None if Anphi_ung is None else Anphi_ung.copy(), AnZ_ung.copy())
            fourier_modes_gauged_A[ntor] = (AnR_gauged.copy(), None if Anphi_gauged is None else Anphi_gauged.copy(), AnZ_gauged.copy())

    bnvac_modes = {}
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
            weight = reconstruction_weight(ntor)
            BR_bnvac += weight * np.real(total_R * phase)
            Bphi_bnvac += weight * np.real(total_phi * phase)
            BZ_bnvac += weight * np.real(total_Z * phase)
            if ntor <= MODE_MAX:
                bnvac_modes[ntor] = (total_R.copy(), total_phi.copy(), total_Z.copy())

    bvac_modes = {}
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

        nphi_bvac = phi_grid_bvac.size
        mode_accum = {
            n: [np.zeros(npts, dtype=complex) for _ in range(3)]
            for n in range(MODE_MAX + 1)
        }
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
            for n in range(MODE_MAX + 1):
                factor = np.exp(-1j * n * phi_sample)
                mode_accum[n][0] += values_R * factor
                mode_accum[n][1] += values_phi * factor
                mode_accum[n][2] += values_Z * factor
        for n in range(MODE_MAX + 1):
            bvac_modes[n] = tuple(acc / nphi_bvac for acc in mode_accum[n])

    direct_cart = compute_direct_biot_savart(axis_eval)
    if direct_cart.shape[0] != npts:
        raise RuntimeError("Direct Biot-Savart solver returned unexpected sample count.")
    cosphi = np.cos(phi_eval)
    sinphi = np.sin(phi_eval)
    BR_direct = direct_cart[:, 0] * cosphi + direct_cart[:, 1] * sinphi
    Bphi_direct = -direct_cart[:, 0] * sinphi + direct_cart[:, 1] * cosphi
    BZ_direct = direct_cart[:, 2]

    # 2D R-Z plane comparison at phi = phi_center
    phi_plane = phi_center
    cos_plane = np.cos(phi_plane)
    sin_plane = np.sin(phi_plane)
    R_mesh_cm_fourier, Z_mesh_cm_fourier = np.meshgrid(grid0.R, grid0.Z, indexing="xy")
    R_mesh_m_fourier = R_mesh_cm_fourier / scale
    Z_mesh_m_fourier = Z_mesh_cm_fourier / scale
    plane_shape_fourier = R_mesh_cm_fourier.shape
    plane_R_flat_cm = R_mesh_cm_fourier.ravel()
    plane_Z_flat_cm = Z_mesh_cm_fourier.ravel()
    plane_R_flat_m = R_mesh_m_fourier.ravel()
    plane_Z_flat_m = Z_mesh_m_fourier.ravel()

    plane_fields: List[Tuple[str, np.ndarray, np.ndarray, np.ndarray]] = []

    try:
        plane_points = np.vstack(
            (
                plane_R_flat_m * cos_plane,
                plane_R_flat_m * sin_plane,
                plane_Z_flat_m,
            )
        )
        direct_plane_cart = compute_direct_biot_savart(plane_points)
        BR_plane_direct = (
            direct_plane_cart[:, 0] * cos_plane + direct_plane_cart[:, 1] * sin_plane
        ).reshape(plane_shape_fourier)
        Bphi_plane_direct = (
            -direct_plane_cart[:, 0] * sin_plane + direct_plane_cart[:, 1] * cos_plane
        ).reshape(plane_shape_fourier)
        BZ_plane_direct = direct_plane_cart[:, 2].reshape(plane_shape_fourier)
        mag_plane_direct = np.sqrt(
            BR_plane_direct**2 + Bphi_plane_direct**2 + BZ_plane_direct**2
        )
        plane_fields.append(("Direct Biot-Savart", R_mesh_m_fourier, Z_mesh_m_fourier, mag_plane_direct))
    except FileNotFoundError as err:
        print(f"Direct Biot-Savart plane evaluation skipped: {err}")

    n_plane = plane_R_flat_cm.size
    BR_plane_ung = np.zeros(n_plane)
    Bphi_plane_ung = np.zeros(n_plane)
    BZ_plane_ung = np.zeros(n_plane)
    BR_plane_gauged = np.zeros(n_plane)
    Bphi_plane_gauged = np.zeros(n_plane)
    BZ_plane_gauged = np.zeros(n_plane)
    for ntor in range(0, nmax + 1):
        grid, AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ = cache[ntor]
        spl_ung = spline_gauged_Anvac(grid, AnR, AnZ, ntor=ntor, Anphi=Anphi)
        BnR_plane_ung, Bnphi_plane_ung, BnZ_plane_ung = evaluate_mode_at_points(
            spl_ung, plane_R_flat_cm, plane_Z_flat_cm, ntor
        )
        weight = reconstruction_weight(ntor)
        phase_plane = np.exp(1j * ntor * phi_plane)
        BR_plane_ung += weight * np.real(BnR_plane_ung * phase_plane)
        Bphi_plane_ung += weight * np.real(Bnphi_plane_ung * phase_plane)
        BZ_plane_ung += weight * np.real(BnZ_plane_ung * phase_plane)
        if ntor == 0:
            BR_plane_gauged += weight * np.real(BnR_plane_ung * phase_plane)
            Bphi_plane_gauged += weight * np.real(Bnphi_plane_ung * phase_plane)
            BZ_plane_gauged += weight * np.real(BnZ_plane_ung * phase_plane)
        else:
            gauged_AnR, gauged_AnZ = gauge_Anvac(
                grid, AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ, ntor=ntor
            )
            spl_g = spline_gauged_Anvac(
                grid,
                gauged_AnR,
                gauged_AnZ,
                ntor=ntor,
                Anphi=None,
            )
            BnR_plane_g, Bnphi_plane_g, BnZ_plane_g = evaluate_mode_at_points(
                spl_g, plane_R_flat_cm, plane_Z_flat_cm, ntor
            )
            BR_plane_gauged += weight * np.real(BnR_plane_g * phase_plane)
            Bphi_plane_gauged += weight * np.real(Bnphi_plane_g * phase_plane)
            BZ_plane_gauged += weight * np.real(BnZ_plane_g * phase_plane)

    BR_plane_ung = BR_plane_ung.reshape(plane_shape_fourier)
    Bphi_plane_ung = Bphi_plane_ung.reshape(plane_shape_fourier)
    BZ_plane_ung = BZ_plane_ung.reshape(plane_shape_fourier)
    BR_plane_gauged = BR_plane_gauged.reshape(plane_shape_fourier)
    Bphi_plane_gauged = Bphi_plane_gauged.reshape(plane_shape_fourier)
    BZ_plane_gauged = BZ_plane_gauged.reshape(plane_shape_fourier)
    mag_plane_ung = np.sqrt(BR_plane_ung**2 + Bphi_plane_ung**2 + BZ_plane_ung**2)
    mag_plane_gauged = np.sqrt(
        BR_plane_gauged**2 + Bphi_plane_gauged**2 + BZ_plane_gauged**2
    )

    plane_fields.append(("Fourier (ungauged)", R_mesh_m_fourier, Z_mesh_m_fourier, mag_plane_ung))
    plane_fields.append(("Fourier (gauged n>0)", R_mesh_m_fourier, Z_mesh_m_fourier, mag_plane_gauged))

    plot_plane_magnitudes(plane_fields, "tilted_coil_RZ_comparison.png")

    # Prepare n=2 datasets
    def to_real_dict(data_dict):
        return {
            n: tuple(np.real(comp) for comp in value)
            for n, value in data_dict.items()
        }

    fourier_modes_ungauged_B_real = to_real_dict(fourier_modes_ungauged)
    fourier_modes_gauged_B_real = to_real_dict(fourier_modes_gauged)
    bnvac_modes_B_real = to_real_dict(bnvac_modes)
    bvac_modes_B_real = to_real_dict(bvac_modes)

    def to_real_dict_A(data_dict):
        return {
            n: tuple(None if comp is None else np.real(comp) for comp in value)
            for n, value in data_dict.items()
        }

    fourier_modes_ungauged_A_real = to_real_dict_A(fourier_modes_ungauged_A)
    fourier_modes_gauged_A_real = to_real_dict_A(fourier_modes_gauged_A)

    def reconstruct_gauged_from_B(data_dict):
        result = {}
        for n, comps in data_dict.items():
            if n == 0:
                continue
            Br, _, Bz = comps
            A_R = 1j * R_eval_cm * Bz / n
            A_Z = -1j * R_eval_cm * Br / n
            result[n] = (
                np.real(A_R),
                np.real(A_Z),
            )
        return result

    An_from_B_fourier = reconstruct_gauged_from_B(fourier_modes_gauged)
    An_from_B_bnvac = reconstruct_gauged_from_B(bnvac_modes)
    An_from_B_bvac = reconstruct_gauged_from_B(bvac_modes)

    axis_analytical = VectorComponents(BR_analytical, Bphi_analytical, BZ_analytical)
    axis_fourier_ung = VectorComponents(BR_fourier_ungauged, Bphi_fourier_ungauged, BZ_fourier_ungauged)
    axis_fourier_g = VectorComponents(BR_fourier_gauged, Bphi_fourier_gauged, BZ_fourier_gauged)
    axis_direct = VectorComponents(BR_direct, Bphi_direct, BZ_direct)
    axis_comparisons: List[Tuple[str, VectorComponents]] = [
        ("Analytical", axis_analytical),
        ("Fourier (ungauged)", axis_fourier_ung),
        ("Fourier (gauged n>0)", axis_fourier_g),
        ("Direct Biot-Savart", axis_direct),
    ]
    if BR_bnvac is not None:
        axis_comparisons.append(("Fourier Bnvac", VectorComponents(BR_bnvac, Bphi_bnvac, BZ_bnvac)))
    if BR_bvac is not None:
        axis_comparisons.append(("Bvac grid", VectorComponents(BR_bvac, Bphi_bvac, BZ_bvac)))

    plot_component_series(
        "tilted_coil_validation.png",
        s_eval,
        axis_comparisons,
        ("$B_R$ (G)", "$B_\\phi$ (G)", "$B_Z$ (G)"),
        title_prefix="Magnetic Field",
    )

    an_total_ung_vec = components_from_totals(An_total_ung)
    an_total_gauged_vec = components_from_totals(An_total_gauged)
    an_comparisons = [
        ("An (ungauged)", an_total_ung_vec),
        ("An (gauged)", an_total_gauged_vec),
    ]
    plot_component_series(
        "tilted_coil_An_total.png",
        s_eval,
        an_comparisons,
        ("Re[A_R] (G·cm)", "Re[A_\\phi] (G·cm)", "Re[A_Z] (G·cm)"),
        title_prefix="Vector Potential",
    )

    harmonic_datasets = [
        ("Fourier (ungauged)", fourier_modes_ungauged_B_real),
        ("Fourier (gauged)", fourier_modes_gauged_B_real),
        ("Fourier Bnvac", bnvac_modes_B_real),
        ("Bvac grid", bvac_modes_B_real),
    ]

    for mode_idx in range(0, MODE_MAX + 1):
        component_sets: List[Tuple[str, VectorComponents]] = []
        if mode_idx == 0:
            component_sets.append(("Analytical", axis_analytical))
        for label, data_dict in harmonic_datasets:
            if mode_idx not in data_dict:
                continue
            component_sets.append((label, components_from_tuple(data_dict[mode_idx])))
        if component_sets:
            plot_mode_series(
                f"tilted_coil_mode{mode_idx}.png",
                mode_idx,
                s_eval,
                component_sets,
                (r"Re[B_R^{(n)}] (G)", r"Re[B_\\phi^{(n)}] (G)", r"Re[B_Z^{(n)}] (G)"),
                title_prefix="Magnetic Field",
            )

    vector_mode_datasets = [
        ("Fourier (ungauged)", fourier_modes_ungauged_A_real),
        ("Fourier (gauged)", fourier_modes_gauged_A_real),
    ]
    vector_recon_sets = [
        ("An from B (Fourier)", An_from_B_fourier),
        ("An from B (Bnvac)", An_from_B_bnvac),
        ("An from B (Bvac)", An_from_B_bvac),
    ]

    for mode_idx in range(0, MODE_MAX + 1):
        component_sets = []
        for label, data_dict in vector_mode_datasets:
            if mode_idx not in data_dict:
                continue
            component_sets.append((label, components_from_tuple(data_dict[mode_idx])))
        for label, data_dict in vector_recon_sets:
            if mode_idx not in data_dict:
                continue
            comp_R, comp_Z = data_dict[mode_idx]
            component_sets.append(
                (label, components_from_radial_vertical(np.asarray(comp_R), np.asarray(comp_Z)))
            )
        if component_sets:
            plot_mode_series(
                f"tilted_coil_An_mode{mode_idx}.png",
                mode_idx,
                s_eval,
                component_sets,
                (r"Re[A_R^{(n)}] (G·cm)", r"Re[A_\\phi^{(n)}] (G·cm)", r"Re[A_Z^{(n)}] (G·cm)"),
                title_prefix="Vector Potential",
            )

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

    comparison_for_errors = axis_comparisons[1:]
    max_err, error_records, xfail_records = summarize_component_errors(
        axis_analytical,
        comparison_for_errors,
        RELATIVE_TOLERANCE_PERCENT,
        XFAIL_LABELS,
    )

    print("\nRelative errors vs analytical solution (mean where |B| > 1e-5 G):")
    for label, err in error_records:
        suffix = " [xfail]" if label in XFAIL_LABELS else ""
        print(
            "  {label:17s}-> BR: {0:.2f}%, Bphi: {1:.2f}%, BZ: {2:.2f}%{suffix}".format(
                *err, label=label, suffix=suffix
            )
        )

    if max_err <= RELATIVE_TOLERANCE_PERCENT:
        if xfail_records:
            print("\nExpected failures (marked xfail):")
            for label, err in xfail_records:
                print(
                    "  {label:17s}-> BR: {0:.2f}%, Bphi: {1:.2f}%, BZ: {2:.2f}%".format(
                        *err, label=label
                    )
                )
        print(
            "\n✓ TEST PASSED: All reconstructions (Fourier ungauged/gauged and direct) match analytical field "
            f"(worst-case error {max_err:.2f}% <= {RELATIVE_TOLERANCE_PERCENT:.1f}%)"
        )
        return 0

    print(
        "\n✗ TEST FAILED: At least one reconstruction deviates more than "
        f"{RELATIVE_TOLERANCE_PERCENT:.1f}% from the analytical field (worst-case error {max_err:.2f}%)"
    )

    report_harmonic_errors(
        "Relative errors for n={n} harmonic (reference = Bvac grid):",
        range(0, MODE_MAX + 1),
        bvac_modes_B_real,
        harmonic_datasets,
        XFAIL_LABELS,
    )

    report_vector_potential_comparisons(
        range(1, MODE_MAX + 1),
        fourier_modes_gauged_A_real,
        vector_recon_sets,
        XFAIL_LABELS,
    )

    return 1

if __name__ == "__main__":
    sys.exit(main())
