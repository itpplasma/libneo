#!/usr/bin/env python3
"""Compare coil_tools Fourier field data against Anvac reconstructions using all modes."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from math import pi
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
from numpy import (
    amax,
    amin,
    array,
    asarray,
    cos,
    linspace,
    log10,
    mean,
    median,
    sin,
    sqrt,
    zeros,
    zeros_like,
)
from numpy.fft import fft
from scipy.interpolate import RectBivariateSpline

from libneo.biotsavart_fourier import (
    field_divfree,
    gauge_Anvac,
    read_Anvac_fourier_all,
    read_Bnvac_fourier_all,
    reconstruct_field_from_modes,
    spline_gauged_Anvac,
)
from _magfie import compute_field_from_gpec_file


FOURIER_REFERENCE_CURRENT = 0.1  # coil_tools stores modes per abampere (=10 A)
DEFAULT_FFT_SAMPLES = 192


@dataclass
class Grid:
    R: np.ndarray
    Z: np.ndarray
    nR: int
    nZ: int
    R_min: float
    R_max: float
    Z_min: float
    Z_max: float


@dataclass
class ModeSet:
    grid: Grid
    mode_numbers: np.ndarray  # (nmodes,)
    BnR: np.ndarray  # (nmodes, ncoil, nR, nZ)
    Bnphi: np.ndarray
    BnZ: np.ndarray

    @property
    def nmodes(self) -> int:
        return self.mode_numbers.size

    @property
    def ncoil(self) -> int:
        return self.BnR.shape[1]


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare coil_tools Fourier, Anvac and direct Biot–Savart results",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("reference", type=Path, help="coil_tools Fourier HDF5 file (Bnvac)")
    parser.add_argument("test", type=Path, help="coil_tools vector-potential NetCDF file (Anvac)")
    parser.add_argument("--currents", type=Path, required=True, help="Text file with coil currents (A)")
    parser.add_argument("--coil-files", type=Path, nargs="+", required=True,
                        help="Coil geometry files in GPEC format")
    parser.add_argument("--ntor", type=int, default=0,
                        help="Toroidal mode to highlight in diagnostic tables")
    parser.add_argument("--fft-samples", type=int, default=DEFAULT_FFT_SAMPLES,
                        help="Number of toroidal samples for direct FFT")
    parser.add_argument("--prefactor", type=float, default=FOURIER_REFERENCE_CURRENT,
                        help="Scale converting SI currents (A) to coil_tools reference units")
    parser.add_argument("--per-coil-output", type=Path,
                        help="Output PNG for per-coil |B| maps (phi=0)")
    parser.add_argument("--sum-output", type=Path,
                        help="Output PNG for total-field |B| comparison (phi=0)")
    parser.add_argument("--deriv-diff-output", type=Path,
                        help="Output PNG comparing stored vs spline dAphi/dR,Z for --ntor")
    parser.add_argument("--axis-output", type=Path,
                        help="Output PNG for on-axis field comparison")
    parser.add_argument("--axis-origin", type=float, nargs=3,
                        help="Axis origin (cm) for validation plot")
    parser.add_argument("--axis-normal", type=float, nargs=3,
                        help="Axis direction vector (not normalised)")
    parser.add_argument("--axis-range", type=float, default=60.0,
                        help="Half-length of axis segment (cm)")
    parser.add_argument("--axis-samples", type=int, default=181,
                        help="Number of samples along the validation axis")
    parser.add_argument("--coil-radius", type=float,
                        help="Optional analytic circular coil radius (cm) for reference curve")
    parser.add_argument("--dpi", type=int, default=150, help="DPI for saved figures")
    parser.add_argument("--show", action="store_true", help="Display figures interactively")
    return parser.parse_args()


def _build_grid(grid_raw) -> Grid:
    return Grid(
        R=array(grid_raw.R, copy=True),
        Z=array(grid_raw.Z, copy=True),
        nR=grid_raw.nR,
        nZ=grid_raw.nZ,
        R_min=float(grid_raw.R_min),
        R_max=float(grid_raw.R_max),
        Z_min=float(grid_raw.Z_min),
        Z_max=float(grid_raw.Z_max),
    )


def _load_fourier_modes(path: Path) -> ModeSet:
    grid_raw, mode_numbers, BnR, Bnphi, BnZ = read_Bnvac_fourier_all(str(path))
    grid = _build_grid(grid_raw)
    return ModeSet(
        grid=grid,
        mode_numbers=mode_numbers,
        BnR=array(BnR, dtype=complex),
        Bnphi=array(Bnphi, dtype=complex),
        BnZ=array(BnZ, dtype=complex),
    )


def _load_anvac_modes(path: Path, target_mode_numbers: np.ndarray) -> Tuple[ModeSet, ModeSet, Dict[str, np.ndarray], ModeSet]:
    grid_raw, mode_numbers, AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ, AnX_raw, AnY_raw, AnZ_raw =         read_Anvac_fourier_all(str(path))
    if not np.array_equal(mode_numbers, target_mode_numbers):
        raise ValueError(
            "Anvac NetCDF file does not contain the same toroidal modes as the Fourier reference"
        )

    grid = _build_grid(grid_raw)
    nmodes, ncoil, nR, nZ = AnR.shape

    BnR_stored = np.empty((nmodes, ncoil, nR, nZ), dtype=complex)
    Bnphi_stored = np.empty_like(BnR_stored)
    BnZ_stored = np.empty_like(BnR_stored)

    BnR_spline = np.empty_like(BnR_stored)
    Bnphi_spline = np.empty_like(BnR_stored)
    BnZ_spline = np.empty_like(BnR_stored)

    dAphi_dR_spline = np.empty_like(Anphi)
    dAphi_dZ_spline = np.empty_like(Anphi)

    for mode_idx, ntor in enumerate(mode_numbers):
        AnR_mode = AnR[mode_idx]
        Anphi_mode = Anphi[mode_idx]
        AnZ_mode = AnZ[mode_idx]
        dAnphi_dR_mode = dAnphi_dR[mode_idx]
        dAnphi_dZ_mode = dAnphi_dZ[mode_idx]

        gauged_AnR, gauged_AnZ = gauge_Anvac(
            grid_raw,
            AnR_mode,
            Anphi_mode,
            AnZ_mode,
            dAnphi_dR_mode,
            dAnphi_dZ_mode,
            ntor=ntor,
        )
        spl = spline_gauged_Anvac(grid_raw, gauged_AnR, gauged_AnZ, ntor=ntor, Anphi=Anphi_mode)
        BnR_tmp, Bnphi_tmp, BnZ_tmp = field_divfree(spl, grid.R, grid.Z, ntor=ntor)
        BnR_stored[mode_idx] = BnR_tmp
        Bnphi_stored[mode_idx] = Bnphi_tmp
        BnZ_stored[mode_idx] = BnZ_tmp

        if ntor == 0:
            for coil_idx in range(ncoil):
                spl_real = RectBivariateSpline(grid.R, grid.Z, Anphi_mode[coil_idx].real)
                spl_imag = RectBivariateSpline(grid.R, grid.Z, Anphi_mode[coil_idx].imag)
                dAphi_dR_spline[mode_idx, coil_idx] = (
                    spl_real(grid.R, grid.Z, dx=1, dy=0, grid=True)
                    + 1j * spl_imag(grid.R, grid.Z, dx=1, dy=0, grid=True)
                )
                dAphi_dZ_spline[mode_idx, coil_idx] = (
                    spl_real(grid.R, grid.Z, dx=0, dy=1, grid=True)
                    + 1j * spl_imag(grid.R, grid.Z, dx=0, dy=1, grid=True)
                )
        else:
            dAphi_dR_spline[mode_idx] = dAnphi_dR_mode
            dAphi_dZ_spline[mode_idx] = dAnphi_dZ_mode

        gauged_AnR_spline, gauged_AnZ_spline = gauge_Anvac(
            grid_raw,
            AnR_mode,
            Anphi_mode,
            AnZ_mode,
            dAphi_dR_spline[mode_idx],
            dAphi_dZ_spline[mode_idx],
            ntor=ntor,
        )
        spl_spline = spline_gauged_Anvac(
            grid_raw,
            gauged_AnR_spline,
            gauged_AnZ_spline,
            ntor=ntor,
            Anphi=Anphi_mode,
        )
        BnR_tmp, Bnphi_tmp, BnZ_tmp = field_divfree(spl_spline, grid.R, grid.Z, ntor=ntor)
        BnR_spline[mode_idx] = BnR_tmp
        Bnphi_spline[mode_idx] = Bnphi_tmp
        BnZ_spline[mode_idx] = BnZ_tmp

    BnR_ungauged, Bnphi_ungauged, BnZ_ungauged = _compute_ungauged_modes(
        mode_numbers,
        AnR,
        Anphi,
        AnZ,
        dAnphi_dR,
        dAnphi_dZ,
        grid_raw.R,
        grid_raw.Z,
    )

    diagnostics = {
        'dAphi_dR_stored': dAnphi_dR,
        'dAphi_dZ_stored': dAnphi_dZ,
        'dAphi_dR_spline': dAphi_dR_spline,
        'dAphi_dZ_spline': dAphi_dZ_spline,
        'AnX_raw': AnX_raw,
        'AnY_raw': AnY_raw,
        'AnZ_raw': AnZ_raw,
    }

    stored = ModeSet(grid=grid, mode_numbers=mode_numbers, BnR=BnR_stored,
                     Bnphi=Bnphi_stored, BnZ=BnZ_stored)
    spline = ModeSet(grid=grid, mode_numbers=mode_numbers, BnR=BnR_spline,
                     Bnphi=Bnphi_spline, BnZ=BnZ_spline)
    ungauged = ModeSet(grid=grid, mode_numbers=mode_numbers, BnR=BnR_ungauged,
                       Bnphi=Bnphi_ungauged, BnZ=BnZ_ungauged)
    return stored, spline, diagnostics, ungauged


def _ensure_same_grid(lhs: Grid, rhs: Grid) -> None:
    if lhs.nR != rhs.nR or lhs.nZ != rhs.nZ:
        raise ValueError("Grid dimensions differ between inputs")
    if not np.allclose(lhs.R, rhs.R) or not np.allclose(lhs.Z, rhs.Z):
        raise ValueError("Reference and test grids do not match")


def _read_currents(path: Path) -> np.ndarray:
    values: List[float] = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            values.extend(float(tok) for tok in stripped.split())
    return asarray(values, dtype=float)


def _load_coil_projections(filenames: Sequence[Path]) -> List[Tuple[np.ndarray, np.ndarray]]:
    projections: List[Tuple[np.ndarray, np.ndarray]] = []
    for filename in filenames:
        with open(filename, "r", encoding="utf-8") as handle:
            header = handle.readline().split()
            if len(header) < 4:
                raise ValueError(f"Malformed GPEC header in {filename}")
            ncoil = int(header[0])
            nseg_plus = int(header[2])
            nseg = nseg_plus - 1
            for _ in range(ncoil):
                coords = []
                for _ in range(nseg):
                    parts = handle.readline().split()
                    if len(parts) != 3:
                        raise ValueError(f"Unexpected end of file in {filename}")
                    coords.append([float(val) for val in parts])
                closing = handle.readline().split()
                if len(closing) != 3:
                    raise ValueError(f"Missing closing point in {filename}")
                coords.append([float(val) for val in closing])
                coords = asarray(coords)
                coords_cm = coords * 100.0
                R_path = np.hypot(coords_cm[:, 0], coords_cm[:, 1])
                Z_path = coords_cm[:, 2]
                projections.append((R_path, Z_path))
    return projections


def _sum_over_coils(weights: np.ndarray, modes: np.ndarray) -> np.ndarray:
    """Sum over coils with given weights; modes has shape (nmodes,ncoil,...)"""
    if modes.shape[1] != weights.size:
        raise ValueError("Weights do not match number of coils")
    return np.tensordot(weights, modes, axes=(0, 1))


def _compute_ungauged_modes(
    mode_numbers: np.ndarray,
    AnR: np.ndarray,
    Anphi: np.ndarray,
    AnZ: np.ndarray,
    dAphi_dR: np.ndarray,
    dAphi_dZ: np.ndarray,
    R_vals: np.ndarray,
    Z_vals: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute curl(A) for ungauged Fourier harmonics in cylindrical coordinates."""
    ntor = mode_numbers.reshape((-1, 1, 1, 1)).astype(float)
    R_grid = R_vals.reshape((1, 1, -1, 1))
    inv_R = np.divide(
        1.0, R_grid,
        out=np.zeros_like(R_grid, dtype=float),
        where=np.abs(R_grid) > 0.0,
    )
    try:
        dAnR_dZ = np.gradient(AnR, Z_vals, axis=3, edge_order=2)
        dAnZ_dR = np.gradient(AnZ, R_vals, axis=2, edge_order=2)
    except ValueError:
        dAnR_dZ = np.gradient(AnR, Z_vals, axis=3, edge_order=1)
        dAnZ_dR = np.gradient(AnZ, R_vals, axis=2, edge_order=1)
    BR = 1j * ntor * AnZ * inv_R - dAphi_dZ
    Bphi = dAnR_dZ - dAnZ_dR
    BZ = Anphi * inv_R + dAphi_dR - 1j * ntor * AnR * inv_R
    return BR, Bphi, BZ


def _reconstruct_field_map_on_grid(
    BnR_modes: np.ndarray,
    Bnphi_modes: np.ndarray,
    BnZ_modes: np.ndarray,
    mode_numbers: np.ndarray,
    phi: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Reconstruct Bx, By, Bz on the cylindrical grid at fixed toroidal angle phi."""
    nmodes, nR, nZ = BnR_modes.shape
    if nmodes != mode_numbers.size:
        raise ValueError("Mode number list does not match mode array shape")

    BR = BnR_modes[0].real.copy()
    Bphi = Bnphi_modes[0].real.copy()
    BZ = BnZ_modes[0].real.copy()

    if nmodes > 1:
        exp_factor = np.exp(1j * mode_numbers[1:] * phi)
        for idx, factor in enumerate(exp_factor, start=1):
            BR += 2.0 * np.real(BnR_modes[idx] * factor)
            Bphi += 2.0 * np.real(Bnphi_modes[idx] * factor)
            BZ += 2.0 * np.real(BnZ_modes[idx] * factor)

    cosphi = cos(phi)
    sinphi = sin(phi)
    Bx = BR * cosphi - Bphi * sinphi
    By = BR * sinphi + Bphi * cosphi
    return Bx, By, BZ


def _field_magnitude(Bx: np.ndarray, By: np.ndarray, Bz: np.ndarray) -> np.ndarray:
    return np.sqrt(Bx**2 + By**2 + Bz**2)


def _coils_field_magnitude(
    mode_set: ModeSet,
    weights: np.ndarray,
    phi: float,
) -> np.ndarray:
    """Return |B| maps (phi fixed) for each individual coil including weights."""
    nmodes, ncoil, _, _ = mode_set.BnR.shape
    results = np.empty((ncoil, mode_set.grid.nR, mode_set.grid.nZ), dtype=float)
    for coil_idx in range(ncoil):
        scaled_BnR = mode_set.BnR[:, coil_idx] * weights[coil_idx]
        scaled_Bnphi = mode_set.Bnphi[:, coil_idx] * weights[coil_idx]
        scaled_BnZ = mode_set.BnZ[:, coil_idx] * weights[coil_idx]
        Bx, By, Bz = _reconstruct_field_map_on_grid(
            scaled_BnR,
            scaled_Bnphi,
            scaled_BnZ,
            mode_set.mode_numbers,
            phi,
        )
        results[coil_idx] = _field_magnitude(Bx, By, Bz)
    return results


def _create_per_coil_plot(
    args: argparse.Namespace,
    fourier: ModeSet,
    anvac_stored: ModeSet,
    anvac_spline: ModeSet,
    weights: np.ndarray,
    coil_projections: List[Tuple[np.ndarray, np.ndarray]] | None,
) -> None:
    ncoil = fourier.ncoil
    phi = 0.0

    maps_fourier = _coils_field_magnitude(fourier, weights, phi)
    maps_anvac = _coils_field_magnitude(anvac_stored, weights, phi)
    maps_spline = _coils_field_magnitude(anvac_spline, weights, phi)

    titles = ("Fourier", "Anvac (stored)", "Anvac (spline)")
    data = (maps_fourier, maps_anvac, maps_spline)

    vmax = max(map_data.max() for map_data in data)
    vmin = min(map_data.min() for map_data in data)
    norm = Normalize(vmin=vmin, vmax=vmax)

    fig, axs = plt.subplots(ncoil, 3, figsize=(12, 3 * ncoil), layout="constrained")
    extent = [fourier.grid.R_min, fourier.grid.R_max, fourier.grid.Z_min, fourier.grid.Z_max]

    for coil_idx in range(ncoil):
        for col, label in enumerate(titles):
            ax = axs[coil_idx, col] if ncoil > 1 else axs[col]
            im = ax.imshow(
                data[col][coil_idx].T,
                origin="lower",
                cmap="magma",
                extent=extent,
                norm=norm,
                interpolation="bilinear",
            )
            ax.set_title(f"Coil {coil_idx + 1}: {label}")
            ax.set_xlabel("R [cm]")
            ax.set_ylabel("Z [cm]")
            if coil_projections:
                for R_path, Z_path in coil_projections:
                    ax.plot(R_path, Z_path, color="black", linewidth=0.8, alpha=0.6)
            ax.set_aspect("equal", adjustable="box")

    fig.colorbar(im, ax=axs, location="bottom", fraction=0.05, pad=0.08, label="|B| [G]")
    fig.savefig(args.per_coil_output, dpi=args.dpi, bbox_inches="tight")
    if not args.show:
        plt.close(fig)


def _create_sum_plot(
    args: argparse.Namespace,
    grid: Grid,
    coil_projections: List[Tuple[np.ndarray, np.ndarray]] | None,
    fourier_map: np.ndarray,
    anvac_map: np.ndarray,
    spline_map: np.ndarray,
    ungauged_map: np.ndarray,
    direct_map: np.ndarray | None,
) -> None:
    panels = []
    titles = []

    if direct_map is not None:
        reference = direct_map
        panels.append(direct_map)
        titles.append("Direct")
    else:
        reference = fourier_map

    panels.extend([fourier_map, anvac_map, spline_map, ungauged_map])
    titles.extend(["Fourier", "Anvac (stored)", "Anvac (spline)", "Anvac (ungauged)"])

    finite_ref = reference[np.isfinite(reference)]
    if finite_ref.size:
        vmin = np.percentile(finite_ref, 1.0)
        vmax = np.percentile(finite_ref, 99.0)
        if vmax <= vmin:
            vmin = finite_ref.min()
            vmax = finite_ref.max()
    else:
        vmin = reference.min()
        vmax = reference.max()

    norm = Normalize(vmin=vmin, vmax=vmax)

    fig, axs = plt.subplots(1, len(panels), figsize=(4 * len(panels), 4.5), layout="constrained")
    if len(panels) == 1:
        axs = [axs]
    extent = [grid.R_min, grid.R_max, grid.Z_min, grid.Z_max]

    for ax, panel, title in zip(axs, panels, titles):
        im = ax.imshow(
            panel.T,
            origin="lower",
            cmap="magma",
            extent=extent,
            norm=norm,
            interpolation="bilinear",
        )
        ax.set_title(title)
        ax.set_xlabel("R [cm]")
        ax.set_ylabel("Z [cm]")
        if coil_projections:
            for R_path, Z_path in coil_projections:
                ax.plot(R_path, Z_path, color="black", linewidth=0.8, alpha=0.6)
        ax.set_aspect("equal", adjustable="box")

    fig.colorbar(im, ax=axs, location="bottom", fraction=0.05, pad=0.08, label="|B(φ=0)| [G]")
    fig.savefig(args.sum_output, dpi=args.dpi, bbox_inches="tight")
    if not args.show:
        plt.close(fig)


def _create_deriv_comparison_plot(
    args: argparse.Namespace,
    diagnostics: Dict[str, np.ndarray],
    grid: Grid,
    weights: np.ndarray,
    mode_numbers: np.ndarray,
) -> None:
    if args.ntor not in mode_numbers:
        raise ValueError(f"Requested ntor={args.ntor} not present in data")
    mode_idx = int(np.where(mode_numbers == args.ntor)[0][0])

    dAphi_dR_stored = diagnostics['dAphi_dR_stored'][mode_idx]
    dAphi_dZ_stored = diagnostics['dAphi_dZ_stored'][mode_idx]
    dAphi_dR_spline = diagnostics['dAphi_dR_spline'][mode_idx]
    dAphi_dZ_spline = diagnostics['dAphi_dZ_spline'][mode_idx]

    weights_expand = weights[:, None, None]
    dR_stored_sum = np.tensordot(weights, dAphi_dR_stored, axes=(0, 0))
    dZ_stored_sum = np.tensordot(weights, dAphi_dZ_stored, axes=(0, 0))
    dR_spline_sum = np.tensordot(weights, dAphi_dR_spline, axes=(0, 0))
    dZ_spline_sum = np.tensordot(weights, dAphi_dZ_spline, axes=(0, 0))

    mag_dR_stored = np.abs(dR_stored_sum)
    mag_dZ_stored = np.abs(dZ_stored_sum)
    mag_dR_spline = np.abs(dR_spline_sum)
    mag_dZ_spline = np.abs(dZ_spline_sum)

    extent = [grid.R_min, grid.R_max, grid.Z_min, grid.Z_max]
    fig, axs = plt.subplots(2, 2, figsize=(12, 10), layout="constrained")

    vmin_dR = min(mag_dR_stored.min(), mag_dR_spline.min())
    vmax_dR = max(mag_dR_stored.max(), mag_dR_spline.max())
    vmin_dZ = min(mag_dZ_stored.min(), mag_dZ_spline.min())
    vmax_dZ = max(mag_dZ_stored.max(), mag_dZ_spline.max())

    im0 = axs[0, 0].imshow(mag_dR_stored.T, origin="lower", cmap="viridis",
                            extent=extent, interpolation="bilinear",
                            vmin=vmin_dR, vmax=vmax_dR)
    axs[0, 0].set_title("|dAphi/dR| (stored)")
    axs[0, 0].set_xlabel("R [cm]")
    axs[0, 0].set_ylabel("Z [cm]")
    fig.colorbar(im0, ax=axs[0, 0])

    im1 = axs[0, 1].imshow(mag_dR_spline.T, origin="lower", cmap="viridis",
                            extent=extent, interpolation="bilinear",
                            vmin=vmin_dR, vmax=vmax_dR)
    axs[0, 1].set_title("|dAphi/dR| (spline)")
    axs[0, 1].set_xlabel("R [cm]")
    axs[0, 1].set_ylabel("Z [cm]")
    fig.colorbar(im1, ax=axs[0, 1])

    im2 = axs[1, 0].imshow(mag_dZ_stored.T, origin="lower", cmap="viridis",
                            extent=extent, interpolation="bilinear",
                            vmin=vmin_dZ, vmax=vmax_dZ)
    axs[1, 0].set_title("|dAphi/dZ| (stored)")
    axs[1, 0].set_xlabel("R [cm]")
    axs[1, 0].set_ylabel("Z [cm]")
    fig.colorbar(im2, ax=axs[1, 0])

    im3 = axs[1, 1].imshow(mag_dZ_spline.T, origin="lower", cmap="viridis",
                            extent=extent, interpolation="bilinear",
                            vmin=vmin_dZ, vmax=vmax_dZ)
    axs[1, 1].set_title("|dAphi/dZ| (spline)")
    axs[1, 1].set_xlabel("R [cm]")
    axs[1, 1].set_ylabel("Z [cm]")
    fig.colorbar(im3, ax=axs[1, 1])

    fig.savefig(args.deriv_diff_output, dpi=args.dpi, bbox_inches="tight")
    if not args.show:
        plt.close(fig)


def _compute_segment_direct_field(
    coil_files: Sequence[Path],
    coil_currents: np.ndarray,
    x_eval: np.ndarray,
    y_eval: np.ndarray,
    z_eval: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    Bx_total = np.zeros_like(x_eval, dtype=float)
    By_total = np.zeros_like(y_eval, dtype=float)
    Bz_total = np.zeros_like(z_eval, dtype=float)

    offset = 0
    for coil_file in coil_files:
        with open(coil_file, "r", encoding="utf-8") as handle:
            header = handle.readline().split()
            if len(header) < 4:
                raise ValueError(f"Malformed GPEC header in {coil_file}")
            ncoil = int(header[0])
        currents_slice = coil_currents[offset:offset + ncoil]
        offset += ncoil

        Bx, By, Bz = compute_field_from_gpec_file(
            str(coil_file),
            currents_slice,
            x_eval.reshape(-1),
            y_eval.reshape(-1),
            z_eval.reshape(-1),
        )
        Bx_total += Bx.reshape(x_eval.shape)
        By_total += By.reshape(y_eval.shape)
        Bz_total += Bz.reshape(z_eval.shape)

    return Bx_total, By_total, Bz_total


def _compute_direct_fourier_modes(
    coil_files: Sequence[Path],
    coil_currents: np.ndarray,
    grid: Grid,
    mode_numbers: np.ndarray,
    n_phi: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    R_grid = grid.R
    Z_grid = grid.Z
    nR, nZ = grid.nR, grid.nZ

    phi_grid = linspace(0.0, 2.0 * pi, n_phi, endpoint=False)
    BR_all = zeros((nR, nZ, n_phi))
    Bphi_all = zeros_like(BR_all)
    BZ_all = zeros_like(BR_all)

    R_mesh, Z_mesh = np.meshgrid(R_grid, Z_grid, indexing="ij")
    base_R = R_mesh.reshape(-1)
    base_Z = Z_mesh.reshape(-1)

    for idx, phi in enumerate(phi_grid):
        cosphi = cos(phi)
        sinphi = sin(phi)
        x_eval = (base_R * cosphi).reshape(nR, nZ)
        y_eval = (base_R * sinphi).reshape(nR, nZ)
        z_eval = base_Z.reshape(nR, nZ)

        Bx, By, Bz = _compute_segment_direct_field(coil_files, coil_currents, x_eval, y_eval, z_eval)

        BR_all[:, :, idx] = Bx * cosphi + By * sinphi
        Bphi_all[:, :, idx] = -Bx * sinphi + By * cosphi
        BZ_all[:, :, idx] = Bz

    fft_norm = 1.0 / n_phi
    fft_BR = fft(BR_all, axis=2) * fft_norm
    fft_Bphi = fft(Bphi_all, axis=2) * fft_norm
    fft_BZ = fft(BZ_all, axis=2) * fft_norm

    nmodes = mode_numbers.size
    BnR = np.empty((nmodes, nR, nZ), dtype=complex)
    Bnphi = np.empty_like(BnR)
    BnZ = np.empty_like(BnR)

    n_phi_available = fft_BR.shape[2]
    for idx, ntor in enumerate(mode_numbers):
        if ntor >= n_phi_available:
            raise ValueError(
                f"Requested mode ntor={ntor} exceeds Nyquist limit for n_phi={n_phi_available}"
            )
        BnR[idx] = fft_BR[:, :, ntor]
        Bnphi[idx] = fft_Bphi[:, :, ntor]
        BnZ[idx] = fft_BZ[:, :, ntor]

    return BnR, Bnphi, BnZ, phi_grid


def _compute_axis_profiles(
    grid: Grid,
    mode_numbers: np.ndarray,
    BnR_modes: np.ndarray,
    Bnphi_modes: np.ndarray,
    BnZ_modes: np.ndarray,
    origin: np.ndarray,
    direction: np.ndarray,
    s_vals: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    points = origin[None, :] + s_vals[:, None] * direction[None, :]
    Bx, By, Bz = reconstruct_field_from_modes(
        BnR_modes,
        Bnphi_modes,
        BnZ_modes,
        mode_numbers,
        grid.R,
        grid.Z,
        points[:, 0],
        points[:, 1],
        points[:, 2],
    )
    B_parallel = Bx * direction[0] + By * direction[1] + Bz * direction[2]
    B_mag = _field_magnitude(Bx, By, Bz)
    return B_parallel, B_mag


def _create_axis_plot(
    args: argparse.Namespace,
    axis_points: np.ndarray,
    direction: np.ndarray,
    s_vals: np.ndarray,
    B_parallel_fourier: np.ndarray,
    B_parallel_anvac: np.ndarray,
    B_parallel_spline: np.ndarray,
    B_parallel_ungauged: np.ndarray,
    B_parallel_direct_modes: np.ndarray,
    B_parallel_direct_segments: np.ndarray,
    analytic_curve: np.ndarray | None,
) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 4.5), layout="constrained")
    ax.plot(s_vals, B_parallel_fourier, label="Fourier", linewidth=1.4)
    ax.plot(s_vals, B_parallel_anvac, label="Anvac (stored)", linewidth=1.2)
    ax.plot(s_vals, B_parallel_spline, label="Anvac (spline)", linewidth=1.2)
    ax.plot(s_vals, B_parallel_ungauged, label="Anvac (ungauged)", linewidth=1.0)
    ax.plot(s_vals, B_parallel_direct_modes, label="Direct FFT", linewidth=1.0)
    ax.plot(s_vals, B_parallel_direct_segments, label="Direct segments", linewidth=1.0, linestyle="--")
    if analytic_curve is not None:
        ax.plot(s_vals, analytic_curve, label="Analytic circle", linewidth=1.0, linestyle=":")
    ax.set_xlabel("Axis coordinate s [cm]")
    ax.set_ylabel("B_parallel [G]")
    ax.legend(loc="best")
    fig.savefig(args.axis_output, dpi=args.dpi, bbox_inches="tight")
    if not args.show:
        plt.close(fig)


def _print_stats(reference_map: np.ndarray, comparison_map: np.ndarray, label: str) -> Tuple[float, float]:
    rel_err = np.abs(comparison_map - reference_map) / np.maximum(reference_map, 1e-15)
    median_err = float(median(rel_err) * 100.0)
    mean_err = float(mean(rel_err) * 100.0)
    print(f"Median relative error ({label}): {median_err:.4f}%")
    print(f"Mean relative error   ({label}): {mean_err:.4f}%")
    return median_err, mean_err


def main() -> None:
    args = _parse_args()

    fourier = _load_fourier_modes(args.reference)
    anvac_stored, anvac_spline, diagnostics, anvac_ungauged = _load_anvac_modes(
        args.test, fourier.mode_numbers
    )
    _ensure_same_grid(fourier.grid, anvac_stored.grid)

    coil_projections = _load_coil_projections(args.coil_files)

    currents = _read_currents(args.currents)
    if currents.size != fourier.ncoil:
        raise ValueError(
            f"Current file provides {currents.size} entries but geometry defines {fourier.ncoil} coils"
        )

    # coil_tools stores Fourier amplitudes in cgs (abampere) units. Convert the
    # supplied SI currents to the reference units for mode summation, but keep
    # the SI values for direct Biot–Savart checks.
    currents_modes = currents * args.prefactor

    # Sum over coils
    BnR_fourier_sum = _sum_over_coils(currents_modes, fourier.BnR)
    Bnphi_fourier_sum = _sum_over_coils(currents_modes, fourier.Bnphi)
    BnZ_fourier_sum = _sum_over_coils(currents_modes, fourier.BnZ)

    BnR_anvac_sum = _sum_over_coils(currents_modes, anvac_stored.BnR)
    Bnphi_anvac_sum = _sum_over_coils(currents_modes, anvac_stored.Bnphi)
    BnZ_anvac_sum = _sum_over_coils(currents_modes, anvac_stored.BnZ)

    BnR_spline_sum = _sum_over_coils(currents_modes, anvac_spline.BnR)
    Bnphi_spline_sum = _sum_over_coils(currents_modes, anvac_spline.Bnphi)
    BnZ_spline_sum = _sum_over_coils(currents_modes, anvac_spline.BnZ)

    BnR_ungauged_sum = _sum_over_coils(currents_modes, anvac_ungauged.BnR)
    Bnphi_ungauged_sum = _sum_over_coils(currents_modes, anvac_ungauged.Bnphi)
    BnZ_ungauged_sum = _sum_over_coils(currents_modes, anvac_ungauged.BnZ)

    # Direct Fourier modes (already weighted by currents)
    BnR_direct_sum, Bnphi_direct_sum, BnZ_direct_sum, _ = _compute_direct_fourier_modes(
        args.coil_files,
        currents,
        fourier.grid,
        fourier.mode_numbers,
        args.fft_samples,
    )

    # Field maps at phi = 0
    Bx_fourier_map, By_fourier_map, Bz_fourier_map = _reconstruct_field_map_on_grid(
        BnR_fourier_sum,
        Bnphi_fourier_sum,
        BnZ_fourier_sum,
        fourier.mode_numbers,
        phi=0.0,
    )
    Bx_anvac_map, By_anvac_map, Bz_anvac_map = _reconstruct_field_map_on_grid(
        BnR_anvac_sum,
        Bnphi_anvac_sum,
        BnZ_anvac_sum,
        fourier.mode_numbers,
        phi=0.0,
    )
    Bx_spline_map, By_spline_map, Bz_spline_map = _reconstruct_field_map_on_grid(
        BnR_spline_sum,
        Bnphi_spline_sum,
        BnZ_spline_sum,
        fourier.mode_numbers,
        phi=0.0,
    )
    Bx_ungauged_map, By_ungauged_map, Bz_ungauged_map = _reconstruct_field_map_on_grid(
        BnR_ungauged_sum,
        Bnphi_ungauged_sum,
        BnZ_ungauged_sum,
        fourier.mode_numbers,
        phi=0.0,
    )
    Bx_direct_map, By_direct_map, Bz_direct_map = _reconstruct_field_map_on_grid(
        BnR_direct_sum,
        Bnphi_direct_sum,
        BnZ_direct_sum,
        fourier.mode_numbers,
        phi=0.0,
    )

    magnitude_fourier = _field_magnitude(Bx_fourier_map, By_fourier_map, Bz_fourier_map)
    magnitude_anvac = _field_magnitude(Bx_anvac_map, By_anvac_map, Bz_anvac_map)
    magnitude_spline = _field_magnitude(Bx_spline_map, By_spline_map, Bz_spline_map)
    magnitude_ungauged = _field_magnitude(Bx_ungauged_map, By_ungauged_map, Bz_ungauged_map)
    magnitude_direct = _field_magnitude(Bx_direct_map, By_direct_map, Bz_direct_map)

    # Plots
    if args.per_coil_output is not None:
        _create_per_coil_plot(args, fourier, anvac_stored, anvac_spline, currents_modes, coil_projections)

    if args.sum_output is not None:
        _create_sum_plot(
            args,
            fourier.grid,
            coil_projections,
            magnitude_fourier,
            magnitude_anvac,
            magnitude_spline,
            magnitude_ungauged,
            magnitude_direct,
        )

    if args.deriv_diff_output is not None:
        _create_deriv_comparison_plot(
            args,
            diagnostics,
            fourier.grid,
            currents_modes,
            fourier.mode_numbers,
        )

    print("Comparison statistics (relative errors wrt Fourier map at φ=0):")
    _print_stats(magnitude_fourier, magnitude_anvac, "Anvac (stored)")
    _print_stats(magnitude_fourier, magnitude_spline, "Anvac (spline)")
    _print_stats(magnitude_fourier, magnitude_ungauged, "Anvac (ungauged)")
    _print_stats(magnitude_fourier, magnitude_direct, "Direct")

    if args.axis_output is not None:
        if args.axis_origin is None or args.axis_normal is None:
            raise ValueError("Axis origin and normal must be provided for axis output")
        origin = asarray(args.axis_origin, dtype=float)
        direction = asarray(args.axis_normal, dtype=float)
        direction /= np.linalg.norm(direction)
        s_vals = linspace(-args.axis_range, args.axis_range, args.axis_samples)
        points = origin[None, :] + s_vals[:, None] * direction[None, :]

        B_parallel_fourier, _ = _compute_axis_profiles(
            fourier.grid,
            fourier.mode_numbers,
            BnR_fourier_sum,
            Bnphi_fourier_sum,
            BnZ_fourier_sum,
            origin,
            direction,
            s_vals,
        )
        B_parallel_anvac, _ = _compute_axis_profiles(
            fourier.grid,
            fourier.mode_numbers,
            BnR_anvac_sum,
            Bnphi_anvac_sum,
            BnZ_anvac_sum,
            origin,
            direction,
            s_vals,
        )
        B_parallel_spline, _ = _compute_axis_profiles(
            fourier.grid,
            fourier.mode_numbers,
            BnR_spline_sum,
            Bnphi_spline_sum,
            BnZ_spline_sum,
            origin,
            direction,
            s_vals,
        )
        B_parallel_ungauged, _ = _compute_axis_profiles(
            fourier.grid,
            fourier.mode_numbers,
            BnR_ungauged_sum,
            Bnphi_ungauged_sum,
            BnZ_ungauged_sum,
            origin,
            direction,
            s_vals,
        )
        B_parallel_direct_modes, _ = _compute_axis_profiles(
            fourier.grid,
            fourier.mode_numbers,
            BnR_direct_sum,
            Bnphi_direct_sum,
            BnZ_direct_sum,
            origin,
            direction,
            s_vals,
        )

        # Direct segment evaluation for reference
        Bx_seg, By_seg, Bz_seg = _compute_segment_direct_field(
            args.coil_files,
            currents,
            points[:, 0].reshape(-1, 1),
            points[:, 1].reshape(-1, 1),
            points[:, 2].reshape(-1, 1),
        )
        B_parallel_segments = (
            Bx_seg[:, 0] * direction[0]
            + By_seg[:, 0] * direction[1]
            + Bz_seg[:, 0] * direction[2]
        )

        analytic_curve = None
        if args.coil_radius is not None:
            radius_m = args.coil_radius / 100.0
            s_vals_m = s_vals / 100.0
            mu0 = 4e-7 * pi
            analytic_curve = (
                mu0
                * currents.sum()
                * radius_m**2
                / (2.0 * (radius_m**2 + s_vals_m**2) ** 1.5)
                * 1.0e4
            )

        _create_axis_plot(
            args,
            points,
            direction,
            s_vals,
            B_parallel_fourier,
            B_parallel_anvac,
            B_parallel_spline,
            B_parallel_ungauged,
            B_parallel_direct_modes,
            B_parallel_segments,
            analytic_curve,
        )


if __name__ == "__main__":
    main()
