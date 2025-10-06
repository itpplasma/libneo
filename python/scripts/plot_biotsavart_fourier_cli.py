#!/usr/bin/env python3
"""Generate comparison plots for coil_tools Biot–Savart Fourier outputs.

This CLI expects a Fourier HDF5 file (Bn) and the matching vector-potential NetCDF
file (An). It produces richer plots than the legacy interactive helper:

* per-coil log10|B|^2 pcolormesh for reference vs. test
* coil-summed log10|B|^2 pcolormesh for ref/test (and optional residual)

The script is intentionally lightweight: it does **not** evaluate the optional
segment-based direct field and therefore runs in headless CI environments.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Dict, Tuple

import matplotlib

if os.environ.get("MPLBACKEND") is None:  # pragma: no cover - CLI usage
    matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from scipy.interpolate import RectBivariateSpline

from libneo.biotsavart_fourier import (  # type: ignore  # noqa: E402
    field_divfree,
    gauged_Anvac_from_Bnvac,
    gauge_Anvac,
    read_Anvac_fourier,
    read_Bnvac_fourier,
    spline_gauged_Anvac,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path, help="Fourier HDF5 file from coil_tools (Bn)")
    parser.add_argument("test", type=Path, help="Vector-potential NetCDF from coil_tools (An)")
    parser.add_argument("--ntor", type=int, default=2, help="Toroidal mode index")
    parser.add_argument("--per-coil-output", type=Path,
                        help="PNG written with per-coil comparison pcolormesh")
    parser.add_argument("--sum-output", type=Path,
                        help="PNG written with coil-summed comparison pcolormesh")
    parser.add_argument("--deriv-diff-output", type=Path,
                        help="PNG comparing spline vs analytic gauge derivatives")
    parser.add_argument("--dpi", type=int, default=150)
    parser.add_argument("--show", action="store_true", help="Display figures interactively")
    parser.add_argument("--currents", type=Path,
                        help="Optional text file with coil currents (A)")
    parser.add_argument("--prefactor", type=float, default=0.1,
                        help="Scaling applied to currents (matches coil_tools tables)")
    parser.add_argument("--coil-files", nargs="*", type=Path,
                        help="Optional coil geometry files (currently unused in this script)")
    return parser.parse_args()


def load_modes(
    reference: Path,
    test: Path,
    ntor: int,
) -> Tuple[
    Tuple[np.ndarray, np.ndarray, np.ndarray],
    Tuple[np.ndarray, np.ndarray, np.ndarray],
    np.ndarray,
    np.ndarray,
    Dict[str, np.ndarray] | None,
]:
    ref_grid, BnR_ref, Bnphi_ref, BnZ_ref = read_Bnvac_fourier(str(reference), ntor=ntor)
    test_grid, AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ = read_Anvac_fourier(str(test), ntor=ntor)

    if ref_grid.nR != test_grid.nR or ref_grid.nZ != test_grid.nZ:
        raise ValueError("Reference and test grids have different dimensions")

    R = np.array(ref_grid.R)
    Z = np.array(ref_grid.Z)

    ref_gauged = gauged_Anvac_from_Bnvac(ref_grid, BnR_ref, BnZ_ref, ntor=ntor)
    ref_spl = spline_gauged_Anvac(ref_grid, *ref_gauged, ntor=ntor)
    ref_field = field_divfree(ref_spl, R, Z, ntor=ntor)

    analytic_gauged_R, analytic_gauged_Z = gauge_Anvac(
        test_grid, AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ, ntor=ntor
    )
    if ntor != 0:
        analytic_gauged_R = analytic_gauged_R - 0.5j * Anphi
    analytic_spl = spline_gauged_Anvac(test_grid, analytic_gauged_R, analytic_gauged_Z, ntor=ntor)
    analytic_field = field_divfree(analytic_spl, R, Z, ntor=ntor)

    if ntor == 0:
        axis_field = _compute_axisymmetric_field(test_grid, AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ)
        return ref_field, axis_field, R, Z, None

    spline_gauged_R, spline_gauged_Z = _gauge_anvac_with_spline(test_grid, AnR, Anphi, AnZ, ntor)
    spline_spl = spline_gauged_Anvac(test_grid, spline_gauged_R, spline_gauged_Z, ntor=ntor)
    spline_field_rect = _field_divfree_rect(spline_spl, R, Z, ntor)
    analytic_field_rect = _field_divfree_rect(analytic_spl, R, Z, ntor)

    diagnostics: Dict[str, np.ndarray] = {
        "gauged_diff_R": spline_gauged_R - analytic_gauged_R,
        "gauged_diff_Z": spline_gauged_Z - analytic_gauged_Z,
        "BnR_spline": spline_field_rect[0],
        "Bnphi_spline": spline_field_rect[1],
        "BnZ_spline": spline_field_rect[2],
        "BnR_analytic": analytic_field_rect[0],
        "Bnphi_analytic": analytic_field_rect[1],
        "BnZ_analytic": analytic_field_rect[2],
    }

    return ref_field, analytic_field, R, Z, diagnostics


def load_currents(path: Path | None, ncoil: int, prefactor: float) -> np.ndarray:
    if path is None:
        return np.ones(ncoil, dtype=float) * prefactor
    values = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            values.extend(float(tok) for tok in line.split())
    if len(values) != ncoil:
        raise ValueError(f"Current file {path} expected {ncoil} values, found {len(values)}")
    return np.asarray(values, dtype=float) * prefactor


def magnitude_squared(BnR: np.ndarray, Bnphi: np.ndarray, BnZ: np.ndarray) -> np.ndarray:
    return (BnR * np.conj(BnR) + Bnphi * np.conj(Bnphi) + BnZ * np.conj(BnZ)).real


def make_per_coil_plot(ref_field, test_field, R: np.ndarray, Z: np.ndarray,
                        output: Path | None, dpi: int, show: bool) -> None:
    if output is None and not show:
        return

    ref_mag = magnitude_squared(*ref_field)
    test_mag = magnitude_squared(*test_field)
    ncoil = ref_mag.shape[0]

    norm = Normalize(vmin=np.nanmin(np.log10(ref_mag)), vmax=np.nanmax(np.log10(ref_mag)))
    fig, axes = plt.subplots(ncoil, 2, figsize=(8, 4 * ncoil), layout='constrained')
    axes = np.atleast_2d(axes)

    R_edges = edges_from_centres(R)
    Z_edges = edges_from_centres(Z)

    pcm = None
    for idx in range(ncoil):
        for col, (title, mag) in enumerate((("Reference", ref_mag), ("Test", test_mag))):
            ax = axes[idx, col]
            pcm = ax.pcolormesh(R_edges, Z_edges, np.log10(mag[idx].T), cmap='magma', norm=norm)
            ax.set_title(f"Coil {idx + 1} – {title}")
            ax.set_xlabel('R [cm]')
            ax.set_ylabel('Z [cm]')

    if pcm is not None:
        fig.colorbar(pcm, ax=axes.ravel().tolist(), shrink=0.5, label=r'$\log_{10}|\vec{B}_n|^2$')

    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, dpi=dpi)
    if show:
        plt.show()
    else:
        plt.close(fig)


def make_sum_plot(ref_field, test_field, R: np.ndarray, Z: np.ndarray,
                   output: Path | None, dpi: int, show: bool) -> None:
    if output is None and not show:
        return

    ref_sum = np.sum(ref_field[0], axis=0), np.sum(ref_field[1], axis=0), np.sum(ref_field[2], axis=0)
    test_sum = np.sum(test_field[0], axis=0), np.sum(test_field[1], axis=0), np.sum(test_field[2], axis=0)

    ref_mag = magnitude_squared(*ref_sum)
    test_mag = magnitude_squared(*test_sum)
    diff_mag = np.abs(test_mag - ref_mag)

    R_edges = edges_from_centres(R)
    Z_edges = edges_from_centres(Z)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), layout='constrained')
    entries = [
        ("Reference", ref_mag, 'magma'),
        ("Test", test_mag, 'magma'),
        ("|Test - Ref|", diff_mag, 'viridis'),
    ]
    for ax, (title, data, cmap) in zip(axes, entries):
        pcm = ax.pcolormesh(R_edges, Z_edges, np.log10(np.maximum(data.T, 1e-30)), cmap=cmap)
        ax.set_title(title)
        ax.set_xlabel('R [cm]')
        ax.set_ylabel('Z [cm]')
        fig.colorbar(pcm, ax=ax, shrink=0.8, label=r'$\log_{10}$ value')

    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, dpi=dpi)
    if show:
        plt.show()
    else:
        plt.close(fig)


def edges_from_centres(values: np.ndarray) -> np.ndarray:
    edges = np.empty(values.size + 1, dtype=values.dtype)
    edges[1:-1] = 0.5 * (values[:-1] + values[1:])
    edges[0] = values[0] - (edges[1] - values[0])
    edges[-1] = values[-1] + (values[-1] - edges[-2])
    return edges


def _gauge_anvac_with_spline(
    grid,
    AnR: np.ndarray,
    Anphi: np.ndarray,
    AnZ: np.ndarray,
    ntor: int,
) -> Tuple[np.ndarray, np.ndarray]:
    ncoil = AnR.shape[0]
    R_vals = np.asarray(grid.R)
    Z_vals = np.asarray(grid.Z)
    gauged_AnR = np.zeros_like(AnR)
    gauged_AnZ = np.zeros_like(AnZ)

    factor = 1j / ntor
    R_matrix = R_vals[:, np.newaxis]

    for k in range(ncoil):
        spl_Aphi_real = RectBivariateSpline(R_vals, Z_vals, Anphi[k].real, kx=5, ky=5)
        spl_Aphi_imag = RectBivariateSpline(R_vals, Z_vals, Anphi[k].imag, kx=5, ky=5)
        Aphi_eval = spl_Aphi_real(R_vals, Z_vals) + 1j * spl_Aphi_imag(R_vals, Z_vals)
        dAphi_dR = spl_Aphi_real(R_vals, Z_vals, dx=1, dy=0) + 1j * spl_Aphi_imag(R_vals, Z_vals, dx=1, dy=0)
        dAphi_dZ = spl_Aphi_real(R_vals, Z_vals, dx=0, dy=1) + 1j * spl_Aphi_imag(R_vals, Z_vals, dx=0, dy=1)

        spl_AnR_real = RectBivariateSpline(R_vals, Z_vals, AnR[k].real, kx=5, ky=5)
        spl_AnR_imag = RectBivariateSpline(R_vals, Z_vals, AnR[k].imag, kx=5, ky=5)
        AnR_eval = spl_AnR_real(R_vals, Z_vals) + 1j * spl_AnR_imag(R_vals, Z_vals)

        spl_AnZ_real = RectBivariateSpline(R_vals, Z_vals, AnZ[k].real, kx=5, ky=5)
        spl_AnZ_imag = RectBivariateSpline(R_vals, Z_vals, AnZ[k].imag, kx=5, ky=5)
        AnZ_eval = spl_AnZ_real(R_vals, Z_vals) + 1j * spl_AnZ_imag(R_vals, Z_vals)

        gauged_AnR[k] = AnR_eval + factor * (R_matrix * dAphi_dR + Aphi_eval)
        gauged_AnZ[k] = AnZ_eval + factor * R_matrix * dAphi_dZ

    return gauged_AnR, gauged_AnZ


def _field_divfree_rect(
    spl: dict,
    R_vals: np.ndarray,
    Z_vals: np.ndarray,
    ntor: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    if ntor == 0:
        raise ValueError("_field_divfree_rect requires ntor != 0")

    ncoil = len(spl['AnR_Re'])
    nR = R_vals.size
    nZ = Z_vals.size
    R_matrix = R_vals[:, np.newaxis]

    BnR = np.zeros((ncoil, nR, nZ), dtype=complex)
    Bnphi = np.zeros_like(BnR)
    BnZ = np.zeros_like(BnR)

    for k in range(ncoil):
        anR = spl['AnR_Re'][k](R_vals, Z_vals) + 1j * spl['AnR_Im'][k](R_vals, Z_vals)
        anZ = spl['AnZ_Re'][k](R_vals, Z_vals) + 1j * spl['AnZ_Im'][k](R_vals, Z_vals)

        dAnR_dZ = spl['AnR_Re'][k](R_vals, Z_vals, dy=1) + 1j * spl['AnR_Im'][k](R_vals, Z_vals, dy=1)
        dAnZ_dR = spl['AnZ_Re'][k](R_vals, Z_vals, dx=1) + 1j * spl['AnZ_Im'][k](R_vals, Z_vals, dx=1)

        BnR[k] = 1j * ntor * anZ / R_matrix
        Bnphi[k] = dAnR_dZ - dAnZ_dR
        BnZ[k] = -1j * ntor * anR / R_matrix

    return BnR, Bnphi, BnZ


def _compute_axisymmetric_field(
    grid,
    AnR: np.ndarray,
    Anphi: np.ndarray,
    AnZ: np.ndarray,
    _dAnphi_dR: np.ndarray,
    _dAnphi_dZ: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    ncoil = AnR.shape[0]
    BnR = np.zeros_like(AnR)
    Bnphi = np.zeros_like(AnR)
    BnZ = np.zeros_like(AnR)

    R_vals = np.asarray(grid.R)
    Z_vals = np.asarray(grid.Z)

    for k in range(ncoil):
        A_R = AnR[k]
        A_phi = Anphi[k]
        A_Z = AnZ[k]

        spline_Aphi_real = RectBivariateSpline(R_vals, Z_vals, A_phi.real, kx=5, ky=5)
        spline_Aphi_imag = RectBivariateSpline(R_vals, Z_vals, A_phi.imag, kx=5, ky=5)
        dAphi_dZ = spline_Aphi_real(R_vals, Z_vals, dy=1) + 1j * spline_Aphi_imag(R_vals, Z_vals, dy=1)
        dAphi_dR = spline_Aphi_real(R_vals, Z_vals, dx=1) + 1j * spline_Aphi_imag(R_vals, Z_vals, dx=1)

        BnR[k] = -dAphi_dZ
        with np.errstate(divide="ignore", invalid="ignore"):
            BnZ[k] = dAphi_dR + A_phi / R_vals[:, np.newaxis]

        spline_AR_real = RectBivariateSpline(R_vals, Z_vals, A_R.real, kx=5, ky=5)
        spline_AR_imag = RectBivariateSpline(R_vals, Z_vals, A_R.imag, kx=5, ky=5)
        dA_R_dZ = spline_AR_real(R_vals, Z_vals, dy=1) + 1j * spline_AR_imag(R_vals, Z_vals, dy=1)

        spline_AZ_real = RectBivariateSpline(R_vals, Z_vals, A_Z.real, kx=5, ky=5)
        spline_AZ_imag = RectBivariateSpline(R_vals, Z_vals, A_Z.imag, kx=5, ky=5)
        dA_Z_dR = spline_AZ_real(R_vals, Z_vals, dx=1) + 1j * spline_AZ_imag(R_vals, Z_vals, dx=1)

        Bnphi[k] = dA_R_dZ - dA_Z_dR

        if R_vals[0] == 0.0 and R_vals.size > 1:
            BnR[k, 0, :] = BnR[k, 1, :]
            Bnphi[k, 0, :] = Bnphi[k, 1, :]
            BnZ[k, 0, :] = BnZ[k, 1, :]

        invalid = ~np.isfinite(BnZ[k])
        if invalid.any():
            finite_vals = BnZ[k][np.isfinite(BnZ[k])]
            if finite_vals.size > 0:
                BnZ[k][invalid] = finite_vals[0]
            else:
                BnZ[k][invalid] = 0.0

    return BnR, Bnphi, BnZ


def _tensordot_currents(weights: np.ndarray, array: np.ndarray) -> np.ndarray:
    return np.tensordot(weights, array, axes=(0, 0))


def make_deriv_diff_plot(
    diagnostics: Dict[str, np.ndarray],
    weights: np.ndarray,
    R: np.ndarray,
    Z: np.ndarray,
    output: Path,
    dpi: int,
    show: bool,
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)

    R_edges = edges_from_centres(R)
    Z_edges = edges_from_centres(Z)

    analytic = (
        _tensordot_currents(weights, diagnostics['BnR_analytic']),
        _tensordot_currents(weights, diagnostics['Bnphi_analytic']),
        _tensordot_currents(weights, diagnostics['BnZ_analytic']),
    )
    spline = (
        _tensordot_currents(weights, diagnostics['BnR_spline']),
        _tensordot_currents(weights, diagnostics['Bnphi_spline']),
        _tensordot_currents(weights, diagnostics['BnZ_spline']),
    )

    diff = tuple(s - a for s, a in zip(spline, analytic))

    fig, axes = plt.subplots(3, 3, figsize=(14, 10), layout='constrained')
    entries = ['B_R', 'B_φ', 'B_Z']

    for row, (label, ana, spl, delta) in enumerate(zip(entries, analytic, spline, diff)):
        data_sets = [ana, spl, delta]
        titles = [f'{label} (analytic)', f'{label} (spline)', f'Δ{label}']
        for col, (data, title) in enumerate(zip(data_sets, titles)):
            ax = axes[row, col]
            if col < 2:
                values = np.maximum(np.abs(data).T, 1e-30)
                pcm = ax.pcolormesh(R_edges, Z_edges, np.log10(values), cmap='magma')
                cbar_label = r'$\log_{10}|B_n|$'
            else:
                values = np.maximum(np.abs(delta).T, 1e-30)
                pcm = ax.pcolormesh(R_edges, Z_edges, np.log10(values), cmap='viridis')
                cbar_label = r'$\log_{10}|ΔB_n|$'
            ax.set_title(title)
            ax.set_xlabel('R [cm]')
            ax.set_ylabel('Z [cm]')
            fig.colorbar(pcm, ax=ax, shrink=0.8, label=cbar_label)

    fig.savefig(output, dpi=dpi)
    if show:
        plt.show()
    else:
        plt.close(fig)


def main() -> None:
    args = parse_args()

    ref_field, test_field, R, Z, diagnostics = load_modes(args.reference, args.test, args.ntor)
    weights = load_currents(args.currents, ref_field[0].shape[0], args.prefactor)

    ref_weighted = tuple(component * weights[:, np.newaxis, np.newaxis] for component in ref_field)
    test_weighted = tuple(component * weights[:, np.newaxis, np.newaxis] for component in test_field)

    make_per_coil_plot(ref_weighted, test_weighted, R, Z, args.per_coil_output, args.dpi, args.show)
    make_sum_plot(ref_weighted, test_weighted, R, Z, args.sum_output, args.dpi, args.show)

    if diagnostics is not None and args.deriv_diff_output is not None:
        max_gauge_R = np.max(np.abs(diagnostics["gauged_diff_R"]))
        max_gauge_Z = np.max(np.abs(diagnostics["gauged_diff_Z"]))
        print(
            f"Gauge diagnostics: max |ΔA_R| = {max_gauge_R:.3e}, "
            f"max |ΔA_Z| = {max_gauge_Z:.3e}"
        )
        make_deriv_diff_plot(
            diagnostics,
            weights,
            R,
            Z,
            args.deriv_diff_output,
            args.dpi,
            args.show,
        )


if __name__ == "__main__":
    main()
