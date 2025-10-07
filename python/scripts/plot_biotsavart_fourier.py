#!/usr/bin/env python3
"""Unified Biot–Savart comparison and plotting tool.

This script combines the functionality of the previous
``compare_superposition.py`` and ``plot_biotsavart_fourier.py`` helpers:

* Loads Fourier-mode (Bn) data and vector-potential (An) data from coil_tools
  outputs and reconstructs the magnetic-field harmonics.
* Optionally evaluates a direct segment-based Biot–Savart reference using the
  GPEC coil geometry and currents.
* Computes statistics (median/mean relative error) between the Fourier,
  vector-potential, and direct fields.
* Generates diagnostic plots for per-coil magnitudes and a four-panel summed
  comparison (Fourier, vector, direct, coil geometry).
* Can validate the magnetic field along a specified axis against an analytic
  circular-coil solution.

Example (matches the AUG CTest invocation):

```
python3 plot_biotsavart_fourier.py \
    test/magfie/test_data/aug_reference.h5 \
    test/magfie/test_data/aug_test.nc \
    --currents test/magfie/test_data/aug_currents.txt \
    --coil-files test/magfie/test_data/aug_bu.dat test/magfie/test_data/aug_bl.dat \
    --ntor 2 \
    --comparison-output test/magfie/test_data/superposition_comparison.png \
    --per-coil-output test/magfie/test_data/biotsavart_fourier.png \
    --sum-output test/magfie/test_data/biotsavart_fourier_sum.png
```
"""

from __future__ import annotations

import argparse
import ctypes
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

from numpy import (
    abs,
    allclose,
    amin,
    amax,
    array,
    asarray,
    conj,
    cos,
    exp,
    full,
    hypot,
    arctan2,
    linspace,
    log10,
    maximum,
    mean,
    median,
    ndarray,
    newaxis,
    outer,
    real,
    sin,
    sqrt,
    searchsorted,
    stack,
    swapaxes,
    tensordot,
    empty,
    zeros,
    zeros_like,
    isnan,
    pi,
    meshgrid,
)

from numpy.fft import fft
from numpy.linalg import norm

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[2]
PYTHON_DIR = REPO_ROOT / "python"
BUILD_DIR = REPO_ROOT / "build"


def _find_site_packages(base: Path) -> Path | None:
    lib_dir = base / "lib"
    if not lib_dir.exists():
        return None
    for subdir in lib_dir.iterdir():
        if not subdir.is_dir():
            continue
        site_packages = subdir / "site-packages"
        if site_packages.exists():
            return site_packages
    return None


for candidate in (PYTHON_DIR, BUILD_DIR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

for venv_candidate in (REPO_ROOT / ".venv", REPO_ROOT.parent / ".venv"):
    if venv_candidate.exists():
        site_packages = _find_site_packages(venv_candidate)
        if site_packages is not None and str(site_packages) not in sys.path:
            sys.path.insert(0, str(site_packages))

from libneo.biotsavart_fourier import (  # type: ignore  # noqa: E402
    read_Bnvac_fourier,
    read_Anvac_fourier,
    gauged_Anvac_from_Bnvac,
    gauge_Anvac,
    spline_gauged_Anvac,
    field_divfree,
)
import _magfie  # type: ignore  # noqa: E402
from _magfie import (
    compute_field_from_gpec_file,
    field_divb0_initialize_from_grid,
    field_divb0_eval,
)


def _set_convex_wall_filename(path: Path) -> None:
    """Write ``path`` into the Fortran ``convexfile`` module variable."""

    lib = ctypes.CDLL(_magfie.__file__)
    buffer = (ctypes.c_char * 1024).in_dll(lib, "__input_files_MOD_convexfile")
    buffer.value = str(path).encode()


FOURIER_REFERENCE_CURRENT = 0.1  # coil_tools Fourier tables assume 0.1 A per coil
DEFAULT_FFT_SAMPLES = 256


@dataclass
class Grid:
    R: ndarray
    Z: ndarray
    nR: int
    nZ: int
    R_min: float
    R_max: float
    Z_min: float
    Z_max: float


@dataclass
class ModeData:
    grid: Grid
    BnR: ndarray  # shape (ncoil, nR, nZ)
    Bnphi: ndarray
    BnZ: ndarray


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare coil_tools Fourier, vector-potential, and direct Biot–Savart outputs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("reference", type=Path, help="coil_tools Fourier HDF5 file (Bn)")
    parser.add_argument("test", type=Path, help="coil_tools vector-potential NetCDF file (An)")
    parser.add_argument("--currents", type=Path, help="Text file with coil currents (ampere)")
    parser.add_argument(
        "--coil-files",
        type=Path,
        nargs="+",
        help="One or more coil geometry files in GPEC format (required for direct comparison)",
    )
    parser.add_argument("--ntor", type=int, default=2, help="Toroidal mode number to analyse")
    parser.add_argument(
        "--fft-samples",
        type=int,
        default=DEFAULT_FFT_SAMPLES,
        help="Number of toroidal samples for direct Biot–Savart FFT",
    )
    parser.add_argument("--prefactor", type=float, default=FOURIER_REFERENCE_CURRENT,
                        help="Reference current used when tabulating the Fourier kernels (A)")
    parser.add_argument("--per-coil-output", type=Path,
                        help="Filename for per-coil magnitude comparison plot")
    parser.add_argument("--sum-output", type=Path,
                        help="Filename for coil-summed magnitude plot (four-panel)")
    parser.add_argument("--deriv-diff-output", type=Path,
                        help="Filename for derivative comparison diagnostic plot")
    parser.add_argument("--axis-output", type=Path,
                        help="Filename for optional axis validation plot")
    parser.add_argument("--dpi", type=int, default=150, help="Figure DPI for saved plots")
    parser.add_argument("--show", action="store_true", help="Display figures interactively")

    axis_group = parser.add_argument_group("Axis validation")
    axis_group.add_argument("--axis-origin", type=float, nargs=3,
                            help="Axis origin (cm) for analytic comparison")
    axis_group.add_argument("--axis-normal", type=float, nargs=3,
                            help="Axis direction vector")
    axis_group.add_argument("--axis-range", type=float, default=100.0,
                            help="Half-length of axis segment (cm)")
    axis_group.add_argument("--axis-samples", type=int, default=200,
                            help="Number of samples along axis")
    axis_group.add_argument("--coil-radius", type=float,
                            help="Analytic circular coil radius (cm)")

    return parser.parse_args()


def _load_mode_from_bnvac(path: Path, ntor: int) -> ModeData:
    grid_raw, BnR, Bnphi, BnZ = read_Bnvac_fourier(str(path), ntor=ntor)
    grid = Grid(
        R=array(grid_raw.R, copy=True),
        Z=array(grid_raw.Z, copy=True),
        nR=grid_raw.nR,
        nZ=grid_raw.nZ,
        R_min=float(grid_raw.R_min),
        R_max=float(grid_raw.R_max),
        Z_min=float(grid_raw.Z_min),
        Z_max=float(grid_raw.Z_max),
    )
    return ModeData(grid=grid, BnR=BnR, Bnphi=Bnphi, BnZ=BnZ)


def _load_mode_from_anvac(path: Path, ntor: int) -> Tuple[ModeData, dict, ModeData, dict]:
    from scipy.interpolate import RectBivariateSpline

    grid_raw, AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ = read_Anvac_fourier(str(path), ntor=ntor)
    grid = Grid(
        R=array(grid_raw.R, copy=True),
        Z=array(grid_raw.Z, copy=True),
        nR=grid_raw.nR,
        nZ=grid_raw.nZ,
        R_min=float(grid_raw.R_min),
        R_max=float(grid_raw.R_max),
        Z_min=float(grid_raw.Z_min),
        Z_max=float(grid_raw.Z_max),
    )

    # Version 1: Use stored derivatives from NetCDF
    gauged_AnR, gauged_AnZ = gauge_Anvac(
        grid_raw, AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ, ntor=ntor
    )
    if ntor != 0:
        gauged_AnR = gauged_AnR - 0.5j * Anphi
    spl = spline_gauged_Anvac(grid_raw, gauged_AnR, gauged_AnZ, ntor=ntor)
    BnR, Bnphi, BnZ = field_divfree(spl, grid.R, grid.Z, ntor=ntor)

    # Version 2: Compute derivatives from RectBivariateSpline of Aphi
    ncoil = Anphi.shape[0]
    dAphi_dR_spline = zeros_like(Anphi)
    dAphi_dZ_spline = zeros_like(Anphi)

    for ic in range(ncoil):
        # Create spline of Aphi for this coil
        spl_Aphi_real = RectBivariateSpline(grid.R, grid.Z, Anphi[ic].real)
        spl_Aphi_imag = RectBivariateSpline(grid.R, grid.Z, Anphi[ic].imag)

        # Evaluate derivatives on grid - need grid=True to return 2D array
        dAphi_dR_spline[ic] = (spl_Aphi_real(grid.R, grid.Z, dx=1, dy=0, grid=True) +
                                1j * spl_Aphi_imag(grid.R, grid.Z, dx=1, dy=0, grid=True))
        dAphi_dZ_spline[ic] = (spl_Aphi_real(grid.R, grid.Z, dx=0, dy=1, grid=True) +
                                1j * spl_Aphi_imag(grid.R, grid.Z, dx=0, dy=1, grid=True))

    gauged_AnR_spline, gauged_AnZ_spline = gauge_Anvac(
        grid_raw, AnR, Anphi, AnZ, dAphi_dR_spline, dAphi_dZ_spline, ntor=ntor
    )
    if ntor != 0:
        gauged_AnR_spline = gauged_AnR_spline - 0.5j * Anphi
    spl_spline = spline_gauged_Anvac(grid_raw, gauged_AnR_spline, gauged_AnZ_spline, ntor=ntor)
    BnR_spline, Bnphi_spline, BnZ_spline = field_divfree(spl_spline, grid.R, grid.Z, ntor=ntor)

    diagnostics = {
        'dAphi_dR_stored': dAphi_dR,
        'dAphi_dZ_stored': dAphi_dZ,
        'dAphi_dR_spline': dAphi_dR_spline,
        'dAphi_dZ_spline': dAphi_dZ_spline,
    }

    mode_stored = ModeData(grid=grid, BnR=BnR, Bnphi=Bnphi, BnZ=BnZ)
    mode_spline = ModeData(grid=grid, BnR=BnR_spline, Bnphi=Bnphi_spline, BnZ=BnZ_spline)

    return mode_stored, spl, mode_spline, diagnostics


def _ensure_same_grid(lhs: Grid, rhs: Grid) -> None:
    if lhs.nR != rhs.nR or lhs.nZ != rhs.nZ:
        raise ValueError("Grid dimensions differ between inputs")
    if not allclose(lhs.R, rhs.R) or not allclose(lhs.Z, rhs.Z):
        raise ValueError("Reference and test grids do not match")


def _read_currents(path: Path) -> ndarray:
    values = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            values.extend(float(tok) for tok in stripped.split())
    return asarray(values, dtype=float)


def _load_coil_projections(filenames: Sequence[Path]) -> List[Tuple[ndarray, ndarray]]:
    projections: List[Tuple[ndarray, ndarray]] = []
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
                R_path = hypot(coords_cm[:, 0], coords_cm[:, 1])
                Z_path = coords_cm[:, 2]
                projections.append((R_path, Z_path))
    return projections


def _read_gpec_header(path: Path) -> int:
    with open(path, "r", encoding="utf-8") as handle:
        header = handle.readline().split()
    if len(header) < 4:
        raise ValueError(f"Malformed GPEC header in {path}")
    return int(header[0])


def _compute_segment_direct_field(
    coil_files: Sequence[Path],
    coil_currents: ndarray,
    x_eval: ndarray,
    y_eval: ndarray,
    z_eval: ndarray,
) -> Tuple[ndarray, ndarray, ndarray]:
    Bx_total = zeros_like(x_eval, dtype=float)
    By_total = zeros_like(y_eval, dtype=float)
    Bz_total = zeros_like(z_eval, dtype=float)

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
            str(coil_file), currents_slice, x_eval, y_eval, z_eval
        )
        Bx_total += Bx
        By_total += By
        Bz_total += Bz

    return Bx_total, By_total, Bz_total


def _compute_direct_fourier_mode(
    coil_files: Sequence[Path],
    coil_currents: ndarray,
    grid: Grid,
    ntor: int,
    n_phi: int,
) -> Tuple[ndarray, ndarray, ndarray, ndarray, ndarray, ndarray, ndarray]:
    R_grid = grid.R
    Z_grid = grid.Z
    nR, nZ = grid.nR, grid.nZ

    phi_grid = linspace(0.0, 2.0 * pi, n_phi, endpoint=False)
    BR_all = zeros((nR, nZ, n_phi))
    Bphi_all = zeros_like(BR_all)
    BZ_all = zeros_like(BR_all)

    R_mesh, Z_mesh = meshgrid(R_grid, Z_grid, indexing="ij")
    base_x = R_mesh.reshape(-1)
    base_z = Z_mesh.reshape(-1)

    for idx, phi in enumerate(phi_grid):
        cosphi = cos(phi)
        sinphi = sin(phi)
        x_eval = base_x * cosphi
        y_eval = base_x * sinphi
        z_eval = base_z

        Bx, By, Bz = _compute_segment_direct_field(
            coil_files, coil_currents, x_eval, y_eval, z_eval
        )

        Bx = Bx.reshape(nR, nZ)
        By = By.reshape(nR, nZ)
        Bz = Bz.reshape(nR, nZ)

        BR_all[:, :, idx] = Bx * cosphi + By * sinphi
        Bphi_all[:, :, idx] = -Bx * sinphi + By * cosphi
        BZ_all[:, :, idx] = Bz

    fft_norm = 1.0 / n_phi
    BnR = fft(BR_all, axis=2)[:, :, ntor] * fft_norm
    Bnphi = fft(Bphi_all, axis=2)[:, :, ntor] * fft_norm
    BnZ = fft(BZ_all, axis=2)[:, :, ntor] * fft_norm

    return BnR, Bnphi, BnZ, BR_all, Bphi_all, BZ_all, phi_grid


def _magnitude_squared(BnR: ndarray, Bnphi: ndarray, BnZ: ndarray) -> ndarray:
    return (BnR * conj(BnR) + Bnphi * conj(Bnphi) + BnZ * conj(BnZ)).real


def _tensordot_currents(currents: ndarray, array: ndarray) -> ndarray:
    if array.shape[0] != currents.size:
        # Single coil case - no coil dimension, just multiply by single current
        return currents[0] * array
    return tensordot(currents, array, axes=(0, 0))


def _summed_magnitude(
    currents: ndarray,
    BnR: ndarray,
    Bnphi: ndarray,
    BnZ: ndarray,
) -> ndarray:
    sum_BnR = _tensordot_currents(currents, BnR)
    sum_Bnphi = _tensordot_currents(currents, Bnphi)
    sum_BnZ = _tensordot_currents(currents, BnZ)
    return _magnitude_squared(sum_BnR, sum_Bnphi, sum_BnZ)


def _bilinear_interpolate(
    R_grid: ndarray,
    Z_grid: ndarray,
    field: ndarray,
    R_val: float,
    Z_val: float,
) -> float:
    if R_val < R_grid[0] or R_val > R_grid[-1] or Z_val < Z_grid[0] or Z_val > Z_grid[-1]:
        return float("nan")

    i_hi = searchsorted(R_grid, R_val)
    if i_hi == 0:
        i0 = i1 = 0
        t = 0.0
    elif i_hi >= R_grid.size:
        i0 = i1 = R_grid.size - 1
        t = 0.0
    else:
        i0, i1 = i_hi - 1, i_hi
        denom = R_grid[i1] - R_grid[i0]
        t = 0.0 if denom == 0.0 else (R_val - R_grid[i0]) / denom

    j_hi = searchsorted(Z_grid, Z_val)
    if j_hi == 0:
        j0 = j1 = 0
        u = 0.0
    elif j_hi >= Z_grid.size:
        j0 = j1 = Z_grid.size - 1
        u = 0.0
    else:
        j0, j1 = j_hi - 1, j_hi
        denom = Z_grid[j1] - Z_grid[j0]
        u = 0.0 if denom == 0.0 else (Z_val - Z_grid[j0]) / denom

    f00 = field[i0, j0]
    f01 = field[i0, j1]
    f10 = field[i1, j0]
    f11 = field[i1, j1]

    return ((1.0 - t) * (1.0 - u) * f00 +
            (1.0 - t) * u * f01 +
            t * (1.0 - u) * f10 +
            t * u * f11)


def _evaluate_modes_at_point(
    BnR_modes: ndarray,
    Bnphi_modes: ndarray,
    BnZ_modes: ndarray,
    ntor: int,
    R_grid: ndarray,
    Z_grid: ndarray,
    x_val: float,
    y_val: float,
    z_val: float,
) -> Tuple[float, float, float]:
    R_val = hypot(x_val, y_val)
    phi_val = arctan2(y_val, x_val)

    BR = _bilinear_interpolate(R_grid, Z_grid, BnR_modes, R_val, z_val)
    Bphi = _bilinear_interpolate(R_grid, Z_grid, Bnphi_modes, R_val, z_val)
    BZ = _bilinear_interpolate(R_grid, Z_grid, BnZ_modes, R_val, z_val)

    if isnan(BR) or isnan(Bphi) or isnan(BZ):
        return float("nan"), float("nan"), float("nan")

    if ntor == 0:
        BR_total = BR.real
        Bphi_total = Bphi.real
        BZ_total = BZ.real
    else:
        factor = exp(1j * ntor * phi_val)
        BR_total = 2.0 * real(BR * factor)
        Bphi_total = 2.0 * real(Bphi * factor)
        BZ_total = 2.0 * real(BZ * factor)

    cosphi = cos(phi_val)
    sinphi = sin(phi_val)
    Bx = BR_total * cosphi - Bphi_total * sinphi
    By = BR_total * sinphi + Bphi_total * cosphi
    return Bx, By, BZ_total


def _axis_validation(
    args: argparse.Namespace,
    coil_files: Sequence[Path],
    coil_currents_amp: ndarray,
    grid: Grid,
    BnR_fourier: ndarray,
    Bnphi_fourier: ndarray,
    BnZ_fourier: ndarray,
    BnR_vector: ndarray,
    Bnphi_vector: ndarray,
    BnZ_vector: ndarray,
    BnR_direct: ndarray | None,
    Bnphi_direct: ndarray | None,
    BnZ_direct: ndarray | None,
) -> None:
    if args.axis_origin is None or args.axis_normal is None or args.coil_radius is None:
        return

    origin = asarray(args.axis_origin, dtype=float)
    normal = asarray(args.axis_normal, dtype=float)
    norm_vec = norm(normal)
    if norm_vec == 0.0:
        raise ValueError("Axis normal vector must be non-zero")
    direction = normal / norm_vec

    s_vals = linspace(-args.axis_range, args.axis_range, args.axis_samples)
    points = origin[newaxis, :] + outer(s_vals, direction)

    x_eval = points[:, 0]
    y_eval = points[:, 1]
    z_eval = points[:, 2]

    B_parallel_fourier = zeros_like(s_vals)
    B_parallel_vector = zeros_like(s_vals)
    B_parallel_direct = zeros_like(s_vals)

    for idx, (x_val, y_val, z_val) in enumerate(zip(x_eval, y_eval, z_eval)):
        Bx_f, By_f, Bz_f = _evaluate_modes_at_point(
            BnR_fourier, Bnphi_fourier, BnZ_fourier,
            args.ntor, grid.R, grid.Z,
            x_val, y_val, z_val,
        )
        B_parallel_fourier[idx] = Bx_f * direction[0] + By_f * direction[1] + Bz_f * direction[2]

        Bx_v, By_v, Bz_v = _evaluate_modes_at_point(
            BnR_vector, Bnphi_vector, BnZ_vector,
            args.ntor, grid.R, grid.Z,
            x_val, y_val, z_val,
        )
        B_parallel_vector[idx] = Bx_v * direction[0] + By_v * direction[1] + Bz_v * direction[2]

        if BnR_direct is not None and Bnphi_direct is not None and BnZ_direct is not None:
            Bx_d, By_d, Bz_d = _evaluate_modes_at_point(
                BnR_direct, Bnphi_direct, BnZ_direct,
                args.ntor, grid.R, grid.Z,
                x_val, y_val, z_val,
            )
            B_parallel_direct[idx] = Bx_d * direction[0] + By_d * direction[1] + Bz_d * direction[2]

    if BnR_direct is None or Bnphi_direct is None or BnZ_direct is None:
        Bx_raw, By_raw, Bz_raw = _compute_segment_direct_field(coil_files, coil_currents_amp, x_eval, y_eval, z_eval)
        B_parallel_direct = Bx_raw * direction[0] + By_raw * direction[1] + Bz_raw * direction[2]

    mu0 = 4e-7 * pi
    radius_m = args.coil_radius / 100.0
    s_vals_m = s_vals / 100.0
    current_ampere = sum(coil_currents_amp)
    B_analytic = mu0 * current_ampere * radius_m**2 / (2.0 * (radius_m**2 + s_vals_m**2)**1.5)
    B_analytic *= 1.0e4  # Tesla -> Gauss

    abs_err_fourier = abs(B_parallel_fourier - B_analytic)
    abs_err_vector = abs(B_parallel_vector - B_analytic)
    abs_err_direct = abs(B_parallel_direct - B_analytic)

    rel_err_fourier = abs_err_fourier / maximum(abs(B_analytic), 1e-20)
    rel_err_vector = abs_err_vector / maximum(abs(B_analytic), 1e-20)
    rel_err_direct = abs_err_direct / maximum(abs(B_analytic), 1e-20)

    if args.axis_output is not None:
        fig, ax = plt.subplots(figsize=(8, 4), layout="constrained")
        ax.plot(s_vals, B_analytic, label="Analytic", linestyle="--", linewidth=1.6)
        ax.plot(s_vals, B_parallel_fourier, label="Fourier", linewidth=1.2)
        ax.plot(s_vals, B_parallel_vector, label="Vector", linewidth=1.2)
        ax.plot(s_vals, B_parallel_direct, label="Direct", linewidth=1.2)
        ax.set_xlabel("Axis coordinate s [cm]")
        ax.set_ylabel("B_parallel [Gauss]")
        ax.legend(loc="best")
        fig.savefig(args.axis_output, dpi=args.dpi, bbox_inches="tight")
        if not args.show:
            plt.close(fig)

    print("Axis validation (max relative errors):")
    print(f"  Fourier vs analytic: {amax(rel_err_fourier) * 100.0:.4f}%")
    print(f"  Vector  vs analytic: {amax(rel_err_vector) * 100.0:.4f}%")
    print(f"  Direct  vs analytic: {amax(rel_err_direct) * 100.0:.4f}%")


def _create_per_coil_plot(
    args: argparse.Namespace,
    mode_fourier: ModeData,
    mode_vector: ModeData,
    mode_vector_spline: ModeData,
    BnR_spline_all: ndarray,
    Bnphi_spline_all: ndarray,
    BnZ_spline_all: ndarray,
    output: Path,
    coil_projections: Sequence[Tuple[ndarray, ndarray]] | None,
) -> None:
    ncoil = mode_fourier.BnR.shape[0]
    grid = mode_fourier.grid
    titles = ("reference", "Anvac (stored)", "Anvac (spline)", "Bnvac")

    log_bn2 = empty((4, ncoil, grid.nR, grid.nZ), dtype=float)
    log_bn2[0] = log10(maximum(_magnitude_squared(mode_fourier.BnR, mode_fourier.Bnphi, mode_fourier.BnZ), 1e-300))
    log_bn2[1] = log10(maximum(_magnitude_squared(mode_vector.BnR, mode_vector.Bnphi, mode_vector.BnZ), 1e-300))
    log_bn2[2] = log10(maximum(_magnitude_squared(mode_vector_spline.BnR, mode_vector_spline.Bnphi, mode_vector_spline.BnZ), 1e-300))
    log_bn2[3] = log10(maximum(_magnitude_squared(BnR_spline_all, Bnphi_spline_all, BnZ_spline_all), 1e-300))

    norm = Normalize(vmin=amin(log_bn2), vmax=amax(log_bn2))
    fig, axs = plt.subplots(ncoil, 4, figsize=(12, 3 * ncoil), layout="constrained")

    for kcoil in range(ncoil):
        for col, label in enumerate(titles):
            ax = axs[kcoil, col] if ncoil > 1 else axs[col]
            im = ax.imshow(
                log_bn2[col][kcoil].T,
                origin="lower",
                cmap="magma",
                extent=[grid.R_min, grid.R_max, grid.Z_min, grid.Z_max],
                norm=norm,
                interpolation="bilinear",
            )
            ax.set_title(f"Coil {kcoil + 1}: {label}")
            ax.set_xlabel("R [cm]")
            ax.set_ylabel("Z [cm]")
            if coil_projections:
                for R_path, Z_path in coil_projections:
                    ax.plot(R_path, Z_path, color="black", linewidth=0.8, alpha=0.6)
            ax.set_aspect("equal", adjustable="box")

    cbar = fig.colorbar(im, ax=axs, location="bottom", fraction=0.05, pad=0.08)
    cbar.set_label("log10 |B_n|^2")
    fig.savefig(output, dpi=args.dpi, bbox_inches="tight")
    if not args.show:
        plt.close(fig)


def _create_sum_plot(
    args: argparse.Namespace,
    magnitude_fourier: ndarray,
    magnitude_vector: ndarray,
    magnitude_spline: ndarray,
    magnitude_field_divb0: ndarray | None,
    magnitude_direct: ndarray | None,
    output: Path,
    grid: Grid,
    coil_projections: Sequence[Tuple[ndarray, ndarray]] | None,
) -> None:
    top_titles = ("reference sum", "Anvac sum", "Bnvac sum", "Direct sum")
    top_panels = [magnitude_fourier, magnitude_vector, magnitude_spline]
    if magnitude_direct is None:
        top_panels.append(magnitude_fourier)
    else:
        top_panels.append(magnitude_direct)

    bottom_titles = ("field_divB0 sum", "", "", "")
    bottom_panel = magnitude_field_divb0

    accum = [log10(maximum(panel, 1e-300)) for panel in top_panels]
    if bottom_panel is not None:
        accum.append(log10(maximum(bottom_panel, 1e-300)))

    norm = Normalize(vmin=min(a.min() for a in accum), vmax=max(a.max() for a in accum))

    fig, axs = plt.subplots(2, 4, figsize=(14, 7.5), layout="constrained")
    extent = [grid.R_min, grid.R_max, grid.Z_min, grid.Z_max]
    plotted_axes = []

    for idx, label in enumerate(top_titles):
        ax = axs[0, idx]
        im = ax.imshow(
            log10(maximum(top_panels[idx], 1e-300)).T,
            origin="lower",
            cmap="magma",
            extent=extent,
            norm=norm,
            interpolation="bilinear",
        )
        ax.set_title(label)
        ax.set_xlabel("R [cm]")
        ax.set_ylabel("Z [cm]")
        if coil_projections:
            for R_path, Z_path in coil_projections:
                ax.plot(R_path, Z_path, color="black", linewidth=0.8, alpha=0.6)
        ax.set_aspect("equal", adjustable="box")
        plotted_axes.append(ax)

    if bottom_panel is not None:
        ax = axs[1, 0]
        im_bottom = ax.imshow(
            log10(maximum(bottom_panel, 1e-300)).T,
            origin="lower",
            cmap="magma",
            extent=extent,
            norm=norm,
            interpolation="bilinear",
        )
        ax.set_title(bottom_titles[0])
        ax.set_xlabel("R [cm]")
        ax.set_ylabel("Z [cm]")
        if coil_projections:
            for R_path, Z_path in coil_projections:
                ax.plot(R_path, Z_path, color="black", linewidth=0.8, alpha=0.6)
        ax.set_aspect("equal", adjustable="box")
        plotted_axes.append(ax)
        im = im_bottom
    else:
        ax = axs[1, 0]
        ax.axis("off")

    for idx in range(1, 4):
        axs[1, idx].axis("off")

    cbar = fig.colorbar(im, ax=plotted_axes, location="bottom", fraction=0.05, pad=0.08)
    cbar.set_label("log10 |Σ B_n|^2")
    fig.savefig(output, dpi=args.dpi, bbox_inches="tight")
    if not args.show:
        plt.close(fig)
def _create_deriv_comparison_plot(
    args: argparse.Namespace,
    diagnostics: dict,
    grid: Grid,
    currents: ndarray,
    output: Path,
) -> None:
    """Plot comparison of stored vs spline-derived dAphi/dR and dAphi/dZ."""
    dAphi_dR_stored = diagnostics['dAphi_dR_stored']
    dAphi_dZ_stored = diagnostics['dAphi_dZ_stored']
    dAphi_dR_spline = diagnostics['dAphi_dR_spline']
    dAphi_dZ_spline = diagnostics['dAphi_dZ_spline']

    # Weight by currents and sum over coils
    weights = currents * args.prefactor
    dAphi_dR_stored_sum = _tensordot_currents(weights, dAphi_dR_stored)
    dAphi_dZ_stored_sum = _tensordot_currents(weights, dAphi_dZ_stored)
    dAphi_dR_spline_sum = _tensordot_currents(weights, dAphi_dR_spline)
    dAphi_dZ_spline_sum = _tensordot_currents(weights, dAphi_dZ_spline)

    # Compute magnitude
    mag_dR_stored = abs(dAphi_dR_stored_sum)
    mag_dZ_stored = abs(dAphi_dZ_stored_sum)
    mag_dR_spline = abs(dAphi_dR_spline_sum)
    mag_dZ_spline = abs(dAphi_dZ_spline_sum)

    # Create 2x2 plot - just the values side by side
    fig, axs = plt.subplots(2, 2, figsize=(12, 10), layout="constrained")
    extent = [grid.R_min, grid.R_max, grid.Z_min, grid.Z_max]

    # Find common colorbar limits for each derivative
    vmin_dR = min(amin(mag_dR_stored), amin(mag_dR_spline))
    vmax_dR = max(amax(mag_dR_stored), amax(mag_dR_spline))
    vmin_dZ = min(amin(mag_dZ_stored), amin(mag_dZ_spline))
    vmax_dZ = max(amax(mag_dZ_stored), amax(mag_dZ_spline))

    # dAphi/dR stored
    im0 = axs[0, 0].imshow(mag_dR_stored.T, origin="lower", cmap="viridis",
                            extent=extent, interpolation="bilinear",
                            vmin=vmin_dR, vmax=vmax_dR)
    axs[0, 0].set_title("dAphi/dR (stored)")
    axs[0, 0].set_xlabel("R [cm]")
    axs[0, 0].set_ylabel("Z [cm]")
    fig.colorbar(im0, ax=axs[0, 0])

    # dAphi/dR spline
    im1 = axs[0, 1].imshow(mag_dR_spline.T, origin="lower", cmap="viridis",
                            extent=extent, interpolation="bilinear",
                            vmin=vmin_dR, vmax=vmax_dR)
    axs[0, 1].set_title("dAphi/dR (spline)")
    axs[0, 1].set_xlabel("R [cm]")
    axs[0, 1].set_ylabel("Z [cm]")
    fig.colorbar(im1, ax=axs[0, 1])

    # dAphi/dZ stored
    im2 = axs[1, 0].imshow(mag_dZ_stored.T, origin="lower", cmap="viridis",
                            extent=extent, interpolation="bilinear",
                            vmin=vmin_dZ, vmax=vmax_dZ)
    axs[1, 0].set_title("dAphi/dZ (stored)")
    axs[1, 0].set_xlabel("R [cm]")
    axs[1, 0].set_ylabel("Z [cm]")
    fig.colorbar(im2, ax=axs[1, 0])

    # dAphi/dZ spline
    im3 = axs[1, 1].imshow(mag_dZ_spline.T, origin="lower", cmap="viridis",
                            extent=extent, interpolation="bilinear",
                            vmin=vmin_dZ, vmax=vmax_dZ)
    axs[1, 1].set_title("dAphi/dZ (spline)")
    axs[1, 1].set_xlabel("R [cm]")
    axs[1, 1].set_ylabel("Z [cm]")
    fig.colorbar(im3, ax=axs[1, 1])

    fig.savefig(output, dpi=args.dpi, bbox_inches="tight")
    if not args.show:
        plt.close(fig)


def _print_stats(
    magnitude_fourier: ndarray,
    magnitude_vector: ndarray,
    magnitude_spline: ndarray,
    magnitude_field_divb0: ndarray | None,
    magnitude_direct: ndarray | None,
) -> None:
    rel_err_vector = abs(magnitude_vector - magnitude_fourier) / (magnitude_fourier + 1e-15)
    rel_err_spline = abs(magnitude_spline - magnitude_fourier) / (magnitude_fourier + 1e-15)

    median_vector = float(median(rel_err_vector) * 100.0)
    mean_vector = float(mean(rel_err_vector) * 100.0)
    median_spline = float(median(rel_err_spline) * 100.0)
    mean_spline = float(mean(rel_err_spline) * 100.0)

    print(f"Median relative error (vector): {median_vector:.4f}%")
    print(f"Mean relative error (vector):   {mean_vector:.4f}%")
    print(f"Median relative error (spline): {median_spline:.4f}%")
    print(f"Mean relative error (spline):   {mean_spline:.4f}%")

    if magnitude_field_divb0 is not None:
        rel_err_field = abs(magnitude_field_divb0 - magnitude_fourier) / (magnitude_fourier + 1e-15)
        median_field = float(median(rel_err_field) * 100.0)
        mean_field = float(mean(rel_err_field) * 100.0)
        print(f"Median relative error (field_divB0): {median_field:.4f}%")
        print(f"Mean relative error (field_divB0):   {mean_field:.4f}%")

    if magnitude_direct is None:
        return

    rel_err_direct = abs(magnitude_direct - magnitude_fourier) / (magnitude_fourier + 1e-15)
    median_direct = float(median(rel_err_direct) * 100.0)
    mean_direct = float(mean(rel_err_direct) * 100.0)

    print("Direct vs Fourier:")
    print(f"  Median relative error: {median_direct:.4f}%")
    print(f"  Mean relative error:   {mean_direct:.4f}%")


def main() -> None:
    args = _parse_args()

    mode_fourier = _load_mode_from_bnvac(args.reference, args.ntor)
    mode_vector, spline_vec, mode_vector_spline, diagnostics = _load_mode_from_anvac(args.test, args.ntor)
    _ensure_same_grid(mode_fourier.grid, mode_vector.grid)

    coil_projections = None
    if args.coil_files:
        coil_projections = _load_coil_projections(args.coil_files)

    if args.currents is None:
        raise ValueError("--currents is required for comparison")
    if not args.coil_files:
        raise ValueError("--coil-files is required for comparison")

    currents = _read_currents(args.currents)
    coil_counts = [_read_gpec_header(path) for path in args.coil_files]
    if sum(coil_counts) != currents.size:
        raise ValueError(
            f"Currents file provides {currents.size} entries but geometry defines {sum(coil_counts)} coils"
        )
    if mode_fourier.BnR.shape[0] != currents.size:
        raise ValueError(
            f"Mismatch between currents ({currents.size}) and Fourier data ({mode_fourier.BnR.shape[0]})"
        )

    BnR_spline_all, Bnphi_spline_all, BnZ_spline_all = field_divfree(
        spline_vec, mode_vector.grid.R, mode_vector.grid.Z, ntor=args.ntor
    )

    weights = currents * args.prefactor
    BnR_fourier_sum = _tensordot_currents(weights, mode_fourier.BnR)
    Bnphi_fourier_sum = _tensordot_currents(weights, mode_fourier.Bnphi)
    BnZ_fourier_sum = _tensordot_currents(weights, mode_fourier.BnZ)

    BnR_vector_sum = _tensordot_currents(weights, mode_vector.BnR)
    Bnphi_vector_sum = _tensordot_currents(weights, mode_vector.Bnphi)
    BnZ_vector_sum = _tensordot_currents(weights, mode_vector.BnZ)

    BnR_spline_sum = _tensordot_currents(weights, BnR_spline_all)
    Bnphi_spline_sum = _tensordot_currents(weights, Bnphi_spline_all)
    BnZ_spline_sum = _tensordot_currents(weights, BnZ_spline_all)

    BnR_direct = Bnphi_direct = BnZ_direct = None
    direct_magnitude = None
    magnitude_field_divb0 = None
    if args.coil_files:
        (BnR_direct, Bnphi_direct, BnZ_direct,
         BR_samples, Bphi_samples, BZ_samples, phi_grid) = _compute_direct_fourier_mode(
            args.coil_files, currents, mode_fourier.grid, args.ntor, args.fft_samples
        )
        direct_magnitude = _magnitude_squared(BnR_direct, Bnphi_direct, BnZ_direct)

        nR = int(mode_fourier.grid.nR)
        nZ = int(mode_fourier.grid.nZ)
        n_phi = phi_grid.size

        br_ordered = swapaxes(BR_samples, 1, 2)
        bp_ordered = swapaxes(Bphi_samples, 1, 2)
        bz_ordered = swapaxes(BZ_samples, 1, 2)

        n_phi_ext = n_phi + 1
        br_ext = zeros((nR, n_phi_ext, nZ), dtype=float, order="F")
        bp_ext = zeros((nR, n_phi_ext, nZ), dtype=float, order="F")
        bz_ext = zeros((nR, n_phi_ext, nZ), dtype=float, order="F")
        br_ext[:, :n_phi, :] = br_ordered
        bp_ext[:, :n_phi, :] = bp_ordered
        bz_ext[:, :n_phi, :] = bz_ordered
        br_ext[:, n_phi, :] = br_ordered[:, 0, :]
        bp_ext[:, n_phi, :] = bp_ordered[:, 0, :]
        bz_ext[:, n_phi, :] = bz_ordered[:, 0, :]

        convex_candidate: Path | None = None
        for coil_file in args.coil_files:
            candidate = coil_file.parent / "convexwall.dat"
            if candidate.exists():
                convex_candidate = candidate
                break

        if convex_candidate is not None:
            _set_convex_wall_filename(convex_candidate)
        else:
            raise FileNotFoundError(
                "convexwall.dat not found alongside coil geometry; required for field_divB0"
            )

        field_divb0_initialize_from_grid(
            args.ntor,
            mode_fourier.grid.R_min,
            mode_fourier.grid.R_max,
            0.0,
            2.0 * pi,
            mode_fourier.grid.Z_min,
            mode_fourier.grid.Z_max,
            br_ext,
            bp_ext,
            bz_ext,
            nR,
            n_phi_ext,
            nZ,
        )

        R_mesh, Z_mesh = meshgrid(mode_fourier.grid.R, mode_fourier.grid.Z, indexing="ij")
        r_flat = R_mesh.reshape(-1)
        z_flat = Z_mesh.reshape(-1)
        n_points = r_flat.size

        BR_field = zeros((nR, nZ, n_phi))
        Bphi_field = zeros_like(BR_field)
        BZ_field = zeros_like(BR_field)

        phi_buffer = zeros(n_points)

        for idx in range(n_phi):
            phi_buffer[:] = phi_grid[idx]
            br_slice, bp_slice, bz_slice = field_divb0_eval(r_flat, phi_buffer, z_flat, n_points)
            BR_field[:, :, idx] = br_slice.reshape(nR, nZ)
            Bphi_field[:, :, idx] = bp_slice.reshape(nR, nZ)
            BZ_field[:, :, idx] = bz_slice.reshape(nR, nZ)

        fft_norm = 1.0 / n_phi
        BnR_field = fft(BR_field, axis=2)[:, :, args.ntor] * fft_norm
        Bnphi_field = fft(Bphi_field, axis=2)[:, :, args.ntor] * fft_norm
        BnZ_field = fft(BZ_field, axis=2)[:, :, args.ntor] * fft_norm

        magnitude_field_divb0 = _magnitude_squared(BnR_field, Bnphi_field, BnZ_field)

    magnitude_fourier = _magnitude_squared(BnR_fourier_sum, Bnphi_fourier_sum, BnZ_fourier_sum)
    magnitude_vector = _magnitude_squared(BnR_vector_sum, Bnphi_vector_sum, BnZ_vector_sum)
    magnitude_spline = _magnitude_squared(BnR_spline_sum, Bnphi_spline_sum, BnZ_spline_sum)

    if args.per_coil_output is not None:
        _create_per_coil_plot(args, mode_fourier, mode_vector, mode_vector_spline,
                              BnR_spline_all, Bnphi_spline_all, BnZ_spline_all,
                              args.per_coil_output, coil_projections)

    if args.deriv_diff_output is not None:
        _create_deriv_comparison_plot(args, diagnostics, mode_fourier.grid, currents,
                                      args.deriv_diff_output)

    if args.sum_output is not None:
        _create_sum_plot(
            args,
            magnitude_fourier,
            magnitude_vector,
            magnitude_spline,
            magnitude_field_divb0,
            direct_magnitude,
            args.sum_output,
            mode_fourier.grid,
            coil_projections,
        )

    _print_stats(
        magnitude_fourier,
        magnitude_vector,
        magnitude_spline,
        magnitude_field_divb0,
        direct_magnitude,
    )

    if args.axis_origin is not None and args.axis_normal is not None and args.coil_radius is not None:
        _axis_validation(
            args,
            args.coil_files,
            currents,
            mode_fourier.grid,
            BnR_fourier_sum,
            Bnphi_fourier_sum,
            BnZ_fourier_sum,
            BnR_vector_sum,
            Bnphi_vector_sum,
            BnZ_vector_sum,
            BnR_direct,
            Bnphi_direct,
            BnZ_direct,
        )


if __name__ == "__main__":
    main()
