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
    isfinite,
    errstate,
    pi,
    meshgrid,
)

from numpy.fft import fft
from numpy.linalg import norm

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from scipy.interpolate import RectBivariateSpline

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
                        help="Filename for diagnostic plot comparing spline vs analytic gauge")
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


def _load_mode_from_anvac(path: Path, ntor: int) -> Tuple[ModeData, dict | None, dict | None]:
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

    if ntor == 0:
        BnR, Bnphi, BnZ = _compute_axisymmetric_field(
            grid, AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ
        )
        return ModeData(grid=grid, BnR=BnR, Bnphi=Bnphi, BnZ=BnZ), None, None

    gauged_spline_R, gauged_spline_Z = _gauge_anvac_with_spline(grid, AnR, Anphi, AnZ, ntor)
    spl_spline = spline_gauged_Anvac(grid_raw, gauged_spline_R, gauged_spline_Z, ntor=ntor)
    BnR_spline, Bnphi_spline, BnZ_spline = _field_divfree_rect(
        spl_spline, grid.R, grid.Z, ntor=ntor
    )

    gauged_analytic_R, gauged_analytic_Z = gauge_Anvac(
        grid_raw, AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ, ntor=ntor
    )
    gauged_analytic_R = gauged_analytic_R - 0.5j * Anphi
    spl_analytic = spline_gauged_Anvac(
        grid_raw, gauged_analytic_R, gauged_analytic_Z, ntor=ntor
    )
    BnR_analytic, Bnphi_analytic, BnZ_analytic = _field_divfree_rect(
        spl_analytic, grid.R, grid.Z, ntor=ntor
    )

    diagnostics = {
        "gauged_diff_R": gauged_spline_R - gauged_analytic_R,
        "gauged_diff_Z": gauged_spline_Z - gauged_analytic_Z,
        "BnR_spline": BnR_spline,
        "Bnphi_spline": Bnphi_spline,
        "BnZ_spline": BnZ_spline,
        "BnR_analytic": BnR_analytic,
        "Bnphi_analytic": Bnphi_analytic,
        "BnZ_analytic": BnZ_analytic,
    }

    return (
        ModeData(grid=grid, BnR=BnR_spline, Bnphi=Bnphi_spline, BnZ=BnZ_spline),
        spl_spline,
        diagnostics,
    )


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
    BnR_analytic: ndarray | None = None,
    Bnphi_analytic: ndarray | None = None,
    BnZ_analytic: ndarray | None = None,
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
    B_parallel_analytic = zeros_like(s_vals) if BnR_analytic is not None else None

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

        if B_parallel_analytic is not None:
            Bx_a, By_a, Bz_a = _evaluate_modes_at_point(
                BnR_analytic, Bnphi_analytic, BnZ_analytic,
                args.ntor, grid.R, grid.Z,
                x_val, y_val, z_val,
            )
            B_parallel_analytic[idx] = Bx_a * direction[0] + By_a * direction[1] + Bz_a * direction[2]

        if BnR_direct is not None and Bnphi_direct is not None and BnZ_direct is not None:
            Bx_d, By_d, Bz_d = _evaluate_modes_at_point(
                BnR_direct, Bnphi_direct, BnZ_direct,
                args.ntor, grid.R, grid.Z,
                x_val, y_val, z_val,
            )
            B_parallel_direct[idx] = Bx_d * direction[0] + By_d * direction[1] + Bz_d * direction[2]

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
        if B_parallel_analytic is not None:
            ax.plot(s_vals, B_parallel_analytic, label="Analytic gauge", linewidth=1.2, linestyle="-." )
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
    if B_parallel_analytic is not None:
        abs_err_analytic = abs(B_parallel_vector - B_parallel_analytic)
        rel_err_vs_analytic = abs_err_analytic / maximum(abs(B_parallel_analytic), 1e-20)
        print(f"  Vector vs analytic gauge (max abs): {amax(abs_err_analytic):.4e}")
        print(f"  Vector vs analytic gauge (max rel): {amax(rel_err_vs_analytic) * 100.0:.4f}%")


def _create_per_coil_plot(
    args: argparse.Namespace,
    mode_fourier: ModeData,
    mode_vector: ModeData,
    BnR_spline_all: ndarray,
    Bnphi_spline_all: ndarray,
    BnZ_spline_all: ndarray,
    output: Path,
    coil_projections: Sequence[Tuple[ndarray, ndarray]] | None,
) -> None:
    ncoil = mode_fourier.BnR.shape[0]
    grid = mode_fourier.grid
    titles = ("reference", "Anvac", "Bnvac")

    log_bn2 = empty((3, ncoil, grid.nR, grid.nZ), dtype=float)
    log_bn2[0] = log10(maximum(_magnitude_squared(mode_fourier.BnR, mode_fourier.Bnphi, mode_fourier.BnZ), 1e-300))
    log_bn2[1] = log10(maximum(_magnitude_squared(mode_vector.BnR, mode_vector.Bnphi, mode_vector.BnZ), 1e-300))
    log_bn2[2] = log10(maximum(_magnitude_squared(BnR_spline_all, Bnphi_spline_all, BnZ_spline_all), 1e-300))

    norm = Normalize(vmin=amin(log_bn2), vmax=amax(log_bn2))
    fig, axs = plt.subplots(ncoil, 3, figsize=(9, 3 * ncoil), layout="constrained")

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


def _compute_rect_spline_field(
    BnR: ndarray,
    Bnphi: ndarray,
    BnZ: ndarray,
    grid: Grid,
    order: int,
) -> Tuple[ndarray, ndarray, ndarray]:
    """Compute field using RectBivariateSpline interpolation of a given order."""
    ncoil = BnR.shape[0]
    BnR_interp = zeros_like(BnR)
    Bnphi_interp = zeros_like(Bnphi)
    BnZ_interp = zeros_like(BnZ)

    for k in range(ncoil):
        spline_R_real = RectBivariateSpline(grid.R, grid.Z, BnR[k].real, kx=order, ky=order)
        spline_R_imag = RectBivariateSpline(grid.R, grid.Z, BnR[k].imag, kx=order, ky=order)
        BnR_interp[k] = spline_R_real(grid.R, grid.Z) + 1j * spline_R_imag(grid.R, grid.Z)

        spline_phi_real = RectBivariateSpline(grid.R, grid.Z, Bnphi[k].real, kx=order, ky=order)
        spline_phi_imag = RectBivariateSpline(grid.R, grid.Z, Bnphi[k].imag, kx=order, ky=order)
        Bnphi_interp[k] = spline_phi_real(grid.R, grid.Z) + 1j * spline_phi_imag(grid.R, grid.Z)

        spline_Z_real = RectBivariateSpline(grid.R, grid.Z, BnZ[k].real, kx=order, ky=order)
        spline_Z_imag = RectBivariateSpline(grid.R, grid.Z, BnZ[k].imag, kx=order, ky=order)
        BnZ_interp[k] = spline_Z_real(grid.R, grid.Z) + 1j * spline_Z_imag(grid.R, grid.Z)

    return BnR_interp, Bnphi_interp, BnZ_interp


def _create_deriv_diff_plot(
    output: Path,
    grid: Grid,
    components: dict,
    coil_projections: Sequence[Tuple[ndarray, ndarray]] | None,
    show: bool,
) -> None:
    rows = ("B_R", "B_phi", "B_Z")
    titles = ("Analytic", "Spline", "|Δ|")
    fig, axs = plt.subplots(len(rows), len(titles), figsize=(12, 9), layout="constrained")
    extent = [grid.R_min, grid.R_max, grid.Z_min, grid.Z_max]

    all_values = []
    for comp in rows:
        data = components[comp]
        all_values.extend(
            [
                log10(maximum(abs(data["analytic"]), 1e-300)),
                log10(maximum(abs(data["spline"]), 1e-300)),
                log10(maximum(abs(data["delta"]), 1e-300)),
            ]
        )
    vmin = min(arr.min() for arr in all_values)
    vmax = max(arr.max() for arr in all_values)

    for i, comp in enumerate(rows):
        data = components[comp]
        arrays = (
            log10(maximum(abs(data["analytic"]), 1e-300)),
            log10(maximum(abs(data["spline"]), 1e-300)),
            log10(maximum(abs(data["delta"]), 1e-300)),
        )
        for j, arr in enumerate(arrays):
            ax = axs[i, j]
            im = ax.imshow(
                arr.T,
                origin="lower",
                cmap="magma",
                extent=extent,
                vmin=vmin,
                vmax=vmax,
                interpolation="bilinear",
            )
            ax.set_title(f"{comp} {titles[j]}")
            ax.set_xlabel("R [cm]")
            ax.set_ylabel("Z [cm]")
            if coil_projections:
                for R_path, Z_path in coil_projections:
                    ax.plot(R_path, Z_path, color="black", linewidth=0.8, alpha=0.6)
            ax.set_aspect("equal", adjustable="box")

    fig.colorbar(im, ax=axs, location="bottom", fraction=0.05, pad=0.08, label="log10 |B_n|")
    fig.savefig(output, dpi=150, bbox_inches="tight")
    if not show:
        plt.close(fig)


def _gauge_anvac_with_spline(
    grid: Grid,
    AnR: ndarray,
    Anphi: ndarray,
    AnZ: ndarray,
    ntor: int,
) -> Tuple[ndarray, ndarray]:
    """Apply the non-axisymmetric gauge using spline-derived derivatives."""

    ncoil = AnR.shape[0]
    R_vals = grid.R
    Z_vals = grid.Z
    gauged_AnR = zeros_like(AnR)
    gauged_AnZ = zeros_like(AnZ)

    factor = 1j / ntor
    R_matrix = R_vals[:, newaxis]

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
    R_vals: ndarray,
    Z_vals: ndarray,
    ntor: int,
) -> Tuple[ndarray, ndarray, ndarray]:
    """Reproduce ``field_divfree`` using RectBivariateSpline evaluation."""

    if ntor == 0:
        raise ValueError("_field_divfree_rect requires ntor != 0; use _compute_axisymmetric_field")

    ncoil = len(spl['AnR_Re'])
    nR = R_vals.size
    nZ = Z_vals.size
    R_matrix = R_vals[:, newaxis]

    BnR = zeros((ncoil, nR, nZ), dtype=complex)
    Bnphi = zeros_like(BnR)
    BnZ = zeros_like(BnR)

    for k in range(ncoil):
        anR_real = spl['AnR_Re'][k](R_vals, Z_vals)
        anR_imag = spl['AnR_Im'][k](R_vals, Z_vals)
        anZ_real = spl['AnZ_Re'][k](R_vals, Z_vals)
        anZ_imag = spl['AnZ_Im'][k](R_vals, Z_vals)

        anR = anR_real + 1j * anR_imag
        anZ = anZ_real + 1j * anZ_imag

        dAnR_dZ = spl['AnR_Re'][k](R_vals, Z_vals, dy=1) + 1j * spl['AnR_Im'][k](R_vals, Z_vals, dy=1)
        dAnZ_dR = spl['AnZ_Re'][k](R_vals, Z_vals, dx=1) + 1j * spl['AnZ_Im'][k](R_vals, Z_vals, dx=1)

        BnR[k] = 1j * ntor * anZ / R_matrix
        Bnphi[k] = dAnR_dZ - dAnZ_dR
        BnZ[k] = -1j * ntor * anR / R_matrix

    return BnR, Bnphi, BnZ


def _compute_axisymmetric_field(
    grid: Grid,
    AnR: ndarray,
    Anphi: ndarray,
    AnZ: ndarray,
    _dAnphi_dR: ndarray,
    _dAnphi_dZ: ndarray,
) -> Tuple[ndarray, ndarray, ndarray]:
    """Compute the axisymmetric (ntor=0) magnetic field directly from vector potential data."""

    ncoil = AnR.shape[0]
    BnR = zeros_like(AnR)
    Bnphi = zeros_like(AnR)
    BnZ = zeros_like(AnR)

    R_vals = grid.R
    Z_vals = grid.Z

    for k in range(ncoil):
        A_R = AnR[k]
        A_phi = Anphi[k]
        A_Z = AnZ[k]

        spline_Aphi_real = RectBivariateSpline(R_vals, Z_vals, A_phi.real, kx=5, ky=5)
        spline_Aphi_imag = RectBivariateSpline(R_vals, Z_vals, A_phi.imag, kx=5, ky=5)
        dAphi_dZ = spline_Aphi_real(R_vals, Z_vals, dy=1) + 1j * spline_Aphi_imag(R_vals, Z_vals, dy=1)
        dAphi_dR = spline_Aphi_real(R_vals, Z_vals, dx=1) + 1j * spline_Aphi_imag(R_vals, Z_vals, dx=1)

        BnR[k] = -dAphi_dZ
        with errstate(divide="ignore", invalid="ignore"):
            BnZ[k] = dAphi_dR + A_phi / R_vals[:, newaxis]

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

        invalid = ~isfinite(BnZ[k])
        if invalid.any():
            finite_vals = BnZ[k][isfinite(BnZ[k])]
            if finite_vals.size > 0:
                BnZ[k][invalid] = finite_vals[0]
            else:
                BnZ[k][invalid] = 0.0

    return BnR, Bnphi, BnZ


def _create_sum_plot(
    args: argparse.Namespace,
    magnitude_fourier: ndarray,
    magnitude_vector: ndarray,
    magnitude_spline: ndarray,
    magnitude_field_divb0: ndarray | None,
    magnitude_rect5: ndarray | None,
    magnitude_rect3: ndarray | None,
    magnitude_rect1: ndarray | None,
    magnitude_direct: ndarray | None,
    output: Path,
    grid: Grid,
    coil_projections: Sequence[Tuple[ndarray, ndarray]] | None,
    components: dict[str, dict[str, ndarray]] | None = None,
) -> None:
    top_titles = ("reference sum", "Anvac sum", "Bnvac sum", "Direct sum")
    top_panels = [magnitude_fourier, magnitude_vector, magnitude_spline]
    if magnitude_direct is None:
        top_panels.append(magnitude_fourier)
    else:
        top_panels.append(magnitude_direct)

    bottom_titles = ("field_divB0 sum", "Rect k=5 sum", "Rect k=3 sum", "Rect k=1 sum")
    bottom_panels = [magnitude_field_divb0, magnitude_rect5, magnitude_rect3, magnitude_rect1]

    accum = [log10(maximum(panel, 1e-300)) for panel in top_panels]
    if magnitude_field_divb0 is not None:
        accum.append(log10(maximum(magnitude_field_divb0, 1e-300)))
    if magnitude_rect5 is not None:
        accum.append(log10(maximum(magnitude_rect5, 1e-300)))
    if magnitude_rect3 is not None:
        accum.append(log10(maximum(magnitude_rect3, 1e-300)))
    if magnitude_rect1 is not None:
        accum.append(log10(maximum(magnitude_rect1, 1e-300)))
    if components is not None:
        for dataset in components.values():
            for value in dataset.values():
                if value is not None:
                    accum.append(log10(maximum(value, 1e-300)))

    norm = Normalize(vmin=min(a.min() for a in accum), vmax=max(a.max() for a in accum))

    n_component_rows = 0 if components is None else len(components)
    nrows = 2 + n_component_rows
    fig, axs = plt.subplots(nrows, 4, figsize=(14, 3.5 * nrows), layout="constrained")
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

    # Plot bottom panels
    for bottom_idx, (panel, title) in enumerate(zip(bottom_panels, bottom_titles)):
        if panel is not None and title:
            ax = axs[1, bottom_idx]
            im_bottom = ax.imshow(
                log10(maximum(panel, 1e-300)).T,
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
            plotted_axes.append(ax)
            im = im_bottom
        else:
            axs[1, bottom_idx].axis("off")

    # Turn off remaining bottom panels (none, ensured by padding)
    for idx in range(len(bottom_panels), 4):
        axs[1, idx].axis("off")

    if components is not None:
        row = 2
        for comp_label, datasets in components.items():
            titles = [f"{comp_label} {name}" for name in ["reference", "Anvac", "Bnvac", "Direct"]]
            panels = [datasets.get(key) for key in ["reference", "anvac", "bnvac", "direct"]]
            # ensure rows exist even if fewer than 4 panels
            for col in range(4):
                ax = axs[row, col]
                if col >= len(panels) or panels[col] is None:
                    ax.axis("off")
                    continue
                panel = panels[col]
                im = ax.imshow(
                    log10(maximum(panel, 1e-300)).T,
                    origin="lower",
                    cmap="magma",
                    extent=extent,
                    norm=norm,
                    interpolation="bilinear",
                )
                ax.set_title(titles[col])
                ax.set_xlabel("R [cm]")
                ax.set_ylabel("Z [cm]")
                if coil_projections:
                    for R_path, Z_path in coil_projections:
                        ax.plot(R_path, Z_path, color="black", linewidth=0.8, alpha=0.6)
                ax.set_aspect("equal", adjustable="box")
                plotted_axes.append(ax)
            row += 1

    cbar = fig.colorbar(im, ax=plotted_axes, location="bottom", fraction=0.05, pad=0.08)
    cbar.set_label("log10 |Σ B_n|^2")
    fig.savefig(output, dpi=args.dpi, bbox_inches="tight")
    if not args.show:
        plt.close(fig)




def _print_stats(
    magnitude_fourier: ndarray,
    magnitude_vector: ndarray,
    magnitude_spline: ndarray,
    magnitude_field_divb0: ndarray | None,
    magnitude_rect5: ndarray | None,
    magnitude_rect3: ndarray | None,
    magnitude_rect1: ndarray | None,
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

    if magnitude_rect5 is not None:
        rel_err_rect5 = abs(magnitude_rect5 - magnitude_fourier) / (magnitude_fourier + 1e-15)
        median_rect5 = float(median(rel_err_rect5) * 100.0)
        mean_rect5 = float(mean(rel_err_rect5) * 100.0)
        print(f"Median relative error (Rect k=5): {median_rect5:.4f}%")
        print(f"Mean relative error (Rect k=5):   {mean_rect5:.4f}%")

    if magnitude_rect3 is not None:
        rel_err_rect3 = abs(magnitude_rect3 - magnitude_fourier) / (magnitude_fourier + 1e-15)
        median_rect3 = float(median(rel_err_rect3) * 100.0)
        mean_rect3 = float(mean(rel_err_rect3) * 100.0)
        print(f"Median relative error (Rect k=3): {median_rect3:.4f}%")
        print(f"Mean relative error (Rect k=3):   {mean_rect3:.4f}%")

    if magnitude_rect1 is not None:
        rel_err_rect1 = abs(magnitude_rect1 - magnitude_fourier) / (magnitude_fourier + 1e-15)
        median_rect1 = float(median(rel_err_rect1) * 100.0)
        mean_rect1 = float(mean(rel_err_rect1) * 100.0)
        print(f"Median relative error (Rect k=1): {median_rect1:.4f}%")
        print(f"Mean relative error (Rect k=1):   {mean_rect1:.4f}%")

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
    mode_vector, spline_vec, diagnostics = _load_mode_from_anvac(args.test, args.ntor)
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

    BnR_spline_all = mode_vector.BnR
    Bnphi_spline_all = mode_vector.Bnphi
    BnZ_spline_all = mode_vector.BnZ

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

    BnR_analytic_sum = Bnphi_analytic_sum = BnZ_analytic_sum = None
    if diagnostics is not None:
        BnR_analytic_sum = _tensordot_currents(weights, diagnostics["BnR_analytic"])
        Bnphi_analytic_sum = _tensordot_currents(weights, diagnostics["Bnphi_analytic"])
        BnZ_analytic_sum = _tensordot_currents(weights, diagnostics["BnZ_analytic"])

        max_gauge_R = amax(abs(diagnostics["gauged_diff_R"]))
        max_gauge_Z = amax(abs(diagnostics["gauged_diff_Z"]))
        print(
            f"Gauge diagnostics: max |ΔA_R| = {max_gauge_R:.3e}, max |ΔA_Z| = {max_gauge_Z:.3e}"
        )

        if args.deriv_diff_output is not None:
            diff_components = {
                "B_R": {
                    "analytic": BnR_analytic_sum,
                    "spline": BnR_spline_sum,
                    "delta": BnR_spline_sum - BnR_analytic_sum,
                },
                "B_phi": {
                    "analytic": Bnphi_analytic_sum,
                    "spline": Bnphi_spline_sum,
                    "delta": Bnphi_spline_sum - Bnphi_analytic_sum,
                },
                "B_Z": {
                    "analytic": BnZ_analytic_sum,
                    "spline": BnZ_spline_sum,
                    "delta": BnZ_spline_sum - BnZ_analytic_sum,
                },
            }
            max_delta = max(
                amax(abs(arr["delta"])) for arr in diff_components.values()
            )
            print(f"Gauge diagnostics: max |ΔB| (coil-weighted) = {max_delta:.3e}")

            _create_deriv_diff_plot(
                args.deriv_diff_output,
                mode_fourier.grid,
                diff_components,
                coil_projections,
                args.show,
            )

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

    # Compute RectBivariateSpline field variants for comparison
    BnR_rect5_all, Bnphi_rect5_all, BnZ_rect5_all = _compute_rect_spline_field(
        BnR_spline_all, Bnphi_spline_all, BnZ_spline_all, mode_fourier.grid, order=5
    )
    BnR_rect5_sum = _tensordot_currents(weights, BnR_rect5_all)
    Bnphi_rect5_sum = _tensordot_currents(weights, Bnphi_rect5_all)
    BnZ_rect5_sum = _tensordot_currents(weights, BnZ_rect5_all)
    magnitude_rect5 = _magnitude_squared(BnR_rect5_sum, Bnphi_rect5_sum, BnZ_rect5_sum)

    BnR_rect3_all, Bnphi_rect3_all, BnZ_rect3_all = _compute_rect_spline_field(
        BnR_spline_all, Bnphi_spline_all, BnZ_spline_all, mode_fourier.grid, order=3
    )
    BnR_rect3_sum = _tensordot_currents(weights, BnR_rect3_all)
    Bnphi_rect3_sum = _tensordot_currents(weights, Bnphi_rect3_all)
    BnZ_rect3_sum = _tensordot_currents(weights, BnZ_rect3_all)
    magnitude_rect3 = _magnitude_squared(BnR_rect3_sum, Bnphi_rect3_sum, BnZ_rect3_sum)

    BnR_rect1_all, Bnphi_rect1_all, BnZ_rect1_all = _compute_rect_spline_field(
        BnR_spline_all, Bnphi_spline_all, BnZ_spline_all, mode_fourier.grid, order=1
    )
    BnR_rect1_sum = _tensordot_currents(weights, BnR_rect1_all)
    Bnphi_rect1_sum = _tensordot_currents(weights, Bnphi_rect1_all)
    BnZ_rect1_sum = _tensordot_currents(weights, BnZ_rect1_all)
    magnitude_rect1 = _magnitude_squared(BnR_rect1_sum, Bnphi_rect1_sum, BnZ_rect1_sum)

    if args.per_coil_output is not None:
        _create_per_coil_plot(args, mode_fourier, mode_vector,
                              BnR_spline_all, Bnphi_spline_all, BnZ_spline_all,
                              args.per_coil_output, coil_projections)

    components = None
    if args.sum_output is not None:
        def _build_component_dict(ref: ndarray, vec: ndarray, spline_arr: ndarray, direct_arr: ndarray | None) -> dict[str, ndarray | None]:
            mapping: dict[str, ndarray | None] = {
                'reference': abs(ref)**2,
                'anvac': abs(vec)**2,
                'bnvac': abs(spline_arr)**2,
                'direct': abs(direct_arr)**2 if direct_arr is not None else None,
            }
            return mapping

        components = {
            'B_R': _build_component_dict(BnR_fourier_sum, BnR_vector_sum, BnR_spline_sum, BnR_direct),
            'B_phi': _build_component_dict(Bnphi_fourier_sum, Bnphi_vector_sum, Bnphi_spline_sum, Bnphi_direct),
            'B_Z': _build_component_dict(BnZ_fourier_sum, BnZ_vector_sum, BnZ_spline_sum, BnZ_direct),
        }

        _create_sum_plot(
            args,
            magnitude_fourier,
            magnitude_vector,
            magnitude_spline,
            magnitude_field_divb0,
            magnitude_rect5,
            magnitude_rect3,
            magnitude_rect1,
            direct_magnitude,
            args.sum_output,
            mode_fourier.grid,
            coil_projections,
            components,
        )

    _print_stats(
        magnitude_fourier,
        magnitude_vector,
        magnitude_spline,
        magnitude_field_divb0,
        magnitude_rect5,
        magnitude_rect3,
        magnitude_rect1,
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
            BnR_analytic_sum,
            Bnphi_analytic_sum,
            BnZ_analytic_sum,
        )


if __name__ == "__main__":
    main()
