#!/usr/bin/env python
"""Cross-check libneo Cerfon-Freidberg solver against ASCOT5's B_GS implementation."""
from __future__ import annotations

import argparse
import ctypes
import math
import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ASCOT_COEFF = np.array([
    2.21808016e-02,
   -1.28841781e-01,
   -4.17718173e-02,
   -6.22680280e-02,
    6.20083978e-03,
   -1.20524711e-03,
   -3.70147050e-05,
    0.00000000e+00,
    0.00000000e+00,
    0.00000000e+00,
    0.00000000e+00,
    0.00000000e+00,
   -1.55000000e-01,
    0.0,
], dtype=np.float64)

R0_DEFAULT = 6.2
Z0_DEFAULT = 0.0
BPHI0_DEFAULT = 5.3
PSI_MULT_DEFAULT = 200.0


class AscotGS(ctypes.Structure):
    _fields_ = [
        ("R0", ctypes.c_double),
        ("z0", ctypes.c_double),
        ("raxis", ctypes.c_double),
        ("zaxis", ctypes.c_double),
        ("B_phi0", ctypes.c_double),
        ("psi0", ctypes.c_double),
        ("psi1", ctypes.c_double),
        ("psi_mult", ctypes.c_double),
        ("psi_coeff", ctypes.c_double * 14),
        ("Nripple", ctypes.c_int),
        ("a0", ctypes.c_double),
        ("alpha0", ctypes.c_double),
        ("delta0", ctypes.c_double),
    ]


def run(cmd: list[str], cwd: Path | None = None, **kwargs) -> None:
    subprocess.run(cmd, cwd=cwd, check=True, **kwargs)


def compile_ascot_library(clone_dir: Path) -> Path:
    compiler = shutil.which("gcc-14") or shutil.which("gcc")
    if compiler is None:
        raise RuntimeError("gcc or gcc-14 not found; required to build ascot5 reference library")

    flags: list[str] = []
    link_flags: list[str] = ["-fopenmp", "-lhdf5_hl", "-lhdf5"]

    homebrew_root = Path("/opt/homebrew")
    if homebrew_root.exists():
        hb_include = homebrew_root / "include"
        if hb_include.exists():
            flags.append(f"-I{hb_include}")
        omp_include = homebrew_root / "opt" / "libomp" / "include"
        if omp_include.exists():
            flags.append(f"-I{omp_include}")
        hb_lib = homebrew_root / "lib"
        if hb_lib.exists():
            link_flags.append(f"-L{hb_lib}")
        omp_lib = homebrew_root / "opt" / "libomp" / "lib"
        if omp_lib.exists():
            link_flags.append(f"-L{omp_lib}")
            link_flags.append("-lomp")

    run(["make", "clean"], cwd=clone_dir)
    run([
        "make",
        f"CC={compiler}",
        f"FLAGS={' '.join(flags)}",
        f"LFLAGS={' '.join(link_flags)}",
        "libascot",
        "-j2",
    ], cwd=clone_dir)

    candidates = [clone_dir / "build" / "libascot.so", clone_dir / "build" / "libascot.dylib"]
    for lib in candidates:
        if lib.exists():
            return lib
    raise FileNotFoundError("ASCOT5 build did not produce libascot shared library")


def load_ascot_functions(lib_path: Path) -> tuple[ctypes.CDLL, AscotGS]:
    lib = ctypes.CDLL(str(lib_path))
    lib.B_GS_init.argtypes = [ctypes.POINTER(AscotGS), ctypes.c_double, ctypes.c_double,
                              ctypes.c_double, ctypes.c_double, ctypes.c_double,
                              ctypes.c_double, ctypes.c_double, ctypes.c_double,
                              ctypes.POINTER(ctypes.c_double), ctypes.c_int,
                              ctypes.c_double, ctypes.c_double, ctypes.c_double]
    lib.B_GS_init.restype = ctypes.c_int
    lib.B_GS_eval_psi.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_double,
                                  ctypes.c_double, ctypes.c_double, ctypes.POINTER(AscotGS)]
    lib.B_GS_eval_psi.restype = ctypes.c_int
    lib.B_GS_get_axis_rz.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(AscotGS)]
    lib.B_GS_get_axis_rz.restype = ctypes.c_int
    lib.B_GS_free.argtypes = [ctypes.POINTER(AscotGS)]
    lib.B_GS_free.restype = None
    return lib, AscotGS()


def evaluate_ascot_flux(lib, data: AscotGS, points: np.ndarray) -> np.ndarray:
    flux = np.empty(points.shape[0], dtype=np.float64)
    psi_val = ctypes.c_double()
    for idx, (r, z) in enumerate(points):
        err = lib.B_GS_eval_psi(ctypes.byref(psi_val), float(r), 0.0, float(z), ctypes.byref(data))
        if err != 0:
            raise RuntimeError(f"B_GS_eval_psi failed at point index {idx}")
        flux[idx] = psi_val.value
    return flux


def read_flux_csv(csv_path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    data = np.loadtxt(csv_path, delimiter=",", skiprows=1)
    R = data[:, 0]
    Z = data[:, 1]
    psi = data[:, 2]
    return R, Z, psi


def reshape_grid(R: np.ndarray, Z: np.ndarray, values: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    r_unique = np.unique(R)
    z_unique = np.unique(Z)
    if r_unique.size * z_unique.size != values.size:
        raise ValueError("Flux data cannot be reshaped to a rectangular grid.")
    order = np.argsort(R * 1e6 + Z)
    values_sorted = values[order]
    R_sorted = R[order]
    Z_sorted = Z[order]
    V_grid = values_sorted.reshape((r_unique.size, z_unique.size))
    R_grid = R_sorted.reshape((r_unique.size, z_unique.size))
    Z_grid = Z_sorted.reshape((r_unique.size, z_unique.size))
    return R_grid, Z_grid, V_grid


def compute_axis_shift(R_grid: np.ndarray, Z_grid: np.ndarray, psi_grid: np.ndarray, R0: float) -> float:
    z_axis_idx = np.argmin(np.abs(Z_grid[0, :]))
    psi_slice = psi_grid[:, z_axis_idx]
    axis_idx = np.argmin(psi_slice)
    return float(R_grid[axis_idx, z_axis_idx] - R0)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path, default=None, help="Path to libneo build directory")
    parser.add_argument("--output-dir", type=Path, default=None, help="Directory for artifacts")
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[2]
    build_dir = args.build_dir or repo_root / "build"
    build_dir = build_dir.resolve()
    if not build_dir.exists():
        raise FileNotFoundError(f"Build directory {build_dir} does not exist")

    fortran_exe = build_dir / "test" / "test_analytical_circular.x"
    if not fortran_exe.exists():
        raise FileNotFoundError(f"Required executable not found: {fortran_exe}")

    artifact_root = (args.output_dir or build_dir / "test" / "ascot5_compare").resolve()
    artifact_root.mkdir(parents=True, exist_ok=True)

    # 1. Run libneo analytical test to dump circular flux data
    libneo_dir = artifact_root / "libneo_flux"
    libneo_dir.mkdir(exist_ok=True)
    run([str(fortran_exe), str(libneo_dir)])
    flux_csv = libneo_dir / "flux_circular.csv"
    if not flux_csv.exists():
        raise FileNotFoundError(f"Expected flux output missing: {flux_csv}")

    R, Z, psi_lib = read_flux_csv(flux_csv)
    R_grid, Z_grid, psi_lib_grid = reshape_grid(R, Z, psi_lib)

    # 2. Clone ascot5 and build the B_GS shared library
    ascot_clone = artifact_root / "ascot5"
    if ascot_clone.exists():
        shutil.rmtree(ascot_clone)
    run(["git", "clone", "--depth", "1", "https://github.com/ascot4fusion/ascot5.git", str(ascot_clone)])

    lib_path = compile_ascot_library(ascot_clone)
    lib, bgs_data = load_ascot_functions(lib_path)

    coeff_array = (ctypes.c_double * len(ASCOT_COEFF))(*ASCOT_COEFF)
    init_code = lib.B_GS_init(
        ctypes.byref(bgs_data),
        R0_DEFAULT,
        Z0_DEFAULT,
        R0_DEFAULT,
        Z0_DEFAULT,
        BPHI0_DEFAULT,
        0.0,
        1.0,
        PSI_MULT_DEFAULT,
        coeff_array,
        0,
        1.0,
        0.0,
        0.0,
    )
    if init_code != 0:
        raise RuntimeError("B_GS_init failed")

    pts = np.column_stack((R, Z))
    psi_asc = evaluate_ascot_flux(lib, bgs_data, pts)
    psi_asc_grid = psi_asc.reshape(R_grid.shape)
    lib.B_GS_free(ctypes.byref(bgs_data))

    # 3. Align fluxes by subtracting their minima
    psi_lib_aligned = psi_lib_grid - psi_lib_grid.min()
    psi_asc_aligned = psi_asc_grid - psi_asc_grid.min()

    diff = psi_asc_aligned - psi_lib_aligned
    rms = float(np.sqrt(np.mean(diff**2)))

    shift_lib = compute_axis_shift(R_grid, Z_grid, psi_lib_aligned, R0_DEFAULT)
    shift_asc = compute_axis_shift(R_grid, Z_grid, psi_asc_aligned, R0_DEFAULT)
    shift_diff = abs(shift_lib - shift_asc)

    print(f"Flux RMS difference: {rms:.6e} Wb")
    print(f"libneo axis shift:   {shift_lib:.6e} m")
    print(f"ascot5 axis shift:   {shift_asc:.6e} m")
    print(f"Shift difference:    {shift_diff:.6e} m")

    rms_tol = 3.0e-2
    shift_tol = 2.0e-3
    if not math.isfinite(rms) or rms > rms_tol:
        raise AssertionError(f"Flux RMS difference {rms:.3e} exceeds tolerance {rms_tol:.3e}")
    if not math.isfinite(shift_diff) or shift_diff > shift_tol:
        raise AssertionError(f"Axis shift mismatch {shift_diff:.3e} exceeds tolerance {shift_tol:.3e}")

    # 4. Generate diagnostic plot
    levels = np.linspace(min(psi_lib_aligned.min(), psi_asc_aligned.min()),
                         max(psi_lib_aligned.max(), psi_asc_aligned.max()), 15)
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.contour(R_grid, Z_grid, psi_lib_aligned, levels=levels, colors='C0', linestyles='solid', linewidths=1.1)
    ax.contour(R_grid, Z_grid, psi_asc_aligned, levels=levels, colors='C1', linestyles='dashed', linewidths=1.1)
    ax.set_xlabel('R [m]')
    ax.set_ylabel('Z [m]')
    ax.set_title('Cerfon-Freidberg vs ASCOT5 B_GS flux surfaces')
    handles = [plt.Line2D([], [], color='C0', label='libneo'), plt.Line2D([], [], color='C1', linestyle='--', label='ascot5')]
    ax.legend(handles=handles, loc='best')
    ax.set_aspect('equal')
    plot_path = artifact_root / 'ascot5_flux_comparison.png'
    fig.tight_layout()
    fig.savefig(plot_path, dpi=150)
    plt.close(fig)

    print(f"Saved flux comparison plot to {plot_path}")


if __name__ == "__main__":
    main()
