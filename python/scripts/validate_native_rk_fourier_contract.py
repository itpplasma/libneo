#!/usr/bin/env python3
"""Validate chartmap native-RK fields against an independent Fourier source.

The script deliberately does not import the boozmn-to-chartmap converter.  It
reconstructs the radial interpolation, centered-difference derivative channel,
and Boozer Fourier sums directly, then compares them with the ``rk_*`` arrays
stored in a chartmap NetCDF file.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import make_interp_spline


def _fourier_sum(coefficients, ixm, ixn, theta, zeta, sine=False):
    angle = (
        ixm[:, None, None] * theta[None, :, None]
        - ixn[:, None, None] * zeta[None, None, :]
    )
    basis = np.sin(angle) if sine else np.cos(angle)
    return np.einsum("sm,mnt->snt", coefficients, basis, optimize=True)


def _extend_half(values):
    extended = np.empty((values.shape[0] + 2, *values.shape[1:]))
    extended[1:-1] = values
    extended[0] = 1.5 * extended[1] - 0.5 * extended[2]
    extended[-1] = 1.5 * extended[-2] - 0.5 * extended[-3]
    return extended


def _radial_value_and_derivative(values, s_half, s_full, s_target):
    extended = _extend_half(values)
    half_with_edges = np.concatenate(([0.0], s_half, [1.0]))
    value = make_interp_spline(half_with_edges, extended, k=3, axis=0)(s_target)
    ds = s_full[1] - s_full[0]
    centered = (extended[2:-1] - extended[1:-2]) / ds
    derivative = make_interp_spline(
        s_full[1:-1], centered, k=3, axis=0
    )(s_target)
    return value, derivative


def _reconstruct(boozmn: Path, chartmap: Path):
    with Dataset(boozmn) as source, Dataset(chartmap) as target:
        ns = int(np.asarray(source.variables["ns_b"][:]))
        nfp = int(np.asarray(source.variables["nfp_b"][:]))
        jlist = np.asarray(source.variables["jlist"][:], dtype=int)
        ixm = np.asarray(source.variables["ixm_b"][:], dtype=int)
        ixn = np.asarray(source.variables["ixn_b"][:], dtype=int)
        bmnc = np.asarray(source.variables["bmnc_b"][:], dtype=float)
        lasym = bool(np.asarray(source.variables["lasym__logical__"][:]))
        bmns = (
            np.asarray(source.variables["bmns_b"][:], dtype=float)
            if lasym
            else None
        )

        s_half = (jlist - 1.5) / (ns - 1)
        s_full = np.linspace(0.0, 1.0, ns)
        s = np.asarray(target.variables["s"][:], dtype=float)
        theta = np.asarray(target.variables["theta"][:], dtype=float)
        zeta = np.asarray(target.variables["zeta"][:], dtype=float)

        bmnc_s, dbmnc_s = _radial_value_and_derivative(
            bmnc, s_half, s_full, s
        )
        bmod = _fourier_sum(bmnc_s, ixm, ixn, theta, zeta)
        dbmod_ds = _fourier_sum(dbmnc_s, ixm, ixn, theta, zeta)
        dbmod_dtheta = _fourier_sum(
            -bmnc_s * ixm[None, :], ixm, ixn, theta, zeta, sine=True
        )
        dbmod_dzeta = _fourier_sum(
            bmnc_s * ixn[None, :], ixm, ixn, theta, zeta, sine=True
        )
        if bmns is not None:
            bmns_s, dbmns_s = _radial_value_and_derivative(
                bmns, s_half, s_full, s
            )
            bmod += _fourier_sum(bmns_s, ixm, ixn, theta, zeta, sine=True)
            dbmod_ds += _fourier_sum(
                dbmns_s, ixm, ixn, theta, zeta, sine=True
            )
            dbmod_dtheta += _fourier_sum(
                bmns_s * ixm[None, :], ixm, ixn, theta, zeta
            )
            dbmod_dzeta += _fourier_sum(
                -bmns_s * ixn[None, :], ixm, ixn, theta, zeta
            )

        reconstructed = {
            "Bmod": bmod * 1.0e4,
            "dBmod_ds": dbmod_ds * 1.0e4,
            "dBmod_dtheta": dbmod_dtheta * 1.0e4,
            "dBmod_dzeta": dbmod_dzeta * 1.0e4,
        }
        actual = {
            name: np.asarray(target.variables[f"rk_{name}"][:], dtype=float)
            .transpose(2, 1, 0)
            for name in reconstructed
        }
    return reconstructed, actual


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("boozmn", type=Path)
    parser.add_argument("chartmap", type=Path)
    parser.add_argument("--rtol", type=float, default=1.0e-13)
    parser.add_argument("--atol", type=float, default=1.0e-8)
    args = parser.parse_args()

    reconstructed, actual = _reconstruct(args.boozmn, args.chartmap)
    for name in reconstructed:
        difference = np.abs(reconstructed[name] - actual[name])
        scale = max(float(np.max(np.abs(actual[name]))), np.finfo(float).tiny)
        print(
            f"{name}: max_abs={np.max(difference):.6e} "
            f"max_rel={np.max(difference) / scale:.6e}"
        )
        np.testing.assert_allclose(
            reconstructed[name], actual[name], rtol=args.rtol, atol=args.atol
        )
    print("native RK Fourier contract: PASS")


if __name__ == "__main__":
    main()
