"""Compose NEO-2 Boozer spectra and write a SIMPLE-compatible chartmap.

The NEO-2 ASDEX-U reader evaluates ``m*theta + n*nper*phi``, while chartmaps
use ``m*theta - n*nfp*zeta``. This converter changes the perturbation mode
sign at that boundary. Axisymmetric modes remain unchanged, and the separate
perturbation file defines the output field periodicity.
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
from types import SimpleNamespace
import tempfile

import numpy as np


_COEFFICIENTS = (
    "rmnc", "rmns", "zmnc", "zmns", "vmnc", "vmns", "bmnc", "bmns"
)
_GEOMETRY_COEFFICIENTS = _COEFFICIENTS[:6]
_PROFILE_NAMES = (
    "s", "iota", "Jpol_divided_by_nper", "Itor", "pprime", "sqrt_g_00"
)
_SOURCE_PHASES = {
    "neo2-plus": "m*theta+n*nper*phi",
    "minus": "m*theta-n*nper*phi",
}
_OUTPUT_PHASE = "m*theta-n*nfp*zeta"


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _require_close(name, left, right, *, rtol=1.0e-8, atol=1.0e-12):
    if not np.allclose(left, right, rtol=rtol, atol=atol):
        raise ValueError(f"axisymmetric and perturbation {name} are incompatible")


def _resample_surfaces(data, surfaces):
    """Copy spectra/profiles to covered surfaces; never extrapolate a field."""
    from scipy.interpolate import CubicSpline

    source = np.asarray(data.s, dtype=float)
    if (source.size < 2 or not np.all(np.isfinite(source))
            or np.any(np.diff(source) <= 0)):
        raise ValueError("radial grids must be finite and strictly increasing")
    if surfaces[0] < source[0] or surfaces[-1] > source[-1]:
        raise ValueError("radial grids have insufficient overlap")
    for surface in range(data.nsurf):
        if (list(zip(data.m[surface], data.n[surface]))
                != list(zip(data.m[0], data.n[0]))):
            raise ValueError("mode ordering changes between surfaces")
    result = SimpleNamespace(**vars(data))
    result.nsurf = len(surfaces)
    result.s = surfaces.copy()
    for name in _PROFILE_NAMES[1:] + _COEFFICIENTS:
        values = np.asarray(getattr(data, name), dtype=float)
        setattr(result, name, CubicSpline(source, values, axis=0,
                                         extrapolate=False)(surfaces))
    result.m = [np.asarray(data.m[0]).copy() for _ in surfaces]
    result.n = [np.asarray(data.n[0]).copy() for _ in surfaces]
    return result


def _align_pair(axis, perturbation):
    # Retain the equilibrium's uniform half-grid on the common support. GPEC
    # spectra may have a nonuniform radial grid and zero (unspecified) metadata.
    axis_s = np.asarray(axis.s, dtype=float)
    pert_s = np.asarray(perturbation.s, dtype=float)
    for grid in (axis_s, pert_s):
        if (grid.size < 2 or not np.all(np.isfinite(grid))
                or np.any(np.diff(grid) <= 0)):
            raise ValueError("radial grids must be finite and strictly increasing")
    if np.array_equal(axis_s, pert_s):
        # Still validate mode ordering and copy both inputs.
        surfaces = axis_s
    else:
        surfaces = axis_s[(axis_s >= pert_s[0]) & (axis_s <= pert_s[-1])]
    if len(surfaces) < 2:
        raise ValueError("radial grids have insufficient overlap")
    return (_resample_surfaces(axis, surfaces),
            _resample_surfaces(perturbation, surfaces))


def _validate_pair(axis, perturbation):
    for name, left, right in (
        ("iota profiles", axis.iota, perturbation.iota),
        ("toroidal-current profiles", axis.Itor, perturbation.Itor),
        ("poloidal-current profiles",
         np.asarray(axis.Jpol_divided_by_nper) * axis.nper,
         np.asarray(perturbation.Jpol_divided_by_nper) * perturbation.nper),
        ("toroidal fluxes", axis.flux, perturbation.flux),
        ("minor radii", axis.a, perturbation.a),
        ("major radii", axis.R, perturbation.R),
    ):
        # A perturbation-only spectrum has no equilibrium profiles/header.
        if np.any(np.asarray(right) != 0):
            _require_close(name, left, right)
    for surface in range(axis.nsurf):
        if np.any(np.asarray(axis.n[surface]) != 0):
            raise ValueError("the axisymmetric .bc file contains n != 0 modes")


def combine_boozer_data(axis, perturbation=None, *, scale=1.0,
                        geometry="axisymmetric", source_phase="neo2-plus"):
    """Return composed Boozer data without modifying either input.

    ``geometry='axisymmetric'`` adds only the perturbation's magnetic-field
    coefficients. ``geometry='full'`` also adds its R, Z, and toroidal phase
    coefficients. The coefficients are appended because the files may use
    disjoint mode sets and different field-period declarations. Differing radial
    grids are matched by cubic interpolation onto equilibrium surfaces within
    their common support. Equilibrium profiles always come from ``axis``;
    zero perturbation profiles/header values mean unspecified metadata.
    """
    if geometry not in ("axisymmetric", "full"):
        raise ValueError("geometry must be 'axisymmetric' or 'full'")
    if source_phase not in _SOURCE_PHASES:
        raise ValueError("source_phase must be 'neo2-plus' or 'minus'")
    if not np.isfinite(scale):
        raise ValueError("perturbation scale must be finite")

    if perturbation is None:
        output_nper = axis.nper
        profile_source = axis
    else:
        axis, perturbation = _align_pair(axis, perturbation)
        _validate_pair(axis, perturbation)
        output_nper = perturbation.nper
        profile_source = axis

    combined = SimpleNamespace(
        comments=list(axis.comments),
        nsurf=axis.nsurf,
        nper=output_nper,
        axis_nper=axis.nper,
        flux=float(axis.flux),
        a=float(axis.a),
        R=float(axis.R),
    )
    for name in _PROFILE_NAMES:
        setattr(combined, name, np.asarray(getattr(profile_source, name)).copy())

    combined.m = []
    combined.n = []
    for name in _COEFFICIENTS:
        setattr(combined, name, [])

    for surface in range(axis.nsurf):
        axis_m = np.asarray(axis.m[surface], dtype=int)
        axis_n = np.asarray(axis.n[surface], dtype=int)
        if perturbation is None:
            perturbation_m = np.empty(0, dtype=int)
            perturbation_n = np.empty(0, dtype=int)
        else:
            perturbation_m = np.asarray(perturbation.m[surface], dtype=int)
            perturbation_n = np.asarray(perturbation.n[surface], dtype=int)
            if source_phase == "neo2-plus":
                perturbation_n = -perturbation_n
        combined.m.append(np.concatenate((axis_m, perturbation_m)))
        combined.n.append(np.concatenate((axis_n, perturbation_n)))

        for name in _COEFFICIENTS:
            axis_values = np.asarray(getattr(axis, name)[surface], dtype=float)
            if name in ("vmnc", "vmns"):
                axis_values = axis_values * output_nper / axis.nper
            if perturbation is None:
                perturbation_values = np.empty(0)
            elif name in _GEOMETRY_COEFFICIENTS and geometry == "axisymmetric":
                perturbation_values = np.zeros(perturbation_m.size)
            else:
                perturbation_values = scale * np.asarray(
                    getattr(perturbation, name)[surface], dtype=float
                )
            getattr(combined, name).append(
                np.concatenate((axis_values, perturbation_values))
            )

    combined.m0b = max(np.max(np.abs(modes)) for modes in combined.m)
    combined.n0b = max(np.max(np.abs(modes)) for modes in combined.n)
    return combined


def combine_neo2_boozer_files(axis, perturbation=None, *, scale=1.0,
                               geometry="axisymmetric",
                               source_phase="neo2-plus"):
    """Read and compose an axisymmetric NEO-2 .bc and optional perturbation."""
    from libneo.boozer import BoozerFile

    axis_data = BoozerFile(str(axis))
    perturbation_data = (
        BoozerFile(str(perturbation)) if perturbation is not None else None
    )
    return combine_boozer_data(
        axis_data,
        perturbation_data,
        scale=scale,
        geometry=geometry,
        source_phase=source_phase,
    )


def convert_neo2bc_to_chartmap(
    axis,
    output,
    *,
    perturbation=None,
    scale=1.0,
    geometry="axisymmetric",
    nrho=50,
    ntheta=48,
    nzeta=96,
    covariant_sign=1,
    source_phase="neo2-plus",
):
    """Write a composed NEO-2 field as an extended Boozer chartmap."""
    from libneo.bc_to_booz_xform import write_boozmn
    from libneo.booz_xform_to_boozer_chartmap import convert_boozmn_to_chartmap

    axis = Path(axis)
    perturbation_path = Path(perturbation) if perturbation is not None else None
    combined = combine_neo2_boozer_files(
        axis,
        perturbation_path,
        scale=scale,
        geometry=geometry,
        source_phase=source_phase,
    )
    attrs = {
        "booz2chartmap_source": str(axis),
        "neo2bc_axis_source": str(axis),
        "neo2bc_axis_sha256": _sha256(axis),
        "neo2bc_fourier_phase": _OUTPUT_PHASE,
        "neo2bc_source_fourier_phase": _SOURCE_PHASES[source_phase],
        "neo2bc_output_fourier_phase": _OUTPUT_PHASE,
        "neo2bc_source_length_unit": "m",
        "neo2bc_source_field_unit": "T",
        "neo2bc_output_unit_system": "CGS-Gaussian",
        "neo2bc_perturbation_scale": float(scale),
        "neo2bc_axis_phase_scale": float(combined.nper / combined.axis_nper),
        "neo2bc_geometry": geometry,
        "neo2bc_s_min": float(combined.s[0]),
        "neo2bc_s_max": float(combined.s[-1]),
        "neo2bc_radial_matching": "equilibrium half-grid on common support; cubic",

    }
    if perturbation_path is not None:
        attrs.update(
            neo2bc_perturbation_source=str(perturbation_path),
            neo2bc_perturbation_sha256=_sha256(perturbation_path),
        )

    with tempfile.TemporaryDirectory(prefix="libneo-neo2bc-") as directory:
        boozmn = Path(directory) / "composed.boozmn.nc"
        write_boozmn(combined, boozmn, source=axis)
        return convert_boozmn_to_chartmap(
            boozmn,
            output,
            nrho=nrho,
            ntheta=ntheta,
            nzeta=nzeta,
            covariant_sign=covariant_sign,
            chartmap_attrs=attrs,
            netcdf_format="NETCDF3_64BIT_OFFSET",
        )


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Convert NEO-2 axisymmetric and perturbation .bc files "
        "to an extended Boozer chartmap"
    )
    parser.add_argument("axis", help="axisymmetric NEO-2 .bc file")
    parser.add_argument("output", help="output chartmap NetCDF")
    parser.add_argument("--perturbation", help="separate NEO-2 perturbation .bc")
    parser.add_argument("--scale", type=float, default=1.0)
    parser.add_argument(
        "--geometry", choices=("axisymmetric", "full"), default="axisymmetric"
    )
    parser.add_argument("--nrho", type=int, default=50)
    parser.add_argument("--ntheta", type=int, default=48)
    parser.add_argument("--nzeta", type=int, default=96)
    parser.add_argument("--covariant-sign", type=int, choices=(-1, 1), default=1)
    parser.add_argument(
        "--source-phase",
        choices=tuple(_SOURCE_PHASES),
        default="neo2-plus",
        help="NEO-2 plus phase is the safe default; use minus only for a "
        "source known to use the chartmap sign",
    )
    args = parser.parse_args(argv)

    convert_neo2bc_to_chartmap(
        args.axis,
        args.output,
        perturbation=args.perturbation,
        scale=args.scale,
        geometry=args.geometry,
        nrho=args.nrho,
        ntheta=args.ntheta,
        nzeta=args.nzeta,
        covariant_sign=args.covariant_sign,
        source_phase=args.source_phase,
    )


if __name__ == "__main__":
    main()
