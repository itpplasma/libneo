"""Checked axisymmetric G-EQDSK export to the asymmetric LHS Boozer format."""
from pathlib import Path
import os
import tempfile

import numpy as np


def convert_eqdsk_to_bc(gfile, output, *, nsurf=64, mpol=32, n_theta=1024,
                        n_q=128, route=1, psimax=None, exp_bpol=-1):
    """Write a physical LHS ``.bc`` file from checked native contour data.

    Surfaces are the midpoints of ``nsurf`` uniform normalized toroidal-flux
    intervals. Routes 1 and 2 construct Boozer angles. ``psimax`` optionally
    selects an axis-relative poloidal flux boundary in CGS, while ``exp_bpol``
    is -1 to infer flux units, 0 for Wb/rad, or 1 for Wb.

    Native equilibrium checks are mandatory. Nonfinite or rejected data never
    replace an existing output. Successful output uses 16 decimal places and
    is replaced atomically in its destination directory. Return its Path.
    """
    from efit_to_boozer.boozer import get_contour_transform_extended
    from .boozer import BoozerFile

    for name, value, minimum in (("nsurf", nsurf, 1), ("mpol", mpol, 0),
                                  ("n_theta", n_theta, 16), ("n_q", n_q, 6)):
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value < minimum:
            raise ValueError(f"{name} must be an integer >= {minimum}")
    if isinstance(route, (bool, np.bool_)) or not isinstance(route, (int, np.integer)) or route not in (1, 2):
        raise ValueError("Boozer .bc export requires route 1 or 2")
    gfile, output = Path(gfile).resolve(), Path(output)
    if output.resolve() == gfile:
        raise ValueError("Output must differ from the input G-EQDSK file")
    s = (np.arange(nsurf, dtype=float) + 0.5) / nsurf
    data = get_contour_transform_extended(gfile, s, mpol, n_theta=n_theta,
                                          n_q=n_q, route=route, psimax=psimax,
                                          exp_bpol=exp_bpol)
    if data["status_code"] != 0 or data["overridden"]:
        raise ValueError(data["message"] or "A checked equilibrium is required for .bc export")
    physical = ("spectra", "iota", "curr_i", "curr_g", "axis", "psitor_max",
                "pp_file", "sqrtg00", "rminor", "sig_theta")
    if any(not np.all(np.isfinite(data[key])) for key in physical):
        raise ValueError("Nonfinite native metadata cannot be written to a .bc file")
    if data["sig_theta"] not in (-1, 1) or data["rminor"] <= 0 or data["psitor_max"] == 0:
        raise ValueError("Invalid native geometry or flux metadata for .bc export")

    result = BoozerFile("")
    result.m0b, result.n0b, result.nsurf, result.nper = mpol, 0, nsurf, 1
    try:
        with np.errstate(over="raise", invalid="raise"):
            result.flux = data["sig_theta"] * data["psitor_max"] * (2*np.pi*1e-8)
            result.a, result.R = data["rminor"] * 0.01, data["axis"][0]
            result.s, result.iota = s, data["iota"]
            result.Jpol_divided_by_nper = -5*data["curr_g"]
            result.Itor = -5*data["curr_i"]
            result.pprime = 0.1*data["pp_file"]
            result.sqrt_g_00 = -np.abs(data["sqrtg00"]*data["psitor_max"])*1e-6*(2*np.pi)**2
            spectra = data["spectra"].copy()
            spectra[:4] *= 0.01
            spectra[6:] *= 1e-4
    except FloatingPointError as error:
        raise ValueError("Nonfinite physical values cannot be written to a .bc file") from error
    converted = (result.flux, result.a, result.R, result.iota,
                 result.Jpol_divided_by_nper, result.Itor, result.pprime,
                 result.sqrt_g_00, spectra)
    if any(not np.all(np.isfinite(value)) for value in converted):
        raise ValueError("Nonfinite physical values cannot be written to a .bc file")
    result.m = np.broadcast_to(np.arange(mpol + 1), (nsurf, mpol + 1))
    result.n = np.zeros_like(result.m)
    for name, coefficients in zip(("rmnc", "rmns", "zmnc", "zmns", "vmnc", "vmns", "bmnc", "bmns"),
                                  spectra):
        setattr(result, name, coefficients.T)

    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=output.parent, prefix=f".{output.name}.",
                                          suffix=".tmp", delete=False) as stream:
            temporary = Path(stream.name)
        result.write(temporary, precision=16)
        os.replace(temporary, output)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return output
