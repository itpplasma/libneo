"""Convert a Strumberger Boozer .bc file to a booz_xform-style boozmn NetCDF.

The output variable layout matches what _read_boozmn in
booz_xform_to_boozer_chartmap.py consumes, so the produced NetCDF can be
passed directly to convert_boozmn_to_chartmap.

Radial staggering convention (matches booz_xform output):
  Full grid (ns points):  iota_b, buco_b, bvco_b, phi_b
  Half grid (ns-2 points): bmnc_b, rmnc_b, zmns_b, pmns_b
  jlist[k] = k+2  (1-based, half-grid surface k sits between full-grid k+1 and k+2)

Unit conversions:
  .bc curr_pol/nper [A] -> bvco_b [T*m] = curr_pol/nper * nfp * mu_0 / (2*pi)
  .bc curr_tor [A]      -> buco_b [T*m] = curr_tor * mu_0 / (2*pi)
  .bc vmns (phase) = (phib-phi)*nper/(2*pi) -> pmns_b [rad] = vmns * 2*pi/nfp
  .bc bmnc [T]          -> bmnc_b [T]  (no conversion needed)
  .bc rmnc, zmns [m]    -> rmnc_b, zmns_b [m] (no conversion needed)
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

MU0 = 4.0e-7 * np.pi  # T*m/A
TWOPI = 2.0 * np.pi


def _bc_to_mode_arrays(bc):
    """Return flat (ixm, ixn) and (nsurf, nmn) coefficient arrays from BoozerFile.

    Mode ordering: all modes from the first surface (assumed identical across
    all surfaces). For axisymmetric equilibria (n0b=0) only n=0 modes appear.
    """
    m0 = np.array(bc.m[0], dtype=int)
    n0 = np.array(bc.n[0], dtype=int)
    nmn = len(m0)
    nsurf = bc.nsurf

    rmnc = np.zeros((nsurf, nmn))
    zmns = np.zeros((nsurf, nmn))
    pmns = np.zeros((nsurf, nmn))
    bmnc = np.zeros((nsurf, nmn))
    lasym = False
    rmns = np.zeros((nsurf, nmn)) if lasym else None
    zmnc = np.zeros((nsurf, nmn)) if lasym else None
    bmns = np.zeros((nsurf, nmn)) if lasym else None

    for k in range(nsurf):
        rmnc[k] = np.array(bc.rmnc[k], dtype=float)
        zmns[k] = np.array(bc.zmns[k], dtype=float)
        # vmns in .bc = (phib-phi)*nfp/(2*pi); pmns_b = phib-phi in radians
        pmns[k] = np.array(bc.vmns[k], dtype=float) * TWOPI / bc.nper
        bmnc[k] = np.array(bc.bmnc[k], dtype=float)

        r_s = np.array(bc.rmns[k], dtype=float)
        z_c = np.array(bc.zmnc[k], dtype=float)
        b_s = np.array(bc.bmns[k], dtype=float)
        if np.any(r_s != 0.0) or np.any(z_c != 0.0) or np.any(b_s != 0.0):
            lasym = True

    if lasym:
        rmns = np.zeros((nsurf, nmn))
        zmnc = np.zeros((nsurf, nmn))
        bmns = np.zeros((nsurf, nmn))
        for k in range(nsurf):
            rmns[k] = np.array(bc.rmns[k], dtype=float)
            zmnc[k] = np.array(bc.zmnc[k], dtype=float)
            bmns[k] = np.array(bc.bmns[k], dtype=float)

    return m0, n0, rmnc, zmns, pmns, bmnc, rmns, zmnc, bmns, lasym


def write_boozmn(bc, output, source=None):
    """Write a boozmn NetCDF from a BoozerFile object.

    Parameters
    ----------
    bc : libneo.boozer.BoozerFile
        Parsed .bc file.
    output : path-like
        Destination NetCDF path.
    source : path-like, optional
        Original .bc path recorded in the bc2boozmn_source provenance attribute.
    """
    import netCDF4

    nsurf = bc.nsurf
    nfp = bc.nper
    ns = nsurf + 2  # axis (j=1) + nsurf half-grid + LCFS (j=ns)

    # Half-grid surface indices (1-based): 2 .. ns-1
    jlist = np.arange(2, ns, dtype=int)  # shape (nsurf,)

    # Build mode arrays from .bc.
    m0, n0, rmnc_h, zmns_h, pmns_h, bmnc_h, rmns_h, zmnc_h, bmns_h, lasym = (
        _bc_to_mode_arrays(bc)
    )
    ixm = m0.astype(np.int32)
    ixn = (n0 * nfp).astype(np.int32)  # booz_xform stores n*nfp in ixn_b

    # Surface functions on half grid from .bc.
    s_half = np.array(bc.s, dtype=float)  # shape (nsurf,)
    iota_h = np.array(bc.iota, dtype=float)
    bvco_h = np.array(bc.Jpol_divided_by_nper, dtype=float) * nfp * MU0 / TWOPI
    buco_h = np.array(bc.Itor, dtype=float) * MU0 / TWOPI

    # Full-grid arrays by linear extrapolation from half grid (axis + LCFS).
    from scipy.interpolate import CubicSpline

    spl_iota = CubicSpline(s_half, iota_h, extrapolate=True)
    spl_bvco = CubicSpline(s_half, bvco_h, extrapolate=True)
    spl_buco = CubicSpline(s_half, buco_h, extrapolate=True)

    # Full-grid s: 0, 1/(ns-1), 2/(ns-1), ..., 1
    s_full = np.linspace(0.0, 1.0, ns)
    iota_full = spl_iota(s_full)
    bvco_full = spl_bvco(s_full)
    buco_full = spl_buco(s_full)
    iota_full[0] = spl_iota(0.0)
    bvco_full[0] = spl_bvco(0.0)
    buco_full[0] = 0.0   # zero poloidal current at axis

    # Toroidal flux: bc.flux is the edge value in T*m^2.
    phi_full = s_full * float(bc.flux)

    with netCDF4.Dataset(str(output), "w", format="NETCDF4") as ds:
        # Dimensions
        ds.createDimension("radius", ns)
        ds.createDimension("mn_mode", len(ixm))
        ds.createDimension("comput_surfs", nsurf)

        def _var(name, dtype, dims, data, **kw):
            v = ds.createVariable(name, dtype, dims)
            v[:] = data
            for k, val in kw.items():
                setattr(v, k, val)
            return v

        _var("ns_b", "i4", (), ns)
        _var("nfp_b", "i4", (), nfp)
        _var("mboz_b", "i4", (), int(np.max(m0)))
        _var("nboz_b", "i4", (), int(np.max(np.abs(n0))))
        _var("mnboz_b", "i4", (), len(ixm))
        _var("lasym__logical__", "i4", (), int(lasym))

        _var("jlist", "i4", ("comput_surfs",), jlist)
        _var("ixm_b", "i4", ("mn_mode",), ixm)
        _var("ixn_b", "i4", ("mn_mode",), ixn)

        _var("iota_b", "f8", ("radius",), iota_full)
        _var("buco_b", "f8", ("radius",), buco_full)
        _var("bvco_b", "f8", ("radius",), bvco_full)
        _var("phi_b", "f8", ("radius",), phi_full)

        _var("bmnc_b", "f8", ("comput_surfs", "mn_mode"), bmnc_h)
        _var("rmnc_b", "f8", ("comput_surfs", "mn_mode"), rmnc_h)
        _var("zmns_b", "f8", ("comput_surfs", "mn_mode"), zmns_h)
        _var("pmns_b", "f8", ("comput_surfs", "mn_mode"), pmns_h)

        if lasym:
            _var("bmns_b", "f8", ("comput_surfs", "mn_mode"), bmns_h)
            _var("rmns_b", "f8", ("comput_surfs", "mn_mode"), rmns_h)
            _var("zmnc_b", "f8", ("comput_surfs", "mn_mode"), zmnc_h)

        ds.bc2boozmn_source = str(source if source is not None else output)


def convert_bc_to_boozmn(bc_file, output):
    """Read a .bc file and write a boozmn NetCDF."""
    from libneo.boozer import BoozerFile

    bc = BoozerFile(str(bc_file))
    write_boozmn(bc, output, source=bc_file)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Convert a Strumberger Boozer .bc file to a booz_xform boozmn NetCDF"
    )
    parser.add_argument("bc_file", help="Input .bc file")
    parser.add_argument("output", help="Output boozmn NetCDF file")
    args = parser.parse_args(argv)

    convert_bc_to_boozmn(args.bc_file, args.output)
    print(f"Done. Output: {args.output}")


if __name__ == "__main__":
    main()
