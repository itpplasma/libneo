"""Convert booz_xform output (boozmn*.nc) to a libneo Boozer chartmap.

Reads the Boozer-coordinate Fourier harmonics produced by the booz_xform
library, splines them to a uniform rho grid, Fourier-sums Bmod/R/Z/phi onto a
tensor angle grid, recovers the toroidal flux, and calls
write_boozer_chartmap to emit the extended chartmap SIMPLE reads at runtime.

The boozmn file alone is sufficient: the toroidal flux is taken from phi_b
when present and otherwise recovered from the surface geometry via the Boozer
identity sqrt(g) B^2 = psi' (G + iota I).

Radial staggering: bmnc_b/rmnc_b/zmns_b/pmns_b live on the half grid
(surfaces jlist), iota_b/buco_b/bvco_b/phi_b on the full grid; index jlist[k]-1
holds the matching half-grid value for surface k.

Half-grid s: s_half[k] = (jlist[k] - 1.5) / (ns_b - 1), the midpoint of the
interval between full-grid surfaces jlist[k]-1 and jlist[k] (1-based).
"""

from __future__ import annotations

import argparse

import numpy as np

from libneo.boozer_chartmap_writer import (
    METER_TO_CM,
    TESLA_METER2_TO_GAUSS_CM2,
    TESLA_METER_TO_GAUSS_CM,
    TESLA_TO_GAUSS,
    write_boozer_chartmap,
)

TWOPI = 2.0 * np.pi


def _read_boozmn(filename):
    """Read a boozmn file into a plain dict of numpy arrays."""
    import netCDF4

    d = {}
    with netCDF4.Dataset(filename) as ds:

        def var(n):
            return np.asarray(ds.variables[n][:])

        d["ns"] = int(var("ns_b"))
        d["nfp"] = int(var("nfp_b"))
        d["lasym"] = bool(var("lasym__logical__"))
        d["jlist"] = var("jlist").astype(int)
        d["ixm"] = var("ixm_b").astype(int)
        d["ixn"] = var("ixn_b").astype(int)
        d["iota"] = var("iota_b")
        d["buco"] = var("buco_b")
        d["bvco"] = var("bvco_b")
        d["phi"] = var("phi_b")
        d["aminor_m"] = (
            float(ds.libneo_vmec_aminor_m)
            if "libneo_vmec_aminor_m" in ds.ncattrs()
            else None
        )
        d["rmajor_m"] = (
            float(ds.libneo_vmec_rmajor_m)
            if "libneo_vmec_rmajor_m" in ds.ncattrs()
            else None
        )
        d["bmnc"] = var("bmnc_b")
        d["rmnc"] = var("rmnc_b")
        d["zmns"] = var("zmns_b")
        d["pmns"] = var("pmns_b")
        if d["lasym"]:
            d["bmns"] = var("bmns_b")
            d["rmns"] = var("rmns_b")
            d["zmnc"] = var("zmnc_b")
            d["pmnc"] = var("pmnc_b")
    d["s_half"] = (d["jlist"] - 1.5) / (d["ns"] - 1)
    return d


def _interp_coeffs(coeffs, ixm, rho_half, rho_out):
    """Interpolate packed Fourier coefficients from the half grid to rho_out.

    Cubic spline in rho per mode; below the innermost half-grid surface the
    poloidal mode number sets the near-axis behaviour c ~ rho^m, so modes are
    continued with a power law instead of the spline's polynomial tail.
    """
    from scipy.interpolate import CubicSpline

    out = CubicSpline(rho_half, coeffs, axis=0)(rho_out)
    inner = rho_out < rho_half[0]
    if np.any(inner):
        ratio = rho_out[inner, None] / rho_half[0]
        m = np.minimum(ixm[None, :], 50)
        out[inner, :] = coeffs[0][None, :] * ratio**m
    return out


def _fourier_eval(coeffs_grid, ixm, ixn, theta, zeta, parity):
    """Sum coeffs(rho, mn) * trig(m*theta - n*zeta) on a tensor angle grid.

    parity is 'cos' or 'sin'. Returns an (n_rho, n_theta, n_zeta) array.
    """
    angle = (
        ixm[:, None, None] * theta[None, :, None]
        - ixn[:, None, None] * zeta[None, None, :]
    )
    basis = np.cos(angle) if parity == "cos" else np.sin(angle)
    return np.tensordot(coeffs_grid, basis, axes=([1], [0]))


def _derive_psi_prime(d):
    """Recover psi' = Phi'(s)/(2*pi) from geometry when phi_b is absent.

    On each surface sqrt(g) B^2 = psi' (G + iota I) with sqrt(g) the Jacobian
    of (s, theta_B, zeta_B), so the angle average of J B^2 / (G + iota I)
    gives psi'. The spread over mid-radius surfaces is a consistency check.
    """
    from scipy.interpolate import CubicSpline

    ixm, ixn = d["ixm"], d["ixn"]
    s_half = d["s_half"]
    mid = np.where((s_half > 0.25) & (s_half < 0.95))[0]
    if len(mid) < 3:
        mid = np.arange(1, len(s_half) - 1)
    nth, nze = 73, 73
    theta = np.linspace(0.0, TWOPI, nth, endpoint=False)
    zeta = np.linspace(0.0, TWOPI / d["nfp"], nze, endpoint=False)
    angle = (
        ixm[:, None, None] * theta[None, :, None]
        - ixn[:, None, None] * zeta[None, None, :]
    )
    cosg, sing = np.cos(angle), np.sin(angle)

    def ev(c, basis):
        return np.tensordot(c, basis, axes=([0], [0]))

    names = ["rmnc", "zmns", "pmns"]
    if d["lasym"]:
        names += ["rmns", "zmnc", "pmnc"]
    ds_coeffs = {
        name: CubicSpline(s_half, d[name], axis=0)(s_half[mid], nu=1)
        for name in names
    }

    psi_primes = []
    for i, k in enumerate(mid):
        R = ev(d["rmnc"][k], cosg)
        B = ev(d["bmnc"][k], cosg)
        R_s = ev(ds_coeffs["rmnc"][i], cosg)
        Z_s = ev(ds_coeffs["zmns"][i], sing)
        phi_s = ev(ds_coeffs["pmns"][i], sing)
        R_t = ev(-ixm * d["rmnc"][k], sing)
        Z_t = ev(ixm * d["zmns"][k], cosg)
        phi_t = ev(ixm * d["pmns"][k], cosg)
        R_z = ev(ixn * d["rmnc"][k], sing)
        Z_z = ev(-ixn * d["zmns"][k], cosg)
        phi_z = 1.0 + ev(-ixn * d["pmns"][k], cosg)
        if d["lasym"]:
            R += ev(d["rmns"][k], sing)
            B += ev(d["bmns"][k], sing)
            R_s += ev(ds_coeffs["rmns"][i], sing)
            Z_s += ev(ds_coeffs["zmnc"][i], cosg)
            phi_s += ev(ds_coeffs["pmnc"][i], cosg)
            R_t += ev(ixm * d["rmns"][k], cosg)
            Z_t += ev(-ixm * d["zmnc"][k], sing)
            phi_t += ev(-ixm * d["pmnc"][k], sing)
            R_z += ev(-ixn * d["rmns"][k], cosg)
            Z_z += ev(ixn * d["zmnc"][k], sing)
            phi_z += ev(ixn * d["pmnc"][k], sing)

        J = R * (
            R_s * (phi_t * Z_z - phi_z * Z_t)
            + phi_s * (Z_t * R_z - Z_z * R_t)
            + Z_s * (R_t * phi_z - R_z * phi_t)
        )

        jidx = d["jlist"][k] - 1
        G, I, iota = d["bvco"][jidx], d["buco"][jidx], d["iota"][jidx]
        psi_primes.append(np.mean(J * B**2) / (G + iota * I))

    psi_primes = np.array(psi_primes)
    spread = np.ptp(psi_primes) / np.abs(np.median(psi_primes))
    if spread > 1.0e-3:
        raise RuntimeError(
            f"psi' from geometry varies by {spread:.2e} across surfaces; "
            "boozmn file is inconsistent"
        )
    return float(np.median(psi_primes))


def convert_boozmn_to_chartmap(
    boozmn,
    output,
    *,
    nrho=50,
    ntheta=48,
    nzeta=96,
    covariant_sign=1,
):
    """Read a boozmn file and write a libneo Boozer chartmap.

    booz_xform/scipy are imported lazily so this module imports without the
    optional dependencies.
    """
    from scipy.interpolate import CubicSpline

    if covariant_sign not in (-1, 1):
        raise ValueError("covariant_sign must be -1 or 1")

    d = _read_boozmn(boozmn)
    nfp = d["nfp"]
    j = d["jlist"] - 1
    s_half = d["s_half"]
    rho_half = np.sqrt(s_half)

    iota_h = d["iota"][j]
    buco_h = d["buco"][j]
    bvco_h = d["bvco"][j]

    if np.any(d["phi"] != 0.0):
        torflux_si = -float(d["phi"][-1]) / TWOPI
    else:
        torflux_si = _derive_psi_prime(d)

    rho_grid = np.linspace(1.0e-3, 1.0, nrho)
    s_grid = rho_grid**2
    s = np.linspace(rho_grid[0] ** 2, 1.0, nrho)
    theta_geom = np.linspace(0.0, TWOPI, ntheta, endpoint=False)
    zeta_geom = np.linspace(0.0, TWOPI / nfp, nzeta, endpoint=False)

    iota_spline = CubicSpline(s_half, iota_h)
    B_theta = covariant_sign * CubicSpline(s_half, buco_h)(s_grid)
    B_phi = covariant_sign * CubicSpline(s_half, bvco_h)(s_grid)
    iota_int = iota_spline.antiderivative()
    A_phi = -torflux_si * (iota_int(s) - iota_int(0.0))

    bmnc_g = _interp_coeffs(d["bmnc"], d["ixm"], rho_half, rho_grid)
    Bmod = _fourier_eval(bmnc_g, d["ixm"], d["ixn"], theta_geom, zeta_geom, "cos")
    if d["lasym"]:
        bmns_g = _interp_coeffs(d["bmns"], d["ixm"], rho_half, rho_grid)
        Bmod += _fourier_eval(
            bmns_g, d["ixm"], d["ixn"], theta_geom, zeta_geom, "sin"
        )

    rmnc_g = _interp_coeffs(d["rmnc"], d["ixm"], rho_half, rho_grid)
    zmns_g = _interp_coeffs(d["zmns"], d["ixm"], rho_half, rho_grid)
    pmns_g = _interp_coeffs(d["pmns"], d["ixm"], rho_half, rho_grid)
    R = _fourier_eval(rmnc_g, d["ixm"], d["ixn"], theta_geom, zeta_geom, "cos")
    Z = _fourier_eval(zmns_g, d["ixm"], d["ixn"], theta_geom, zeta_geom, "sin")
    p = _fourier_eval(pmns_g, d["ixm"], d["ixn"], theta_geom, zeta_geom, "sin")
    if d["lasym"]:
        rmns_g = _interp_coeffs(d["rmns"], d["ixm"], rho_half, rho_grid)
        zmnc_g = _interp_coeffs(d["zmnc"], d["ixm"], rho_half, rho_grid)
        pmnc_g = _interp_coeffs(d["pmnc"], d["ixm"], rho_half, rho_grid)
        R += _fourier_eval(rmns_g, d["ixm"], d["ixn"], theta_geom, zeta_geom, "sin")
        Z += _fourier_eval(zmnc_g, d["ixm"], d["ixn"], theta_geom, zeta_geom, "cos")
        p += _fourier_eval(pmnc_g, d["ixm"], d["ixn"], theta_geom, zeta_geom, "cos")

    phi_cyl = zeta_geom[None, None, :] + p
    X = R * np.cos(phi_cyl)
    Y = R * np.sin(phi_cyl)

    Bmod = Bmod * TESLA_TO_GAUSS
    B_theta = B_theta * TESLA_METER_TO_GAUSS_CM
    B_phi = B_phi * TESLA_METER_TO_GAUSS_CM
    A_phi = A_phi * TESLA_METER2_TO_GAUSS_CM2
    torflux = torflux_si * TESLA_METER2_TO_GAUSS_CM2
    X, Y, Z = X * METER_TO_CM, Y * METER_TO_CM, Z * METER_TO_CM

    attrs = {
        "booz2chartmap_source": str(boozmn),
        "booz2chartmap_covariant_sign": np.int32(covariant_sign),
    }
    if d["aminor_m"] is not None:
        attrs["aminor_m"] = d["aminor_m"]
    if d["rmajor_m"] is not None:
        attrs["vmec_rmajor_m"] = d["rmajor_m"]

    write_boozer_chartmap(
        output,
        rho=rho_grid,
        s=s,
        theta=theta_geom,
        zeta=zeta_geom,
        x=X,
        y=Y,
        z=Z,
        A_phi=A_phi,
        B_theta=B_theta,
        B_phi=B_phi,
        Bmod=Bmod,
        num_field_periods=nfp,
        torflux=torflux,
        **attrs,
    )
    return torflux


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Convert booz_xform boozmn output to a libneo Boozer "
        "chartmap NetCDF file"
    )
    parser.add_argument("boozmn", help="booz_xform output file (boozmn*.nc)")
    parser.add_argument("output", help="Output Boozer chartmap .nc file")
    parser.add_argument("--nrho", type=int, default=50)
    parser.add_argument(
        "--ntheta",
        type=int,
        default=48,
        help="poloidal points, endpoint-excluded geometry grid",
    )
    parser.add_argument(
        "--nzeta",
        type=int,
        default=96,
        help="toroidal points per field period, endpoint-excluded",
    )
    parser.add_argument(
        "--covariant-sign",
        type=int,
        choices=(-1, 1),
        default=1,
        help=(
            "explicit common sign applied to B_theta and B_phi; use -1 "
            "when mapping native VMEC/booz_xform covariants to a target "
            "convention with the opposite current signs"
        ),
    )
    args = parser.parse_args(argv)

    torflux = convert_boozmn_to_chartmap(
        args.boozmn,
        args.output,
        nrho=args.nrho,
        ntheta=args.ntheta,
        nzeta=args.nzeta,
        covariant_sign=args.covariant_sign,
    )
    print(f"Done. torflux={torflux:.6e} G cm^2")


if __name__ == "__main__":
    main()
