"""Convert a GVEC state to a libneo Boozer chartmap.

GVEC evaluates the field and geometry directly in straight-field-line Boozer
coordinates (``state.evaluate_sfl(sfl="boozer")``); this front end repackages
that output, converts SI -> CGS, applies the left-handed flip SIMPLE expects,
and calls write_boozer_chartmap to emit the extended chartmap. The Boozer
transform itself lives in the gvec library, so GVEC reaches the libneo hub via
the chartmap format (the same as booz_xform).

gvec is imported lazily inside convert_gvec_to_chartmap so this module imports
without the optional dependency.
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


def assemble_chartmap(
    output,
    *,
    rho,
    s,
    theta,
    zeta,
    X,
    Y,
    Z,
    Bmod,
    B_theta,
    B_phi,
    A_phi,
    torflux,
    num_field_periods,
    flip,
    **attrs,
):
    """Convert SI Boozer field/geometry to CGS, flip handedness, and write.

    Inputs are in SI (T, m, T*m, T*m^2) on a (n_rho, n_theta, n_zeta) grid.
    ``flip`` is "tor" or "pol": GVEC is right-handed, SIMPLE left-handed, so a
    toroidal flip negates the toroidal surface functions (B_phi, A_phi) and a
    poloidal flip negates the poloidal one (B_theta) and the toroidal flux.
    The angle grid is assumed already built against the flipped angle (as in
    GVEC's evaluate_sfl call). gvec is not needed by this function.
    """
    B_theta = np.asarray(B_theta, dtype=float)
    B_phi = np.asarray(B_phi, dtype=float)
    A_phi = np.asarray(A_phi, dtype=float)
    torflux = float(torflux)

    if flip == "tor":
        B_phi = -B_phi
        A_phi = -A_phi
    elif flip == "pol":
        B_theta = -B_theta
        torflux = -torflux
    else:
        raise ValueError(f"invalid flip option: {flip!r} (expected 'tor' or 'pol')")

    write_boozer_chartmap(
        output,
        rho=rho,
        s=s,
        theta=theta,
        zeta=zeta,
        x=np.asarray(X, dtype=float) * METER_TO_CM,
        y=np.asarray(Y, dtype=float) * METER_TO_CM,
        z=np.asarray(Z, dtype=float) * METER_TO_CM,
        A_phi=A_phi * TESLA_METER2_TO_GAUSS_CM2,
        B_theta=B_theta * TESLA_METER_TO_GAUSS_CM,
        B_phi=B_phi * TESLA_METER_TO_GAUSS_CM,
        Bmod=np.asarray(Bmod, dtype=float) * TESLA_TO_GAUSS,
        num_field_periods=num_field_periods,
        torflux=torflux * TESLA_METER2_TO_GAUSS_CM2,
        gvec2chartmap_flip=flip,
        **attrs,
    )
    return torflux * TESLA_METER2_TO_GAUSS_CM2


def convert_gvec_to_chartmap(
    paramfile,
    statefile,
    output,
    *,
    nrho=50,
    ntheta=36,
    nphi=81,
    boozer_factor=1,
    Bcov="avg",
):
    """Load a GVEC state and write a libneo Boozer chartmap.

    gvec is imported lazily so this module imports without the dependency.
    The geometry grid is built against the toroidally flipped angle, matching
    SIMPLE's GVEC converter (left-handed coordinates).
    """
    import gvec

    flip = "tor"
    state = gvec.load_state(paramfile, statefile)
    nfp = int(state.nfp)
    phi_period = TWOPI / nfp

    rho_grid = np.linspace(1.0e-3, 1.0, nrho)
    s = np.linspace(rho_grid[0] ** 2, 1.0, nrho)
    rho_profile = np.sqrt(s)
    theta_geom = np.linspace(0.0, TWOPI, ntheta, endpoint=False)
    zeta_geom = np.linspace(0.0, phi_period, nphi, endpoint=False)

    ev = state.evaluate_sfl(
        "mod_B",
        "B_theta_B",
        "B_zeta_B",
        "Phi_edge",
        "pos",
        rho=rho_grid,
        theta=-theta_geom,
        zeta=-zeta_geom,
        sfl="boozer",
        MNfactor=boozer_factor,
    )
    ev = ev.transpose("xyz", "rad", "pol", "tor")
    Bmod = ev.mod_B.values

    if Bcov == "avg":
        ev_bcov = state.evaluate("B_theta_avg", "B_zeta_avg", rho=rho_grid)
        B_theta = ev_bcov.B_theta_avg.values
        B_phi = ev_bcov.B_zeta_avg.values
    elif Bcov == "boozer-avg":
        B_theta = ev.B_theta_B.mean(["pol", "tor"]).values
        B_phi = ev.B_zeta_B.mean(["pol", "tor"]).values
    elif Bcov == "boozer-0":
        B_theta = ev.B_theta_B.values[:, 0, 0]
        B_phi = ev.B_zeta_B.values[:, 0, 0]
    else:
        raise ValueError(f"invalid Bcov option: {Bcov!r}")

    ev_aphi = state.evaluate_sfl(
        "chi",
        rho=rho_profile,
        theta=np.array([0.0]),
        zeta=np.array([0.0]),
        sfl="boozer",
        MNfactor=boozer_factor,
    )
    A_phi = -np.asarray(ev_aphi.chi.values).reshape(nrho, -1)[:, 0]

    torflux = ev.Phi_edge.item()
    pos = ev.pos.values
    X, Y, Z = pos[0], pos[1], pos[2]

    return assemble_chartmap(
        output,
        rho=rho_grid,
        s=s,
        theta=theta_geom,
        zeta=zeta_geom,
        X=X,
        Y=Y,
        Z=Z,
        Bmod=Bmod,
        B_theta=B_theta,
        B_phi=B_phi,
        A_phi=A_phi,
        torflux=torflux,
        num_field_periods=nfp,
        flip=flip,
        gvec2chartmap_boozer_factor=int(boozer_factor),
        gvec2chartmap_Bcov_method=Bcov,
    )


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Convert a GVEC state to a libneo Boozer chartmap NetCDF file"
    )
    parser.add_argument("paramfile", help="GVEC parameter file (.ini)")
    parser.add_argument("statefile", help="GVEC state file (.dat)")
    parser.add_argument("output", help="Output Boozer chartmap .nc file")
    parser.add_argument("--nrho", type=int, default=50)
    parser.add_argument("--ntheta", type=int, default=36)
    parser.add_argument("--nphi", type=int, default=81)
    parser.add_argument("--boozer-factor", type=int, default=1)
    parser.add_argument(
        "--Bcov",
        choices=["avg", "boozer-avg", "boozer-0"],
        default="avg",
        help="method for the B_theta/B_phi surface functions",
    )
    args = parser.parse_args(argv)

    torflux = convert_gvec_to_chartmap(
        args.paramfile,
        args.statefile,
        args.output,
        nrho=args.nrho,
        ntheta=args.ntheta,
        nphi=args.nphi,
        boozer_factor=args.boozer_factor,
        Bcov=args.Bcov,
    )
    print(f"Done. torflux={torflux:.6e} G cm^2")


if __name__ == "__main__":
    main()
