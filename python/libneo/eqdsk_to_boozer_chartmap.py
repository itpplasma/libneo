"""Convert an EQDSK/g-file axisymmetric equilibrium to a libneo Boozer chartmap.

Uses the efit_to_boozer Fortran/f2py module (_efit_to_boozer) that is already
built at runtime (see build/_efit_to_boozer*.so).  The module reads the
efit_to_boozer.inp config file from the current directory, which in turn points
at the g-file.  This front end writes that .inp file from the arguments and
then calls the module.

Output convention
-----------------
Axisymmetric: num_field_periods=1, nzeta=1.  Geometry varies with
phi_cyl = zeta - G_of_thb(theta_B); Bmod is independent of zeta.
Surface functions:
  B_phi (G*cm) = C_norm (covariant toroidal B from magdata, CGS)
  B_theta (G*cm) = C_norm / q
  A_phi(s) = -torflux * integral_0^s iota(s') ds'
torflux = -psitor_max (negative: left-handed SIMPLE convention)
"""

from __future__ import annotations

import argparse
import os
import tempfile
from pathlib import Path

import numpy as np

from libneo.boozer_chartmap_writer import (
    METER_TO_CM,
    TESLA_METER2_TO_GAUSS_CM2,
    TESLA_METER_TO_GAUSS_CM,
    TESLA_TO_GAUSS,
    write_boozer_chartmap,
)

TWOPI = 2.0 * np.pi

_INP_TEMPLATE = """\
{nstep}      nstep    - number of integration steps
{nlabel}       nlabel   - grid size over radial variable
{ntheta_int}       ntheta   - grid size over poloidal angle
{nsurfmax}     nsurfmax - number of starting points for separatrix search
{nsurf}        nsurf    - number of flux surfaces in Boozer file
{mpol}         mpol     - number of poloidal modes in Boozer file
{psimax}
"""


def _write_inp(path, gfile, *, nstep=3600, nlabel=500, ntheta_int=500,
               nsurfmax=10000, nsurf=2000, mpol=12, psimax=1e10):
    """Write efit_to_boozer.inp pointing at gfile.

    psimax: poloidal flux boundary value in CGS (G*cm^2); controls when the
    separatrix scan terminates.  For a limiter discharge pass
    abs(psi_edge_Wb) * 1e8 * 1.01 (10% above edge flux) so the scan stops
    at the LCFS without needing to reach an X-point.
    """
    path = Path(path)
    path.write_text(
        _INP_TEMPLATE.format(
            nstep=nstep, nlabel=nlabel, ntheta_int=ntheta_int,
            nsurfmax=nsurfmax, nsurf=nsurf, mpol=mpol, psimax=psimax,
        )
    )


def _write_field_divB0_inp(path, gfile, *, convexfile="'unused'"):
    """Write minimal field_divB0.inp so the field reader locates the g-file.

    Matches the 13-line format read by read_field_input in field_divB0.f90.
    The Fortran reader uses err= not end= for the last line, so all 13 lines
    must be present.
    """
    if not convexfile.startswith("'"):
        convexfile = f"'{convexfile}'"
    path = Path(path)
    path.write_text(
        f"0\n"              # ipert: 0=eq only
        f"1\n"              # iequil: 1=with equilibrium
        f"1.00\n"           # ampl
        f"72\n"             # ntor
        f"0.99\n"           # cutoff
        f"4\n"              # icftype
        f"'{gfile}'\n"      # gfile
        f"'unused'\n"       # pfile
        f"{convexfile}\n"   # convexfile
        f"'unused'\n"       # fluxdatapath
        f"0\n"              # nwindow_r
        f"0\n"              # nwindow_z
        f"1\n"              # ieqfile: 1=EFIT
    )


def _convex_wall_from_lcfs(lcfs_R, lcfs_Z, margin=1.1):
    """Return (R_cm, Z_cm) wall points enclosing the LCFS with a margin."""
    import math
    R0 = float(np.mean(lcfs_R))
    Z0 = float(np.mean(lcfs_Z))
    r = np.sqrt((lcfs_R - R0) ** 2 + (lcfs_Z - Z0) ** 2)
    a = float(np.max(r)) * margin
    n = max(len(lcfs_R), 100)
    phi = np.linspace(0.0, 2.0 * math.pi, n, endpoint=False)
    return (R0 + a * np.cos(phi)) * 100.0, (Z0 + a * np.sin(phi)) * 100.0


def _write_convex_wall_from_lcfs(path, lcfs_R, lcfs_Z, margin=1.1):
    """Write a convex wall file (cm) enclosing the LCFS."""
    Rs, Zs = _convex_wall_from_lcfs(lcfs_R, lcfs_Z, margin)
    path = Path(path)
    with path.open("w") as f:
        for R_cm, Z_cm in zip(Rs, Zs):
            f.write(f"  {R_cm:.16e}  {Z_cm:.16e}\n")


def convert_eqdsk_to_chartmap(
    gfile,
    output,
    *,
    nrho=50,
    ntheta=48,
    nstep=3600,
    nlabel=500,
    ntheta_int=500,
    nsurfmax=10000,
    nsurf=2000,
    convexwall=None,
    psimax=None,
):
    """Read an EQDSK g-file and write a libneo Boozer chartmap.

    The efit_to_boozer f2py module is loaded lazily and run inside a temporary
    directory so its inp/output files do not pollute the caller's cwd.

    Parameters
    ----------
    nsurfmax : int
        Number of starting points for the separatrix search in efit_to_boozer.
        Smaller values run faster; 200 is sufficient for most EFIT files.
    nsurf : int
        Number of flux surfaces in the Boozer file.  Must be >= nlabel.
    convexwall : path-like or None
        Path to a convex wall file for stretch_coords.  If None a default
        circular wall is generated from the g-file boundary.
    psimax : float or None
        Poloidal flux at the plasma boundary in CGS (G*cm^2/rad), absolute
        value as interpolated by field_eq from the g-file psi map.  The
        flux-surface scan stops at the first surface with psi beyond this
        value.  If None, the boundary is located purely by the field-line
        integration leaving the computational box; that fails for synthetic
        equilibria without an X-point, whose psi keeps rising to the box
        edge (the scan then runs to nsurfmax and the last box-confined
        surface, not the LCFS, becomes s = 1).
    """
    import _efit_to_boozer as _etb
    from _efit_to_boozer import efit_to_boozer_mod
    from efit_to_boozer.boozer import get_boozer_transform, get_magnetic_axis
    from scipy.interpolate import CubicSpline

    gfile = str(Path(gfile).resolve())
    if convexwall is not None:
        convexwall = str(Path(convexwall).resolve())

    orig_dir = os.getcwd()
    with tempfile.TemporaryDirectory() as tmpdir:
        os.chdir(tmpdir)
        try:
            # Without an explicit psimax the boundary is located purely by
            # the field-line integration (efit_to_boozer's historical 1e10
            # placeholder disables the flux-based stop).
            inp_kwargs = dict(
                nstep=nstep,
                nlabel=nlabel,
                ntheta_int=ntheta_int,
                nsurfmax=nsurfmax,
                nsurf=nsurf,
            )
            if psimax is not None:
                inp_kwargs["psimax"] = psimax
            _write_inp("efit_to_boozer.inp", gfile, **inp_kwargs)
            if convexwall is None:
                from libneo.eqdsk_base import read_eqdsk
                eq = read_eqdsk(gfile)
                lcfs = eq["Lcfs"]
                _write_convex_wall_from_lcfs(
                    "convexwall.dat", lcfs[:, 0], lcfs[:, 1]
                )
                wall_path = "convexwall.dat"
            else:
                wall_path = convexwall
            _write_field_divB0_inp("field_divB0.inp", gfile,
                                   convexfile=wall_path)

            # Uniform s grid: skip s=0 (axis singularity) and s=1 (LCFS,
            # field-line integration may not converge exactly at the separatrix).
            s_grid = np.linspace(1.0 / nrho, 1.0 - 0.5 / nrho, nrho)
            rho_grid = np.sqrt(s_grid)

            # Boozer transform per surface: splines r(theta_B), dth(theta_B),
            # G(theta_B), q; stored on stor[0..nrho-1] (stor[-1] is not used).
            # get_boozer_transform calls efit_to_boozer.init() internally from
            # the current directory (tmpdir), so do not call init() separately.
            num_theta_boozer = 64
            stor = np.append(s_grid, 1.0 + 1.0 / nrho)  # sentinel
            (r_of_thb, dth_of_thb, G_of_thb,
             _theta_b, _theta_g, _theta_sf, B0_of_thb) = get_boozer_transform(
                stor, num_theta_boozer
            )

            # psitor_max available after init (called inside get_boozer_transform).
            # SIMPLE/boozer_chartmap convention: torflux < 0 (left-handed).
            psitor_max_cgs = float(efit_to_boozer_mod.psitor_max)
            torflux_si = -psitor_max_cgs * 1.0e-8  # G*cm^2 -> T*m^2, negated

            R_axis_si, Z_axis_si = get_magnetic_axis()  # in SI (m)

            # Per-surface surface functions C_norm (= B_phi Boozer, CGS G*cm)
            # and q from magdata at theta=0.
            psi_dum = np.float64(0.0)
            C_norm_arr = np.zeros(nrho)
            q_arr = np.zeros(nrho)
            for k, sk in enumerate(s_grid):
                s_f = np.float64(sk)
                psi_f = np.float64(0.0)
                th_f = np.float64(0.0)
                (q, _dq, Cn, *_rest) = _etb.efit_to_boozer.magdata(
                    1, s_f, psi_f, th_f
                )
                C_norm_arr[k] = Cn
                q_arr[k] = q

            # Covariant B surface functions (CGS G*cm).
            # B_phi = G(s) = R*B_tor = C_norm from magdata.
            # B_theta = I(s), the Boozer covariant poloidal component, from
            # Ampere's law: 2*pi*I = loop integral of B_pol along the flux
            # surface (only the enclosed toroidal current contributes).
            # NOTE: I is NOT G/q. The field-line pitch iota is the ratio of
            # the CONTRAvariant components; the covariant ratio I/G measures
            # enclosed current. For a tokamak I << G/q, and writing G/q here
            # inflates the field-line length factor (G + iota*I) and slows
            # every traced orbit's parallel dynamics by O(1/q^2).
            B_phi_cgs = C_norm_arr          # G*cm
            # B_pol = |grad psi| / R from the g-file psi map itself: taking
            # sqrt(B^2 - B_tor^2) instead would difference two nearly equal
            # interpolants and lose all accuracy on inner surfaces.
            from libneo.eqdsk_base import read_eqdsk
            from scipy.interpolate import RectBivariateSpline
            eq_amp = read_eqdsk(gfile)
            psi_spl = RectBivariateSpline(
                np.asarray(eq_amp["R"]) * METER_TO_CM,
                np.asarray(eq_amp["Z"]) * METER_TO_CM,
                np.asarray(eq_amp["PsiVs"]).T * 1.0e8,  # Wb/rad -> G*cm^2/rad
            )
            B_theta_cgs = np.zeros(nrho)    # G*cm
            ntheta_amp = 2048
            thb_amp = np.linspace(0.0, TWOPI, ntheta_amp, endpoint=False)
            for k in range(nrho):
                r_amp = np.asarray(r_of_thb[k](thb_amp), dtype=float)
                dth_amp = np.asarray(dth_of_thb[k](thb_amp), dtype=float)
                th_geom = thb_amp + dth_amp
                R_amp = (R_axis_si + r_amp * np.cos(th_geom)) * METER_TO_CM
                Z_amp = (Z_axis_si + r_amp * np.sin(th_geom)) * METER_TO_CM
                dR = (np.roll(R_amp, -1) - np.roll(R_amp, 1)) / 2.0
                dZ = (np.roll(Z_amp, -1) - np.roll(Z_amp, 1)) / 2.0
                dl = np.hypot(dR, dZ)  # cm
                dpsi_dR = psi_spl(R_amp, Z_amp, dx=1, grid=False)
                dpsi_dZ = psi_spl(R_amp, Z_amp, dy=1, grid=False)
                # Signed circulation of B_pol = (-psi_Z, psi_R)/R along
                # increasing theta, so the sign of I follows the enclosed
                # current direction instead of being forced positive.
                B_R = -dpsi_dZ / R_amp  # G
                B_Z = dpsi_dR / R_amp   # G
                B_theta_cgs[k] = np.sum(B_R * dR + B_Z * dZ) / TWOPI  # G*cm

            # A_phi(s) = -torflux * integral_0^s iota(s') ds'
            # iota = 1/q.  Use a cubic spline on s_grid.
            iota_spl = CubicSpline(s_grid, 1.0 / q_arr)
            iota_int = iota_spl.antiderivative()
            A_phi_si = -torflux_si * (iota_int(s_grid) - iota_int(0.0))

            # Geometry and Bmod on (nrho, ntheta, 1) grid (nzeta=1 for axisymm).
            theta_B = np.linspace(0.0, TWOPI, ntheta, endpoint=False)
            zeta_B = np.array([0.0])  # one zeta slice; geometry repeats

            Bmod = np.zeros((nrho, ntheta, 1))
            X = np.zeros((nrho, ntheta, 1))
            Y = np.zeros((nrho, ntheta, 1))
            Z = np.zeros((nrho, ntheta, 1))

            for k in range(nrho):
                r_spl = r_of_thb[k]
                dth_spl = dth_of_thb[k]
                G_spl = G_of_thb[k]
                B0_spl = B0_of_thb[k]

                for j, thb in enumerate(theta_B):
                    r_min = float(r_spl(thb))      # minor radius (m)
                    dth = float(dth_spl(thb))      # theta_geom - theta_B
                    G_val = float(G_spl(thb))       # local phi-shift
                    B0_val = float(B0_spl(thb))     # |B| in CGS (G)

                    theta_geom = thb + dth
                    # phi_cyl = zeta_B - G_val (from boozer.py convention)
                    phi_cyl = 0.0 - G_val

                    R_cyl = R_axis_si + r_min * np.cos(theta_geom)  # m
                    Z_cyl = Z_axis_si + r_min * np.sin(theta_geom)  # m

                    Bmod[k, j, 0] = B0_val  # already in CGS (G)
                    X[k, j, 0] = R_cyl * np.cos(phi_cyl)  # m
                    Y[k, j, 0] = R_cyl * np.sin(phi_cyl)  # m
                    Z[k, j, 0] = Z_cyl                     # m

            _etb.efit_to_boozer.deinit()
        finally:
            os.chdir(orig_dir)

    # Convert to CGS.
    X_cgs = X * METER_TO_CM
    Y_cgs = Y * METER_TO_CM
    Z_cgs = Z * METER_TO_CM
    A_phi_cgs = A_phi_si * TESLA_METER2_TO_GAUSS_CM2
    torflux_cgs = torflux_si * TESLA_METER2_TO_GAUSS_CM2
    # B_phi and B_theta are already in CGS G*cm from magdata.

    write_boozer_chartmap(
        output,
        rho=rho_grid,
        s=s_grid,
        theta=theta_B,
        zeta=zeta_B,
        x=X_cgs,
        y=Y_cgs,
        z=Z_cgs,
        A_phi=A_phi_cgs,
        B_theta=B_theta_cgs,
        B_phi=B_phi_cgs,
        Bmod=Bmod,
        num_field_periods=1,
        torflux=torflux_cgs,
        eqdsk2chartmap_source=str(gfile),
    )
    return torflux_cgs


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Convert an EQDSK g-file to a libneo Boozer chartmap NetCDF file"
    )
    parser.add_argument("gfile", help="EQDSK equilibrium g-file")
    parser.add_argument("output", help="Output Boozer chartmap .nc file")
    parser.add_argument("--nrho", type=int, default=50)
    parser.add_argument(
        "--ntheta",
        type=int,
        default=48,
        help="poloidal Boozer angle points (endpoint-excluded)",
    )
    parser.add_argument("--nstep", type=int, default=3600)
    parser.add_argument("--nlabel", type=int, default=500)
    parser.add_argument("--ntheta-int", type=int, default=500,
                        dest="ntheta_int",
                        help="poloidal grid in efit_to_boozer integration")
    parser.add_argument("--nsurfmax", type=int, default=10000,
                        help="starting points for separatrix search")
    parser.add_argument("--nsurf", type=int, default=2000,
                        help="flux surfaces in Boozer file")
    args = parser.parse_args(argv)

    torflux = convert_eqdsk_to_chartmap(
        args.gfile,
        args.output,
        nrho=args.nrho,
        ntheta=args.ntheta,
        nstep=args.nstep,
        nlabel=args.nlabel,
        ntheta_int=args.ntheta_int,
        nsurfmax=args.nsurfmax,
        nsurf=args.nsurf,
    )
    print(f"Done. torflux={torflux:.6e} G cm^2")


if __name__ == "__main__":
    main()
