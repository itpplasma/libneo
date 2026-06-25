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
        f"0\n"              # ntor
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
    nsurfmax=200,
    nsurf=500,
    convexwall=None,
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
            from libneo.eqdsk_base import read_eqdsk
            eq = read_eqdsk(gfile)
            # psimax in CGS (G*cm^2): stop the surface scan ~1% INSIDE the LCFS.
            # Scanning to or past psi_edge can place the outermost surface in the
            # flat-psi SOL of a limiter equilibrium, where B_pol -> 0 and the
            # field-line transit time -> infinity, stalling the ODE integrator
            # (Maximum steps exceeded). Staying inside keeps every scanned surface
            # on a well-defined flux surface. The magnitude matches psif inside
            # field_divB0, which normalises psi to zero at the axis.
            psi_edge_cgs = abs(eq["PsiedgeVs"] - eq["PsiaxisVs"]) * 0.99e8
            _write_inp(
                "efit_to_boozer.inp",
                gfile,
                nstep=nstep,
                nlabel=nlabel,
                ntheta_int=ntheta_int,
                nsurfmax=nsurfmax,
                nsurf=nsurf,
                psimax=psi_edge_cgs,
            )
            if convexwall is None:
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

            # Covariant B surface functions in SI.
            # B_phi (G) = C_norm [G*cm] -> SI: C_norm * 1e-6 T*m
            # B_theta (I) = C_norm / q [G*cm] -> SI
            B_phi_cgs = C_norm_arr          # G*cm
            B_theta_cgs = C_norm_arr / q_arr  # G*cm; iota=I/G => I=G/q

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
    parser.add_argument("--nsurfmax", type=int, default=200,
                        help="starting points for separatrix search")
    parser.add_argument("--nsurf", type=int, default=500,
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
