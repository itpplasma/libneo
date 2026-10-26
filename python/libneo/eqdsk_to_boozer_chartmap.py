"""Convert an EQDSK/g-file axisymmetric equilibrium to a libneo Boozer chartmap.

Uses the native EFIT contour construction with hard equilibrium-consistency
checks. The optional override requires a reason, which is logged by Fortran.

Output convention
-----------------
Axisymmetric: num_field_periods=1. Cartesian geometry covers the full torus,
even though the field strength and surface functions are axisymmetric.
Geometry varies with
phi_cyl = zeta - G_of_thb(theta_B); Bmod is independent of zeta.
Surface functions:
  B_phi (G*cm) = native contour covariant toroidal field
  B_theta (G*cm) = signed poloidal circulation / (2*pi)
  A_phi(s) = -torflux * integral_0^s iota(s') ds'
torflux = -sign(F)*abs(psitor_max) (left-handed SIMPLE convention)
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
{construction}         construction - 1 contours, 2 legacy field lines
"""


def _write_inp(path, gfile, *, nstep=3600, nlabel=500, ntheta_int=500,
               nsurfmax=10000, nsurf=2000, mpol=12, psimax=1e10, construction=1):
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
            construction=construction,
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
    nzeta=48,
    nstep=3600,
    nlabel=500,
    ntheta_int=500,
    nsurfmax=10000,
    nsurf=2000,
    convexwall=None,
    psimax=None,
    construction=1,
    n_q=128,
    allow_poor_equilibrium=False,
    override_reason="",
    exp_bpol=-1,
):
    """Read an EQDSK g-file and write a libneo Boozer chartmap.

    The contour construction uses the file LCFS by default. ``psimax`` is an
    optional axis-relative poloidal flux boundary in CGS. ``ntheta_int`` and
    ``n_q`` set the contour and toroidal-flux integration grids. Historical
    field-line scan arguments are retained for call compatibility.

    ``construction=2`` explicitly selects historical field-line integration.
    ``exp_bpol`` selects flux units: -1 to infer, 0 for Wb/rad, 1 for Wb.
    Poor equilibria are rejected unless ``allow_poor_equilibrium=True`` and
    a nonempty ``override_reason`` is supplied. The override is logged.
    """
    if construction not in (1, 2):
        raise ValueError("Chartmap construction must be 1 (contour) or 2 (legacy)")
    if not isinstance(nzeta, int) or nzeta < 8:
        raise ValueError("nzeta must be an integer >=8 for periodic Cartesian splines")
    from efit_to_boozer.boozer import get_boozer_transform
    from scipy.interpolate import CubicSpline

    gfile = str(Path(gfile).resolve())
    from libneo.eqdsk_base import read_eqdsk
    eq = read_eqdsk(gfile)
    # A flux interval and a toroidal field are both needed for flux labels.
    if eq["fprof"][0] * (eq["PsiedgeVs"]-eq["PsiaxisVs"]) == 0:
        raise ValueError("Distinct boundary flux and nonzero toroidal field required")
    if convexwall is not None:
        convexwall = str(Path(convexwall).resolve())

    orig_dir = os.getcwd()
    with tempfile.TemporaryDirectory() as tmpdir:
        os.chdir(tmpdir)
        try:
            # Retain the historical deck files for compatibility. Contours
            # use the file LCFS unless psimax explicitly selects an inner edge.
            inp_kwargs = dict(
                nstep=nstep,
                nlabel=nlabel,
                ntheta_int=ntheta_int,
                nsurfmax=nsurfmax,
                nsurf=nsurf, construction=construction,
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

            # Native chartmap readers require a uniform rho grid, reaching
            # the normalized outer surface and beginning close to the axis.
            # The outer surface is the file LCFS or the explicit psimax contour.
            rho_grid = np.linspace(0.01, 1.0, nrho)
            surface_s = rho_grid**2
            # A_phi has its own independent uniform s abscissa. Its entries
            # are not paired by index with the geometry's uniform rho grid.
            s_grid = np.linspace(surface_s[0], 1.0, nrho)

            # The geometry and I/G use the same native contour quadrature.
            num_theta_boozer = 64
            stor = np.append(surface_s, 1.0 + 1.0 / nrho)  # sentinel
            if construction == 2:
                transform, metadata = _legacy_transform(stor, num_theta_boozer, eq)
            else:
                transform, metadata = get_boozer_transform(
                    stor, num_theta_boozer, gfile=gfile,
                    n_theta=ntheta_int + ntheta_int % 2, n_q=n_q, psimax=psimax,
                    allow_poor_equilibrium=allow_poor_equilibrium,
                    override_reason=override_reason, exp_bpol=exp_bpol, return_metadata=True,
                )
            (r_of_thb, dth_of_thb, G_of_thb,
             _theta_b, _theta_g, _theta_sf, B0_of_thb) = transform
            torflux_si = -np.copysign(abs(metadata["psitor_max"]),
                                      eq["fprof"][0]) * 1.0e-8
            R_axis_si, Z_axis_si = metadata["axis"]
            q_arr = 1.0 / metadata["iota"]
            B_theta_cgs = metadata["curr_i"]
            B_phi_cgs = metadata["curr_g"]

            # A_phi(s) = -torflux * integral_0^s iota(s') ds'
            # iota = 1/q.  Use a cubic spline on s_grid.
            iota_spl = CubicSpline(surface_s, 1.0 / q_arr)
            iota_int = iota_spl.antiderivative()
            A_phi_si = -torflux_si * (iota_int(s_grid) - iota_int(0.0))

            # Cartesian geometry rotates with zeta even in an axisymmetric field.
            theta_B = np.linspace(0.0, TWOPI, ntheta, endpoint=False)
            zeta_B = np.linspace(0.0, TWOPI, nzeta, endpoint=False)

            Bmod = np.zeros((nrho, ntheta, nzeta))
            X = np.zeros((nrho, ntheta, nzeta))
            Y = np.zeros((nrho, ntheta, nzeta))
            Z = np.zeros((nrho, ntheta, nzeta))

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
                    phi_cyl = zeta_B - G_val

                    R_cyl = R_axis_si + r_min * np.cos(theta_geom)  # m
                    Z_cyl = Z_axis_si + r_min * np.sin(theta_geom)  # m

                    Bmod[k, j, :] = B0_val  # already in CGS (G)
                    X[k, j, :] = R_cyl * np.cos(phi_cyl)  # m
                    Y[k, j, :] = R_cyl * np.sin(phi_cyl)  # m
                    Z[k, j, :] = Z_cyl                     # m

        finally:
            os.chdir(orig_dir)

    # Convert to CGS.
    X_cgs = X * METER_TO_CM
    Y_cgs = Y * METER_TO_CM
    Z_cgs = Z * METER_TO_CM
    A_phi_cgs = A_phi_si * TESLA_METER2_TO_GAUSS_CM2
    torflux_cgs = torflux_si * TESLA_METER2_TO_GAUSS_CM2
    # B_phi and B_theta are already in CGS G*cm from the native contours.

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


def _legacy_transform(stor, num_theta, eq):
    """Retain historical field-line geometry and circulation for construction 2."""
    import _efit_to_boozer as native
    from efit_to_boozer.boozer import get_boozer_transform, get_magnetic_axis
    from scipy.interpolate import RectBivariateSpline

    transform = get_boozer_transform(stor, num_theta, 2)
    radius, dtheta = transform[:2]
    axis = get_magnetic_axis(method="legacy")
    surfaces = stor[:-1]
    values = np.array([native.efit_to_boozer.magdata(
        1, np.float64(s), np.float64(0.0), np.float64(0.0))[:3] for s in surfaces])
    q = np.sign(eq["fprof"][0]*(eq["PsiedgeVs"]-eq["PsiaxisVs"]))*np.abs(values[:, 0])
    psi = RectBivariateSpline(np.asarray(eq["R"])*METER_TO_CM,
                              np.asarray(eq["Z"])*METER_TO_CM,
                              np.asarray(eq["PsiVs"]).T*1e8)
    theta = np.linspace(0, TWOPI, 2048, endpoint=False)
    current_i = np.zeros(len(surfaces))
    for k in range(len(surfaces)):
        angle = theta+dtheta[k](theta)
        r = (axis[0]+radius[k](theta)*np.cos(angle))*METER_TO_CM
        z = (axis[1]+radius[k](theta)*np.sin(angle))*METER_TO_CM
        dr = (np.roll(r, -1)-np.roll(r, 1))/2
        dz = (np.roll(z, -1)-np.roll(z, 1))/2
        current_i[k] = np.sum(-psi(r, z, dy=1, grid=False)/r*dr
                              + psi(r, z, dx=1, grid=False)/r*dz)/TWOPI
    metadata = dict(axis=axis, iota=1/q, curr_i=current_i,
                    curr_g=values[:, 2]-current_i/q,
                    psitor_max=float(native.efit_to_boozer_mod.psitor_max))
    native.efit_to_boozer.deinit()
    return transform, metadata


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
    parser.add_argument("--construction", type=int, choices=[1, 2], default=1)
    parser.add_argument("--exp-bpol", type=int, choices=[-1, 0, 1], default=-1,
                        help="poloidal-flux units: -1 infer, 0 Wb/rad, 1 Wb")
    parser.add_argument("--allow-poor-equilibrium", action="store_true")
    parser.add_argument("--override-reason", default="")
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
        construction=args.construction, exp_bpol=args.exp_bpol,
        allow_poor_equilibrium=args.allow_poor_equilibrium,
        override_reason=args.override_reason,
    )
    print(f"Done. torflux={torflux:.6e} G cm^2")


if __name__ == "__main__":
    main()
