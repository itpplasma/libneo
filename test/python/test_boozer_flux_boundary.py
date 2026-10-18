"""A finite psimax is a flux boundary, independently of scan spacing."""
import os
from pathlib import Path
import subprocess

import numpy as np
import pytest

from libneo.boozer import BoozerFile
from libneo.eqdsk_to_boozer_chartmap import (
    _write_inp, _write_field_divB0_inp, _write_convex_wall_from_lcfs,
)


@pytest.mark.parametrize("scan", [80, 160])
@pytest.mark.parametrize("polarity", [-1, 1])
def test_prescribed_circular_flux_boundary(tmp_path, scan, polarity):
    binary = os.environ.get("EFIT_TO_BOOZER_BINARY")
    if not binary:
        pytest.skip("set EFIT_TO_BOOZER_BINARY to the native converter")
    # Exactly representable quadratic psi; q=F/sqrt(R0^2-r^2).
    r0, f, radius, nr = 4., 3., .8, 65
    r, z = np.linspace(2.5, 5.5, nr), np.linspace(-1.5, 1.5, nr)
    R, Z = np.meshgrid(r, z)
    psi = polarity*((R-r0)**2+Z**2-radius**2)/2
    header = [3, 3, r0, 2.5, 0, r0, 0, -polarity*radius**2/2, 0, f/r0,
              polarity, -polarity*radius**2/2, 0, r0, 0, 0, 0, 0, 0, 0]
    q = polarity*f/np.sqrt(r0*r0-radius**2*np.linspace(0, 1, nr))
    gfile = tmp_path/"circular.g"
    with gfile.open("w") as out:
        out.write(f"{'exact circular; COCOS 3':48s}{0:4d}{nr:4d}{nr:4d}\n")
        for a in [header, np.full(nr, f), np.zeros(nr), np.zeros(nr),
                  np.zeros(nr), psi.ravel(), q]:
            for i in range(0, len(a), 5):
                out.write("".join(f"{v:16.9E}" for v in a[i:i+5])+"\n")
        theta = np.linspace(0, 2*np.pi, 65)
        boundary = np.column_stack((r0+radius*np.cos(theta), radius*np.sin(theta)))
        out.write("   65   65\n")
        for a in [boundary.ravel(), boundary.ravel()]:
            for i in range(0, len(a), 5):
                out.write("".join(f"{v:16.9E}" for v in a[i:i+5])+"\n")
    # A wall encloses the requested smooth surface, not a separatrix.
    _write_convex_wall_from_lcfs(tmp_path/"convexwall.dat", boundary[:, 0], boundary[:, 1])
    _write_inp(tmp_path/"efit_to_boozer.inp", str(gfile), nlabel=512,
               ntheta_int=256, nsurfmax=scan, nsurf=2000, mpol=12,
               psimax=polarity*radius**2/2*1e8)
    _write_field_divB0_inp(tmp_path/"field_divB0.inp", str(gfile),
                          convexfile="convexwall.dat")
    proc = subprocess.run([str(Path(binary).resolve())], cwd=tmp_path,
                          capture_output=True, text=True, timeout=60,
                          env={**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"})
    (tmp_path/"converter.log").write_text(proc.stdout+proc.stderr)
    assert proc.returncode == 0, proc.stdout+proc.stderr
    bc = BoozerFile(str(tmp_path/"fromefit_neo_lhs.bc"))
    flux = 2*np.pi*f*(r0-np.sqrt(r0*r0-radius**2))
    # The existing six-digit header has a separate serialization floor.
    assert bc.flux == pytest.approx(flux, rel=5e-6)
    # Exact q at fixed normalized toroidal flux, not at a fitted radius.
    q_exact = polarity*f/(r0-np.asarray(bc.s)*flux/(2*np.pi*f))
    np.testing.assert_allclose(1/np.asarray(bc.iota), q_exact, rtol=2e-7)
    native_flux = np.loadtxt(tmp_path/"flux_functions.dat")[-1, 5]*2*np.pi*1e-8
    # The intermediate contour integral follows the poloidal field; the
    # delivered left-handed file restores physical toroidal-flux orientation.
    assert bc.flux == pytest.approx(polarity*native_flux, rel=2e-14)

    # Independent circular oracle for the Boozer metric: J=-Phi*R^2/(2*pi*F).
    # Differentiating densely sampled geometry exposes coefficient-rounding
    # errors that scalar B/q/flux readback cannot see.
    from scipy.interpolate import CubicSpline
    angle = np.arange(256)*2*np.pi/256
    modes = np.asarray(bc.m[0])[:, None]
    co, si = np.cos(modes*angle), np.sin(modes*angle)
    samples = np.array([.05, .4, .9])
    geometry = {}
    for name in ("rmnc", "rmns", "zmnc", "zmns"):
        spline = CubicSpline(bc.s, np.asarray(getattr(bc, name)), axis=0)
        geometry[name] = (spline(samples), spline(samples, 1))
    rc, rs, zc, zs = (geometry[name] for name in
                      ("rmnc", "rmns", "zmnc", "zmns"))
    R = rc[0]@co+rs[0]@si
    R_s = rc[1]@co+rs[1]@si
    Z_s = zc[1]@co+zs[1]@si
    R_t = -(rc[0]*modes.T)@si+(rs[0]*modes.T)@co
    Z_t = -(zc[0]*modes.T)@si+(zs[0]*modes.T)@co
    jacobian = R*(R_t*Z_s-R_s*Z_t)
    expected = -flux*R**2/(2*np.pi*f)
    np.testing.assert_allclose(jacobian, expected, rtol=1e-6)
