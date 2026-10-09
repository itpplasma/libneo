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
def test_prescribed_circular_flux_boundary(tmp_path, scan):
    binary = os.environ.get("EFIT_TO_BOOZER_BINARY")
    if not binary:
        pytest.skip("set EFIT_TO_BOOZER_BINARY to the native converter")
    # Exactly representable quadratic psi; q=F/sqrt(R0^2-r^2).
    r0, f, radius, nr = 4., 3., .8, 65
    r, z = np.linspace(2.5, 5.5, nr), np.linspace(-1.5, 1.5, nr)
    R, Z = np.meshgrid(r, z)
    psi = ((R-r0)**2+Z**2-radius**2)/2
    header = [3, 3, r0, 2.5, 0, r0, 0, -radius**2/2, 0, f/r0,
              1, -radius**2/2, 0, r0, 0, 0, 0, 0, 0, 0]
    q = f/np.sqrt(r0*r0-radius**2*np.linspace(0, 1, nr))
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
    _write_inp(tmp_path/"efit_to_boozer.inp", str(gfile), nlabel=128,
               ntheta_int=256, nsurfmax=scan, nsurf=40, mpol=12,
               psimax=radius**2/2*1e8)
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
    q_exact = f/(r0-np.asarray(bc.s)*flux/(2*np.pi*f))
    np.testing.assert_allclose(1/np.asarray(bc.iota), q_exact, rtol=2e-7)
