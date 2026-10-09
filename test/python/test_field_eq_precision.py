"""The f2py field reader must preserve the double-precision Fortran ABI."""
import importlib.util
import os
from pathlib import Path

import numpy as np
import pytest


def test_field_eq_quadratic_roundtrip(tmp_path, monkeypatch):
    # An explicit build path avoids an editable install's import hook selecting
    # an old extension. Run this test in its own process (Fortran global state).
    build = os.environ.get("LIBNEO_EFIT_EXTENSION")
    if build:
        spec = importlib.util.spec_from_file_location("_efit_to_boozer", build)
        ext = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(ext)
    else:
        ext = pytest.importorskip("_efit_to_boozer")
    monkeypatch.chdir(tmp_path)
    r, z = np.linspace(3, 5, 9), np.linspace(-1, 1, 9)
    R, Z = np.meshgrid(r, z)
    psi = ((R-4)**2+Z**2-1)/2
    h = [2., 2., 4., 3., 0., 4., 0., -.5, 0., .75,
         1., -.5, 0., 4., 0., 0., 0., 0., 0., 0.]
    with Path("field.g").open("w") as stream:
        stream.write(f"{'quadratic':48s}{0:4d}{9:4d}{9:4d}\n")
        for values in (h, [3.]*9, [0.]*9, [0.]*9, [0.]*9, psi.ravel(), [2.]*9):
            for i in range(0, len(values), 5):
                stream.write("".join(f"{v:16.9E}" for v in values[i:i+5])+"\n")
        stream.write("    0    0\n\n\n")
    Path("field_divB0.inp").write_text(
        "0\n1\n1.0\n72\n0.99\n4\n'field.g'\n'unused'\n'unused'\n'unused'\n0\n0\n1\n")
    ext.field_sub.read_field_input("field_divB0.inp")
    ext.field_eq_mod.use_fpol = True
    # Off-grid points exercise the full reader→spline→field wrapper path.
    for R, Z in ((4.125, .375), (3.625, -.125), (4.5, -.5)):
        br, bf, bz, *_ = ext.field_sub.field_eq(100*R, 0., 100*Z)
        np.testing.assert_allclose([br*1e-4, bf*1e-4, bz*1e-4],
                                   [-Z/R, 3/R, (R-4)/R], rtol=0, atol=2e-13)
        np.testing.assert_allclose(ext.field_sub.psif*1e-8,
                                   ((R-4)**2+Z**2)/2, rtol=0, atol=2e-13)
