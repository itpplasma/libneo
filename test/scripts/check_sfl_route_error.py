#!/usr/bin/env python3
"""
boozer_spectrum_ad.x with route = 4 (sfl_straight_current, Hamada II) on a
wout without pressure must stop with a nonzero exit code and name the cause
(p' = 0) and the alternative routes, while route = 3 (Hamada I) runs and
writes the route diagnostics.
"""

from __future__ import annotations

import argparse
import subprocess
import tempfile
from pathlib import Path

from netCDF4 import Dataset


def write_input(path: Path, route: int) -> None:
    with Dataset(path, "w") as ds:
        ds.createDimension("s", 2)
        ds.createVariable("s", "f8", ("s",))[:] = [0.25, 0.5]
        for key, value in (("mboz", 6), ("nboz", 4), ("grid_refinement", 0),
                           ("route", route), ("diagnostics", 1)):
            ds.createVariable(key, "i4", ())[...] = value


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("exe")
    parser.add_argument("wout")
    args = parser.parse_args()
    with tempfile.TemporaryDirectory() as tmp:
        fin, fout = Path(tmp) / "in.nc", Path(tmp) / "out.nc"
        write_input(fin, 4)
        run = subprocess.run([args.exe, "primal", args.wout, str(fin), str(fout)],
                             capture_output=True, text=True)
        text = run.stdout + run.stderr
        if run.returncode == 0:
            print("FAIL: route 4 without pressure did not stop")
            return 1
        for needed in ("sfl_straight_current", "p' /= 0", "sfl_jacobian_const"):
            if needed not in text:
                print(f"FAIL: message lacks '{needed}':\n{text}")
                return 1
        print("ok   route 4 without pressure stops:", text.splitlines()[0][:120])
        write_input(fin, 3)
        run = subprocess.run([args.exe, "primal", args.wout, str(fin), str(fout)],
                             capture_output=True, text=True)
        if run.returncode != 0:
            print("FAIL: route 3 did not run:\n" + run.stdout + run.stderr)
            return 1
        with Dataset(fout) as ds:
            if ds.getncattr("route_name") != "sfl_jacobian_const":
                print("FAIL: route_name attribute")
                return 1
            d_jac = ds["diag_d_jac"][:]
            if not (d_jac < 1e-3).all():
                print(f"FAIL: Hamada I Jacobian variation {d_jac}")
                return 1
        print("ok   route 3 runs, d_jac =", list(map(float, d_jac)))
    print("check_sfl_route_error PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
