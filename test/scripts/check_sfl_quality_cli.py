#!/usr/bin/env python3
"""boozer_spectrum_ad.x on a coarse wout (li383, ns 49): without allow_poor_equilibrium the
equilibrium-quality check stops with an actionable message; with it the run completes and the
override, its reason and the quality record are written to out.nc."""
from __future__ import annotations

import argparse
import subprocess
import tempfile
from pathlib import Path

from netCDF4 import Dataset

REASON = "check_sfl_quality_cli: coarse public test wout"


def write_input(path: Path, allow: int) -> None:
    with Dataset(path, "w") as ds:
        ds.createDimension("s", 2)
        ds.createVariable("s", "f8", ("s",))[:] = [0.25, 0.5]
        for key, value in (("mboz", 6), ("nboz", 4), ("grid_refinement", 0),
                           ("allow_poor_equilibrium", allow)):
            ds.createVariable(key, "i4", ())[...] = value
        ds.override_reason = REASON


def fail(msg: str) -> int:
    print("FAIL:", msg)
    return 1


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("exe")
    parser.add_argument("wout")
    args = parser.parse_args()
    with tempfile.TemporaryDirectory() as tmp:
        fin, fout = Path(tmp) / "in.nc", Path(tmp) / "out.nc"
        write_input(fin, 0)
        run = subprocess.run([args.exe, "primal", args.wout, str(fin), str(fout)],
                             capture_output=True, text=True)
        text = run.stdout + run.stderr
        if run.returncode == 0:
            return fail("coarse wout passed the quality check without override")
        for needed in ("inconsistent equilibrium", "eps_rad", "Raise ns",
                       "allow_poor_equilibrium", "nearest resonance"):
            if needed not in text:
                return fail(f"message lacks '{needed}':\n{text}")
        print("ok   stops without override:", text.splitlines()[0][:110])
        write_input(fin, 1)
        run = subprocess.run([args.exe, "primal", args.wout, str(fin), str(fout)],
                             capture_output=True, text=True)
        if run.returncode != 0:
            return fail(f"override run failed:\n{run.stderr[-800:]}")
        if "sfl: OVERRIDE (" + REASON + ")" not in run.stderr:
            return fail("override not logged")
        with Dataset(fout) as ds:
            if ds.sfl_allow_poor_equilibrium != 1 or ds.sfl_override_reason.strip() != REASON:
                return fail("override not written to out.nc")
            if ds.sfl_quality_code == 0 or "eps_rad" not in ds.sfl_quality_message:
                return fail("quality status not written")
            if ds["quality_value"][0] <= ds["quality_tol"][0]:
                return fail("eps_rad below tolerance for a coarse file")
            if "bmnc_b" not in ds.variables:
                return fail("spectrum missing after override")
        print("ok   override logged and written")
    print("check_sfl_quality_cli PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
