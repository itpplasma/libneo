#!/usr/bin/env python3
"""
Download a small VMEC wout NetCDF fixture for libneo tests.

The file is placed in the specified output directory as ``wout.nc``.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from netCDF4 import Dataset

from test_data import (PINNED, WOUT_CIRCULAR_TOKAMAK_URL, WOUT_LI383_URL,
                       WOUT_QA_URL, WOUT_QH_URL, WOUT_UPDOWN_TOKAMAK_URL, download)

# Landreman-Paul QA and QH (nfp = 2, 4); NCSX li383 with net toroidal current
# (Boozer I != 0); circular and up-down asymmetric (lasym = T) tokamaks.
FIXTURES = {
    "qa": ("wout.nc", WOUT_QA_URL),
    "li383": ("wout_li383.nc", WOUT_LI383_URL),
    "qh": ("wout_qh.nc", WOUT_QH_URL),
    "circular_tokamak": ("wout_circular_tokamak.nc", WOUT_CIRCULAR_TOKAMAK_URL),
    "updown_tokamak": ("wout_updown_tokamak.nc", WOUT_UPDOWN_TOKAMAK_URL),
}


def validate_wout(path: Path) -> None:
    with Dataset(path, "r") as ds:
        if "nfp" not in ds.variables:
            raise RuntimeError("wout fixture missing nfp variable")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="setup_vmec_wout")
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory in which to place wout.nc",
    )
    parser.add_argument("--case", choices=sorted(FIXTURES), default="qa")
    args = parser.parse_args(argv)

    name, url = FIXTURES[args.case]
    outfile = download(url, Path(args.output_dir).resolve() / name, PINNED[url])
    validate_wout(outfile)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

