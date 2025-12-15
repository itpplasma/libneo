#!/usr/bin/env python3
"""
Download a small VMEC wout NetCDF fixture for libneo tests.

The file is placed in the specified output directory as ``wout.nc``.
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import urllib.request

from netCDF4 import Dataset

WOUT_URL = (
    "https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/"
    "wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc"
)
WOUT_SHA256 = "333ab8cef64c4b4f406e76209266fc3bf4d7976e9c28c162983c2120e112e771"


def sha256sum(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def validate_wout(path: Path) -> None:
    with Dataset(path, "r") as ds:
        if "nfp" not in ds.variables:
            raise RuntimeError("wout fixture missing nfp variable")


def download_wout(path: Path) -> None:
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    try:
        with urllib.request.urlopen(WOUT_URL, timeout=120) as response, tmp_path.open(
            "wb"
        ) as handle:
            handle.write(response.read())
        got = sha256sum(tmp_path)
        if got != WOUT_SHA256:
            raise RuntimeError(f"sha256 mismatch: expected {WOUT_SHA256}, got {got}")
        tmp_path.replace(path)
    finally:
        tmp_path.unlink(missing_ok=True)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="setup_vmec_wout")
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory in which to place wout.nc",
    )
    args = parser.parse_args(argv)

    out_dir = Path(args.output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    outfile = out_dir / "wout.nc"

    if outfile.exists() and outfile.stat().st_size > 0:
        got = sha256sum(outfile)
        if got == WOUT_SHA256:
            validate_wout(outfile)
            return 0
        outfile.unlink()

    download_wout(outfile)
    validate_wout(outfile)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

