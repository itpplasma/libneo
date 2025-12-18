#!/usr/bin/env python3
"""
Generate a VMEC-boundary chartmap file using the public libneo Python API.

The output NetCDF file is placed in the specified output directory as
``wout_vmec_python.chartmap.nc`` and must satisfy libneo_coordinates
chartmap validation rules.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


def _import_libneo_from_repo() -> None:
    repo_root = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(repo_root / "python"))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="setup_vmec_chartmap_python")
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args(argv)

    out_dir = Path(args.output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    wout = out_dir / "wout.nc"
    outfile_api = out_dir / "wout_vmec_python.chartmap.nc"
    outfile_cli = out_dir / "wout_vmec_cli.chartmap.nc"
    outfile_api_padded = out_dir / "wout_vmec_python_padded.chartmap.nc"

    if not wout.exists() or wout.stat().st_size == 0:
        raise RuntimeError(f"missing wout fixture: {wout}")

    _import_libneo_from_repo()
    from libneo.chartmap import write_chartmap_from_vmec_boundary
    from libneo.chartmap_cli import main as chartmap_cli_main
    from netCDF4 import Dataset

    write_chartmap_from_vmec_boundary(
        wout,
        outfile_api,
        nrho=17,
        ntheta=33,
        nzeta=17,
        M=12,
        Nt=128,
        Ng=(128, 128),
    )

    write_chartmap_from_vmec_boundary(
        wout,
        outfile_api_padded,
        nrho=9,
        ntheta=17,
        nzeta=9,
        boundary_offset=0.5,
        M=8,
        Nt=64,
        Ng=(64, 64),
    )

    chartmap_cli_main(
        [
            "from-vmec",
            str(wout),
            str(outfile_cli),
            "--nrho",
            "17",
            "--ntheta",
            "33",
            "--nzeta",
            "17",
            "--boundary-offset",
            "0.0",
            "--M",
            "12",
            "--Nt",
            "128",
            "--Ng",
            "128",
            "128",
        ]
    )

    for path in (outfile_api, outfile_cli, outfile_api_padded):
        if not path.exists() or path.stat().st_size == 0:
            raise RuntimeError(f"failed to create {path}")
        with Dataset(path, "r") as ds:
            if ds.getncattr("zeta_convention") != "cyl":
                raise RuntimeError("expected zeta_convention=cyl")
            if ds.getncattr("rho_convention") != "unknown":
                raise RuntimeError("expected rho_convention=unknown")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
