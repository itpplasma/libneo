from __future__ import annotations

import argparse
from pathlib import Path

from .nemec import read_nemec, write_vmec_wout


def _parse_args(argv: list[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(prog="nemec2vmec")
    parser.add_argument("input", type=Path, help="NEMEC text wout file")
    parser.add_argument("output", type=Path, help="Output VMEC wout NetCDF file")
    parser.add_argument("--force", action="store_true", help="Overwrite output file")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    if args.output.exists() and not args.force:
        raise SystemExit(f"{args.output} exists; pass --force to overwrite")
    data = read_nemec(args.input)
    write_vmec_wout(data, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
