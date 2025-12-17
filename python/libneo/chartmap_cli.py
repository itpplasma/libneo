from __future__ import annotations

import argparse
from pathlib import Path

from .chartmap import write_chartmap_from_stl, write_chartmap_from_vmec_boundary


def _parse_args(argv: list[str] | None) -> argparse.Namespace:
    p = argparse.ArgumentParser(prog="libneo-write-chartmap")
    sub = p.add_subparsers(dest="cmd", required=True)

    p_vmec = sub.add_parser("from-vmec", help="Generate chartmap from VMEC boundary")
    p_vmec.add_argument("wout", type=Path, help="Path to VMEC wout NetCDF")
    p_vmec.add_argument("out", type=Path, help="Output chartmap NetCDF path")
    p_vmec.add_argument("--nrho", type=int, default=33)
    p_vmec.add_argument("--ntheta", type=int, default=65)
    p_vmec.add_argument("--nzeta", type=int, default=33)
    p_vmec.add_argument("--num-field-periods", type=int)
    p_vmec.add_argument("--M", type=int, default=16)
    p_vmec.add_argument("--Nt", type=int, default=256)
    p_vmec.add_argument("--Ng", type=int, nargs=2, default=(256, 256))

    p_stl = sub.add_parser("from-stl", help="Generate chartmap from an STL boundary")
    p_stl.add_argument("stl", type=Path, help="Path to STL file")
    p_stl.add_argument("out", type=Path, help="Output chartmap NetCDF path")
    p_stl.add_argument("--nrho", type=int, default=33)
    p_stl.add_argument("--ntheta", type=int, default=65)
    p_stl.add_argument("--nphi", type=int, default=32)
    p_stl.add_argument("--axis-x", type=float)
    p_stl.add_argument("--axis-y", type=float)
    p_stl.add_argument("--num-field-periods", type=int, default=1)
    p_stl.add_argument("--n-boundary-points", type=int, default=512)
    p_stl.add_argument("--stitch-tol", type=float, default=1.0e-6)
    p_stl.add_argument("--M", type=int, default=16)
    p_stl.add_argument("--Nt", type=int, default=256)
    p_stl.add_argument("--Ng", type=int, nargs=2, default=(256, 256))

    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)

    if args.cmd == "from-vmec":
        write_chartmap_from_vmec_boundary(
            args.wout,
            args.out,
            nrho=int(args.nrho),
            ntheta=int(args.ntheta),
            nzeta=int(args.nzeta),
            num_field_periods=None if args.num_field_periods is None else int(args.num_field_periods),
            M=int(args.M),
            Nt=int(args.Nt),
            Ng=(int(args.Ng[0]), int(args.Ng[1])),
        )
        return 0

    if args.cmd == "from-stl":
        axis_xy = None
        if args.axis_x is not None or args.axis_y is not None:
            if args.axis_x is None or args.axis_y is None:
                raise SystemExit("--axis-x and --axis-y must be provided together")
            axis_xy = (float(args.axis_x), float(args.axis_y))

        write_chartmap_from_stl(
            args.stl,
            args.out,
            nrho=int(args.nrho),
            ntheta=int(args.ntheta),
            nphi=int(args.nphi),
            axis_xy=axis_xy,
            num_field_periods=int(args.num_field_periods),
            n_boundary_points=int(args.n_boundary_points),
            stitch_tol=float(args.stitch_tol),
            M=int(args.M),
            Nt=int(args.Nt),
            Ng=(int(args.Ng[0]), int(args.Ng[1])),
        )
        return 0

    raise RuntimeError("unreachable")


if __name__ == "__main__":
    raise SystemExit(main())

