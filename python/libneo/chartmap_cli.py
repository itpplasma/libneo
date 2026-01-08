from __future__ import annotations

import argparse
from pathlib import Path

from .chartmap import (
    write_chartmap_from_stl,
    write_chartmap_from_vmec_and_stl,
    write_chartmap_from_vmec_boundary,
    write_chartmap_from_vmec_extended,
)


def _parse_args(argv: list[str] | None) -> argparse.Namespace:
    p = argparse.ArgumentParser(prog="libneo-write-chartmap")
    sub = p.add_subparsers(dest="cmd", required=True)

    p_vmec = sub.add_parser("from-vmec", help="Generate chartmap from VMEC boundary")
    p_vmec.add_argument("wout", type=Path, help="Path to VMEC wout NetCDF")
    p_vmec.add_argument("out", type=Path, help="Output chartmap NetCDF path")
    p_vmec.add_argument("--nrho", type=int, default=33)
    p_vmec.add_argument("--ntheta", type=int, default=65)
    p_vmec.add_argument("--nzeta", type=int, default=33)
    p_vmec.add_argument("--s-boundary", type=float, default=1.0)
    p_vmec.add_argument("--boundary-offset", type=float, default=0.0)
    p_vmec.add_argument(
        "--boundary-param",
        type=str,
        default="arc",
        choices=("arc", "theta"),
        help="Boundary parameterization for map2disc (arc or theta).",
    )
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

    p_ext = sub.add_parser(
        "from-vmec-extended",
        help="Generate chartmap from VMEC with cubic Hermite extension beyond LCFS",
    )
    p_ext.add_argument("wout", type=Path, help="Path to VMEC wout NetCDF")
    p_ext.add_argument("out", type=Path, help="Output chartmap NetCDF path")
    p_ext.add_argument("--nrho", type=int, default=33)
    p_ext.add_argument("--ntheta", type=int, default=65)
    p_ext.add_argument("--nzeta", type=int, default=33)
    p_ext.add_argument(
        "--nrho-vmec",
        type=int,
        help="Number of rho points inside LCFS (alternative to --rho-lcfs)",
    )
    p_ext.add_argument(
        "--rho-lcfs",
        type=float,
        help="Radial location of LCFS in [0,1] (alternative to --nrho-vmec)",
    )
    p_ext.add_argument(
        "--boundary-offset",
        type=float,
        default=0.1,
        help="Extension distance beyond LCFS in meters",
    )
    p_ext.add_argument("--num-field-periods", type=int)

    p_wall = sub.add_parser(
        "from-vmec-to-wall",
        help="Generate chartmap from VMEC extended to STL wall boundary",
    )
    p_wall.add_argument("wout", type=Path, help="Path to VMEC wout NetCDF")
    p_wall.add_argument("stl", type=Path, help="Path to STL wall geometry file")
    p_wall.add_argument("out", type=Path, help="Output chartmap NetCDF path")
    p_wall.add_argument("--nrho", type=int, default=33)
    p_wall.add_argument("--ntheta", type=int, default=65)
    p_wall.add_argument("--nzeta", type=int, default=33)
    p_wall.add_argument(
        "--rho-lcfs",
        type=float,
        default=0.8,
        help="Radial location of LCFS in (0,1)",
    )
    p_wall.add_argument(
        "--n-boundary-points",
        type=int,
        default=512,
        help="Points for STL boundary extraction",
    )
    p_wall.add_argument(
        "--stitch-tol",
        type=float,
        default=1.0e-6,
        help="Tolerance for stitching open contours",
    )
    p_wall.add_argument("--num-field-periods", type=int)

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
            s_boundary=float(args.s_boundary),
            boundary_offset=float(args.boundary_offset),
            boundary_param=str(args.boundary_param),
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

    if args.cmd == "from-vmec-extended":
        write_chartmap_from_vmec_extended(
            args.wout,
            args.out,
            nrho=int(args.nrho),
            ntheta=int(args.ntheta),
            nzeta=int(args.nzeta),
            nrho_vmec=None if args.nrho_vmec is None else int(args.nrho_vmec),
            rho_lcfs=None if args.rho_lcfs is None else float(args.rho_lcfs),
            boundary_offset=float(args.boundary_offset),
            num_field_periods=None if args.num_field_periods is None else int(args.num_field_periods),
        )
        return 0

    if args.cmd == "from-vmec-to-wall":
        write_chartmap_from_vmec_and_stl(
            args.wout,
            args.stl,
            args.out,
            nrho=int(args.nrho),
            ntheta=int(args.ntheta),
            nzeta=int(args.nzeta),
            rho_lcfs=float(args.rho_lcfs),
            n_boundary_points=int(args.n_boundary_points),
            stitch_tol=float(args.stitch_tol),
            num_field_periods=None if args.num_field_periods is None else int(args.num_field_periods),
        )
        return 0

    raise RuntimeError("unreachable")


if __name__ == "__main__":
    raise SystemExit(main())
