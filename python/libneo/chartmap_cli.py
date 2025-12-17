from __future__ import annotations

import argparse
from pathlib import Path

from .chartmap import (
    write_chartmap_from_coils_offset_surface,
    write_chartmap_from_stl,
    write_chartmap_from_vmec_boundary,
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

    p_coils = sub.add_parser(
        "from-coils-offset",
        help="Generate chartmap from SIMPLE coils by offset surface + map2disc",
    )
    p_coils.add_argument("coils", type=Path, help="Path to coils file (SIMPLE format)")
    p_coils.add_argument("out", type=Path, help="Output chartmap NetCDF path")
    p_coils.add_argument("--offset-cm", type=float, required=True)
    p_coils.add_argument("--nrho", type=int, default=33)
    p_coils.add_argument("--ntheta", type=int, default=65)
    p_coils.add_argument("--nphi", type=int, default=32)
    p_coils.add_argument("--num-field-periods", type=int, default=1)
    p_coils.add_argument("--grid-shape", type=int, nargs=3, default=(72, 72, 72))
    p_coils.add_argument("--padding-cm", type=float, default=20.0)
    p_coils.add_argument("--sample-step-cm", type=float)
    p_coils.add_argument("--axis-x", type=float)
    p_coils.add_argument("--axis-y", type=float)
    p_coils.add_argument("--seed-R", type=float)
    p_coils.add_argument("--seed-Z", type=float)
    p_coils.add_argument("--n-boundary-points", type=int, default=512)
    p_coils.add_argument("--stitch-tol", type=float, default=1.0e-6)
    p_coils.add_argument("--M", type=int, default=16)
    p_coils.add_argument("--Nt", type=int, default=256)
    p_coils.add_argument("--Ng", type=int, nargs=2, default=(256, 256))

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

    if args.cmd == "from-coils-offset":
        grid_shape = (int(args.grid_shape[0]), int(args.grid_shape[1]), int(args.grid_shape[2]))
        axis_xy = None
        if args.axis_x is not None or args.axis_y is not None:
            if args.axis_x is None or args.axis_y is None:
                raise SystemExit("--axis-x and --axis-y must be provided together")
            axis_xy = (float(args.axis_x), float(args.axis_y))
        seed_rz = None
        if args.seed_R is not None or args.seed_Z is not None:
            if args.seed_R is None or args.seed_Z is None:
                raise SystemExit("--seed-R and --seed-Z must be provided together")
            seed_rz = (float(args.seed_R), float(args.seed_Z))
        write_chartmap_from_coils_offset_surface(
            args.coils,
            args.out,
            offset_cm=float(args.offset_cm),
            nrho=int(args.nrho),
            ntheta=int(args.ntheta),
            nphi=int(args.nphi),
            num_field_periods=int(args.num_field_periods),
            grid_shape=grid_shape,
            padding_cm=float(args.padding_cm),
            sample_step_cm=None if args.sample_step_cm is None else float(args.sample_step_cm),
            axis_xy=axis_xy,
            seed_rz=seed_rz,
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
