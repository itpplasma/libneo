#!/usr/bin/env python3
"""
Create a PNG artifact comparing a chartmap boundary against the VMEC boundary used to generate it.

This is intended for CI artifact generation for map2disc chartmap tests, including
outward-padded boundaries.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np
from netCDF4 import Dataset


def _import_libneo_from_repo() -> None:
    repo_root = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(repo_root / "python"))


def _load_chartmap_boundary(chartmap_path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with Dataset(chartmap_path, "r") as ds:
        rho = np.array(ds.variables["rho"][:], dtype=float)
        theta = np.array(ds.variables["theta"][:], dtype=float)
        zeta = np.array(ds.variables["zeta"][:], dtype=float)
        x = np.array(ds.variables["x"][:], dtype=float)
        y = np.array(ds.variables["y"][:], dtype=float)
        z = np.array(ds.variables["z"][:], dtype=float)

        if ds.getncattr("zeta_convention") != "cyl":
            raise RuntimeError("expected zeta_convention=cyl")
        if ds.getncattr("rho_convention") != "unknown":
            raise RuntimeError("expected rho_convention=unknown")

    if rho.ndim != 1 or theta.ndim != 1 or zeta.ndim != 1:
        raise RuntimeError("expected rho/theta/zeta to be 1D")
    if rho.size < 2 or theta.size < 2 or zeta.size < 2:
        raise RuntimeError("expected rho/theta/zeta lengths >= 2")

    if x.shape != (zeta.size, theta.size, rho.size):
        raise RuntimeError(f"unexpected x shape {x.shape}")

    ir = int(rho.size - 1)
    R_chart = np.sqrt(x[:, :, ir] ** 2 + y[:, :, ir] ** 2) / 100.0
    Z_chart = z[:, :, ir] / 100.0
    return theta, zeta, np.stack([R_chart, Z_chart], axis=-1)


def _select_zeta_indices(nzeta: int) -> list[int]:
    if nzeta < 2:
        return [0]
    candidates = [0, 1, nzeta // 2, nzeta - 2, nzeta - 1]
    out: list[int] = []
    for i in candidates:
        ii = int(max(0, min(nzeta - 1, i)))
        if ii not in out:
            out.append(ii)
    return out


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="plot_vmec_map2disc_chartmap_boundary")
    p.add_argument("--wout", type=Path, required=True)
    p.add_argument("--chartmap", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--s-boundary", type=float, default=1.0)
    p.add_argument("--boundary-scale", type=float, default=1.0)
    p.add_argument("--boundary-padding", type=float, default=0.0)
    args = p.parse_args(argv)

    wout = args.wout.resolve()
    chartmap = args.chartmap.resolve()
    out = args.out.resolve()

    s_boundary = float(args.s_boundary)
    if not (0.0 < s_boundary <= 1.0):
        raise ValueError("s_boundary must be in (0, 1]")
    boundary_scale = float(args.boundary_scale)
    if boundary_scale <= 0.0:
        raise ValueError("boundary_scale must be > 0")
    boundary_padding = float(args.boundary_padding)
    if boundary_padding < 0.0:
        raise ValueError("boundary_padding must be >= 0")

    if not wout.exists() or wout.stat().st_size == 0:
        raise RuntimeError(f"missing wout file: {wout}")
    if not chartmap.exists() or chartmap.stat().st_size == 0:
        raise RuntimeError(f"missing chartmap file: {chartmap}")

    _import_libneo_from_repo()
    from libneo.vmec import VMECGeometry

    theta, zeta, rz_chart = _load_chartmap_boundary(chartmap)
    geom = VMECGeometry.from_file(str(wout))

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # noqa: E402  (after Agg selection)

    iz_list = _select_zeta_indices(int(zeta.size))
    nrows = int(len(iz_list))
    fig, axes = plt.subplots(nrows, 1, figsize=(6, 2.8 * nrows), constrained_layout=True)
    if nrows == 1:
        axes = [axes]

    dist_min = None
    dist_max = None
    for ax, iz in zip(axes, iz_list):
        zeta_val = float(zeta[iz])
        R_vmec, Z_vmec, _ = geom.boundary_rz(s_boundary, theta, zeta_val, use_asym=True)
        R_vmec_adj, Z_vmec_adj, _ = geom.boundary_rz(
            s_boundary,
            theta,
            zeta_val,
            boundary_scale=boundary_scale,
            boundary_padding=boundary_padding,
            use_asym=True,
        )
        R_axis, Z_axis, _ = geom.coords_s(0.0, np.array([0.0]), zeta_val, use_asym=True)
        R0 = float(R_axis[0])
        Z0 = float(Z_axis[0])

        R_chart = rz_chart[iz, :, 0]
        Z_chart = rz_chart[iz, :, 1]

        ax.plot(R_vmec, Z_vmec, ":", lw=1.0, label="VMEC boundary (unpadded)")
        ax.plot(R_vmec_adj, Z_vmec_adj, "-", lw=1.5, label="VMEC boundary (used)")
        ax.plot(R_chart, Z_chart, "--", lw=1.0, label="chartmap rho=1")
        ax.plot([R0], [Z0], "o", ms=3, label="axis")
        ax.set_title(
            f"zeta={zeta_val:.6g} rad  s={s_boundary:.6g}  scale={boundary_scale:.6g}  pad={boundary_padding:.6g} m"
        )
        ax.set_xlabel("R [m]")
        ax.set_ylabel("Z [m]")
        ax.axis("equal")
        ax.legend(loc="best")

        dist = np.sqrt((R_vmec_adj - R_vmec) ** 2 + (Z_vmec_adj - Z_vmec) ** 2)
        dmin = float(np.min(dist))
        dmax = float(np.max(dist))
        dist_min = dmin if dist_min is None else min(dist_min, dmin)
        dist_max = dmax if dist_max is None else max(dist_max, dmax)

    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    plt.close(fig)

    if not out.exists() or out.stat().st_size == 0:
        raise RuntimeError(f"failed to write plot: {out}")

    if dist_min is not None and dist_max is not None:
        print(
            f"INFO: boundary adjustment |d| range over plotted slices: "
            f"[{dist_min:.6g}, {dist_max:.6g}] m"
        )
    print(f"PASS: wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
