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


def _load_chartmap_boundary(chartmap_path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, object]]:
    with Dataset(chartmap_path, "r") as ds:
        attrs = {str(k): ds.getncattr(k) for k in ds.ncattrs()}
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
    return theta, zeta, np.stack([R_chart, Z_chart], axis=-1), attrs


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
    p.add_argument("--s-boundary", type=float, default=None)
    p.add_argument("--boundary-offset", type=float, default=None)
    args = p.parse_args(argv)

    wout = args.wout.resolve()
    chartmap = args.chartmap.resolve()
    out = args.out.resolve()

    _import_libneo_from_repo()
    from libneo.vmec import VMECGeometry
    from shapely.geometry import LineString, Point

    theta, zeta, rz_chart, attrs = _load_chartmap_boundary(chartmap)

    if args.s_boundary is None:
        s_boundary = float(attrs.get("s_boundary", 1.0))
    else:
        s_boundary = float(args.s_boundary)
    if not (0.0 < s_boundary <= 1.0):
        raise ValueError("s_boundary must be in (0, 1]")
    if args.boundary_offset is None:
        boundary_offset = float(attrs.get("boundary_offset_m", 0.0))
    else:
        boundary_offset = float(args.boundary_offset)
    if boundary_offset < 0.0:
        raise ValueError("boundary_offset must be >= 0")

    if not wout.exists() or wout.stat().st_size == 0:
        raise RuntimeError(f"missing wout file: {wout}")
    if not chartmap.exists() or chartmap.stat().st_size == 0:
        raise RuntimeError(f"missing chartmap file: {chartmap}")

    geom = VMECGeometry.from_file(str(wout))
    fourier_M = attrs.get("map2disc_M")
    if fourier_M is not None:
        fourier_M = int(fourier_M)

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
    max_chart_vs_vmec = None
    for ax, iz in zip(axes, iz_list):
        zeta_val = float(zeta[iz])
        theta_ref = np.linspace(0.0, 2.0 * np.pi, 4096, endpoint=False, dtype=float)
        R_vmec, Z_vmec, _ = geom.coords_s(s_boundary, theta_ref, zeta_val, use_asym=True)
        R_vmec_adj, Z_vmec_adj, _ = geom.boundary_rz(
            s_boundary,
            theta,
            zeta_val,
            boundary_offset=boundary_offset,
            fourier_M=fourier_M,
            use_asym=True,
        )
        R_axis, Z_axis, _ = geom.coords_s(0.0, np.array([0.0]), zeta_val, use_asym=True)
        R0 = float(R_axis[0])
        Z0 = float(Z_axis[0])

        R_chart = rz_chart[iz, :, 0]
        Z_chart = rz_chart[iz, :, 1]

        ax.plot(R_vmec, Z_vmec, ":", lw=1.0, label="VMEC boundary (LCFS)")
        ax.plot(R_vmec_adj, Z_vmec_adj, "-", lw=1.5, label="VMEC boundary (offset used)")
        ax.plot(R_chart, Z_chart, "--", lw=1.0, label="chartmap rho=1")
        ax.plot([R0], [Z0], "o", ms=3, label="axis")
        ax.set_title(
            f"zeta={zeta_val:.6g} rad  s={s_boundary:.6g}  normal_offset={boundary_offset:.6g} m"
        )
        ax.set_xlabel("R [m]")
        ax.set_ylabel("Z [m]")
        ax.axis("equal")
        ax.legend(loc="best")

        base = LineString(np.column_stack([R_vmec, Z_vmec]))
        dist = np.array([base.distance(Point(float(r), float(z))) for r, z in zip(R_vmec_adj, Z_vmec_adj)])
        dmin = float(np.min(dist))
        dmax = float(np.max(dist))
        dist_min = dmin if dist_min is None else min(dist_min, dmin)
        dist_max = dmax if dist_max is None else max(dist_max, dmax)

        d_chart = np.sqrt((R_chart - R_vmec_adj) ** 2 + (Z_chart - Z_vmec_adj) ** 2)
        dmax_chart = float(np.max(d_chart))
        max_chart_vs_vmec = dmax_chart if max_chart_vs_vmec is None else max(max_chart_vs_vmec, dmax_chart)

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
    if max_chart_vs_vmec is not None:
        print(f"INFO: chartmap rho=1 vs VMEC offset boundary max|d|={max_chart_vs_vmec:.3e} m")
    print(f"PASS: wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
