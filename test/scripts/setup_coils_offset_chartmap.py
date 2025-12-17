#!/usr/bin/env python3
"""
Generate a chartmap from a SIMPLE coils file by offset-surface + map2disc.

Writes:
- coils.simple
- coils_offset.chartmap.nc
- coils_offset_3d.png (coils + surface)
- coils_offset_slices.png (RZ slices)
- coils_offset_rho_contours.png (interior rho curves in RZ slices)
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np


def _import_libneo_from_repo() -> None:
    repo_root = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(repo_root / "python"))


def _write_circular_coil_simple(path: Path, *, R0: float, z0: float, n: int) -> None:
    t = np.linspace(0.0, 2.0 * np.pi, int(n), endpoint=True)
    x = R0 * np.cos(t)
    y = R0 * np.sin(t)
    z = z0 + 0.0 * t

    # SIMPLE format: current=0 marks coil separator (use last point as separator)
    cur = np.ones_like(t)
    cur[-1] = 0.0

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write(f"{t.size}\n")
        for xi, yi, zi, ii in zip(x, y, z, cur):
            handle.write(f"{xi:.14E} {yi:.14E} {zi:.14E} {ii:.14E}\n")


def _write_stellarator_like_coils_simple(
    path: Path,
    *,
    n_coils: int,
    n_points_per_coil: int,
    R0: float,
    a: float,
    n_waves: int,
    z_scale: float,
) -> None:
    t = np.linspace(0.0, 2.0 * np.pi, int(n_points_per_coil), endpoint=True)
    coils_xyz: list[np.ndarray] = []
    currents: list[np.ndarray] = []

    for k in range(int(n_coils)):
        phi = 2.0 * np.pi * float(k) / float(n_coils)
        R = R0 + a * np.cos(t)
        Z = z_scale * a * np.sin(t)
        xx = R * np.cos(phi)
        yy = R * np.sin(phi)
        zz = Z
        xyz = np.column_stack((xx, yy, zz))
        coils_xyz.append(xyz)

        cur = np.ones((t.size,), dtype=float)
        cur[-1] = 0.0
        currents.append(cur)

    n_total = sum(arr.shape[0] for arr in coils_xyz)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write(f"{n_total}\n")
        for xyz, cur in zip(coils_xyz, currents):
            for (xi, yi, zi), ii in zip(xyz, cur):
                handle.write(f"{xi:.14E} {yi:.14E} {zi:.14E} {ii:.14E}\n")


def _plot_coils_and_surface(out_dir: Path, coils_xyz: np.ndarray, mesh: object) -> None:
    import matplotlib.pyplot as plt
    import trimesh

    if not isinstance(mesh, trimesh.Trimesh):
        raise TypeError("mesh must be a trimesh.Trimesh")

    fig = plt.figure(figsize=(7.0, 6.0))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot(coils_xyz[:, 0], coils_xyz[:, 1], coils_xyz[:, 2], color="C3", linewidth=1.2)

    faces = mesh.faces
    if faces.shape[0] > 8000:
        idx = np.linspace(0, faces.shape[0] - 1, 8000).astype(int)
        faces = faces[idx]
    tri = mesh.vertices[faces]
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    col = Poly3DCollection(tri, linewidths=0.05, alpha=0.25)
    col.set_facecolor((0.4, 0.55, 0.8, 0.25))
    col.set_edgecolor((0.1, 0.1, 0.1, 0.25))
    ax.add_collection3d(col)

    ax.set_xlabel("X [m]")
    ax.set_ylabel("Y [m]")
    ax.set_zlabel("Z [m]")

    mins = np.minimum(mesh.vertices.min(axis=0), coils_xyz.min(axis=0))
    maxs = np.maximum(mesh.vertices.max(axis=0), coils_xyz.max(axis=0))
    center = 0.5 * (mins + maxs)
    span = float(np.max(maxs - mins))
    span = max(span, 1.0e-12)
    half = 0.6 * span
    ax.set_xlim(center[0] - half, center[0] + half)
    ax.set_ylim(center[1] - half, center[1] + half)
    ax.set_zlim(center[2] - half, center[2] + half)
    ax.set_box_aspect((1.0, 1.0, 1.0))
    ax.view_init(elev=18.0, azim=35.0)
    fig.tight_layout()
    fig.savefig(out_dir / "coils_offset_3d.png", dpi=160)
    plt.close(fig)


def _plot_rz_slices(out_dir: Path, slices: list[object], rho_curves: list[np.ndarray]) -> None:
    import matplotlib.pyplot as plt

    n_phi = len(slices)
    ncols = min(4, max(1, n_phi))
    nrows = int(np.ceil(n_phi / ncols))
    def make_fig() -> tuple[object, np.ndarray]:
        fig, axes = plt.subplots(nrows, ncols, figsize=(4.0 * ncols, 4.0 * nrows))
        return fig, np.atleast_1d(axes).reshape(nrows, ncols)

    fig, axes_arr = make_fig()
    for i, s in enumerate(slices):
        ax = axes_arr[i // ncols, i % ncols]
        ax.plot(s.outer_filled[:, 0], s.outer_filled[:, 1], color="k", linewidth=1.2)
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(f"phi={s.phi:.3f}")
        ax.set_xlabel("R [m]")
        ax.set_ylabel("Z [m]")
    for j in range(n_phi, nrows * ncols):
        axes_arr[j // ncols, j % ncols].axis("off")
    fig.tight_layout()
    fig.savefig(out_dir / "coils_offset_slices.png", dpi=160)
    plt.close(fig)

    fig, axes_arr = make_fig()
    for i, s in enumerate(slices):
        ax = axes_arr[i // ncols, i % ncols]
        ax.plot(s.outer_filled[:, 0], s.outer_filled[:, 1], color="k", linewidth=1.2)
        for rr, zz in rho_curves[i]:
            ax.plot(rr, zz, color="C0", linewidth=0.6, alpha=0.8)
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(f"phi={s.phi:.3f}")
        ax.set_xlabel("R [m]")
        ax.set_ylabel("Z [m]")
    for j in range(n_phi, nrows * ncols):
        axes_arr[j // ncols, j % ncols].axis("off")
    fig.tight_layout()
    fig.savefig(out_dir / "coils_offset_rho_contours.png", dpi=160)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="setup_coils_offset_chartmap")
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args(argv)

    out_dir = Path(args.output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    coils_file = out_dir / "coils.simple"
    chartmap_file = out_dir / "coils_offset.chartmap.nc"

    _write_stellarator_like_coils_simple(
        coils_file,
        n_coils=64,
        n_points_per_coil=181,
        R0=1.7,
        a=0.45,
        n_waves=1,
        z_scale=1.0,
    )

    _import_libneo_from_repo()
    from libneo.chartmap import write_chartmap_from_coils_offset_surface
    from libneo.coils_simple import read_simple_coils
    from libneo.coils_simple import coils_to_segments
    from libneo.coils_offset_surface import build_inner_offset_surface_from_segments
    from libneo.stl_boundary import extract_boundary_slices_from_mesh

    coils = read_simple_coils(coils_file)
    a, b = coils_to_segments(coils)
    pts = np.vstack(coils.coils)
    axis_xy = (float(np.mean(pts[:, 0])), float(np.mean(pts[:, 1])))

    a = a.copy()
    b = b.copy()
    a[:, 0] -= axis_xy[0]
    a[:, 1] -= axis_xy[1]
    b[:, 0] -= axis_xy[0]
    b[:, 1] -= axis_xy[1]

    pts_raw = np.vstack(coils.coils)
    rr_raw = np.sqrt((pts_raw[:, 0] - axis_xy[0]) ** 2 + (pts_raw[:, 1] - axis_xy[1]) ** 2)
    seed = np.array([float(np.mean(rr_raw)), 0.0, float(np.mean(pts_raw[:, 2]))], dtype=float)

    surf = build_inner_offset_surface_from_segments(
        a,
        b,
        offset_m=0.10,
        seed=seed,
        grid_shape=(72, 72, 72),
        padding_m=0.20,
    )

    import trimesh

    mesh = trimesh.Trimesh(vertices=surf.vertices, faces=surf.faces, process=True)
    mesh.merge_vertices()
    mesh.remove_unreferenced_vertices()
    trimesh.repair.fix_normals(mesh)
    trimesh.repair.fill_holes(mesh)
    trimesh.smoothing.filter_taubin(mesh, lamb=0.5, nu=-0.53, iterations=12)
    if not mesh.is_watertight:
        raise RuntimeError("offset surface mesh is not watertight")

    coils_xyz = np.vstack(coils.coils).copy()
    coils_xyz[:, 0] -= axis_xy[0]
    coils_xyz[:, 1] -= axis_xy[1]
    _plot_coils_and_surface(out_dir, coils_xyz, mesh)

    slices = extract_boundary_slices_from_mesh(
        mesh,
        n_phi=8,
        n_boundary_points=512,
        stitch_tol=1.0e-6,
    )

    # Generate chartmap (full nphi) for validation + rho curves
    write_chartmap_from_coils_offset_surface(
        coils_file,
        chartmap_file,
        offset_cm=10.0,
        nrho=33,
        ntheta=65,
        nphi=32,
        num_field_periods=1,
        grid_shape=(72, 72, 72),
        padding_cm=20.0,
        axis_xy=axis_xy,
        seed_rz=(float(seed[0]), float(seed[2])),
        M=12,
        Nt=128,
        Ng=(128, 128),
    )

    if not chartmap_file.exists() or chartmap_file.stat().st_size == 0:
        raise RuntimeError(f"failed to create {chartmap_file}")

    # Load chartmap and plot some interior rho curves for each shown slice
    from netCDF4 import Dataset

    with Dataset(chartmap_file, "r") as ds:
        rho = np.array(ds.variables["rho"][:], dtype=float)
        theta = np.array(ds.variables["theta"][:], dtype=float)
        zeta = np.array(ds.variables["zeta"][:], dtype=float)
        x = np.array(ds.variables["x"][:], dtype=float)
        y = np.array(ds.variables["y"][:], dtype=float)
        z = np.array(ds.variables["z"][:], dtype=float)

    def rz_at(iz: int, ir: int) -> tuple[np.ndarray, np.ndarray]:
        rr = np.sqrt(x[iz, :, ir] ** 2 + y[iz, :, ir] ** 2) / 100.0
        zz = z[iz, :, ir] / 100.0
        return np.r_[rr, rr[0]], np.r_[zz, zz[0]]

    # Match slices to nearest zeta index in chartmap
    rho_idx = [
        int(0.25 * (rho.size - 1)),
        int(0.5 * (rho.size - 1)),
        int(0.75 * (rho.size - 1)),
    ]
    rho_curves: list[np.ndarray] = []
    for s in slices:
        iz = int(np.argmin(np.abs(zeta - float(s.phi))))
        curves = []
        for ir in rho_idx:
            curves.append(rz_at(iz, ir))
        rho_curves.append(np.array(curves, dtype=object))

    _plot_rz_slices(out_dir, slices, rho_curves)

    for name in (
        "coils_offset_3d.png",
        "coils_offset_slices.png",
        "coils_offset_rho_contours.png",
    ):
        p = out_dir / name
        if not p.exists() or p.stat().st_size == 0:
            raise RuntimeError(f"missing plot artifact: {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
