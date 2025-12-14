from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np


def _unit(vec: np.ndarray) -> np.ndarray:
    nrm = float(np.linalg.norm(vec))
    if nrm == 0.0:
        return np.zeros_like(vec)
    return vec / nrm


def _polyline_lengths(points: np.ndarray) -> np.ndarray:
    if points.shape[0] < 2:
        return np.zeros((0,), dtype=float)
    return np.linalg.norm(np.diff(points, axis=0), axis=1)


def _is_closed(points: np.ndarray, tol: float) -> bool:
    if points.shape[0] < 3:
        return False
    return float(np.linalg.norm(points[0] - points[-1])) <= tol


def _polygon_area(points: np.ndarray) -> float:
    if points.shape[0] < 3:
        return 0.0
    xy = points
    if _is_closed(xy, tol=1.0e-12):
        xy = xy[:-1]
    x = xy[:, 0]
    y = xy[:, 1]
    return 0.5 * float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def _split_by_mask(points: np.ndarray, mask: np.ndarray) -> list[np.ndarray]:
    segments: list[np.ndarray] = []
    if points.shape[0] == 0:
        return segments

    start = None
    for i, keep in enumerate(mask):
        if keep and start is None:
            start = i
        if (not keep or i == mask.size - 1) and start is not None:
            end = i if not keep else i + 1
            seg = points[start:end]
            if seg.shape[0] >= 2:
                segments.append(seg)
            start = None
    return segments


def _stitch_segments_rz(segments: list[np.ndarray], tol: float) -> list[np.ndarray]:
    remaining = [seg.copy() for seg in segments if seg.shape[0] >= 2]
    stitched: list[np.ndarray] = []

    def endpoints(seg: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        return seg[0], seg[-1]

    while remaining:
        cur = remaining.pop()
        changed = True
        while changed and remaining:
            changed = False
            cur_start, cur_end = endpoints(cur)

            best_k = None
            best_mode = None
            best_dist = None
            for k, seg in enumerate(remaining):
                s_start, s_end = endpoints(seg)
                candidates = [
                    ("append", float(np.linalg.norm(cur_end - s_start))),
                    ("append_rev", float(np.linalg.norm(cur_end - s_end))),
                    ("prepend", float(np.linalg.norm(cur_start - s_end))),
                    ("prepend_rev", float(np.linalg.norm(cur_start - s_start))),
                ]
                mode, dist = min(candidates, key=lambda x: x[1])
                if dist <= tol and (best_dist is None or dist < best_dist):
                    best_k = k
                    best_mode = mode
                    best_dist = dist

            if best_k is None or best_mode is None:
                continue

            seg = remaining.pop(best_k)
            if best_mode == "append":
                cur = np.vstack((cur, seg[1:]))
            elif best_mode == "append_rev":
                cur = np.vstack((cur, seg[::-1][1:]))
            elif best_mode == "prepend":
                cur = np.vstack((seg[:-1], cur))
            elif best_mode == "prepend_rev":
                cur = np.vstack((seg[::-1][:-1], cur))
            else:
                raise RuntimeError("unreachable")
            changed = True

        stitched.append(cur)

    return stitched


def slice_mesh_rz(
    mesh: object,
    phi: float,
    *,
    axis_xy: tuple[float, float] = (0.0, 0.0),
    stitch_tol: float = 1.0e-6,
) -> list[np.ndarray]:
    """
    Slice a triangle mesh with a constant-toroidal-angle plane and return RZ
    polylines in the half-plane corresponding to the given phi.
    """
    import trimesh

    if not isinstance(mesh, trimesh.Trimesh):
        raise TypeError("mesh must be a trimesh.Trimesh")

    axis_origin = np.array([axis_xy[0], axis_xy[1], 0.0], dtype=float)
    e_r = np.array([np.cos(phi), np.sin(phi), 0.0], dtype=float)
    plane_normal = np.array([-np.sin(phi), np.cos(phi), 0.0], dtype=float)
    plane_origin = axis_origin

    section = mesh.section(plane_origin=plane_origin, plane_normal=plane_normal)
    if section is None:
        return []

    segments_rz: list[np.ndarray] = []
    polylines_3d = list(section.discrete)
    if not polylines_3d:
        for ent in section.entities:
            pts = getattr(ent, "points", None)
            if pts is None:
                continue
            idx = np.asarray(pts, dtype=int)
            if idx.size < 2:
                continue
            polylines_3d.append(section.vertices[idx])

    for poly3 in polylines_3d:
        poly3_rel = poly3 - axis_origin
        r_signed = poly3_rel @ e_r
        z = poly3_rel[:, 2]
        keep = r_signed >= 0.0
        segs = _split_by_mask(np.column_stack((r_signed, z)), keep)
        segments_rz.extend(segs)

    if not segments_rz:
        return []

    return _stitch_segments_rz(segments_rz, tol=stitch_tol)


def classify_rz_contours(
    contours: list[np.ndarray], *, closed_tol: float = 1.0e-6
) -> tuple[np.ndarray, list[np.ndarray]]:
    """
    Pick an outer contour and return remaining closed contours as holes.

    If multiple closed contours exist, the outer contour is the one with the
    largest absolute signed area.
    """
    if not contours:
        raise ValueError("no contours to classify")

    closed: list[np.ndarray] = []
    open_: list[np.ndarray] = []
    for c in contours:
        if _is_closed(c, tol=closed_tol):
            closed.append(c)
        else:
            open_.append(c)

    if closed:
        areas = [abs(_polygon_area(c)) for c in closed]
        i_outer = int(np.argmax(areas))
        outer = closed[i_outer]
        holes = [c for i, c in enumerate(closed) if i != i_outer]
        return outer, holes

    outer = max(open_, key=lambda c: c.shape[0])
    return outer, []


def fill_open_contour(contour: np.ndarray, *, closed_tol: float = 1.0e-10) -> np.ndarray:
    """
    Close an open RZ polyline by adding a C1 bridge using a cubic Bezier curve.
    """
    if contour.ndim != 2 or contour.shape[1] != 2:
        raise ValueError("contour must have shape (n, 2)")
    if contour.shape[0] < 3:
        raise ValueError("contour must have at least 3 points")

    if _is_closed(contour, tol=closed_tol):
        if float(np.linalg.norm(contour[0] - contour[-1])) == 0.0:
            return contour
        return np.vstack((contour, contour[0]))

    start = contour[0]
    end = contour[-1]

    t_start = _unit(contour[1] - contour[0])
    t_end = _unit(contour[-1] - contour[-2])

    gap = float(np.linalg.norm(start - end))
    seg_lengths = _polyline_lengths(contour)
    ds = float(np.median(seg_lengths)) if seg_lengths.size else gap / 16.0
    ds = max(ds, gap / 128.0)
    n_bridge = int(max(16, np.ceil(gap / max(ds, 1.0e-12))))

    alpha = 0.35
    c1 = end + alpha * gap * t_end
    c2 = start - alpha * gap * t_start

    u = np.linspace(0.0, 1.0, n_bridge + 1)
    u = u[1:-1]
    if u.size:
        uu = u[:, None]
        bridge = (
            (1.0 - uu) ** 3 * end
            + 3.0 * (1.0 - uu) ** 2 * uu * c1
            + 3.0 * (1.0 - uu) * uu**2 * c2
            + uu**3 * start
        )
        closed = np.vstack((contour, bridge, start))
    else:
        closed = np.vstack((contour, start))
    return closed


def _resample_closed_contour(contour: np.ndarray, n: int) -> np.ndarray:
    if n < 4:
        raise ValueError("n must be >= 4")
    if not _is_closed(contour, tol=1.0e-8):
        raise ValueError("contour must be closed for resampling")

    pts = contour
    if float(np.linalg.norm(pts[0] - pts[-1])) == 0.0:
        pts = pts[:-1]
    ds = _polyline_lengths(np.vstack((pts, pts[0])))
    s = np.concatenate(([0.0], np.cumsum(ds)))
    total = float(s[-1])
    if total == 0.0:
        raise ValueError("degenerate contour")

    s_targets = np.linspace(0.0, total, n, endpoint=False)
    x = np.interp(s_targets, s, np.r_[pts[:, 0], pts[0, 0]])
    y = np.interp(s_targets, s, np.r_[pts[:, 1], pts[0, 1]])
    resampled = np.column_stack((x, y))
    return np.vstack((resampled, resampled[0]))


def _nan_gaps_from_original(
    filled: np.ndarray, original: np.ndarray, *, dist_factor: float = 3.0
) -> tuple[np.ndarray, list[tuple[int, int]]]:
    """
    Return an array shaped like filled with NaNs where original is missing.

    Indices in the returned list are 1-based and inclusive.
    """
    from scipy.spatial import cKDTree

    pts = filled[:-1]
    tree = cKDTree(original)

    seg_lengths = _polyline_lengths(original)
    ds = float(np.median(seg_lengths)) if seg_lengths.size else 0.0
    ds = max(ds, 1.0e-12)
    thresh = dist_factor * ds

    dist, _ = tree.query(pts, k=1)
    keep = dist <= thresh

    out = filled.copy()
    out[:-1][~keep] = np.nan

    holes: list[tuple[int, int]] = []
    in_gap = False
    start = 0
    for i, ok in enumerate(keep, start=1):
        if (not ok) and (not in_gap):
            in_gap = True
            start = i
        if ok and in_gap:
            holes.append((start, i - 1))
            in_gap = False
    if in_gap:
        holes.append((start, keep.size))

    return out, holes


@dataclass(frozen=True)
class BoundarySlice:
    axis_xy: tuple[float, float]
    phi: float
    outer_filled: np.ndarray
    outer_original_on_filled: np.ndarray
    hole_spans: list[tuple[int, int]]
    hole_contours: list[np.ndarray]


def extract_boundary_slices(
    stl_path: Path,
    *,
    n_phi: int = 32,
    n_boundary_points: int = 512,
    axis_xy: tuple[float, float] | None = None,
    stitch_tol: float = 1.0e-6,
) -> list[BoundarySlice]:
    import trimesh

    mesh = trimesh.load_mesh(stl_path, process=False)
    if not isinstance(mesh, trimesh.Trimesh):
        raise TypeError("expected a single STL mesh")

    if axis_xy is None:
        axis_xy = (float(np.mean(mesh.vertices[:, 0])), float(np.mean(mesh.vertices[:, 1])))

    phi_vals = np.linspace(0.0, 2.0 * np.pi, n_phi, endpoint=False)
    slices: list[BoundarySlice] = []
    for phi in phi_vals:
        contours = slice_mesh_rz(
            mesh, float(phi), axis_xy=axis_xy, stitch_tol=stitch_tol
        )
        if not contours:
            raise RuntimeError(f"no section contours found at phi={phi}")

        outer, holes = classify_rz_contours(contours)
        outer_filled = fill_open_contour(outer)
        outer_filled = _resample_closed_contour(outer_filled, n_boundary_points)

        if _is_closed(outer, tol=1.0e-6):
            outer_on_filled = outer_filled.copy()
            hole_spans: list[tuple[int, int]] = []
        else:
            outer_on_filled, hole_spans = _nan_gaps_from_original(outer_filled, outer)

        slices.append(
            BoundarySlice(
                axis_xy=axis_xy,
                phi=float(phi),
                outer_filled=outer_filled,
                outer_original_on_filled=outer_on_filled,
                hole_spans=hole_spans,
                hole_contours=holes,
            )
        )

    return slices


def write_boundaries_netcdf(slices: list[BoundarySlice], out_path: Path) -> None:
    from netCDF4 import Dataset

    n_phi = len(slices)
    if n_phi == 0:
        raise ValueError("no slices to write")

    n_points = slices[0].outer_filled.shape[0] - 1
    if any((s.outer_filled.shape[0] - 1) != n_points for s in slices):
        raise ValueError("all slices must have the same number of boundary points")

    max_holes = max((len(s.hole_spans) for s in slices), default=0)
    max_holes = max(max_holes, 1)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with Dataset(out_path, "w") as ds:
        ds.createDimension("n_phi", n_phi)
        ds.createDimension("n_points", n_points)
        ds.createDimension("max_holes", max_holes)

        v_phi = ds.createVariable("phi", "f8", ("n_phi",))
        v_axis_x = ds.createVariable("axis_x", "f8")
        v_axis_y = ds.createVariable("axis_y", "f8")
        v_ro = ds.createVariable("R_original", "f8", ("n_phi", "n_points"))
        v_zo = ds.createVariable("Z_original", "f8", ("n_phi", "n_points"))
        v_rf = ds.createVariable("R_filled", "f8", ("n_phi", "n_points"))
        v_zf = ds.createVariable("Z_filled", "f8", ("n_phi", "n_points"))
        v_nholes = ds.createVariable("n_holes", "i4", ("n_phi",))
        v_hs = ds.createVariable("hole_start", "i4", ("n_phi", "max_holes"))
        v_he = ds.createVariable("hole_end", "i4", ("n_phi", "max_holes"))

        v_hs.index_base = 1
        v_he.index_base = 1

        v_phi[:] = [s.phi for s in slices]
        v_axis_x.assignValue(float(slices[0].axis_xy[0]))
        v_axis_y.assignValue(float(slices[0].axis_xy[1]))

        holes_start = np.zeros((n_phi, max_holes), dtype=np.int32)
        holes_end = np.zeros((n_phi, max_holes), dtype=np.int32)

        for i, s in enumerate(slices):
            filled = s.outer_filled[:-1]
            original = s.outer_original_on_filled[:-1]

            v_rf[i, :] = filled[:, 0]
            v_zf[i, :] = filled[:, 1]
            v_ro[i, :] = original[:, 0]
            v_zo[i, :] = original[:, 1]

            v_nholes[i] = len(s.hole_spans)
            for j, (a, b) in enumerate(s.hole_spans[:max_holes]):
                holes_start[i, j] = int(a)
                holes_end[i, j] = int(b)

        v_hs[:, :] = holes_start
        v_he[:, :] = holes_end


def _curve_from_contour(contour: np.ndarray) -> Callable[[np.ndarray], np.ndarray]:
    pts = contour
    if float(np.linalg.norm(pts[0] - pts[-1])) == 0.0:
        pts = pts[:-1]
    ds = _polyline_lengths(np.vstack((pts, pts[0])))
    s = np.concatenate(([0.0], np.cumsum(ds)))
    total = float(s[-1])
    if total == 0.0:
        raise ValueError("degenerate contour")

    def curve(t: np.ndarray) -> np.ndarray:
        tt = np.asarray(t, dtype=float)
        tt = np.mod(tt, 2.0 * np.pi)
        s_target = tt / (2.0 * np.pi) * total
        r = np.interp(s_target, s, np.r_[pts[:, 0], pts[0, 0]])
        z = np.interp(s_target, s, np.r_[pts[:, 1], pts[0, 1]])
        return np.array([r, z])

    return curve


def plot_rz_slices(slices: list[BoundarySlice], out_path: Path) -> None:
    import matplotlib.pyplot as plt

    n_phi = len(slices)
    ncols = min(4, n_phi)
    nrows = int(np.ceil(n_phi / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.0 * ncols, 4.0 * nrows))
    axes_arr = np.atleast_1d(axes).reshape(nrows, ncols)

    for i, s in enumerate(slices):
        ax = axes_arr[i // ncols, i % ncols]
        ax.plot(
            s.outer_original_on_filled[:-1, 0],
            s.outer_original_on_filled[:-1, 1],
            linestyle="--",
            linewidth=1.0,
            label="original",
        )
        for h in s.hole_contours:
            ax.plot(h[:, 0], h[:, 1], color="red", linewidth=1.0)
        ax.plot(
            s.outer_filled[:-1, 0],
            s.outer_filled[:-1, 1],
            linewidth=1.5,
            label="filled",
        )
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(f"phi={s.phi:.3f}")
        ax.set_xlabel("R")
        ax.set_ylabel("Z")

    for j in range(n_phi, nrows * ncols):
        axes_arr[j // ncols, j % ncols].axis("off")

    handles, labels = axes_arr[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper right")

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def plot_grid(slice_: BoundarySlice, out_path: Path) -> None:
    import matplotlib.pyplot as plt
    from map2disc import map as m2d

    curve = _curve_from_contour(slice_.outer_filled)
    bcm = m2d.BoundaryConformingMapping(curve=curve, M=16, Nt=256, Ng=(256, 256))
    bcm.solve_domain2disk()
    bcm.solve_disk2domain()

    rho = np.linspace(0.0, 1.0, 33)
    theta = np.linspace(0.0, 2.0 * np.pi, 64, endpoint=False)
    xy = bcm.eval_rt_1d(rho, theta)
    rr = xy[0]
    zz = xy[1]

    fig, ax = plt.subplots(figsize=(6.0, 6.0))
    for i in range(rr.shape[0]):
        ax.plot(
            np.r_[rr[i, :], rr[i, 0]],
            np.r_[zz[i, :], zz[i, 0]],
            color="C0",
            linewidth=0.6,
            alpha=0.8,
        )
    for j in range(rr.shape[1]):
        ax.plot(rr[:, j], zz[:, j], color="C1", linewidth=0.6, alpha=0.8)

    ax.plot(slice_.outer_filled[:, 0], slice_.outer_filled[:, 1], color="k", linewidth=1.3)
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(f"map2disc grid at phi={slice_.phi:.3f}")
    ax.set_xlabel("R")
    ax.set_ylabel("Z")
    fig.tight_layout()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def plot_mesh_3d(
    slices: list[BoundarySlice],
    stl_path: Path,
    out_path: Path,
    *,
    max_faces: int | None = 8000,
) -> None:
    import matplotlib.pyplot as plt
    import trimesh
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    mesh = trimesh.load_mesh(stl_path, process=False)
    if not isinstance(mesh, trimesh.Trimesh):
        raise TypeError("expected a single STL mesh")

    faces = mesh.faces
    if max_faces is not None and faces.shape[0] > max_faces:
        idx = np.linspace(0, faces.shape[0] - 1, max_faces).astype(int)
        faces = faces[idx]

    tri = mesh.vertices[faces]

    fig = plt.figure(figsize=(8.0, 6.0))
    ax = fig.add_subplot(111, projection="3d")
    col = Poly3DCollection(tri, alpha=0.15, linewidths=0.0)
    col.set_facecolor((0.3, 0.3, 0.7))
    ax.add_collection3d(col)

    for s in slices[:: max(1, len(slices) // 8)]:
        rr = s.outer_filled[:-1, 0]
        zz = s.outer_filled[:-1, 1]
        xx = s.axis_xy[0] + rr * np.cos(s.phi)
        yy = s.axis_xy[1] + rr * np.sin(s.phi)
        ax.plot(xx, yy, zz, linewidth=1.0)

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.autoscale_view()
    fig.tight_layout()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def _parse_args(argv: list[str] | None) -> argparse.Namespace:
    p = argparse.ArgumentParser(prog="libneo.stl_boundary")
    p.add_argument("stl", type=Path)
    p.add_argument("--n-phi", type=int, default=32)
    p.add_argument("--n-boundary-points", type=int, default=512)
    p.add_argument("--axis-x", type=float)
    p.add_argument("--axis-y", type=float)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--plot-3d", type=Path)
    p.add_argument(
        "--plot-3d-max-faces",
        type=int,
        default=8000,
        help="Max triangles to draw in --plot-3d; 0 means draw all",
    )
    p.add_argument("--plot-rz", type=Path)
    p.add_argument("--plot-grid", type=Path)
    p.add_argument("--plot-summary", type=Path)
    p.add_argument("--show", action="store_true")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)

    axis_xy = None
    if args.axis_x is not None or args.axis_y is not None:
        if args.axis_x is None or args.axis_y is None:
            raise SystemExit("--axis-x and --axis-y must be provided together")
        axis_xy = (float(args.axis_x), float(args.axis_y))

    slices = extract_boundary_slices(
        args.stl,
        n_phi=args.n_phi,
        n_boundary_points=args.n_boundary_points,
        axis_xy=axis_xy,
    )
    write_boundaries_netcdf(slices, args.output)

    if args.plot_rz is not None:
        plot_rz_slices(slices, args.plot_rz)
    if args.plot_grid is not None:
        plot_grid(slices[0], args.plot_grid)
    if args.plot_3d is not None:
        max_faces = int(args.plot_3d_max_faces)
        plot_mesh_3d(
            slices,
            args.stl,
            args.plot_3d,
            max_faces=None if max_faces == 0 else max_faces,
        )
    if args.plot_summary is not None:
        plot_rz_slices(slices, args.plot_summary)

    if args.show:
        import matplotlib.pyplot as plt

        plt.show()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
