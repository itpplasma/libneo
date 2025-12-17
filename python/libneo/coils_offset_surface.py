from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class OffsetSurfaceResult:
    vertices: np.ndarray
    faces: np.ndarray
    origin: np.ndarray
    spacing: np.ndarray
    grid_shape: tuple[int, int, int]


def _sample_segments(a: np.ndarray, b: np.ndarray, step: float) -> np.ndarray:
    if step <= 0.0:
        raise ValueError("step must be > 0")
    ab = b - a
    seg_len = np.linalg.norm(ab, axis=1)
    n = np.maximum(1, np.ceil(seg_len / step).astype(int))

    points: list[np.ndarray] = []
    for i in range(a.shape[0]):
        ni = int(n[i])
        t = np.linspace(0.0, 1.0, ni + 1, dtype=float, endpoint=True)
        pts = a[i][None, :] + t[:, None] * ab[i][None, :]
        points.append(pts[:-1, :])
    return np.vstack(points)


def build_offset_surface_from_segments(
    a: np.ndarray,
    b: np.ndarray,
    *,
    offset_m: float,
    grid_shape: tuple[int, int, int] = (64, 64, 64),
    padding_m: float = 0.20,
    sample_step_m: float | None = None,
) -> OffsetSurfaceResult:
    """
    Build an offset surface dist(x) = offset_m around polyline segments.

    Uses a sampled point cloud + KDTree to approximate distance to the polyline.
    """
    from scipy.spatial import cKDTree
    from skimage.measure import marching_cubes

    if offset_m <= 0.0:
        raise ValueError("offset_m must be > 0")
    if padding_m < 0.0:
        raise ValueError("padding_m must be >= 0")
    if len(grid_shape) != 3 or min(grid_shape) < 16:
        raise ValueError("grid_shape must be 3 ints, each >= 16")

    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.shape != b.shape or a.ndim != 2 or a.shape[1] != 3:
        raise ValueError("a and b must be arrays with shape (nseg, 3)")
    if a.shape[0] < 4:
        raise ValueError("need at least 4 segments")

    mins = np.minimum(a.min(axis=0), b.min(axis=0))
    maxs = np.maximum(a.max(axis=0), b.max(axis=0))
    origin = mins - float(padding_m) - float(offset_m)
    top = maxs + float(padding_m) + float(offset_m)

    nx, ny, nz = (int(grid_shape[0]), int(grid_shape[1]), int(grid_shape[2]))
    span = top - origin
    spacing = span / np.array([nx - 1, ny - 1, nz - 1], dtype=float)

    if sample_step_m is None:
        sample_step_m = 0.5 * float(np.min(spacing))
    sample_step_m = float(sample_step_m)

    sample_points = _sample_segments(a, b, step=sample_step_m)
    if sample_points.shape[0] < 32:
        raise ValueError("insufficient sampled coil points")
    tree = cKDTree(sample_points)

    xs = origin[0] + spacing[0] * np.arange(nx, dtype=float)
    ys = origin[1] + spacing[1] * np.arange(ny, dtype=float)
    zs = origin[2] + spacing[2] * np.arange(nz, dtype=float)
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    pts = np.column_stack((X.ravel(), Y.ravel(), Z.ravel()))

    dist, _ = tree.query(pts, k=1, workers=-1)
    field = dist.reshape((nx, ny, nz))

    verts_zyx, faces, _, _ = marching_cubes(
        field.astype(np.float32),
        level=float(offset_m),
        spacing=(float(spacing[0]), float(spacing[1]), float(spacing[2])),
    )
    verts = verts_zyx + origin[None, :]

    return OffsetSurfaceResult(
        vertices=np.asarray(verts, dtype=float),
        faces=np.asarray(faces, dtype=np.int64),
        origin=np.asarray(origin, dtype=float),
        spacing=np.asarray(spacing, dtype=float),
        grid_shape=(nx, ny, nz),
    )

