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
    seed: np.ndarray


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


def build_inner_offset_surface_from_segments(
    a: np.ndarray,
    b: np.ndarray,
    *,
    offset_m: float,
    seed: np.ndarray,
    grid_shape: tuple[int, int, int] = (64, 64, 64),
    padding_m: float = 0.20,
    sample_step_m: float | None = None,
    window_r_quantiles: tuple[float, float] = (0.0, 1.0),
    window_z_quantiles: tuple[float, float] = (0.0, 1.0),
) -> OffsetSurfaceResult:
    """
    Build a closed surface inside coils as the boundary of the region
    that is at least offset_m away from the coil filaments.

    The surface corresponds to dist(x) = offset_m on the connected component
    of the set {dist(x) >= offset_m} that contains the provided seed point.

    Uses a sampled point cloud + KDTree to approximate distance to the polyline.
    """
    from scipy.spatial import cKDTree
    from scipy import ndimage
    from skimage.measure import marching_cubes

    if offset_m <= 0.0:
        raise ValueError("offset_m must be > 0")
    if padding_m < 0.0:
        raise ValueError("padding_m must be >= 0")
    if len(grid_shape) != 3 or min(grid_shape) < 16:
        raise ValueError("grid_shape must be 3 ints, each >= 16")
    if (
        len(window_r_quantiles) != 2
        or len(window_z_quantiles) != 2
    ):
        raise ValueError("window quantiles must be (lo, hi)")

    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.shape != b.shape or a.ndim != 2 or a.shape[1] != 3:
        raise ValueError("a and b must be arrays with shape (nseg, 3)")
    if a.shape[0] < 4:
        raise ValueError("need at least 4 segments")

    seed = np.asarray(seed, dtype=float)
    if seed.shape != (3,):
        raise ValueError("seed must have shape (3,)")

    mins = np.minimum(a.min(axis=0), b.min(axis=0))
    maxs = np.maximum(a.max(axis=0), b.max(axis=0))
    origin = mins - float(padding_m) - float(offset_m)
    top = maxs + float(padding_m) + float(offset_m)

    nx, ny, nz = (int(grid_shape[0]), int(grid_shape[1]), int(grid_shape[2]))
    span = top - origin
    spacing = span / np.array([nx - 1, ny - 1, nz - 1], dtype=float)

    if sample_step_m is None:
        sample_step_m = 0.35 * float(np.min(spacing))
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

    seed_idx = np.rint((seed - origin) / spacing).astype(int)
    if np.any(seed_idx < 0) or seed_idx[0] >= nx or seed_idx[1] >= ny or seed_idx[2] >= nz:
        raise RuntimeError("seed point outside distance grid bounds")

    rr_coil = np.sqrt(np.r_[a[:, 0], b[:, 0]] ** 2 + np.r_[a[:, 1], b[:, 1]] ** 2)
    zz_coil = np.r_[a[:, 2], b[:, 2]]

    r_lo, r_hi = (float(window_r_quantiles[0]), float(window_r_quantiles[1]))
    z_lo, z_hi = (float(window_z_quantiles[0]), float(window_z_quantiles[1]))
    if not (0.0 <= r_lo < r_hi <= 1.0):
        raise ValueError("window_r_quantiles must satisfy 0 <= lo < hi <= 1")
    if not (0.0 <= z_lo < z_hi <= 1.0):
        raise ValueError("window_z_quantiles must satisfy 0 <= lo < hi <= 1")

    if r_lo == 0.0:
        r_min_base = float(np.min(rr_coil))
    else:
        r_min_base = float(np.quantile(rr_coil, r_lo))
    if r_hi == 1.0:
        r_max_base = float(np.max(rr_coil))
    else:
        r_max_base = float(np.quantile(rr_coil, r_hi))

    if z_lo == 0.0:
        z_min_base = float(np.min(zz_coil))
    else:
        z_min_base = float(np.quantile(zz_coil, z_lo))
    if z_hi == 1.0:
        z_max_base = float(np.max(zz_coil))
    else:
        z_max_base = float(np.quantile(zz_coil, z_hi))

    r_min = max(0.0, r_min_base + float(offset_m))
    r_max = r_max_base - float(offset_m)
    z_min = z_min_base + float(offset_m)
    z_max = z_max_base - float(offset_m)

    if r_max <= r_min or z_max <= z_min:
        raise RuntimeError("could not form interior window for requested offset")

    R = np.sqrt(X**2 + Y**2)
    window = (R >= r_min) & (R <= r_max) & (Z >= z_min) & (Z <= z_max)
    mask = (field >= float(offset_m)) & window
    labels, nlab = ndimage.label(mask)
    if nlab < 1:
        raise RuntimeError("could not form any interior region at requested offset")

    lab = int(labels[int(seed_idx[0]), int(seed_idx[1]), int(seed_idx[2])])
    if lab == 0:
        r = 8
        i0, j0, k0 = int(seed_idx[0]), int(seed_idx[1]), int(seed_idx[2])
        i1 = max(0, i0 - r)
        i2 = min(nx, i0 + r + 1)
        j1 = max(0, j0 - r)
        j2 = min(ny, j0 + r + 1)
        k1 = max(0, k0 - r)
        k2 = min(nz, k0 + r + 1)

        sub_mask = mask[i1:i2, j1:j2, k1:k2]
        if not np.any(sub_mask):
            raise RuntimeError(
                "seed point is not inside the interior region at requested offset"
            )
        sub_field = field[i1:i2, j1:j2, k1:k2]
        flat = int(np.argmax(np.where(sub_mask, sub_field, -np.inf)))
        di, dj, dk = np.unravel_index(flat, sub_mask.shape)
        seed_idx = np.array([i1 + di, j1 + dj, k1 + dk], dtype=int)
        lab = int(labels[int(seed_idx[0]), int(seed_idx[1]), int(seed_idx[2])])
        if lab == 0:
            raise RuntimeError(
                "seed point is not inside the interior region at requested offset"
            )

    keep = labels == lab
    if (
        np.any(keep[0, :, :])
        or np.any(keep[-1, :, :])
        or np.any(keep[:, 0, :])
        or np.any(keep[:, -1, :])
        or np.any(keep[:, :, 0])
        or np.any(keep[:, :, -1])
    ):
        raise RuntimeError(
            "interior region is not enclosed (touches grid boundary); "
            "increase offset/padding or adjust seed"
        )
    signed = field - float(offset_m)
    signed[~keep] = -np.abs(signed[~keep])

    sigma = 0.75
    signed = ndimage.gaussian_filter(signed, sigma=float(sigma), mode="nearest")

    verts_zyx, faces, _, _ = marching_cubes(
        signed.astype(np.float32),
        level=0.0,
        spacing=(float(spacing[0]), float(spacing[1]), float(spacing[2])),
    )
    verts = verts_zyx + origin[None, :]

    return OffsetSurfaceResult(
        vertices=np.asarray(verts, dtype=float),
        faces=np.asarray(faces, dtype=np.int64),
        origin=np.asarray(origin, dtype=float),
        spacing=np.asarray(spacing, dtype=float),
        grid_shape=(nx, ny, nz),
        seed=np.asarray(seed, dtype=float),
    )
