from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class SimpleCoils:
    coils: list[np.ndarray]


def read_simple_coils(path: str | Path) -> SimpleCoils:
    """
    Read SIMPLE coils format:
      n_points
      x y z current
      ...

    Convention:
    - Coils are polygonal chains separated by rows with current == 0.
    - Points are interpreted as meters.
    """
    p = Path(path)
    with p.open("r", encoding="utf-8") as handle:
        header = handle.readline()
        if not header:
            raise ValueError("empty coils file")
        try:
            n_points = int(header.strip().split()[0])
        except Exception as exc:
            raise ValueError("invalid coils header (expected integer n_points)") from exc

        data = np.loadtxt(handle, dtype=float)
    if data.ndim != 2 or data.shape[1] < 4:
        raise ValueError("expected columns: x y z current")
    if data.shape[0] != n_points:
        raise ValueError(f"expected {n_points} points, got {data.shape[0]}")

    xyz = np.asarray(data[:, 0:3], dtype=float)
    cur = np.asarray(data[:, 3], dtype=float)

    coils: list[np.ndarray] = []
    start = 0
    for i in range(cur.size):
        if cur[i] == 0.0:
            if i - start >= 2:
                coils.append(xyz[start:i].copy())
            start = i + 1
    if cur.size - start >= 2:
        coils.append(xyz[start:].copy())

    if not coils:
        raise ValueError("no coils found (need at least 2 points with nonzero current)")
    return SimpleCoils(coils=coils)


def coils_to_segments(coils: SimpleCoils) -> tuple[np.ndarray, np.ndarray]:
    """
    Convert coils into segment endpoints arrays (A, B) with shape (nseg, 3).
    """
    a_list: list[np.ndarray] = []
    b_list: list[np.ndarray] = []
    for pts in coils.coils:
        if pts.shape[0] < 2:
            continue
        a_list.append(pts[:-1, :])
        b_list.append(pts[1:, :])
    if not a_list:
        raise ValueError("no segments found")
    return np.vstack(a_list), np.vstack(b_list)

