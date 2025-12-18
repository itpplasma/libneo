from __future__ import annotations

from pathlib import Path
from urllib.request import urlopen

import numpy as np
import pytest


def _download(url: str, out_path: Path) -> None:
    with urlopen(url, timeout=60) as resp, open(out_path, "wb") as f:
        f.write(resp.read())


def _has_self_intersection(poly: np.ndarray) -> bool:
    if poly.ndim != 2 or poly.shape[1] != 2 or poly.shape[0] < 4:
        raise ValueError("expected polyline (n,2) with n>=4")

    def _seg_intersect(a, b, c, d) -> bool:
        if max(a[0], b[0]) < min(c[0], d[0]) or max(c[0], d[0]) < min(a[0], b[0]):
            return False
        if max(a[1], b[1]) < min(c[1], d[1]) or max(c[1], d[1]) < min(a[1], b[1]):
            return False

        def orient(p, q, r) -> float:
            return (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0])

        def onseg(p, q, r) -> bool:
            return (
                min(p[0], r[0]) <= q[0] <= max(p[0], r[0])
                and min(p[1], r[1]) <= q[1] <= max(p[1], r[1])
            )

        o1 = orient(a, b, c)
        o2 = orient(a, b, d)
        o3 = orient(c, d, a)
        o4 = orient(c, d, b)
        if (o1 > 0) != (o2 > 0) and (o3 > 0) != (o4 > 0):
            return True

        eps = 1.0e-12
        if abs(o1) < eps and onseg(a, c, b):
            return True
        if abs(o2) < eps and onseg(a, d, b):
            return True
        if abs(o3) < eps and onseg(c, a, d):
            return True
        if abs(o4) < eps and onseg(c, b, d):
            return True
        return False

    pts = np.asarray(poly, dtype=float)
    if np.linalg.norm(pts[0] - pts[-1]) > 1.0e-12:
        pts = np.vstack([pts, pts[0]])

    n = pts.shape[0] - 1
    for i in range(n):
        a = pts[i]
        b = pts[i + 1]
        for j in range(i + 2, n):
            if i == 0 and j == n - 1:
                continue
            c = pts[j]
            d = pts[j + 1]
            if _seg_intersect(a, b, c, d):
                return True
    return False


@pytest.mark.network
def test_vmec_boundary_normal_offset_is_normal_and_non_intersecting(tmp_path) -> None:
    pytest.importorskip("shapely")
    from libneo.vmec import VMECGeometry
    from shapely.geometry import LineString, Point

    wout_url = "https://princetonuniversity.github.io/STELLOPT/examples/wout_ncsx_c09r00_fixed.nc"
    wout = tmp_path / "wout_ncsx.nc"
    _download(wout_url, wout)

    geom = VMECGeometry.from_file(str(wout))
    theta = np.linspace(0.0, 2.0 * np.pi, 512, endpoint=False, dtype=float)
    zeta = 0.0

    Rn, Zn, _ = geom.boundary_rz(1.0, theta, zeta, boundary_offset=0.1, use_asym=True)

    theta_base = np.linspace(0.0, 2.0 * np.pi, 4096, endpoint=False, dtype=float)
    R_base, Z_base, _ = geom.coords_s(1.0, theta_base, zeta, use_asym=True)
    base = LineString(np.column_stack([R_base, Z_base]))

    dist = np.array([base.distance(Point(float(r), float(z))) for r, z in zip(Rn, Zn)])
    assert float(np.max(np.abs(dist - 0.1))) < 5.0e-3

    poly = np.column_stack([Rn, Zn])
    assert _has_self_intersection(poly) is False
