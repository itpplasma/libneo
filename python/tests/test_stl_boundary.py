from __future__ import annotations

import numpy as np


def _is_closed(polyline: np.ndarray, tol: float = 1.0e-10) -> bool:
    if polyline.ndim != 2 or polyline.shape[1] != 2 or polyline.shape[0] < 3:
        return False
    return np.linalg.norm(polyline[0] - polyline[-1]) <= tol


def test_fill_open_contour_closes_and_preserves_circle() -> None:
    from libneo.stl_boundary import fill_open_contour

    r0 = 1.7
    a = 0.35

    t = np.linspace(0.2 * np.pi, 1.8 * np.pi, 200)
    open_circle = np.column_stack((r0 + a * np.cos(t), a * np.sin(t)))

    filled = fill_open_contour(open_circle)

    assert _is_closed(filled, tol=1.0e-8)

    d = np.sqrt((filled[:, 0] - r0) ** 2 + filled[:, 1] ** 2)
    assert np.max(np.abs(d - a)) < 5.0e-3


def test_slice_mesh_rz_torus_returns_closed_outer_contour() -> None:
    import trimesh

    from libneo.stl_boundary import classify_rz_contours, slice_mesh_rz

    r0 = 1.7
    a = 0.35

    mesh = trimesh.creation.torus(
        major_radius=r0,
        minor_radius=a,
        major_sections=128,
        minor_sections=96,
    )

    contours = slice_mesh_rz(mesh, phi=0.0)
    outer, holes = classify_rz_contours(contours)

    assert holes == []
    assert _is_closed(outer, tol=1.0e-6)

    d = np.sqrt((outer[:, 0] - r0) ** 2 + outer[:, 1] ** 2)
    assert float(np.std(d)) < 5.0e-3


def test_write_boundaries_netcdf_format(tmp_path) -> None:
    from netCDF4 import Dataset

    from libneo.stl_boundary import BoundarySlice, write_boundaries_netcdf

    axis_xy = (0.1, -0.2)
    filled = np.array(
        [
            [1.0, 0.0],
            [0.0, 1.0],
            [-1.0, 0.0],
            [0.0, -1.0],
            [1.0, 0.0],
        ],
        dtype=float,
    )
    original = filled.copy()
    original[1, :] = np.nan

    slices = [
        BoundarySlice(
            axis_xy=axis_xy,
            phi=0.0,
            outer_filled=filled,
            outer_original_on_filled=original,
            hole_spans=[(2, 2)],
            hole_contours=[],
        ),
        BoundarySlice(
            axis_xy=axis_xy,
            phi=np.pi,
            outer_filled=filled,
            outer_original_on_filled=filled,
            hole_spans=[],
            hole_contours=[],
        ),
    ]

    out = tmp_path / "boundaries.nc"
    write_boundaries_netcdf(slices, out)

    with Dataset(out, "r") as ds:
        assert set(ds.dimensions) == {"n_phi", "n_points", "max_holes"}
        assert ds.dimensions["n_phi"].size == 2
        assert ds.dimensions["n_points"].size == 4

        for name in [
            "phi",
            "axis_x",
            "axis_y",
            "R_original",
            "Z_original",
            "R_filled",
            "Z_filled",
            "n_holes",
            "hole_start",
            "hole_end",
        ]:
            assert name in ds.variables

        assert float(ds.variables["axis_x"].getValue()) == axis_xy[0]
        assert float(ds.variables["axis_y"].getValue()) == axis_xy[1]

        assert int(ds.variables["n_holes"][0]) == 1
        assert int(ds.variables["hole_start"][0, 0]) == 2
        assert int(ds.variables["hole_end"][0, 0]) == 2
