import subprocess
import sys

import numpy as np
import pytest
from netCDF4 import Dataset

from libneo.nemec import read_nemec, write_vmec_wout
from libneo.vmec import VMECGeometry


def _mode_order(mpol, ntor):
    modes = []
    for m in range(mpol):
        n_min = 0 if m == 0 else -ntor
        for n in range(n_min, ntor + 1):
            modes.append((m, n))
    return modes


def _write_tiny_nemec(path):
    nfp = 2
    nrho = 3
    mpol = 2
    ntor = 1
    modes = _mode_order(mpol, ntor)
    records = np.zeros((nrho, len(modes), 16), dtype=float)

    for j in range(nrho):
        s_scale = 1.0 + 0.2 * j
        for k, (m, n) in enumerate(modes):
            base = s_scale * (10.0 + 2.0 * m + 0.5 * n)
            records[j, k, 0] = base
            records[j, k, 1] = 0.1 * base
            records[j, k, 2] = -0.03 * base
            records[j, k, 3] = 0.07 * base
            records[j, k, 4] = 0.02 * base
            records[j, k, 5] = 2.0 if (m, n) == (0, 0) else 0.01 * base
            records[j, k, 6] = -0.01 * base
            records[j, k, 7] = 0.005 * base
            records[j, k, 8] = -0.02 * base
            records[j, k, 9] = 0.04 * base
            records[j, k, 10] = 0.2 * base
            records[j, k, 11] = 0.3 * base
            records[j, k, 12] = 0.4 * base
            records[j, k, 13] = -0.5 * base
            records[j, k, 14] = 0.06 * base
            records[j, k, 15] = -0.7 * base

    profiles = np.array(
        [
            [0.31, 1.0, 2.0, 0.5, -3.0, 4.0, -0.25, 9.0, 0.1, 0.2, 0.3, 0.4],
            [0.41, 1.1, 2.1, 0.5, -3.1, 4.1, -0.50, 9.1, 0.2, 0.3, 0.4, 0.5],
        ],
        dtype=float,
    )

    values = [0.0, nfp, nrho, mpol, ntor, len(modes), 1, -0.5]
    values.extend(records.ravel())
    values.extend(profiles.ravel())
    path.write_text("\n".join(f"{value:.17e}" for value in values) + "\n")


def _evaluate(parsed, surface, theta, zeta):
    out = np.zeros_like(theta, dtype=float)
    for mode, coeff in enumerate(parsed.rmnc[surface]):
        angle = parsed.xm[mode] * theta - parsed.xn[mode] * zeta
        out += coeff * np.cos(angle) + parsed.rmns[surface, mode] * np.sin(angle)
    return out


def _evaluate_z(parsed, surface, theta, zeta):
    out = np.zeros_like(theta, dtype=float)
    for mode, coeff in enumerate(parsed.zmns[surface]):
        angle = parsed.xm[mode] * theta - parsed.xn[mode] * zeta
        out += coeff * np.sin(angle) + parsed.zmnc[surface, mode] * np.cos(angle)
    return out


def _evaluate_pair(cos_coeff, sin_coeff, xm, xn, theta, zeta):
    out = np.zeros_like(theta, dtype=float)
    for mode, coeff in enumerate(cos_coeff):
        angle = xm[mode] * theta - xn[mode] * zeta
        out += coeff * np.cos(angle) + sin_coeff[mode] * np.sin(angle)
    return out


def _boundary_samples(parsed):
    theta = np.linspace(0.0, 2.0 * np.pi, 21, endpoint=False)
    zeta = np.linspace(0.0, 2.0 * np.pi / parsed.nfp, 21, endpoint=False)
    theta_grid, zeta_grid = np.meshgrid(theta, zeta, indexing="ij")
    r = _evaluate(parsed, -1, theta_grid.ravel(), zeta_grid.ravel())
    z = _evaluate_z(parsed, -1, theta_grid.ravel(), zeta_grid.ravel())
    return r, z


def test_read_nemec_mode_order_and_profiles(tmp_path):
    src = tmp_path / "wout.nemec"
    _write_tiny_nemec(src)

    parsed = read_nemec(src)

    np.testing.assert_array_equal(parsed.xm, np.array([0.0, 0.0, 1.0, 1.0, 1.0]))
    np.testing.assert_array_equal(parsed.xn, np.array([0.0, 2.0, -2.0, 0.0, 2.0]))
    assert parsed.nfp == 2
    assert parsed.ns == 3
    assert parsed.phiedge == pytest.approx(-0.5)
    np.testing.assert_allclose(parsed.phi, np.array([0.0, -0.25, -0.50]))
    np.testing.assert_allclose(parsed.phipf, np.array([-np.pi, -np.pi, -np.pi]))


def test_write_vmec_wout_preserves_geometry_signs(tmp_path):
    src = tmp_path / "wout.nemec"
    out = tmp_path / "wout.nc"
    _write_tiny_nemec(src)

    parsed = read_nemec(src)
    write_vmec_wout(parsed, out)

    geom = VMECGeometry.from_file(str(out))
    theta = np.array([0.2, 1.3, 2.7])
    zeta = 0.37
    surface = 1

    R, Z, zeta_out = geom.coords(surface, theta, zeta)
    assert zeta_out == pytest.approx(zeta)
    np.testing.assert_allclose(R, _evaluate(parsed, surface, theta, zeta))
    np.testing.assert_allclose(Z, _evaluate_z(parsed, surface, theta, zeta))

    with Dataset(out) as ds:
        assert int(ds.variables["signgs"][...]) == -1
        assert int(ds.variables["lasym__logical__"][...]) == 1
        raxis_cc = np.array(ds.variables["raxis_cc"][:])
        zaxis_cs = np.array(ds.variables["zaxis_cs"][:])
        zaxis_cc = np.array(ds.variables["zaxis_cc"][:])
        bmnc = np.array(ds.variables["bmnc"][:])
        bmns = np.array(ds.variables["bmns"][:])
        gmnc = np.array(ds.variables["gmnc"][:])
        gmns = np.array(ds.variables["gmns"][:])
        rmax_surf = float(ds.variables["rmax_surf"][...])
        rmin_surf = float(ds.variables["rmin_surf"][...])
        zmax_surf = float(ds.variables["zmax_surf"][...])
        assert "m*theta - n*zeta*nfp" in ds.source_conventions
    np.testing.assert_allclose(raxis_cc, parsed.rmnc[0, :2])
    np.testing.assert_allclose(zaxis_cs, parsed.zmns[0, :2])

    r_boundary, z_boundary = _boundary_samples(parsed)
    assert rmax_surf >= float(np.max(r_boundary))
    assert rmin_surf <= float(np.min(r_boundary))
    assert zmax_surf >= float(np.max(z_boundary))

    axis_xn = np.array([0.0, 2.0])
    stored_axis_z = np.sum(
        zaxis_cs * np.sin(-axis_xn * zeta)
        + zaxis_cc * np.cos(-axis_xn * zeta)
    )
    expected_axis_z = np.sum(
        parsed.zmns[0, :2] * np.sin(-axis_xn * zeta)
        + parsed.zmnc[0, :2] * np.cos(-axis_xn * zeta)
    )
    assert stored_axis_z == pytest.approx(expected_axis_z)

    theta_grid = np.linspace(0.1, 1.4, 4)
    zeta_grid = np.linspace(0.2, 0.9, 4)
    bsupu = _evaluate_pair(
        parsed.bsupumnc[surface], parsed.bsupumns[surface], parsed.xm, parsed.xn, theta_grid, zeta_grid
    )
    bsupv = _evaluate_pair(
        parsed.bsupvmnc[surface], parsed.bsupvmns[surface], parsed.xm, parsed.xn, theta_grid, zeta_grid
    )
    bsubu = _evaluate_pair(
        parsed.bsubumnc[surface], parsed.bsubumns[surface], parsed.xm, parsed.xn, theta_grid, zeta_grid
    )
    bsubv = _evaluate_pair(
        parsed.bsubvmnc[surface], parsed.bsubvmns[surface], parsed.xm, parsed.xn, theta_grid, zeta_grid
    )
    bmod_expected = np.sqrt(np.maximum(bsupu * bsubu + bsupv * bsubv, 0.0))
    bmod_stored = _evaluate_pair(
        bmnc[surface], bmns[surface], parsed.xm, parsed.xn, theta_grid, zeta_grid
    )
    np.testing.assert_allclose(bmod_stored, bmod_expected, atol=0.8)

    lam_theta = _evaluate_pair(
        parsed.xm * parsed.lmns[surface],
        -parsed.xm * parsed.lmnc[surface],
        parsed.xm,
        parsed.xn,
        theta_grid,
        zeta_grid,
    )
    sqrtg_expected = (1.0 + lam_theta) * parsed.phipf[surface] / bsupv
    sqrtg_stored = _evaluate_pair(
        gmnc[surface], gmns[surface], parsed.xm, parsed.xn, theta_grid, zeta_grid
    )
    np.testing.assert_allclose(sqrtg_stored, sqrtg_expected, atol=0.4)


def test_nemec2vmec_cli_refuses_overwrite(tmp_path):
    src = tmp_path / "wout.nemec"
    out = tmp_path / "wout.nc"
    _write_tiny_nemec(src)
    out.write_text("existing")

    result = subprocess.run(
        [sys.executable, "-m", "libneo.nemec2vmec", str(src), str(out)],
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode != 0
    assert "exists" in result.stderr
