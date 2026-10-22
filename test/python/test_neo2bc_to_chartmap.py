"""Independent manufactured-mode tests against the installed libneo package.

Install the branch being tested in the test environment. Modifying a package's
__path__ can mix worktree modules with modules from an existing editable install.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

SURFACES = (0.25, 5.0 / 12.0, 7.0 / 12.0, 0.75, 11.0 / 12.0)


def _write_bc(path: Path, *, perturbation: bool, surfaces=SURFACES) -> None:
    nper = 2 if perturbation else 1
    with path.open("w", encoding="utf-8") as stream:
        stream.write("CC Boozer-coordinate data file\n")
        stream.write("CC manufactured independent oracle\n")
        stream.write(
            " m0b   n0b  nsurf  nper    flux [Tm^2]        a [m]          R [m]\n"
        )
        stream.write(f"  1  0  {len(surfaces)}  {nper}  -2.0  0.5  1.6\n")
        for surface in surfaces:
            stream.write(
                "        s               iota           Jpol/nper          Itor"
                "            pprime         sqrt g(0,0)\n"
            )
            stream.write(
                "                                          [A]           [A]"
                "             [Pa]         (dV/ds)/nper\n"
            )
            stream.write(
                f" {surface:.8e}  4.0e-1  {5.0e5 / nper:.8e}  2.0e4"
                "  0.0  1.0\n"
            )
            stream.write(
                "    m    n      rmnc [m]         rmns [m]         zmnc [m]"
                "         zmns [m]         vmnc [ ]         vmns [ ]"
                "         bmnc [T]         bmns [T]\n"
            )
            if perturbation:
                stream.write(
                    " -1  1  2.0e-3  -3.0e-3  4.0e-3  -5.0e-3"
                    "  6.0e-4  -7.0e-4  2.0e-2  -3.0e-2\n"
                )
                stream.write(
                    "  1  1  -8.0e-3  9.0e-3  -1.0e-2  1.1e-2"
                    "  -1.2e-3  1.3e-3  4.0e-2  5.0e-2\n"
                )
            else:
                stream.write(
                    "  0  0  1.6  0.0  0.0  0.0  0.0  0.0  2.0  0.0\n"
                )
                stream.write(
                    "  1  0  1.0e-1  0.0  0.0  1.0e-1"
                    "  0.0  2.0e-2  1.0e-1  0.0\n"
                )


@pytest.fixture()
def manufactured_bc(tmp_path):
    axis = tmp_path / "axis.bc"
    perturbation = tmp_path / "pert.bc"
    _write_bc(axis, perturbation=False)
    _write_bc(perturbation, perturbation=True)
    return axis, perturbation


def _minus_phase_sum(m, n, cos_coeff, sin_coeff, theta, zeta, nper):
    angle = np.asarray(m) * theta - np.asarray(n) * nper * zeta
    return np.sum(np.asarray(cos_coeff) * np.cos(angle)) + np.sum(
        np.asarray(sin_coeff) * np.sin(angle)
    )


def test_composition_matches_direct_fourier_oracle(manufactured_bc):
    from libneo.neo2bc_to_chartmap import combine_neo2_boozer_files

    axis, perturbation = manufactured_bc
    combined = combine_neo2_boozer_files(axis, perturbation, scale=0.25)

    assert combined.nper == 2
    np.testing.assert_array_equal(combined.m[0], [0, 1, -1, 1])
    np.testing.assert_array_equal(combined.n[0], [0, 0, -1, -1])

    theta = 0.37
    zeta = 0.21
    actual = _minus_phase_sum(
        combined.m[2],
        combined.n[2],
        combined.bmnc[2],
        combined.bmns[2],
        theta,
        zeta,
        combined.nper,
    )
    expected = 2.0 + 0.1 * np.cos(theta)
    expected += 0.25 * (
        0.02 * np.cos(-theta + 2.0 * zeta)
        - 0.03 * np.sin(-theta + 2.0 * zeta)
        + 0.04 * np.cos(theta + 2.0 * zeta)
        + 0.05 * np.sin(theta + 2.0 * zeta)
    )
    assert actual == pytest.approx(expected)


def test_bmod_only_keeps_axisymmetric_geometry(manufactured_bc):
    from libneo.neo2bc_to_chartmap import combine_neo2_boozer_files

    axis, perturbation = manufactured_bc
    combined = combine_neo2_boozer_files(
        axis, perturbation, scale=0.5, geometry="axisymmetric"
    )

    np.testing.assert_allclose(combined.rmnc[1][2:], 0.0)
    np.testing.assert_allclose(combined.rmns[1][2:], 0.0)
    np.testing.assert_allclose(combined.zmnc[1][2:], 0.0)
    np.testing.assert_allclose(combined.zmns[1][2:], 0.0)
    np.testing.assert_allclose(combined.vmnc[1][2:], 0.0)
    np.testing.assert_allclose(combined.vmns[1][2:], 0.0)
    np.testing.assert_allclose(combined.bmnc[1][2:], [0.01, 0.02])
    np.testing.assert_allclose(combined.bmns[1][2:], [-0.015, 0.025])


def test_full_geometry_scales_every_perturbation_coefficient(manufactured_bc):
    from libneo.neo2bc_to_chartmap import combine_neo2_boozer_files

    axis, perturbation = manufactured_bc
    combined = combine_neo2_boozer_files(
        axis, perturbation, scale=0.5, geometry="full"
    )

    np.testing.assert_allclose(combined.rmnc[1][2:], [0.001, -0.004])
    np.testing.assert_allclose(combined.rmns[1][2:], [-0.0015, 0.0045])
    np.testing.assert_allclose(combined.zmnc[1][2:], [0.002, -0.005])
    np.testing.assert_allclose(combined.zmns[1][2:], [-0.0025, 0.0055])
    np.testing.assert_allclose(combined.vmnc[1][2:], [0.0003, -0.0006])
    np.testing.assert_allclose(combined.vmns[1][2:], [-0.00035, 0.00065])


def test_converter_writes_units_sign_and_hash_provenance(
    manufactured_bc, tmp_path
):
    netCDF4 = pytest.importorskip("netCDF4")
    pytest.importorskip("scipy")
    from libneo.neo2bc_to_chartmap import convert_neo2bc_to_chartmap

    axis, perturbation = manufactured_bc
    output = tmp_path / "field.chartmap.nc"
    convert_neo2bc_to_chartmap(
        axis,
        output,
        perturbation=perturbation,
        scale=0.25,
        geometry="axisymmetric",
        nrho=8,
        ntheta=16,
        nzeta=12,
    )

    with netCDF4.Dataset(output) as dataset:
        assert dataset.data_model == "NETCDF3_64BIT_OFFSET"
        assert dataset.booz2chartmap_source == str(axis)
        assert dataset.neo2bc_source_fourier_phase == "m*theta+n*nper*phi"
        assert dataset.neo2bc_output_fourier_phase == "m*theta-n*nfp*zeta"
        assert dataset.neo2bc_source_length_unit == "m"
        assert dataset.neo2bc_source_field_unit == "T"
        assert dataset.neo2bc_output_unit_system == "CGS-Gaussian"
        assert dataset.neo2bc_perturbation_scale == pytest.approx(0.25)
        assert dataset.neo2bc_axis_phase_scale == pytest.approx(2.0)
        assert dataset.neo2bc_geometry == "axisymmetric"
        assert len(dataset.neo2bc_axis_sha256) == 64
        assert len(dataset.neo2bc_perturbation_sha256) == 64
        assert int(dataset.variables["num_field_periods"][:]) == 2

        theta = np.asarray(dataset.variables["theta"][:])
        zeta = np.asarray(dataset.variables["zeta"][:])
        bmod = np.asarray(dataset.variables["Bmod"][:])
        itheta = 3
        izeta = 4
        expected_tesla = 2.0 + 0.1 * np.cos(theta[itheta])
        expected_tesla += 0.25 * (
            0.02 * np.cos(-theta[itheta] + 2.0 * zeta[izeta])
            - 0.03 * np.sin(-theta[itheta] + 2.0 * zeta[izeta])
            + 0.04 * np.cos(theta[itheta] + 2.0 * zeta[izeta])
            + 0.05 * np.sin(theta[itheta] + 2.0 * zeta[izeta])
        )
        assert bmod[izeta, itheta, -1] == pytest.approx(expected_tesla * 1.0e4)
        assert np.all(np.isfinite(bmod))


def test_axis_cylindrical_angle_is_invariant_when_output_nper_changes(
    manufactured_bc, tmp_path
):
    netCDF4 = pytest.importorskip("netCDF4")
    pytest.importorskip("scipy")
    from libneo.neo2bc_to_chartmap import (
        combine_neo2_boozer_files,
        convert_neo2bc_to_chartmap,
    )

    axis, perturbation = manufactured_bc
    combined = combine_neo2_boozer_files(axis, perturbation)
    assert combined.nper == 2
    assert combined.vmns[0][1] == pytest.approx(4.0e-2)

    output = tmp_path / "angle.chartmap.nc"
    convert_neo2bc_to_chartmap(
        axis,
        output,
        perturbation=perturbation,
        geometry="axisymmetric",
        nrho=8,
        ntheta=16,
        nzeta=12,
    )

    with netCDF4.Dataset(output) as dataset:
        theta = np.asarray(dataset.variables["theta"][:])
        zeta = np.asarray(dataset.variables["zeta"][:])
        x = np.asarray(dataset.variables["x"][:])
        y = np.asarray(dataset.variables["y"][:])
        itheta = 3
        izeta = 4
        radius_m = 1.6 + 0.1 * np.cos(theta[itheta])
        axis_phase = 2.0e-2 * 2.0 * np.pi * np.sin(theta[itheta])
        phi_cyl = zeta[izeta] + axis_phase
        assert x[izeta, itheta, -1] == pytest.approx(
            100.0 * radius_m * np.cos(phi_cyl)
        )
        assert y[izeta, itheta, -1] == pytest.approx(
            100.0 * radius_m * np.sin(phi_cyl)
        )


def test_rejects_incompatible_radial_grids(manufactured_bc):
    from libneo.boozer import BoozerFile
    from libneo.neo2bc_to_chartmap import combine_boozer_data

    axis_path, perturbation_path = manufactured_bc
    axis = BoozerFile(str(axis_path))
    perturbation = BoozerFile(str(perturbation_path))
    perturbation.s[1] += 1.0e-3

    with pytest.raises(ValueError, match="radial grids"):
        combine_boozer_data(axis, perturbation, scale=1.0)


def test_minus_phase_source_preserves_perturbation_mode_sign(manufactured_bc):
    from libneo.neo2bc_to_chartmap import combine_neo2_boozer_files

    axis, perturbation = manufactured_bc
    combined = combine_neo2_boozer_files(
        axis, perturbation, source_phase="minus"
    )

    np.testing.assert_array_equal(combined.n[0], [0, 0, 1, 1])


def test_converter_accepts_first_half_grid_surface(tmp_path):
    pytest.importorskip("netCDF4")
    pytest.importorskip("scipy")
    from libneo.neo2bc_to_chartmap import convert_neo2bc_to_chartmap

    surfaces = (0.1, 0.3, 0.5, 0.7, 0.9)
    axis = tmp_path / "axis-half.bc"
    perturbation = tmp_path / "pert-half.bc"
    output = tmp_path / "half.chartmap.nc"
    _write_bc(axis, perturbation=False, surfaces=surfaces)
    _write_bc(perturbation, perturbation=True, surfaces=surfaces)

    convert_neo2bc_to_chartmap(
        axis,
        output,
        perturbation=perturbation,
        nrho=8,
        ntheta=8,
        nzeta=8,
    )
    assert output.is_file()


def test_nonuniform_surfaces_keep_radial_field_and_currents(tmp_path):
    """Nonuniform source labels must survive the intermediate boozmn file."""
    netCDF4 = pytest.importorskip("netCDF4")
    from libneo.neo2bc_to_chartmap import convert_neo2bc_to_chartmap

    surfaces = (0.02, 0.07, 0.19, 0.36, 0.61, 0.83, 0.98)
    axis = tmp_path / "nonuniform.bc"
    _write_bc(axis, perturbation=False, surfaces=surfaces)
    text = axis.read_text()
    for s in surfaces:
        text = text.replace(
            f" {s:.8e}  4.0e-1  5.00000000e+05  2.0e4",
            f" {s:.8e} {0.4 + 0.1*s:.12e} {5e5 + 1e5*s:.12e} 2.0e4",
        )
    lines = text.splitlines()
    index = 0
    for j, line in enumerate(lines):
        if line.startswith("  0  0 "):
            lines[j] = line.replace("  2.0  0.0", f"  {2 + 0.2*surfaces[index]**2:.12e}  0.0")
            index += 1
    axis.write_text("\n".join(lines) + "\n")
    output = tmp_path / "chartmap.nc"
    convert_neo2bc_to_chartmap(axis, output, nrho=19, ntheta=12, nzeta=4)
    with netCDF4.Dataset(output) as ds:
        rho = np.asarray(ds["rho"][:])
        theta = np.asarray(ds["theta"][:])
        interior = (rho**2 > surfaces[1]) & (rho**2 < surfaces[-2])
        expected = 2 + 0.2*rho[interior]**4 + 0.1*np.cos(theta[:, None])
        np.testing.assert_allclose(ds["Bmod"][0, :, interior]/1e4, expected, rtol=2e-10)
        s = np.asarray(ds["s"][:])
        np.testing.assert_allclose(ds["B_phi"][:]/1e6, (5e5 + 1e5*rho**2)*2e-7, rtol=2e-10)
        torflux = float(ds.torflux)
        expected_a = -torflux*(0.4*s + 0.05*s*s)
        np.testing.assert_allclose(ds["A_phi"][:], expected_a, rtol=2e-10)


def test_converter_retains_perturbation_knots_on_distinct_grids(tmp_path):
    netCDF4 = pytest.importorskip("netCDF4")
    from scipy.interpolate import CubicSpline
    from libneo.boozer import BoozerFile
    from libneo.neo2bc_to_chartmap import convert_neo2bc_to_chartmap

    axis, pert = tmp_path/"axis.bc", tmp_path/"pert.bc"
    coarse = (.02, .07, .19, .36, .61, .83, .98)
    fine = (.02, .07, .19, .33, .45, .49, .5, .51, .57, .61, .83, .98)
    _write_bc(axis, perturbation=False, surfaces=coarse)
    _write_bc(pert, perturbation=True, surfaces=fine)
    lines = pert.read_text().splitlines()
    index = 0
    for j, line in enumerate(lines):
        if line.startswith(" -1  1 "):
            amplitude = .02*(1 + np.exp(-((fine[index] - .5)/.03)**2))
            lines[j] = line.replace("  2.0e-2  -3.0e-2", f"  {amplitude:.14e}  -3.0e-2")
            index += 1
    pert.write_text("\n".join(lines)+"\n")
    output = tmp_path/"field.nc"
    convert_neo2bc_to_chartmap(axis, output, perturbation=pert, nrho=25, ntheta=12, nzeta=6)
    source = BoozerFile(str(pert))
    with netCDF4.Dataset(output) as ds:
        s = np.asarray(ds["rho"][:])**2
        interior = (s > .19) & (s < .83)
        s = s[interior]
        theta, zeta = float(ds["theta"][2]), float(ds["zeta"][1])
        phase = np.asarray(source.m[0])*theta + np.asarray(source.n[0])*source.nper*zeta
        cos = CubicSpline(source.s, np.asarray(source.bmnc), axis=0)(s)
        sin = CubicSpline(source.s, np.asarray(source.bmns), axis=0)(s)
        expected = 2 + .1*np.cos(theta) + cos@np.cos(phase) + sin@np.sin(phase)
        np.testing.assert_allclose(ds["Bmod"][1, 2, interior]/1e4, expected, atol=2e-10)
