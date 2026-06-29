import importlib.util

import numpy as np
import pytest

netCDF4 = pytest.importorskip("netCDF4")

from libneo.boozer_chartmap_writer import (
    METER_TO_CM,
    TESLA_METER2_TO_GAUSS_CM2,
    TESLA_METER_TO_GAUSS_CM,
    TESLA_TO_GAUSS,
    apply_left_handed_flip,
    write_boozer_chartmap,
)

TWOPI = 2.0 * np.pi


def _synthetic_chartmap():
    """Build a small analytic circular-torus chartmap in CGS units."""
    nfp = 3
    n_rho, n_theta, n_zeta = 6, 8, 5
    r_major = 150.0  # cm
    a_minor = 40.0  # cm

    rho = np.linspace(1.0e-3, 1.0, n_rho)
    s = rho**2
    theta = np.linspace(0.0, TWOPI, n_theta, endpoint=False)
    zeta = np.linspace(0.0, TWOPI / nfp, n_zeta, endpoint=False)

    minor = a_minor * rho[:, None, None]
    th = theta[None, :, None]
    ze = zeta[None, None, :]
    R = r_major + minor * np.cos(th)
    Z = minor * np.sin(th) * np.ones_like(ze)
    x = R * np.cos(ze)
    y = R * np.sin(ze)
    z = Z * np.ones((n_rho, n_theta, n_zeta))

    Bmod = 2.0e4 * (1.0 - 0.1 * (minor / r_major) * np.cos(th)) * np.ones_like(ze)

    A_phi = -1.0e8 * s
    B_theta = 1.0e6 * rho**2
    B_phi = 5.0e6 * np.ones(n_rho)
    torflux = -2.5e8

    return dict(
        rho=rho,
        s=s,
        theta=theta,
        zeta=zeta,
        x=x,
        y=y,
        z=z,
        A_phi=A_phi,
        B_theta=B_theta,
        B_phi=B_phi,
        Bmod=Bmod,
        num_field_periods=nfp,
        torflux=torflux,
    )


def test_unit_constants():
    assert TESLA_TO_GAUSS == 1.0e4
    assert METER_TO_CM == 1.0e2
    assert TESLA_METER_TO_GAUSS_CM == 1.0e6
    assert TESLA_METER2_TO_GAUSS_CM2 == 1.0e8


def test_left_handed_flip():
    bphi = np.array([1.0, 2.0, 3.0])
    aphi = np.array([-4.0, 5.0, -6.0])
    fb, fa = apply_left_handed_flip(bphi, aphi)
    np.testing.assert_array_equal(fb, -bphi)
    np.testing.assert_array_equal(fa, -aphi)


def test_write_roundtrip(tmp_path):
    data = _synthetic_chartmap()
    out = tmp_path / "chartmap.nc"

    write_boozer_chartmap(
        out,
        rho=data["rho"],
        s=data["s"],
        theta=data["theta"],
        zeta=data["zeta"],
        x=data["x"],
        y=data["y"],
        z=data["z"],
        A_phi=data["A_phi"],
        B_theta=data["B_theta"],
        B_phi=data["B_phi"],
        Bmod=data["Bmod"],
        num_field_periods=data["num_field_periods"],
        torflux=data["torflux"],
        booz2chartmap_source="synthetic",
    )

    assert out.is_file()

    with netCDF4.Dataset(out) as ds:
        n_rho = data["rho"].size
        n_theta = data["theta"].size
        n_zeta = data["zeta"].size

        assert ds.dimensions["rho"].size == n_rho
        assert ds.dimensions["s"].size == n_rho
        assert ds.dimensions["theta"].size == n_theta
        assert ds.dimensions["zeta"].size == n_zeta

        for name in ("rho", "s", "theta", "zeta", "x", "y", "z",
                     "A_phi", "B_theta", "B_phi", "Bmod",
                     "num_field_periods"):
            assert name in ds.variables, f"missing variable {name}"

        # Global attributes defined by the libneo chartmap format.
        assert ds.getncattr("rho_convention") == "rho_tor"
        assert ds.getncattr("zeta_convention") == "boozer"
        assert ds.getncattr("boozer_field") == 1
        assert float(ds.getncattr("rho_lcfs")) == pytest.approx(data["rho"][-1])
        assert float(ds.getncattr("torflux")) == pytest.approx(data["torflux"])
        assert ds.getncattr("booz2chartmap_source") == "synthetic"

        np.testing.assert_allclose(ds.variables["rho"][:], data["rho"])
        np.testing.assert_allclose(ds.variables["s"][:], data["s"])
        np.testing.assert_allclose(ds.variables["theta"][:], data["theta"])
        np.testing.assert_allclose(ds.variables["zeta"][:], data["zeta"])

        np.testing.assert_allclose(ds.variables["A_phi"][:], data["A_phi"])
        np.testing.assert_allclose(ds.variables["B_theta"][:], data["B_theta"])
        np.testing.assert_allclose(ds.variables["B_phi"][:], data["B_phi"])
        assert int(ds.variables["num_field_periods"][...]) == data[
            "num_field_periods"
        ]

        # Geometry/field stored as (zeta, theta, rho); transpose back to read.
        for name, arr in [("x", data["x"]), ("y", data["y"]),
                          ("z", data["z"]), ("Bmod", data["Bmod"])]:
            stored = np.asarray(ds.variables[name][:])
            assert stored.shape == (n_zeta, n_theta, n_rho)
            np.testing.assert_allclose(
                np.transpose(stored, (2, 1, 0)), arr, rtol=1e-12, atol=1e-12
            )

        assert ds.variables["x"].units == "cm"
        assert ds.variables["A_phi"].radial_abscissa == "s"


def test_shape_validation(tmp_path):
    data = _synthetic_chartmap()
    out = tmp_path / "bad.nc"
    with pytest.raises(ValueError):
        write_boozer_chartmap(
            out,
            rho=data["rho"],
            s=data["s"],
            theta=data["theta"],
            zeta=data["zeta"],
            x=data["x"][:-1],  # wrong shape
            y=data["y"],
            z=data["z"],
            A_phi=data["A_phi"],
            B_theta=data["B_theta"],
            B_phi=data["B_phi"],
            Bmod=data["Bmod"],
            num_field_periods=data["num_field_periods"],
            torflux=data["torflux"],
        )


@pytest.mark.skipif(
    importlib.util.find_spec("booz_xform") is None,
    reason="booz_xform not installed",
)
def test_booz_xform_front_end(tmp_path):
    pytest.importorskip("scipy")
    import booz_xform as bx

    from libneo.booz_xform_to_boozer_chartmap import convert_boozmn_to_chartmap

    # Minimal axisymmetric tokamak via booz_xform's own example if available.
    b = bx.Booz_xform()
    if not hasattr(b, "read_wout"):
        pytest.skip("booz_xform API lacks read_wout for a tiny case")
    pytest.skip("no bundled VMEC wout fixture available for a tiny real case")
