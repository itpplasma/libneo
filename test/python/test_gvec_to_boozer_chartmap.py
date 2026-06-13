import importlib.util

import numpy as np
import pytest

netCDF4 = pytest.importorskip("netCDF4")

from libneo.gvec_to_boozer_chartmap import assemble_chartmap

TWOPI = 2.0 * np.pi


def _synthetic_si():
    """Small analytic circular-torus Boozer field/geometry in SI units."""
    nfp = 2
    n_rho, n_theta, n_zeta = 5, 8, 6
    r_major = 1.5  # m
    a_minor = 0.4  # m

    rho = np.linspace(1.0e-3, 1.0, n_rho)
    s = rho**2
    theta = np.linspace(0.0, TWOPI, n_theta, endpoint=False)
    zeta = np.linspace(0.0, TWOPI / nfp, n_zeta, endpoint=False)

    minor = a_minor * rho[:, None, None]
    th = theta[None, :, None]
    ze = zeta[None, None, :]
    R = r_major + minor * np.cos(th)
    X = R * np.cos(ze)
    Y = R * np.sin(ze)
    Z = minor * np.sin(th) * np.ones_like(ze)

    Bmod = 2.0 * (1.0 - 0.1 * (minor / r_major) * np.cos(th)) * np.ones_like(ze)
    B_theta = 0.2 * rho**2
    B_phi = 5.0 * np.ones(n_rho)
    A_phi = -0.3 * s
    torflux = -1.2

    return dict(
        rho=rho, s=s, theta=theta, zeta=zeta, X=X, Y=Y, Z=Z, Bmod=Bmod,
        B_theta=B_theta, B_phi=B_phi, A_phi=A_phi, torflux=torflux,
        num_field_periods=nfp,
    )


@pytest.mark.parametrize("flip", ["tor", "pol"])
def test_assemble_units_and_flip(tmp_path, flip):
    d = _synthetic_si()
    out = tmp_path / f"gvec_{flip}.nc"

    returned = assemble_chartmap(
        out, rho=d["rho"], s=d["s"], theta=d["theta"], zeta=d["zeta"],
        X=d["X"], Y=d["Y"], Z=d["Z"], Bmod=d["Bmod"], B_theta=d["B_theta"],
        B_phi=d["B_phi"], A_phi=d["A_phi"], torflux=d["torflux"],
        num_field_periods=d["num_field_periods"], flip=flip,
    )

    # Expected signs after the handedness flip.
    if flip == "tor":
        exp_B_theta = d["B_theta"]
        exp_B_phi = -d["B_phi"]
        exp_A_phi = -d["A_phi"]
        exp_torflux = d["torflux"]
    else:  # pol
        exp_B_theta = -d["B_theta"]
        exp_B_phi = d["B_phi"]
        exp_A_phi = d["A_phi"]
        exp_torflux = -d["torflux"]

    with netCDF4.Dataset(out) as ds:
        # SI -> CGS conversions applied by assemble_chartmap.
        np.testing.assert_allclose(ds.variables["Bmod"][:].max(),
                                   (d["Bmod"] * 1.0e4).max(), rtol=1e-12)
        np.testing.assert_allclose(ds.variables["B_theta"][:],
                                   exp_B_theta * 1.0e6, rtol=1e-12)
        np.testing.assert_allclose(ds.variables["B_phi"][:],
                                   exp_B_phi * 1.0e6, rtol=1e-12)
        np.testing.assert_allclose(ds.variables["A_phi"][:],
                                   exp_A_phi * 1.0e8, rtol=1e-12)
        assert float(ds.getncattr("torflux")) == pytest.approx(
            exp_torflux * 1.0e8)
        assert returned == pytest.approx(exp_torflux * 1.0e8)

        # Geometry m -> cm; stored (zeta, theta, rho).
        stored_x = np.transpose(np.asarray(ds.variables["x"][:]), (2, 1, 0))
        np.testing.assert_allclose(stored_x, d["X"] * 1.0e2, rtol=1e-12)

        assert ds.getncattr("zeta_convention") == "boozer"
        assert ds.getncattr("gvec2chartmap_flip") == flip


def test_assemble_rejects_bad_flip(tmp_path):
    d = _synthetic_si()
    with pytest.raises(ValueError):
        assemble_chartmap(
            tmp_path / "bad.nc", rho=d["rho"], s=d["s"], theta=d["theta"],
            zeta=d["zeta"], X=d["X"], Y=d["Y"], Z=d["Z"], Bmod=d["Bmod"],
            B_theta=d["B_theta"], B_phi=d["B_phi"], A_phi=d["A_phi"],
            torflux=d["torflux"], num_field_periods=d["num_field_periods"],
            flip="nonsense",
        )


@pytest.mark.skipif(
    importlib.util.find_spec("gvec") is None,
    reason="gvec library not installed",
)
def test_gvec_front_end_importable():
    # The real conversion needs a GVEC state file; here we only assert the
    # front end and its gvec dependency import together where gvec is present.
    from libneo.gvec_to_boozer_chartmap import convert_gvec_to_chartmap

    assert callable(convert_gvec_to_chartmap)
