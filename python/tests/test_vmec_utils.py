import numpy as np
import pytest

from libneo import vmec_utils


def synthetic_eqdsk():
    theta = 2.0 * np.pi * np.linspace(0.0, 1.0, 73, endpoint=False) ** 1.4
    nprof = 32
    return {
        "Lcfs": np.column_stack((6.2 + 2.0 * np.cos(theta),
                                  0.2 + 1.4 * np.sin(theta))),
        "PsiaxisVs": 0.0,
        "PsiedgeVs": 1.0,
        "fprof": np.ones(nprof),
        "ptotprof": np.linspace(2.0e5, 0.0, nprof),
        "qprof": np.linspace(1.0, 3.0, nprof),
        "dpressdpsiprof": np.linspace(-2.0e5, -1.0e4, nprof),
        "fdfdpsiprof": np.linspace(5.0e4, 1.0e3, nprof),
        "Ip": 15.0e6,
        "Rpsi0": 6.2,
        "Zpsi0": 0.2,
    }


def test_boundary_fourier_resamples_nonuniform_geometric_angle():
    theta = 2.0 * np.pi * np.linspace(0.0, 1.0, 97, endpoint=False) ** 1.7
    major_radius = 6.2
    minor_radius_r = 2.0
    minor_radius_z = 1.4
    r = major_radius + minor_radius_r * np.cos(theta)
    z = 0.25 + minor_radius_z * np.sin(theta)
    data = {}

    vmec_utils.compute_boundary_fourier(
        r, z, data, mpol=4, ntor=0, theta=theta
    )

    assert data["rbc"][0, 0] == pytest.approx(major_radius, abs=1.0e-6)
    assert data["rbc"][1, 0] == pytest.approx(minor_radius_r, abs=2.0e-3)
    assert data["zbc"][0, 0] == pytest.approx(0.25, abs=1.0e-6)
    assert data["zbs"][1, 0] == pytest.approx(minor_radius_z, abs=2.0e-3)
    assert np.max(np.abs(data["rbs"])) < 2.0e-3
    assert np.max(np.abs(data["zbc"][1:])) < 2.0e-3


def test_eqdsk_iota_polynomial_closes_on_iota_endpoint(monkeypatch):
    monkeypatch.setattr(vmec_utils.eqdsk_base, "read_eqdsk",
                        lambda _: synthetic_eqdsk())

    converted = vmec_utils.eqdsk2vmec_gfile("synthetic.eqdsk")

    assert np.sum(converted["ai"]) == pytest.approx(
        converted["ai_aux_f"][-1], abs=1.0e-12
    )
    assert np.sum(converted["ai"]) != pytest.approx(
        converted["ac_aux_f"][-1]
    )


def test_eqdsk_conversion_uses_discrete_profile_quadrature(monkeypatch):
    monkeypatch.setattr(vmec_utils.eqdsk_base, "read_eqdsk",
                        lambda _: synthetic_eqdsk())
    monkeypatch.setattr(
        vmec_utils,
        "quad",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("adaptive quadrature must not be used")
        ),
        raising=False,
    )

    converted = vmec_utils.eqdsk2vmec_gfile("synthetic.eqdsk")

    assert converted["phiedge"] == pytest.approx(4.0 * np.pi)


def test_eqdsk_default_retains_twenty_boundary_harmonics(monkeypatch):
    monkeypatch.setattr(vmec_utils.eqdsk_base, "read_eqdsk",
                        lambda _: synthetic_eqdsk())

    converted = vmec_utils.eqdsk2vmec_gfile("synthetic.eqdsk")
    vmec_input = vmec_utils.default_vmec_input(converted)

    assert converted["rbc"].shape == (21, 1)
    assert vmec_input["mpol"] == 21
    assert vmec_input["ntheta"] >= 48
