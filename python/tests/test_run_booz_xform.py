import sys
from types import SimpleNamespace

import netCDF4
import numpy as np
import pytest

from libneo.run_booz_xform import run_booz_xform


class FakeBoozXform:
    latest = None

    def __init__(self):
        self.ns_in = 5
        self.compute_surfs = np.arange(self.ns_in)
        self.mboz = 0
        self.nboz = 0
        self.ran = False
        FakeBoozXform.latest = self

    def read_wout(self, filename, flux=False):
        self.input = filename
        self.flux = flux

    def run(self):
        self.ran = True

    def write_boozmn(self, filename):
        with netCDF4.Dataset(filename, "w"):
            pass


def install_fake(monkeypatch):
    monkeypatch.setitem(
        sys.modules,
        "booz_xform",
        SimpleNamespace(Booz_xform=FakeBoozXform, __version__="test-1.0"),
    )


def test_run_booz_xform_sets_resolution_surfaces_and_provenance(
    tmp_path, monkeypatch
):
    install_fake(monkeypatch)
    wout = tmp_path / "wout_case.nc"
    with netCDF4.Dataset(wout, "w") as dataset:
        aminor = dataset.createVariable("Aminor_p", "f8")
        aminor.assignValue(2.6)
        rmajor = dataset.createVariable("Rmajor_p", "f8")
        rmajor.assignValue(6.1)
    output = tmp_path / "boozmn_case.nc"

    result = run_booz_xform(
        wout, output, mboz=64, nboz=0, surfaces=[0, 2, 4]
    )

    transform = FakeBoozXform.latest
    assert result == output
    assert transform.mboz == 64
    assert transform.nboz == 0
    assert transform.flux is True
    assert np.array_equal(transform.compute_surfs, [0, 2, 4])
    assert transform.ran
    with netCDF4.Dataset(output) as ds:
        assert ds.libneo_booz_xform_source == str(wout.resolve())
        assert len(ds.libneo_booz_xform_source_sha256) == 64
        assert ds.libneo_booz_xform_version == "test-1.0"
        assert ds.libneo_booz_xform_mboz == 64
        assert ds.libneo_booz_xform_nboz == 0
        assert ds.libneo_booz_xform_surfaces == "0,2,4"
        assert ds.libneo_vmec_aminor_m == pytest.approx(2.6)
        assert ds.libneo_vmec_rmajor_m == pytest.approx(6.1)


@pytest.mark.parametrize(
    "kwargs,match",
    [
        ({"mboz": 0}, "mboz"),
        ({"nboz": -1}, "nboz"),
        ({"surfaces": []}, "surfaces"),
        ({"surfaces": [5]}, "outside"),
    ],
)
def test_run_booz_xform_rejects_invalid_controls(
    tmp_path, monkeypatch, kwargs, match
):
    install_fake(monkeypatch)
    wout = tmp_path / "wout_case.nc"
    with netCDF4.Dataset(wout, "w") as dataset:
        aminor = dataset.createVariable("Aminor_p", "f8")
        aminor.assignValue(0.46)
        rmajor = dataset.createVariable("Rmajor_p", "f8")
        rmajor.assignValue(1.64)

    with pytest.raises(ValueError, match=match):
        run_booz_xform(wout, tmp_path / "boozmn.nc", **kwargs)
