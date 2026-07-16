"""Tests for the plain EQDSK reader/writer in libneo.eqdsk_base."""

import numpy as np
import pytest

from libneo.eqdsk_base import read_eqdsk, write_eqdsk


def _synthetic_eqdata(nrgr: int, nzgr: int) -> dict:
    boundary = np.column_stack(
        [2.0 + 0.5*np.cos(np.linspace(0.0, 2.0*np.pi, 33)),
         0.5*np.sin(np.linspace(0.0, 2.0*np.pi, 33))]
    )
    return {
        "header": "synthetic eqdsk for header format test          ",
        "nrgr": nrgr,
        "nzgr": nzgr,
        "rboxlength": 1.6,
        "zboxlength": 1.6,
        "R0": 2.0,
        "rboxleft": 1.2,
        "zboxmid": 0.0,
        "Rpsi0": 2.0,
        "Zpsi0": 0.0,
        "PsiaxisVs": 0.0,
        "PsiedgeVs": 0.15,
        "Btor_at_R0": 2.0,
        "Ip": 5.0e5,
        "fprof": np.full(nrgr, 4.0),
        "ptotprof": np.linspace(5.0e3, 0.0, nrgr),
        "fdfdpsiprof": np.zeros(nrgr),
        "dpressdpsiprof": np.full(nrgr, -3.0e4),
        "PsiVs": np.linspace(0.0, 0.15, nrgr*nzgr).reshape(nzgr, nrgr),
        "qprof": np.linspace(1.0, 4.0, nrgr),
        "npbound": boundary.shape[0],
        "nplimiter": boundary.shape[0],
        "Lcfs": boundary,
        "Limiter": 1.1*boundary,
    }


@pytest.mark.parametrize("nrgr,nzgr", [(65, 65), (129, 129), (257, 129)])
def test_write_eqdsk_header_matches_3i4_columns(tmp_path, nrgr, nzgr):
    """The g-file header ends in (3i4) integer fields at columns 49-60.

    Fortran readers (libneo read_eqfile1/read_eqfile2 and EFIT-compatible
    codes) parse the grid dimensions with format (6a8,3i4), so the writer
    must place idum, nw, nh in fixed four-character columns. Free-form
    spacing breaks every grid dimension with three or more digits.
    """
    path = tmp_path / "test.eqdsk"
    write_eqdsk(path, _synthetic_eqdata(nrgr, nzgr))

    header = path.read_text().splitlines()[0]
    assert int(header[48:52]) == 0
    assert int(header[52:56]) == nrgr
    assert int(header[56:60]) == nzgr


@pytest.mark.parametrize("nrgr,nzgr", [(65, 65), (129, 129)])
def test_write_eqdsk_read_eqdsk_roundtrip(tmp_path, nrgr, nzgr):
    path = tmp_path / "test.eqdsk"
    eqdata = _synthetic_eqdata(nrgr, nzgr)
    write_eqdsk(path, eqdata)
    back = read_eqdsk(path)

    assert back["nrgr"] == nrgr
    assert back["nzgr"] == nzgr
    np.testing.assert_allclose(back["PsiVs"], eqdata["PsiVs"],
                               atol=1e-6*eqdata["PsiedgeVs"])
    np.testing.assert_allclose(back["qprof"], eqdata["qprof"], atol=1e-6)
