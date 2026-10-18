"""
Tests for eqdsk.py
"""

import pytest
from numpy.testing import assert_allclose

import json
import gzip
import os.path

from libneo import eqdsk_base

test_files = [
    # Local
    "test/resources/input_efit_file.dat",
]


def test_read_eqdsk():
    for test_file in test_files:
        _ = eqdsk_base.read_eqdsk(test_file)

def test_write_eqdsk(tmp_path):
    for test_file in test_files:
        tmp_file = (tmp_path / "output_efit_file.dat").as_posix()
        eqdsk = eqdsk_base.read_eqdsk(test_file)
        eqdsk_base.write_eqdsk(tmp_file, eqdsk)
        eqdsk2 = eqdsk_base.read_eqdsk(tmp_file)
        for key in eqdsk:
            # Check if string
            if isinstance(eqdsk[key], str):
                assert eqdsk[key] == eqdsk2[key]
            else:
                assert_allclose(eqdsk[key], eqdsk2[key])


def test_read_short_header(tmp_path):
    source = test_files[0]
    expected = eqdsk_base.read_eqdsk(source)
    with open(source) as stream:
        stream.readline()
        body = stream.read()
    target = tmp_path / "short-header.geqdsk"
    target.write_text(f"EFIT 3 {expected['nrgr']} {expected['nzgr']}\n" + body)
    actual = eqdsk_base.read_eqdsk(target)
    for key in ["nrgr", "nzgr", "PsiVs", "fprof", "Lcfs", "Limiter"]:
        assert_allclose(actual[key], expected[key])

if __name__ == "__main__":
    pytest.main([__file__, "-s"])
