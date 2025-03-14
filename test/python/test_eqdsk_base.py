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

if __name__ == "__main__":
    pytest.main([__file__, "-s"])
