"""
Tests for eqdsk.py
"""

import pytest

import json
import gzip
import os.path

from libneo import eqdsk

# Directory where the golden records are stored for regression tests
golden_record_dir = "/proj/plasma/DATA/TESTS/libneo/eqdsk"


test_files = [
    # Local
    "tests/resources/input_efit_file.dat",
]


def test_eqdsk_read():
    for test_file in test_files:
        _ = eqdsk.eqdsk_file(test_file)