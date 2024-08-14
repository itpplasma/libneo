"""
Tests for eqdsk.py
"""

import pytest

import json
import gzip
import os.path

from libneo import eqdsk

test_files = [
    # Local
    "test/resources/input_efit_file.dat",
]


def test_eqdsk_read():
    for test_file in test_files:
        _ = eqdsk.eqdsk_file(test_file)
