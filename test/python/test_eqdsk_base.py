"""
Tests for eqdsk.py
"""

import pytest

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


if __name__ == "main":
