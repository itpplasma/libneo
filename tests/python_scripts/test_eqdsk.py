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
    "../resources/input_efit_file.dat",

    # PROCESS code
    "/proj/plasma/DATA/DEMO/Equil_2021_PMI_QH_mode_betap_1d04_li_1d02_Ip_18d27MA_SOF.eqdsk",

    # Standardized
    "/proj/plasma/DATA/DEMO/Equil_2021_PMI_QH_mode_betap_1d04_li_1d02_Ip_18d27MA_SOF_std.eqdsk",

    # TODO: CHEASE
    # "/proj/plasma/DATA/DEMO/teams/Equilibrium_DEMO2019_CHEASE/MOD_Qprof_Test/EQDSK_DEMO2019_q1_COCOS_02.OUT"
]


def test_eqdsk_read():
    for test_file in test_files:
        _ = eqdsk.eqdsk_file(test_file)


def test_eqdsk_golden_records():
    for test_file in test_files:
        eqdsk_object = eqdsk.eqdsk_file(test_file)

        data = eqdsk_object.__dict__
        replace_array_members_by_lists(data)

        assert are_dicts_equal(data, get_golden_record(test_file))


def get_golden_record(filename):
    """
    Returns the golden record for the given test file.
    """

    file_base, _ = os.path.splitext(os.path.basename(filename))
    with gzip.open(f"{golden_record_dir}/{file_base}.json.gz", "rb") as f:
        compressed_data = f.read()
        json_data = compressed_data.decode("utf-8")
        return json.loads(json_data)


def store_golden_records():
    """
    Stores reference data for regression tests. Call only manually after fixing
    a bug and checking that the new results are correct.
    """

    for file in test_files:
        basedir = os.path.dirname(file)
        file_base, ext = os.path.splitext(os.path.basename(file))
        data = eqdsk.eqdsk_file(f"{basedir}/{file_base}{ext}").__dict__
        replace_array_members_by_lists(data)

        outfile = f"{golden_record_dir}/{file_base}.json.gz"

        if os.path.isfile(outfile):
            raise RuntimeError(f"Golden record file {outfile} already exists.")

        with gzip.open(outfile, "wb") as f:
            json_data = json.dumps(data, indent=4)
            json_bytes = json_data.encode("utf-8")
            f.write(json_bytes)


def replace_array_members_by_lists(data):
    """
    Replaces numpy arrays by lists in a dictionary. This is needed because
    numpy arrays are not JSON serializable.
    """
    import numpy as np

    for key, value in data.items():
        if isinstance(value, dict):
            replace_array_members_by_lists(value)
        elif isinstance(value, list):
            for i, item in enumerate(value):
                if isinstance(item, np.ndarray):
                    value[i] = item.tolist()
        elif isinstance(value, np.ndarray):
            data[key] = value.tolist()


def are_dicts_equal(dict1, dict2, tolerance=1e-9):
    """
    Compares two dictionaries. Returns True if they match, False otherwise.
    For float members, the comparison is done using math.isclose() with
    the given tolerance.
    """
    import math

    for key, value1 in dict1.items():
        if key not in dict2:
            return False
        value2 = dict2[key]
        if isinstance(value1, list) and isinstance(value2, list):
            if not are_float_lists_equal(value1, value2, tolerance):
                return False
        elif isinstance(value1, float) and isinstance(value2, float):
            if not math.isclose(value1, value2, rel_tol=tolerance):
                return False
        elif value1 != value2:
            return False
    return True


def are_float_lists_equal(list1, list2, tolerance=1e-9):
    """
    Compares two lists of floats. Returns True if they match, False otherwise.
    The comparison is done using math.isclose() with the given tolerance.
    """
    import math

    if len(list1) != len(list2):
        return False
    for a, b in zip(list1, list2):
        if isinstance(a, float) and isinstance(b, float):
            if not math.isclose(a, b, rel_tol=tolerance):
                return False
        elif isinstance(a, list) and isinstance(b, list):
            if not are_float_lists_equal(a, b, tolerance):
                return False
    return True
