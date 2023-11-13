"""
Tests for eqdsk.py
"""

import pytest

from libneo import eqdsk

@pytest.mark.skip(reason="TODO: fix. Fails at q profile.")
def test_eqdsk_read():
    e = eqdsk.eqdsk_file('../resources/input_efit_file.dat')


def test_eqdsk_read_demo():
    e = eqdsk.eqdsk_file('/proj/plasma/DATA/DEMO/Equil_2021_PMI_QH_mode_betap_1d04_li_1d02_Ip_18d27MA_SOF.eqdsk')
