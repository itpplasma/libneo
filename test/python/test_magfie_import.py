"""Regression tests for the f2py-generated VMEC wrappers."""

import importlib


def test_magfie_module_exposes_vmec_wrappers():
    mod = importlib.import_module("_magfie")

    wrappers = getattr(mod, "f2py_vmec_wrappers", None)
    assert wrappers is not None, "_magfie must expose f2py_vmec_wrappers attribute"

    expected = (
        "splint_vmec_data_wrapper",
        "vmec_field_wrapper",
        "vmec_field_cylindrical_wrapper",
    )
    for name in expected:
        assert hasattr(wrappers, name), f"VMEC wrapper '{name}' missing"
