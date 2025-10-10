# TODO

## ASCOT5 support
- [x] Diagnose `_magfie` import failure (`_f2pyinitspline_vmec_sub_` unresolved):
    - [x] Inspect `_magfie.cpython-*.so` link deps (`otool -L`) to confirm `libspline_vmec_sub` objects missing.
    - [x] Review `CMakeLists.txt`/f2py build recipe to ensure `spline_vmec_sub.f90` is included in extension sources.
    - [x] Rebuild `_magfie` after adjusting link order; validate import via `python - <<'PY' import _magfie; PY`.
    - [x] Add regression check ensuring `_magfie` imports under pytest (see `test_magfie_import`).
- [x] Audit f2py-exposed VMEC wrappers once the module loads:
    - [x] Confirm the key wrappers are exposed to Python (guarded by `test_magfie_import`).
    - [x] Document argument order and unit conventions (see `doc/ascot5.md`).
- [x] Implement and exercise VMEC/mgrid conversions and B₃DS writer:
    - [x] Conversion helpers live in `python/libneo/ascot5/__init__.py` with end-to-end pytest coverage.
- [x] Diagnostic PNGs are generated per-test and stored under `pytest_artifacts/` for review.
- [x] Documentation & developer notes:
    - [x] Added `doc/ascot5.md` for workflow, units, and troubleshooting.
    - [x] Extended README with an ASCOT5 helper summary and pointers to the docs.
