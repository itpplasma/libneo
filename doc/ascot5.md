# ASCOT5 Field Export Support

The ASCOT5 helpers live in `python/libneo/ascot5/__init__.py`. They bridge the
VMEC/`_magfie` Fortran wrappers and ASCOT5's `B_3DS` magnetic-field files.

## Units

- The underlying VMEC Fortran routines continue to operate in CGS
  (centimetres/Gauss).  Python converts to SI only at the final export stage.
- Radii and heights in the returned `B3DSField` dataclass are expressed in
  metres; magnetic-field components in the in-memory object are still in Gauss
  until `write_b3ds_hdf5` multiplies by `GAUSS_TO_TESLA`.

## Workflow

1. Load VMEC metadata via `_magfie.f2py_vmec_wrappers.splint_vmec_data_wrapper`.
2. Sample the last closed flux surface to establish a cylindrical bounding box.
3. Walk the structured grid, solving for `(s, \vartheta)` and evaluating
   `_magfie.f2py_vmec_wrappers.vmec_field_cylindrical_wrapper`.
4. Store the field in a `B3DSField` dataclass and, if needed, write it to HDF5
   using `write_b3ds_hdf5`.

`field_from_mgrid` mirrors this process for ASCOT5 mgrid NetCDF inputs.

## Testing

- `pytest -q` exercises VMEC and mgrid conversion via
  `test/python/test_ascot5_vmec.py` and `test/python/test_ascot5_mgrid.py`.  The
  tests also leave diagnostic PNGs in their respective temporary directories. If
  the environment variable `PYTEST_ARTIFACTS` is set, the tests copy both the
  generated HDF5 and PNG files into that directory for easy inspection.
- `test/python/test_magfie_import.py` provides a regression guard that the
  `_magfie` extension and its key VMEC wrappers remain importable.

## Troubleshooting

If `_magfie` fails to import, rebuild the shared object by running `make` (or
  `cmake --build --preset default`).  Confirm that the active Python interpreter
  matches the virtual environment expected by CMake; `bash -lc 'which python'`
  should resolve to the project `/.venv` directory.

If the exported field appears hollow it usually means the Newton solver could
not converge for large portions of the padded R/Z bounding box.  Increase
`nr`/`nz`/`nphi`, tighten the padding, or raise `max_iter`/`tol` when calling
`field_from_vmec` to obtain a denser valid region.  Points very close to the
magnetic axis fall back to an on-axis field evaluation so the B field remains
finite there.
