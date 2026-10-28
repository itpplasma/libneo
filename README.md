# libneo

libneo is the shared Fortran library and Python package of the plasma theory
group at ITPcp, TU Graz. It collects equilibrium and magnetic-field readers,
coordinate conversions and numerical utilities used by NEO-2, NEO-RT, SIMPLE,
KAMEL, MEPHIT and related codes.

Status: research infrastructure under active development; interfaces may
change.

## Contents

- Magnetic-field representations and readers: G-EQDSK, VMEC `wout`, Boozer
  (`.bc`, `boozmn`), chartmap, mgrid, SPECTRE, JOREK, coil files and
  Biot-Savart fields (`src/field`, `src/magfie`).
- Coordinate systems and conversions, including VMEC, chartmap and geometric
  flux coordinates (`src/coordinates`).
- Converters from EFIT/G-EQDSK and Hamada coordinates to Boozer coordinates
  (`efit_to_boozer.x`, `hamada_to_boozer.x`) and from VMEC to chartmaps
  (`vmec_to_chartmap.x`).
- Splines and batch interpolation, ODE integration, collision frequencies,
  Poincare plots of field lines.
- HDF5 and NetCDF access through [Fortio](https://github.com/lazy-fortran/fortio)
  (`hdf5_tools`, `nctools`), without linking the system HDF5 or NetCDF libraries.
- Python package `libneo` with chartmap generators, VMEC/NEMEC/EQDSK tools,
  coil conversion and ASCOT5 field export; MATLAB classes in `matlab/`.

## Building

Requirements: CMake 3.22 or newer, Ninja, GCC and gfortran, BLAS/LAPACK.
CMake fetches Fortio and [Fortnum](https://github.com/lazy-fortran/fortnum) at
configure time (revisions pinned in `fpm.toml`). On Debian or Ubuntu,
`setup/debian.sh` installs the system packages; for computing clusters, load the
modules in `setup/<machine>.sh`.

    make          # configure and build in build/ with CMake and Ninja
    make test     # run the CTest suite

Codes that use libneo usually fetch it through CMake; NEO-2, NEO-RT and KAMEL
accept `LIBNEO_REF=<branch|tag|sha>` to select a revision.

## Python package

Inside a virtual environment, from the repository root:

    pip install -e .

The Fortran wrappers are built only if `numpy` and `f90wrap` are available.
Optional extras: `.[chartmap]` for the chartmap generators, `.[dev]` for tests.
Example:

    libneo-write-chartmap from-vmec wout.nc wout.chartmap.nc --nrho 33 --ntheta 65 --nzeta 33

## Documentation

- [doc/python.md](doc/python.md): Python package, chartmaps, VMEC and NEMEC
  conventions, coil and ASCOT5 tools.
- [doc/efit_to_boozer/](doc/efit_to_boozer/README.md): EFIT to Boozer
  converter.
- [doc/Poincare.md](doc/Poincare.md): Poincare plot module.
- [doc/ascot5.md](doc/ascot5.md): ASCOT5 field export.
- [doc/matlab.md](doc/matlab.md): MATLAB library.
- `doc/UserDoc/`: LaTeX user documentation of older modules.
- `docs/`: Sphinx documentation of the Python interface, built in CI.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for dependency pins and the merge
workflow, and [CODING_STANDARD.md](CODING_STANDARD.md) for conventions.

## License

Original libneo code is MIT-licensed, see [LICENSE](LICENSE). Vendored
third-party code and linked libraries carry their own licenses;
[THIRD_PARTY.md](THIRD_PARTY.md) lists them.
