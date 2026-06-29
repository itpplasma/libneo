# Third-party code and dependencies

libneo's own code is MIT-licensed (see `LICENSE`). The MIT text covers only
code written for this project. This file lists everything else: code vendored
into the tree, libraries linked into binaries, and tools used at build or test
time, each with its license and the obligations it carries.

## Vendored code

### `src/transport/gen_laguerre_rule.f90` (GNU LGPL)

John Burkardt's Fortran 90 version of IQPACK (Sylvan Elhay and Jaroslav
Kautsky, ACM TOMS Algorithm 655), distributed under the GNU LGPL per its
header. It is compiled into the `transport` static library. The LGPL, not the
MIT license, governs this file. Distributing a binary that contains `transport`
obliges the distributor to let recipients replace and relink the LGPL part
(LGPL section 4); distributing libneo in source form satisfies this.
Replacement with an MIT-licensed implementation is tracked in
[#288](https://github.com/itpplasma/libneo/issues/288). Keep the file's license
header intact until then.

### `python/scripts/ninja_syntax.py` (Apache-2.0)

Google's Ninja build-file writer, Apache License 2.0. Build tooling only; it is
not part of the installed Python wheel (`[tool.scikit-build.wheel.packages]`
ships `python/libneo` and `python/efit_to_boozer` only). Apache-2.0 is
MIT-compatible; the header notice must stay.

### `cmake/Modules/findFFTW/` (BSD-3-Clause)

CMake find module by Wenzel Jakob and Patrick Bos, BSD 3-Clause, with its own
`LICENSE` file in the directory. Build tooling only.

## Linked libraries

Licenses of libraries that end up inside libneo binaries or are loaded by them:

| Library | License | Linked into |
|---|---|---|
| FFTW3 | GPL-2.0-or-later | `magfie` |
| HDF5 | BSD-3-Clause-style | `hdf5_tools` |
| NetCDF-C, NetCDF-Fortran | BSD-3-Clause | `magfie`, field I/O |
| OpenBLAS / reference LAPACK | BSD-3-Clause / modified BSD | linear algebra throughout |
| Intel MKL (optional alternative) | Intel Simplified Software License | linear algebra, user-supplied |
| GCC runtimes (libgfortran, libgomp, libquadmath) | GPL-3.0 with GCC Runtime Library Exception | all gfortran binaries |

FFTW is the one entry that conflicts with MIT labeling: a conveyed binary
containing `magfie` combines MIT code with GPL FFTW, so the GPL terms govern
that distribution. Removal is tracked in
[#287](https://github.com/itpplasma/libneo/issues/287); the GSL link with the
same problem was removed in #286.

The GCC runtime libraries are GPL but their Runtime Library Exception
explicitly permits linking them into programs under any license; they impose no
copyleft obligation here.

## Test-only: GSL (GPL-3.0-or-later)

`test/collisions/test_incomplete_gamma_gsl_oracle.f90` checks the in-tree
incomplete gamma function against GSL. The test source is original libneo code
under MIT, which is GPL-compatible. CMake builds the test executable only when
it finds GSL (`find_package(GSL QUIET)`); the executable links GSL and is
therefore a combined work under GPL terms. GPL obligations attach when a binary
is conveyed to someone else, not when it is built or run. This test executable
lives in the build tree, is never installed, and is never distributed, so no
obligation arises, and nothing GSL-related reaches the installed libraries or
wheels. Do not install or package this executable; if that is ever wanted,
running GSL in a separate process instead of linking it removes the combined
work entirely.

## Python dependencies

The Python interface depends on packages installed as separate distributions
(numpy, scipy, matplotlib, netCDF4, h5py, shapely, trimesh, f90wrap, optionally
map2disc). They are not vendored; each is governed by its own license and none
of their code is copied into this repository. f90wrap is LGPL-3.0: it generates
the wrapper sources at build time and the generated Python modules import the
separate `f90wrap` package at runtime, which the LGPL permits without
conditions on libneo's license.
