# spectre_test.h5

SPECTRE output used as fixture by `test_spectre_reader`.

- Generator: SPECTRE (https://gitlab.com/spectre-eq/spectre), commit
  `f51638520d2dbdee9b0b4e4174ace5b2349b7460`.
- Case: `G3V3L3Fi` test input (Igeometry=3, Nfp=5, Nvol=3, Mpol=4, Ntor=2,
  Lrad=8,4,4; stellarator symmetric).
- Generation scripts: https://gitlab.tugraz.at/plasma/proj/spectre-orbits.

Reference values hardcoded in `test/spectre/test_spectre_reader.f90` were
extracted from this file with h5py at double precision.
