# libneo MATLAB library

The `matlab/` directory contains MATLAB classes, functions and code interfaces.
Instructions for running the balance workflow are in
[matlab/README.md](../matlab/README.md).

- `EFIT`: read, write and modify equilibrium g-files (G-EQDSK, often called
  EFIT files).
- `GPEC_interface`: minimal interface to run GPEC together with DCON.
- `InputFile`: read, create and modify input files containing Fortran
  namelists.
- `KIN`: write `.kin` files (kinetic profiles stored column-wise).
- `KiLCA_interface`: prepare, run and post-process KiLCA runs. Requires a
  compiled [KiLCA](https://github.com/itpplasma/KiLCA).
- `Shared`: classes and functions shared by the other parts.
- `Utility`: general-purpose classes and functions.
