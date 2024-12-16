# libneo
Common code for plasma codes of ITPcp, in particular for different Neo-2 versions.


## Getting started

### Prerequisites
`libneo` requires the following libraries to be available with
development headers: GSL, FGSL, BLAS, LAPACK, SuiteSparse and SuperLU.

For your own Debian or Ubuntu system, run `setup/debian.sh`

### Build
`libneo` is built using `cmake`. The build process is as follows:

    mkdir build
    cd build
    cmake ..
    make

### Install Python interface
The Python interface is located in the `python` subdirectory. Follow the
`README.md` in that directory for installation instructions.

## doc
Documentation, this includes user documentation and interface
documentation. The latter is generated via doxygen.

## matlab
Matlab library that contains various functions, classes and code interfaces for general and specific uses.

### BALANCE
Interface for the balance code. Contains Matlab code that is used to run the code, make pre- and post processing. Contains Fortran source files that are used in the interface.

### EFIT
A class to read/write/modify equilibrium g-files (commonly wrongly called efit files).

### GPEC_interface
A minimalist interface to run the code GPEC together with DCON.

### InputFile
A class that can be used to read/create/modify inputfiles that contain Fortran Namelists.

### KIN
A class that can write .kin files (kinetic profiles saved column wise)

### KiLCA_interface
Interface to prepare/run/postprocess the code KiLCA. A compiled version of KiLCA is needed (see [here](https://github.com/itpplasma/KiLCA))

### Shared
Shared classes and functions of the Matlab library.

### Utility
Classes and functions that have a general use for many application.

### mnDAT
Class to write mnDAT files.

## python
Python interface package and scripts/functions.

## src
Fortran source files for the library.
Subfolders contain source for additional tools, also compiled into
libneo.

### hdf5_tools
Interface to hdf5, to simplify calls.

### magfie

### MC

### contrib

### MyMPILib
Interface to MPI, so no actual mpi routines need to be called in
programs.

## tests
Contains resources for tests of the library.
