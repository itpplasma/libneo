# libneo
Common code for plasma codes of ITPcp, in particular for different NEO-2 versions.

## Getting started

### Prerequisites
`libneo` requires CMake, Ninja, GCC+GFortran, OpenMPI, 
a BLAS/LAPACK (OpenBLAS or MKL) and FFTW, NetCDF+NetCDF-Fortran 
and HDF5 including development headers.

For your own Debian or Ubuntu system, run `setup/debian.sh`.
For usage on common computing clusters, load the modules in `setup/<machine>.sh`.

### Build
For convenience, the build process can be automatically started by running

    make

directly in the `libneo` directory. This will create a `build` subdirectory
and run `cmake` with `ninja` internally. 

### Install Python interface
The Python interface is located in the `python` subdirectory. This interface
is only built if `python`, `numpy` and `f90wrap` are available. Please activate
a suitable `venv` first, and then run

    pip install -e .

in the `libneo` directory for an editable build.

## doc
Documentation, this includes user documentation and interface
documentation. The latter is generated via doxygen.

## matlab
Matlab library that contains various functions, classes and code interfaces for general and specific uses.

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

## python
Python interface package and scripts/functions.

### Coils format converter
Convert STELLOPT coils format to simple biotsavart format:

    python -m libneo.convert_coils_to_simple input.coils output.simple

The simple format is compatible with libneo's `neo_biotsavart` module and SIMPLE code.

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

### poincare
Poincare plot generation for magnetic field lines. Computes field line trajectories 
and their intersections with toroidal cross-sections to visualize magnetic islands, 
chaos regions, and flux surface structure.

## tests
Contains resources for tests of the library.
