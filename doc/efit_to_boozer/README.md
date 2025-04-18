This program converts a geqdsk file into a Boozer file in `.bc` text format.

# Installation

## Prerequisites
`efit_to_boozer` requires a checked-out copy of https://github.com/itpplasma/libneo on the same directory level.

## Main executable
Run

    make

This will yield `build/efit_to_boozer.x`

## Python package
    
    pip install -e .

This places an _efit_to_boozer compile module in `venv/lib/python../site-packages`


# Usage
Run

    build/efit_to_boozer.x

Input files:

    efit_to_boozer.inp
    field_divB0.inp
    convexwall.dat

The file 'efit_to_boozer.inp' contains calculation parameters.
The file 'field_divB0.inp' contains the name of the input-data files and
some switches regarding the calculation.

The file `convexwall.dat` describes the wall of the device. So as long as
the geometry remains the same, this needs not to be switched.

Output files:

    box_size_axis.dat
    btor_rbig.dat
    flux_functions.dat
    fromefit.bc
    out.06

The file `fromefit.bc` will contain the boozer data.
