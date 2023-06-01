# libneo Python interface and scripts

## Getting started
First you require a build of `liblibneo.so` in some build directory. Follow
the `README.md` in the top-level `libneo` directory for this.

The Python package is installed via

    LIBNEO_BUILD_PATH=/path/to/libneo/build pip install -e .

(you may need to use `pip3` instead of `pip` on certain systems)

Currently the interface supports the magnetic field module in its VMEC variant.
It is imported with

    from libneo import magfie

Be sure that you have the build directory directly under libneo, and add it to LD_LIBRARY_PATH via

    export LD_LIBRARY_PATH=/path/to/libneo/build:$LD_LIBRARY_PATH

## Subdirectories
### libneo
This directory just contains the package `__init__.py` file.

### scripts
Here, ad-hoc scripts, notebooks and functions are located together with
some example inputs. Needs cleanup.
