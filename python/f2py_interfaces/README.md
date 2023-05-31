The package is now installed via

    cd python/f2py_interfaces/
    pip[3] install -e .

and imported with

    from libneo import magfie

Be sure that you have the build directory directly under libneo, and add it to LD_LIBRARY_PATH via

    export LD_LIBRARY_PATH=/path/to/libneo/build:$LD_LIBRARY_PATH

If you used another name for the build folder you have to adapt the
above line and setup.py.
