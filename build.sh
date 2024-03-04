#!/bin/sh

cmake --preset default
cmake --build --preset default

if [ ! -d "tools/h5merge/build" ] ; then
    echo "Creating 'build' directory..."
    mkdir tools/h5merge/build
fi

cd tools/h5merge/build
cmake ..
make
cd -
