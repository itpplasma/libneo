#!/bin/bash

cmake --preset default
cmake --build --preset default

if [ ! -d "tools/h5merge/build" ] ; then
    echo "Creating 'build' directory..."
    mkdir -p tools/h5merge/build
fi

pushd tools/h5merge/build
cmake ..
make
popd
