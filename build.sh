#!/bin/sh

if [ ! -d "build" ] ; then
    echo "Creating 'build' directory..."
    mkdir build
fi

cd build
cmake ..
make
cd -

if [ ! -d "tools/h5merge/build" ] ; then
    echo "Creating 'build' directory..."
    mkdir tools/h5merge/build
fi

cd tools/h5merge/build
cmake ..
make
cd -