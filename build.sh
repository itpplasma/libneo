#!/bin/sh

if [ ! -d "build" ] ; then
    echo "Creating 'build' directory..."
    mkdir build
fi

cd build
cmake ..
make
cd -
