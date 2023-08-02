#!/bin/sh

mkdir build
cd build
cmake ..
make -j24
cd -
