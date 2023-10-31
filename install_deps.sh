#!/bin/sh

# Install dependencies on Debian or Ubuntu systems
apt-get update && apt-get install -y wget git cmake make gcc g++ gfortran \
    libnetcdf-dev libnetcdff-dev openmpi-bin libopenmpi-dev libfftw3-dev \
    python3 python3-numpy python3-scipy
