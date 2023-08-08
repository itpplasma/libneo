#!/bin/sh

# Install dependencies on Debian or Ubuntu systems
apt-get update && apt-get install -y wget git cmake make python3 \
    gcc g++ gfortran libnetcdf-dev libnetcdff-dev openmpi-bin libopenmpi-dev
