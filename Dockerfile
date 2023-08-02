FROM debian:bookworm-slim

RUN apt-get update && apt-get install -y \
    wget \
    git \
    cmake \
    g++ \
    make \
    python3 \
    gfortran \
    libnetcdf-dev \
    libnetcdff-dev