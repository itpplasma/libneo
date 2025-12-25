# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

libneo is a Fortran library containing common code for plasma physics codes at ITPcp, particularly for NEO-2 versions. It provides magnetic field representations, transport calculations, collision operators, and coordinate transformations for fusion plasma physics.

## Development Approach

### Test-Driven Development (STRICT REQUIREMENT)
**ALWAYS follow the RED-GREEN-REFACTOR cycle:**

1. **RED**: Write a failing test first - never write code without a failing test
2. **GREEN**: Write minimal code to make the test pass - only enough to pass
3. **REFACTOR**: Clean up code while keeping tests green

**Rules:**
- Write tests before any implementation code
- Run tests after each step
- Keep the cycle short and focused

### Coding Standards
See `CODING_STANDARD.md` for complete details. Key points:
- Code is simple and self-explaining without comments
- Break large constructs into small, single-responsibility units
- Use `pure` procedures without side effects when possible
- Prefer structs (`type`) over classes for data bundling
- Use abstract types with factory pattern for polymorphism
- Arrays indexed starting from 1
- 88-character line limit, 4-space indentation
- `implicit none` in all modules and programs

## Build Commands

### Standard Build
```bash
# Simple build with CMake/Ninja (recommended)
make

# Alternative: using fpm
fpm build
```

### Run Tests
```bash
# Run all tests
make test

# Or directly with ctest
cd build && ctest
```

### Install Python Interface
```bash
# Activate your virtual environment first
pip install -e .

# For chartmap tests (requires map2disc from MPCDF GitLab):
pip install -e ".[chartmap]"

# For full development (includes map2disc and other tools):
pip install -e ".[dev]"
```

### Clean Build
```bash
make clean
```

## Architecture Overview

### Core Structure
The library follows an object-oriented Fortran design with abstract types and interfaces:

- **Abstract Field Interface**: `field_t` in `src/field/field_base.f90` defines the interface for all magnetic field implementations
- **Modular Components**: Each major functionality (transport, collisions, species) is in its own module
- **Python Bindings**: F2PY and f90wrap interfaces in `src/f2py_interfaces/` and `python/`

### Key Components

1. **Magnetic Fields** (`src/field/`): Various field representations (VMEC, EFIT, mesh-based)
2. **MAGFIE** (`src/magfie/`): Magnetic field interpolation engine
3. **Transport** (`src/transport/`): Neoclassical transport coefficient calculations
4. **Collisions** (`src/collisions/`): Collision frequency calculations
5. **Coordinate Transformations** (`src/efit_to_boozer/`): EFIT to Boozer coordinate conversion
6. **HDF5 Tools** (`src/hdf5_tools/`): Simplified HDF5 I/O interface

### Dependencies
- CMake 3.18+, Ninja, GCC/GFortran or NVHPC (nvfortran)
- MPI (OpenMPI for gfortran, HPC-X for nvfortran)
- BLAS/LAPACK (OpenBLAS or MKL)
- FFTW3, NetCDF, HDF5
- Python with numpy and f90wrap (for Python interface)
- CUDA toolkit (for GPU offloading with OpenACC)

### GCC16 Setup (faepcr* machines)

On ITPcp cluster machines (hostname pattern `faepcr*`), GCC 16 with nvptx OpenACC offloading is available at `/temp/AG-plasma/opt/gcc16`. Use it for GPU-accelerated builds:

```bash
cmake -S . -B build_gcc16 -G Ninja \
  -DCMAKE_Fortran_COMPILER=/temp/AG-plasma/opt/gcc16/bin/gfortran \
  -DCMAKE_C_COMPILER=/temp/AG-plasma/opt/gcc16/bin/gcc \
  -DCMAKE_BUILD_TYPE=Release
cmake --build build_gcc16 -j
LD_LIBRARY_PATH=/temp/AG-plasma/opt/gcc16/lib64 ctest --test-dir build_gcc16
```

### NVHPC Compiler Setup

When using NVHPC (nvfortran) with OpenACC GPU offloading:

**Quick start:**
```bash
make nvfortran        # Configure and build with nvfortran
make test-nvfortran   # Run tests
```

**Manual setup (if make target does not work):**

1. **Source HPC-X environment** (REQUIRED - the stub wrapper in comm_libs/hpcx/bin does NOT work):
   ```bash
   # Use the FULL HPC-X from comm_libs/13.0 (or similar version subdir)
   source /opt/nvidia/hpc_sdk/Linux_x86_64/VERSION/comm_libs/13.0/hpcx/latest/hpcx-mt-init.sh
   hpcx_load
   ```

2. **Set CUDA home** (auto-detected if CUDA is at /opt/cuda or /usr/local/cuda):
   ```bash
   export NVHPC_CUDA_HOME=/opt/cuda
   ```

3. **Configure with CMake:**
   ```bash
   cmake -S . -B build_nvfortran -G Ninja \
     -DCMAKE_Fortran_COMPILER=/opt/nvidia/hpc_sdk/Linux_x86_64/VERSION/comm_libs/13.0/hpcx/latest/ompi/bin/mpifort \
     -DCMAKE_C_COMPILER=/opt/nvidia/hpc_sdk/Linux_x86_64/VERSION/compilers/bin/nvc \
     -DCMAKE_BUILD_TYPE=Release \
     -DFORCE_FETCH_DEPS=ON
   ```

4. **Build:**
   ```bash
   cmake --build build_nvfortran -j
   ```

**Important notes:**
- NVHPC auto-fetches HDF5, NetCDF, FFTW from source (ABI incompatibility with system libs)
- The stub mpifort at `comm_libs/hpcx/bin/mpifort` has strict CUDA version checks that may fail
- Use the full HPC-X installation at `comm_libs/13.0/hpcx/latest/ompi/bin/mpifort` instead
- CMake auto-detects CUDA at `/opt/cuda` or `/usr/local/cuda` via `cmake/SetupNVHPCCuda.cmake`

### Testing Structure
Tests are organized by component in `test/`:
- Unit tests for each module
- Integration tests for field computations
- Python tests for bindings
- Test utilities in `test/util_for_test/`

**Remember**: Always write tests BEFORE implementation code (TDD)
- You run tests with `make test`