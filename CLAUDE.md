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
- CMake 3.18+, Ninja, GCC/GFortran
- MPI (OpenMPI)
- BLAS/LAPACK (OpenBLAS or MKL)
- FFTW3, NetCDF, HDF5
- Python with numpy and f90wrap (for Python interface)

### Testing Structure
Tests are organized by component in `test/`:
- Unit tests for each module
- Integration tests for field computations
- Python tests for bindings
- Test utilities in `test/util_for_test/`

**Remember**: Always write tests BEFORE implementation code (TDD)