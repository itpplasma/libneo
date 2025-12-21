# DESIGN.md - libneo Architecture and Design

## Overview
This document captures the architectural decisions, design patterns, and implementation strategies for the libneo Fortran library.

## Performance Optimizations

### Trampoline Elimination (Inner Subroutines)
**Problem**: Inner subroutines that access variables from their parent scope create trampolines, which cause:
- Performance overhead from indirect function calls
- Poor compiler optimization opportunities
- Security issues on some systems (executable stack required)
- Cache misses and instruction pipeline stalls

**Solution**: Refactor all inner subroutines to module-level procedures:
1. Move inner subroutines to module scope
2. Pass all required data explicitly as arguments
3. Use derived types to bundle related parameters when argument lists become long

**Implementation Strategy**:
1. **Phase 1 - Discovery**: Identify all inner subroutines accessing outer variables
2. **Phase 2 - Tests**: Ensure adequate tests exist before refactoring
3. **Phase 3 - Refactoring**: Move inner subroutines to module level
4. **Phase 4 - Validation**: Verify performance improvements and correctness

### Identified Cases

#### poincare.f90
- **Location**: `integrate_RZ_along_fieldline` contains inner subroutine `fieldline_derivative`
- **Issue**: `fieldline_derivative` accesses `field` parameter from outer scope
- **Solution**: Move to module level, pass `field` as argument or use context parameter
- **Tests**: Exists in `test/poincare/test_poincare.f90`
- **Priority**: HIGH - Used in field line integration (performance critical)

### Implementation Tasks

#### Task 1: Refactor poincare.f90
**Issue**: #119
**Approach**:
1. The inner subroutine `fieldline_derivative` is already using the context pattern correctly
2. It receives `field` through the `context` parameter in the interface
3. However, it's still an inner subroutine which can cause compiler issues
4. Solution: Keep using context pattern but move subroutine to module level

**Refactoring Steps**:
1. Move `fieldline_derivative` from inside `integrate_RZ_along_fieldline` to module level
2. Rename to `poincare_fieldline_derivative` to avoid naming conflicts
3. Update the call to `odeint_allroutines` to use the module-level procedure
4. Verify tests still pass

#### Task 2: Codebase Audit (COMPLETED)
**Issue**: #119
**Results**:
- Comprehensive scan completed: Only ONE instance found
- Location: `src/poincare.f90` - `fieldline_derivative` inner subroutine
- No other inner subroutines with contains blocks in the codebase
- Many files have multiple contains blocks but they are type definitions or module-level (both OK)

#### Task 3: Performance Benchmarking
**Issue**: #119
**Approach**:
1. Create benchmark for poincare integration before refactoring
2. Measure cache misses using perf tools
3. Compare performance after refactoring
4. Document improvements in issue

## Architecture Patterns

### Field Interface Pattern
The library uses an abstract field interface (`field_t`) with concrete implementations for different field types. This allows polymorphic behavior while maintaining performance.

### Context Parameter Pattern
For ODE integration and similar callbacks, we use a context parameter pattern that allows passing arbitrary data to callback functions without global state.

## Performance Guidelines

### Memory Layout
- Use column-major ordering for arrays (Fortran default)
- Prefer Structure of Arrays (SoA) over Array of Structures (AoS)
- Allocate large arrays as allocatable for heap storage

### Function Design
- Avoid inner subroutines accessing outer variables
- Use `pure` and `elemental` where possible
- Pass data explicitly rather than through module variables
- Keep functions under 100 lines (target 50)

## Testing Strategy

### Testing Requirements
- All refactored code must have tests before modification
- Tests must verify correctness and key performance characteristics where applicable
- Use the existing test framework in `test/`

## Build System
- Primary: CMake with Ninja
- Secondary: fpm (Fortran Package Manager)
- All changes must pass CI pipeline

## GPU Offload (OpenACC / OpenMP target)

### Goal
Support GPU-ready inner loops without forcing a specific backend on downstream codes,
and without accidental device-host copies inside tight iteration loops.

### Current production shape
The batch spline layer now includes many-point evaluation APIs alongside the existing
single-point batch-over-quantities routines:

- 1D: `evaluate_batch_splines_1d_many(spl, x(:), y(:, :))`
- 2D: `evaluate_batch_splines_2d_many(spl, x(2, :), y(:, :))`
- 3D: `evaluate_batch_splines_3d_many(spl, x(3, :), y(:, :))`

For each dimension there is also a resident variant intended for use inside a
downstream-managed device data region:

- 1D: `evaluate_batch_splines_1d_many_resident`
- 2D: `evaluate_batch_splines_2d_many_resident`
- 3D: `evaluate_batch_splines_3d_many_resident`

The resident variants contain `!$acc` kernels and use `present(...)` so they can
run without implicit transfers when the caller has already placed the arrays on
device.

### Construction and data residency
The default `construct_batch_splines_{1,2,3}d(...)` entry points remain host-safe and
produce coefficients that match the existing scalar spline construction used in the
test suite. When built with OpenACC enabled, these constructors also place the
coefficient arrays on the OpenACC device via `!$acc enter data copyin(...)` so the
subsequent many-point evaluation wrappers can offload without additional setup.

For build-on-device performance experiments and for downstream code that wants to
keep construction on the GPU, the explicit device constructor entry points remain
available:

- `construct_batch_splines_2d_resident_device`
- `construct_batch_splines_3d_resident_device`

### How downstream should use OpenACC with no copies
Downstream must keep the data resident across the whole accelerated algorithm:

1. Enter a data region that places the coefficient arrays and inputs on the device.
2. Call the resident routine inside the region.
3. Only update back to host at the boundary where host-side code needs the results.

The same compiler and OpenACC runtime must be used for both libneo and downstream.
For Fortran, this effectively means compiling both projects with the same compiler
because module files are not ABI compatible across compilers.

### OpenACC directive pitfalls
- If a kernel uses `present(...)`, calling it without an active device data region
  will fail at runtime.
- Prefer explicit resident variants and keep host-safe variants separate to avoid
  silent per-call transfers.
- Avoid implicit copies in hot loops: create persistent data regions in the
  downstream algorithm and keep arrays present across iterations.

### OpenMP target offload lessons
OpenMP target offload can be competitive, but libgomp NVPTX may clamp launch geometry
based on kernel metadata and resource usage, which can reduce occupancy and hurt
performance versus OpenACC for the same math kernel. Control data lifetime with
explicit `target enter data` / `target exit data` and avoid implicit maps inside
tight loops.

## Code Standards
See `CODING_STANDARD.md` for detailed coding standards.
