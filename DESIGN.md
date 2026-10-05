# DESIGN.md - libneo Architecture and Design

## Overview
This document captures the architectural decisions, design patterns, and implementation strategies for the libneo Fortran library.

### Ecosystem boundary: interchange rather than equilibrium ownership

libneo owns generic plasma-code interchange, coordinate/convention conversion,
reusable magnetic-field/geometry utilities and compatibility adapters that have
multiple consumers. It is intentionally below KIN6D and TIAGO in the
application stack.

- **KIN6D** owns general plasma equations, stationary/evolution solves,
  differentiability, numerical certification and model-reduction error.
- **TIAGO** owns diagnostic observation/inference, likelihoods, nuisance
  parameters, priors and posterior UQ.
- **libneo** supplies shared readers/writers/converters and mature field
  utilities to both.

EQDSK/gEQDSK, VMEC wout, boozmn/chartmap, SPECTRE/SPEC, JOREK and similar
formats are adapter inputs. A format-specific coordinate model must not leak
into a consumer as a universal physical assumption. In particular, VMEC input
does not imply that KIN6D requires nested flux surfaces.

When KIN6D re-solves or certifies an imported state, the equilibrium residual,
solver and certificate remain KIN6D responsibilities. When TIAGO reconstructs
from measurements, the likelihood/posterior remain TIAGO responsibilities.
Avoid adding wrappers here that duplicate either layer.

New adapters should retain source format/version, units, orientation,
toroidal-angle/field-period/flux conventions and any transformations needed to
obtain the normalized physical representation. Prefer independently checked
round trips or cross-code oracles for convention-sensitive conversions.

### Boy-scout modernization driven by real consumers

libneo is historically grown and is not scheduled for a wholesale rewrite.
Instead, every new production consumer is an opportunity to harden the exact
path it exercises. KIN6D in particular should trigger this process as soon as
it depends on a reader, converter, field evaluator or geometry utility.

For each touched path:

1. freeze the KIN6D/downstream reproducer and expected physical semantics;
2. add a libneo-level analytic, manufactured, round-trip or cross-code oracle;
3. fix generic correctness/convention/accuracy defects in libneo, not in a
   downstream shadow implementation;
4. make the smallest refactor needed for reentrancy, explicit state,
   performance, derivative access or independent checking;
5. preserve existing public behavior when correct, or provide a bounded
   compatibility/migration path when a bug fix changes it;
6. run libneo tests and the relevant downstream regression before promotion.

The purpose is ecosystem compounding: fixes discovered while developing KIN6D
should improve NEO-2, SIMPLE, MEPHIT, TIAGO and future consumers as they update
their libneo revision.

### Differentiability and error capabilities

The current `field_t` value interface is intentionally small. Do not turn it
into a mandatory monster interface that every historical backend must
implement. Add capability-specific interfaces/composition only when real
consumers require them, for example:

- spatial field Jacobian/JVP evaluation;
- parameter-independent interpolation derivatives;
- interpolation/truncation diagnostics;
- interval/ball evaluation or other representation-level enclosures.

A value-only backend remains valid. A consumer that needs stronger semantics
may either select a backend with those capabilities or use libneo to parse and
normalize the source data, then rebuild it in a native differentiable/certified
representation.

Ownership remains strict:

- libneo owns derivatives and error information of its generic
  representation/conversion/evaluation mathematics;
- KIN6D owns derivatives of the physical stationary/evolution equations with
  respect to physical parameters, implicit equilibrium adjoints, validated
  solution existence/error and physical model-reduction bounds;
- TIAGO owns inference/posterior UQ.

This split allows rigorous KIN6D results without requiring every historical
libneo routine to become interval-arithmetic or AD-enabled.

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
