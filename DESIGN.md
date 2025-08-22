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
2. **Phase 2 - Test Coverage**: Ensure adequate test coverage before refactoring
3. **Phase 3 - Refactoring**: Move inner subroutines to module level
4. **Phase 4 - Validation**: Verify performance improvements and correctness

### Identified Cases

#### poincare.f90
- **Location**: `integrate_RZ_along_fieldline` contains inner subroutine `fieldline_derivative`
- **Issue**: `fieldline_derivative` accesses `field` parameter from outer scope
- **Solution**: Move to module level, pass `field` as argument or use context parameter
- **Test Coverage**: Exists in `test/poincare/test_poincare.f90`
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

### Test Coverage Requirements
- All refactored code must have tests before modification
- Tests must verify both correctness and performance characteristics
- Use the existing test framework in `test/`

## Build System
- Primary: CMake with Ninja
- Secondary: fpm (Fortran Package Manager)
- All changes must pass CI pipeline

## Code Standards
See `CODING_STANDARD.md` for detailed coding standards.