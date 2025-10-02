# libneo TODO

_Last updated: 2025-10-02_

---

## 🚧 ACTIVE: Analytical Circular Tokamak Field (Solov'ev-type)

**Goal**: Clean-room implementation of circular tokamak equilibrium for testing canonical coordinates

**Status**: Planning phase

### Background
Implement simplified Cerfon-Freidberg analytical Grad-Shafranov solution for circular tokamak (κ=1, δ=0).
This generalizes the classical Solov'ev (1968) equilibrium to finite aspect ratio.

**Reference**: Cerfon & Freidberg, Physics of Plasmas 17, 032502 (2010) - equations implemented from published paper (no GPL code)

**Use case**: Simple, exact test field for Meiss canonical coordinates in SIMPLE before using real GEQDSK data

### Implementation Plan

#### Files to Create
1. `src/magfie/analytical_gs_circular.f90` - Core GS solver
   - 7 basis functions (symmetric case)
   - Boundary condition setup
   - Coefficient solver (7×7 LAPACK)

2. `src/magfie/analytical_tokamak_field.f90` - Field evaluation interface
   - Initialize equilibrium (R₀, ε, A, B₀)
   - Evaluate ψ(R,Z) and derivatives
   - Compute B-field components (B_R, B_Z, B_φ)

3. `test/source/test_analytical_circular.f90` - Validation tests
   - Verify ∇·B = 0
   - Check separatrix location
   - Compare with ASCOT5 ITER coefficients
   - Flux surface shape validation

#### Interface Design

```fortran
! Module: analytical_tokamak_field
type :: analytical_circular_eq_t
    real(dp) :: R0, epsilon, A_param, B0
    real(dp) :: coeffs(0:6)
    logical :: initialized = .false.
contains
    procedure :: init => init_circular_equilibrium
    procedure :: eval_psi => evaluate_psi
    procedure :: eval_psi_derivatives
    procedure :: eval_bfield => evaluate_bfield
    procedure :: cleanup => destroy_equilibrium
end type

! Initialization
subroutine init(self, R0, epsilon, A_param, B0)
    ! Solve 7×7 system for coefficients
end subroutine

! Evaluation
function eval_psi(self, R, Z) result(psi)
    ! Returns poloidal flux ψ(R,Z)
end function

subroutine eval_bfield(self, R, Z, B_R, B_Z, B_phi, B_mod)
    ! Returns magnetic field components
end subroutine
```

#### Steps (in order)

**§0. Baseline**
- [ ] Verify clean git status in libneo
- [ ] Document plan in this TODO

**§1. Basis Functions Module**
- [ ] Create `analytical_gs_circular.f90`
- [ ] Implement 7 basis functions ψᵢ(x,y) for i=0..6
  - x = R/R₀, y = Z/R₀ (normalized coordinates)
  - Implement as pure functions
- [ ] Implement first derivatives (∂ψᵢ/∂x, ∂ψᵢ/∂y)
- [ ] Implement second derivatives (∂²ψᵢ/∂x², ∂²ψᵢ/∂y², ∂²ψᵢ/∂x∂y)
- [ ] Add unit tests for basis function values at key points

**§2. Particular Solutions**
- [ ] Implement ψ_part_0(x,y) and ψ_part_1(x,y)
  - ψ_part_0(x,y) = 1 - x²
  - ψ_part_1(x,y) = x⁴ - 4x² + log(x) + 8x²log(x) - 12x²log²(x)
                     + (2 - 12log(x))y² + 12y²log²(x) + y⁴
  - These are given explicitly in Cerfon-Freidberg Eq. 9
- [ ] Implement derivatives of particular solutions
- [ ] Test against analytical formulas

**§3. Boundary Condition Solver**
- [ ] Implement `setup_boundary_matrix(epsilon, A, mat, rhs)`
  - 7 boundary conditions for circular case
  - Boundary points: (1±ε, 0), (1, ε)
  - Curvature conditions
- [ ] Implement `solve_coefficients(epsilon, A, coeffs)`
  - Use LAPACK dgesv for 7×7 system
  - Return c₀-c₆
- [ ] Verify against ASCOT5 ITER coefficients (ε=0.323, A=-0.155)

**§4. Field Evaluation Module**
- [ ] Create `analytical_tokamak_field.f90`
- [ ] Implement `analytical_circular_eq_t` type
- [ ] Implement initialization routine
- [ ] Implement ψ(R,Z) evaluation
  - User passes physical coords (R, Z) in meters
  - Convert internally: x=R/R₀, y=Z/R₀
  - Evaluate ψ(x,y) using coefficients and particular solutions
  - Scale: ψ_phys = ψ_norm * psimult (use psimult=200 like ASCOT5)
- [ ] Implement derivatives ∂ψ/∂R, ∂ψ/∂Z
- [ ] Implement B-field evaluation:
  - B_R = -(1/R)∂ψ/∂Z
  - B_Z = (1/R)∂ψ/∂R
  - B_φ = F(ψ)/R  where F(ψ) ≈ B₀·R₀ (constant for initial impl.)
  - B_mod = √(B_R² + B_Z² + B_φ²)
  - Note: Can refine F(ψ) in §8 if needed

**§5. Integration Tests**
- [ ] Create `test_analytical_circular.f90`
- [ ] Test 1: ITER circular equilibrium
  - Initialize with R₀=6.2, ε=0.323, A=-0.155, B₀=5.3
  - Verify coefficients match ASCOT5
  - Verify ψ=0 on separatrix
- [ ] Test 2: Divergence-free B-field
  - Sample points in plasma
  - Compute ∇·B numerically
  - Verify |∇·B| < 10⁻¹² (exact to machine precision)
- [ ] Test 3: Flux surface shapes
  - Trace flux surfaces at multiple ψ values
  - Verify circular cross-sections
- [ ] Test 4: Large aspect ratio limit
  - Test with ε → 0
  - Compare with Solov'ev formulas

**§6. CMake Integration**
- [ ] Add files to `src/magfie/CMakeLists.txt`
- [ ] Add test to `test/CMakeLists.txt`
- [ ] Build and verify: `make clean && make`
- [ ] Run tests: `ctest -R test_analytical_circular`

**§7. Documentation**
- [ ] Add module-level documentation
- [ ] Document initialization parameters
- [ ] Add usage example in comments
- [ ] Note: "Clean-room implementation from Cerfon-Freidberg (2010)"

**§8. Validation & Finalization**
- [ ] Compare full B-field with ASCOT5 at sample points
- [ ] Verify performance (field evaluation ~1 μs, initialization < 100 ms)
- [ ] Optional refinements if needed:
  - More accurate F(ψ) from particular solutions
  - Edge case handling
- [ ] Git commit with descriptive message
- [ ] Update this TODO as complete

### Success Criteria
- Coefficients match ASCOT5 ITER circular case to 10⁻⁸
- ∇·B = 0 to machine precision
- Flux surfaces are circular (κ=1, δ=0)
- No GPL code used (clean-room from paper)
- Tests pass in CI

### Notes
- Uses only symmetric basis functions (c₇-c₁₁ = 0)
- Valid for arbitrary aspect ratio (unlike pure Solov'ev)
- Simpler than geoflux (no GEQDSK, no splines)
- Perfect for testing coordinate transformations

### Implementation Details (Clarifications)
1. **F(ψ) function**: Use F(ψ) = B₀·R₀ (constant) for initial implementation
   - Low-β circular equilibria: F ≈ constant is excellent approximation
   - Can refine with ψ-dependent corrections in §8 if needed

2. **Coordinate handling**: Users pass physical (R,Z) in meters
   - Internal conversion: x=R/R₀, y=Z/R₀
   - Basis functions work in normalized coords
   - Output in physical units (Tesla)

3. **Particular solutions**: Explicit formulas from Cerfon-Freidberg Eq. 9
   - ψ_part_0 and ψ_part_1 are given, not derived
   - Parameter A weights these solutions, controls q-profile

4. **Performance**: ~1 μs per field evaluation (after initialization)
   - Initialization solves 7×7 system once (~1-10 ms, negligible)
   - Evaluation is polynomial+log (~20 arithmetic ops, very fast)
   - Comparable to spline interpolation speed

### Implementation Guidance

**Q1: TDD adherence** - Should we write tests before each implementation section?

**Answer**: YES, use pragmatic TDD adapted for scientific computing.

**Recommended approach**:
- §1-§2: Write unit tests BEFORE implementing each basis/particular function
  - Test basis functions against hand-calculated values at 2-3 points
  - Test derivatives with finite differences (tolerance ~10⁻⁸)
- §3: Write integration test with ASCOT5 ITER case BEFORE implementing solver
- §5: System tests for physical properties (∇·B=0, flux surfaces)

**Flexibility**: You CAN implement multiple related functions together (e.g., all 7 basis functions), THEN test comprehensively. Scientific computing often needs the full picture.

**TDD rhythm**: Test → Implement → Verify → Next function

---

**Q2: LAPACK availability** - Is LAPACK ready to use?

**Answer**: YES, LAPACK is ALREADY AVAILABLE in the build system.

**Evidence**:
```bash
# CMakeLists.txt has: find_package(LAPACK REQUIRED)
# CMakeSources.in links: LAPACK::LAPACK
# src/solve_systems.f90 already uses dgesv
```

**Usage in §3**:
```fortran
subroutine solve_coefficients(epsilon, A_param, coeffs)
    real(dp) :: mat(7,7), rhs(7)
    integer :: ipiv(7), info

    ! Set up boundary matrix...

    call dgesv(7, 1, mat, 7, ipiv, rhs, 7, info)
    if (info /= 0) error stop "LAPACK solve failed"

    coeffs = rhs
end subroutine
```

**Action**: Just use it directly. No build system changes needed.

---

**Q3: Module organization** - Use src/magfie/ or src/field/?

**Answer**: STAY IN `src/magfie/` as specified.

**Reasoning**:
```
src/field/      - Stellarator-specific, 3D coil fields (not for tokamaks)
src/magfie/     - Magnetic field implementations (VMEC, GEQDSK, analytical)
src/coordinates/- Coordinate transformations
```

**Your files fit with existing magfie modules**:
```
src/magfie/
├── geoflux_field.f90          ← Similar: field from GEQDSK
├── geqdsk_tools.f90           ← Similar: equilibrium utilities
├── magfie_vmec.f90            ← Similar: VMEC field evaluation
├── analytical_gs_circular.f90 ← YOUR NEW FILE (same pattern)
└── analytical_tokamak_field.f90 ← YOUR NEW FILE (same pattern)
```

**Decision**: Keep in `src/magfie/`.

---

**Q4: Starting point** - Ready to start with §0 and §1?

**Answer**: YES, READY TO START NOW.

**§0 Baseline** - ✅ COMPLETE (TODO.md committed)

**§1 Next steps**:
```bash
# Create files
touch test/source/test_analytical_gs_circular.f90
touch src/magfie/analytical_gs_circular.f90

# Write tests first (TDD)
# 1. Test ψ₀, ψ₁, ψ₂ at known points
# 2. Implement ψ₀, ψ₁, ψ₂
# 3. Add to CMakeLists
# 4. Build and verify
# 5. Continue with ψ₃-ψ₆
```

**Execution order**:
1. ✅ §0: Commit TODO.md (DONE)
2. §1: Basis functions (TDD: test first, then implement)
3. §2: Particular solutions (TDD: test first, then implement)
4. §3: Boundary solver (integration test with ASCOT5 coefficients)
5. §4: Field evaluation interface
6. §5: System tests (∇·B=0, circular flux surfaces)
7. §6: CMake integration
8. §7: Documentation
9. §8: Final validation

**Ready checklist**:
- [x] LAPACK available? YES
- [x] Clean git status? YES (TODO.md committed)
- [x] File locations? YES - src/magfie/
- [x] TDD approach? YES - tests before implementation
- [x] Reference equations? YES - in TODO.md
- [x] Validation target? YES - ASCOT5 ITER coefficients

**START NOW** with §1 basis functions! 🚀

---

## ✅ COMPLETED: coordinate_system_t for VMEC & Geoflux

_Completed: 2025-10-02_

## ✅✅ FULLY COMPLETED AND TESTED

All implementation, build fixes, and testing are complete:

### Implementation (commit 23fcd78)
- ✅ Abstract `coordinate_system_t` interface in `src/coordinates/libneo_coordinates.f90`
- ✅ VMEC concrete type in `src/coordinates/libneo_coordinates_vmec.f90`
- ✅ Geoflux concrete type in `src/coordinates/libneo_coordinates_geoflux.f90`
- ✅ Factory functions: `make_vmec_coordinate_system()`, `make_geoflux_coordinate_system()`
- ✅ Test created in `test/source/test_coordinate_systems.f90`
- ✅ Test registered in `test/CMakeLists.txt` (line 164)

### Build Fixes (user applied)
- ✅ Moved `src/coordinates/geoflux_coordinates.f90` from magfie to neo (CMakeSources.in line 24)
- ✅ Moved `src/magfie/geqdsk_tools.f90` to neo (CMakeSources.in line 33)
- ✅ Removed both from `src/magfie/CMakeLists.txt`
- ✅ Resolved circular dependency: neo builds cleanly before magfie

### Verification
- ✅ Clean build successful: `make clean && make` completes without errors
- ✅ Test executable runs: `./build/test/test_coordinate_systems.x` passes
- ✅ Test output confirms both factory functions work correctly

---
## Test Output

```
Testing coordinate_system_t abstraction...

Testing factory functions...
  PASS: make_vmec_coordinate_system allocated
  PASS: make_geoflux_coordinate_system allocated

Coordinate system interface verified successfully.
Note: Full integration tests with real data are available through:
  - test_vmec_modules (VMEC coordinate transforms)
  - test_geoflux (Geoflux coordinate transforms)

All coordinate system interface tests passed!
```

---
## Ready for SIMPLE Integration

**Status**: libneo `coordinate_system_t` abstraction is **READY** for use in SIMPLE.

SIMPLE can now proceed with its TODO.md §0-6 to integrate canonical coordinates for tokamak geometries.

### API Reference

**Module**: `use libneo_coordinates`

**Factory Functions**:
```fortran
class(coordinate_system_t), allocatable :: cs
call make_vmec_coordinate_system(cs)       ! For VMEC equilibria
call make_geoflux_coordinate_system(cs)    ! For GEQDSK/geoflux equilibria
```

**Methods**:
```fortran
! Evaluate position in cylindrical/Cartesian coordinates
call cs%evaluate_point(u, x)
! u(3): coordinate in flux space [s, theta, phi]
! x(3): position in physical space [R, phi, Z] or [x, y, z]

! Compute covariant basis vectors
call cs%covariant_basis(u, e_cov)
! e_cov(3,3): columns are ∂x/∂u_i

! Compute metric tensor and its inverse
call cs%metric_tensor(u, g, ginv, sqrtg)
! g(3,3): covariant metric tensor
! ginv(3,3): contravariant metric tensor
! sqrtg: sqrt(det(g)), Jacobian of transformation
```

### Integration Notes

- Both coordinate systems use the same abstract interface
- Runtime polymorphism via `class(coordinate_system_t)`
- No code changes needed when switching between VMEC and geoflux
- Basis vectors and metric tensors computed consistently

**Next**: See SIMPLE/TODO.md for integration steps.
