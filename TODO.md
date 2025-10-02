# libneo TODO

_Last updated: 2025-10-02_

---

## 🚧 ACTIVE: Cerfon-Freidberg vs ASCOT5 Regression

**Goal**: Cross-validate libneo's Cerfon-Freidberg analytical equilibrium against ASCOT5's reference B_GS implementation.

**Status**: ⚙️ **IN PROGRESS** – Solovev code paths have been removed from libneo; a new regression test now clones ASCOT5, builds the B_GS shared library, and overlays flux surfaces plus Shafranov shifts against libneo output. Plot generation and tighter tolerances are the remaining polish.

### Current Implementation Status
- ✅ Removed legacy Solovev Fortran modules and scrubbed build/test references.
- ✅ Restored the committed Cerfon-Freidberg analytical solver sources (`analytical_gs_solver.f90`, `analytical_tokamak_field.f90`) in magfie.
- ✅ `test_analytical_circular.x` focuses on Cerfon-Freidberg checks and emits CSV flux data for downstream comparisons.
- ✅ Added `test/scripts/test_ascot5_compare.py` which
  - runs the libneo analytical test to export circular flux data,
  - shallow-clones `ascot5` and builds a minimal shared library exposing `B_GS`,
  - compares flux contours/norms and saves an overlay plot under `build/test/ascot5_compare`.
- 🚧 Fine-tune RMS/axis-shift tolerances once the comparison is stable in CI.

### Motivation
ASCOT5's B_GS implementation is a de facto community reference for the Cerfon-Freidberg solution. Maintaining a regression against it ensures libneo stays numerically aligned.

### References to Double/Triple Check

#### 1. ASCOT5 analytical template
- `a5py/templates/analyticalinputs.py` – provides the ITER-like coefficient set used in the test.
- `src/Bfield/B_GS.c` – compiled into a shared object and called through `ctypes`.

#### 2. Cerfon & Freidberg paper (2010)
**Reference**: Physics of Plasmas 17, 032502 (2010), DOI: 10.1063/1.3328818
- Equation (8) lists the basis functions that both codes implement.
- Shafranov shift Δ ≈ 0.409 m (ITER parameters) provides the quantitative target.

### Implementation Plan

**§0. Reference data** – ✅ Use ASCOT5 template coefficients (no SciPy dependency).

**§1. Regression test** – ✅ Add Python test that
- runs libneo's analytical solver to dump flux surfaces,
- builds ASCOT5's `B_GS` into a shared library via gcc,
- compares flux RMS and Shafranov shift, and
- stores an overlay plot in `build/test/ascot5_compare`.

**§2. Tolerance tuning** – 🚧 Inspect RMS/axis logs and tighten thresholds once the comparison stabilises.


## ✅ COMPLETED: Analytical Tokamak Field (Cerfon-Freidberg) - FULL CASE

**Goal**: Implement complete Cerfon-Freidberg "One size fits all" analytical equilibrium solver with elongation and triangularity

**Status**: ✅ **IMPLEMENTATION COMPLETE** - Full general solver with κ (elongation) and δ (triangularity) support

**Completed**: 2025-10-02

### Summary
- ✅ All 7 Cerfon-Freidberg basis functions implemented with verified derivatives
- ✅ Particular solution implemented and verified
- ✅ General 7×7 boundary condition solver using LAPACK dgesv
- ✅ Supports arbitrary ε, κ, δ parameters (circular and shaped plasmas)
- ✅ All 56 tests passing (boundary conditions, ∇·B=0 verification)
- ✅ PNG flux surface visualization for both circular and shaped cases
- ✅ Derivative verification tool: `tools/verify_gs_derivatives.py`
- ✅ 2 derivative errors found and fixed in ψ₇ via SymPy verification

### References
1. **Primary**: Cerfon & Freidberg, "One size fits all" analytic solutions to the Grad-Shafranov equation,
   Physics of Plasmas 17, 032502 (2010), DOI: 10.1063/1.3328818
2. **Implementation source**: Verena Eslbauer, "Two analytical solutions to the Grad-Shafranov equation
   using Solovev pressure and poloidal current profiles", Bachelor thesis, TU Graz, November 20, 2017
3. **MATLAB code**: `/Users/ert/Dropbox/proj/stud/Bacc_Verena_Eslbauer/One_size_fits_all.m`

### Background
Fortran port of Eslbauer's COMPLETE MATLAB implementation supporting arbitrary shapes.
Uses 7 basis functions to solve linearized Grad-Shafranov equation with Solovev profiles.
**Supports arbitrary ε (aspect ratio), κ (elongation), and δ (triangularity)**.

### Implementation Plan

### Mathematical Formulation (from Eslbauer thesis, Section 3.2)

**Normalized coordinates**: x = R/R₀, y = Z/R₀

**Linearized Grad-Shafranov equation**:
```
Δ*ψ = x·∂/∂x(1/x·∂ψ/∂x) + ∂²ψ/∂y² = (1-A)x² + A
```

**Solution decomposition**:
```
ψ(x,y) = ψ_p(x,y) + Σᵢ₌₁⁷ cᵢ·ψᵢ(x,y)
```

**Basis functions** (Cerfon-Freidberg Eq. 8, up-down symmetric):
```
ψ₁ = 1
ψ₂ = x²
ψ₃ = y² - x²ln(x)
ψ₄ = x⁴ - 4x²y²
ψ₅ = 2y⁴ - 9x²y² + [3x⁴ - 12x²y²]ln(x)
ψ₆ = x⁶ - 12x⁴y² + 8x²y⁴
ψ₇ = 8y⁶ - 140x²y⁴ + 75x⁴y² + [180x⁴y² - 120x²y⁴ - 15x⁶]ln(x)
```

**Particular solution** (Cerfon-Freidberg Eq. 9):
```
ψ_p = x⁴/8 + A(x²ln(x)/2 - x⁴/8)
```

**Parametric plasma boundary** (Eslbauer Eq. 3.14):
```
x(τ) = 1 + ε·cos(τ + δ·sin(τ))
y(τ) = ε·κ·sin(τ)
τ ∈ [0, 2π]
```
where ε = inverse aspect ratio, κ = elongation, δ = triangularity

**Boundary conditions** (7 equations for GENERAL shaped plasma):
1. ψ(1+ε, 0) = 0  (outer equatorial, τ=0)
2. ψ(1-ε, 0) = 0  (inner equatorial, τ=π)
3. ψ(1-δε, κε) = 0  (high point, τ=π/2)
4. ∂ψ/∂x(1-δε, κε) = 0  (high point is extremum)
5. ∂²ψ/∂y²(1+ε, 0) = -N₁·∂ψ/∂x(1+ε, 0)  (outer curvature)
6. ∂²ψ/∂y²(1-ε, 0) = -N₂·∂ψ/∂x(1-ε, 0)  (inner curvature)
7. ∂²ψ/∂x²(1-δε, κε) = -N₃·∂ψ/∂y(1-δε, κε)  (high point curvature)

**Curvature coefficients** (Eslbauer page 19, general case):
```
α = arcsin(δ)  (triangularity angle)
N₁ = -(1+α)²/(ε·κ²)
N₂ = (1-α)²/(ε·κ²)
N₃ = -κ/(ε·cos²(α))
```

#### Files to Update/Create
1. **RENAME**: `analytical_gs_circular.f90` → `analytical_gs_solver.f90`
   - 7 basis functions ψᵢ(x,y) with first and second derivatives
   - Particular solution ψ_p(x,y, A)
   - General boundary solver supporting ε, κ, δ
   - LAPACK dgesv for 7×7 system
   - **DELETE ALL HARDCODED COEFFICIENTS**

2. `src/magfie/analytical_tokamak_field.f90` - **UPDATE EXISTING**
   - Update to support κ, δ parameters (not just circular!)
   - Change init signature: `init(R0, epsilon, kappa, delta, A_param, B0)`
   - Update coefficient indexing (1:7 not 0:6)
   - Add proper Cerfon-Freidberg + Eslbauer attribution

3. `test/source/test_analytical_circular.f90` - **UPDATE EXISTING**
   - Test ITER shaped case (ε=0.32, κ=1.7, δ=0.33, A=-0.142)
   - Test circular limit (κ=1, δ=0)
   - Remove all ASCOT5 references
   - Verify ∇·B=0 for shaped plasmas

#### Interface Design

```fortran
! Module: analytical_tokamak_field
type :: analytical_tokamak_eq_t  ! Renamed from circular!
    real(dp) :: R0, epsilon, kappa, delta, A_param, B0
    real(dp) :: coeffs(7)  ! 1:7 indexing
    real(dp) :: psimult
    logical :: initialized = .false.
contains
    procedure :: init => init_tokamak_equilibrium
    procedure :: eval_psi => evaluate_psi
    procedure :: eval_psi_derivatives
    procedure :: eval_bfield => evaluate_bfield
    procedure :: cleanup => destroy_equilibrium
end type

! Initialization - NOW WITH KAPPA AND DELTA
subroutine init(self, R0, epsilon, kappa, delta, A_param, B0, psimult_in)
    ! Solves 7×7 boundary condition system for general shape
    ! kappa = elongation (κ)
    ! delta = triangularity (δ)
end subroutine

! Evaluation (unchanged interface)
function eval_psi(self, R, Z) result(psi)
subroutine eval_bfield(self, R, Z, B_R, B_Z, B_phi, B_mod)
```

#### Detailed Implementation Steps

**§0. Baseline** ✅
- [x] Document complete mathematical formulation from Eslbauer thesis
- [x] Update TODO.md with all equations and references
- [x] Note: Need to REPLACE hardcoded implementation

**§1. Rewrite Basis Functions** ✅
- [x] **CRITICAL**: Basis functions are indexed 1-7 (not 0-6)
- [x] Implement 7 basis functions from Cerfon-Freidberg Eq. 8
- [x] Implement first derivatives for all 7 functions
- [x] Implement second derivatives (∂²/∂x², ∂²/∂y²) needed for boundary conditions
- [x] Remove old psi_0 through psi_6 functions
- [x] **VERIFIED**: All derivatives checked with SymPy (2 errors found and fixed in psi_7)

**§2. Fix Particular Solutions** ✅
- [x] Implement correct particular solution (single function):
  ```fortran
  ψ_p(x,y) = x⁴/8 + A·(x²ln(x)/2 - x⁴/8)
  ```
- [x] This replaces the current split into psi_part_0 and psi_part_1
- [x] Implement derivatives ∂ψ_p/∂x, ∂ψ_p/∂y
- [x] Implement second derivatives ∂²ψ_p/∂x², ∂²ψ_p/∂y²
- [x] **VERIFIED**: All derivatives checked with SymPy

**§3. Implement GENERAL Boundary Condition Solver** ✅
- [x] **DELETE** hardcoded coefficient section (lines 207-218)
- [x] Implement 7×7 matrix setup for GENERAL shaped case:
  ```fortran
  subroutine solve_coefficients(epsilon, kappa, delta, A_param, coeffs)
    real(dp), intent(in) :: epsilon, kappa, delta, A_param
    real(dp), intent(out) :: coeffs(7)
    real(dp) :: mat(7,7), rhs(7)
    real(dp) :: x_out, x_in, x_high, y_high
    real(dp) :: alpha, N1, N2, N3
    integer :: ipiv(7), info

    ! Boundary points from parametric representation
    x_out = 1.0_dp + epsilon     ! τ=0
    x_in = 1.0_dp - epsilon      ! τ=π
    x_high = 1.0_dp - delta*epsilon  ! τ=π/2
    y_high = kappa*epsilon

    ! Curvature coefficients
    alpha = asin(delta)
    N1 = -(1.0_dp + alpha)**2 / (epsilon * kappa**2)
    N2 = (1.0_dp - alpha)**2 / (epsilon * kappa**2)
    N3 = -kappa / (epsilon * cos(alpha)**2)

    ! Row 1: ψ(x_out, 0) = 0
    mat(1,:) = [eval basis at (x_out, 0)]
    rhs(1) = -psi_p(x_out, 0.0_dp, A_param)

    ! Row 2: ψ(x_in, 0) = 0
    mat(2,:) = [eval basis at (x_in, 0)]
    rhs(2) = -psi_p(x_in, 0.0_dp, A_param)

    ! Row 3: ψ(x_high, y_high) = 0
    mat(3,:) = [eval basis at (x_high, y_high)]
    rhs(3) = -psi_p(x_high, y_high, A_param)

    ! Row 4: ∂ψ/∂x(x_high, y_high) = 0
    mat(4,:) = [eval ∂ψᵢ/∂x at (x_high, y_high)]
    rhs(4) = -dpsi_p_dx(x_high, y_high, A_param)

    ! Row 5: ∂²ψ/∂y²(x_out, 0) + N₁·∂ψ/∂x(x_out, 0) = 0
    ! Row 6: ∂²ψ/∂y²(x_in, 0) + N₂·∂ψ/∂x(x_in, 0) = 0
    ! Row 7: ∂²ψ/∂x²(x_high, y_high) + N₃·∂ψ/∂y(x_high, y_high) = 0

    call dgesv(7, 1, mat, 7, ipiv, rhs, 7, info)
    if (info /= 0) error stop "Boundary condition matrix singular"
    coeffs = rhs
  end subroutine
  ```
- [x] Test with ITER shaped: ε=0.32, κ=1.7, δ=0.33, A=-0.142
- [x] Test circular limit: κ=1, δ=0

**§4. Update Field Evaluation Module** ✅
- [x] Add kappa, delta fields to type (kept name for backward compat)
- [x] Update init signature: add kappa, delta optional parameters
- [x] Update coefficient array: `coeffs(0:6)` → `coeffs(7)`
- [x] Update solve_coefficients call to pass kappa, delta
- [x] Update ψ evaluation loop (1:7 indexing)

**§5. Update Tests** ✅
- [x] Test ITER shaped case: R₀=6.2, ε=0.32, κ=1.7, δ=0.33, A=-0.142, B₀=5.3
- [x] Test circular limit: same but κ=1, δ=0
- [x] Test ∇·B=0 for both cases (56 tests total, all passing)
- [x] Added PNG visualization output (flux_circular.png, flux_shaped.png in build/test/)
- [x] Proper attribution to Cerfon-Freidberg and Eslbauer

**§6. Build & Verify** ✅
- [x] Build: `make clean && make` - SUCCESS
- [x] Run tests: `ctest -R test_analytical_circular` - ALL PASS (56/56 tests)
- [x] PNG artifacts generated successfully

**§7. Documentation & Final Cleanup** ✅
- [x] Proper attribution in all files
- [x] Module-level documentation updated
- [x] Derivative verification tool created: `tools/verify_gs_derivatives.py`

**§8. Validation & Commit** - IN PROGRESS
- [x] All analytical tokamak tests pass (56/56)
- [x] Boundary conditions satisfied to machine precision (tol=1e-10)
- [x] ∇·B=0 verified for both circular and shaped plasmas
- [x] PNG flux surface plots generated
- [ ] Commit and push changes

### Success Criteria
- General solver works for arbitrary ε and A (not just ITER case)
- ∇·B = 0 to machine precision everywhere in plasma
- Flux surfaces are circular (κ=1, δ=0) for circular case
- Boundary conditions satisfied at separatrix
- All ASCOT5 references removed
- Proper attribution to Cerfon-Freidberg and Eslbauer
- Tests pass with improved accuracy

### Key Implementation Notes from Eslbauer Thesis

**From Section 3.2, pages 17-22**:

1. **Basis function indexing**: ψᵢ for i=1..7 (MATLAB c_1 to c_7, Fortran coeffs(1:7))

2. **Particular solution**: ψ_p = x⁴/8 + A(x²ln(x)/2 - x⁴/8) [Eq. 3.12-3.13]

3. **General boundary points** (parametric, works for ALL shapes):
   - Outer equatorial: x=1+ε, y=0 (τ=0)
   - Inner equatorial: x=1-ε, y=0 (τ=π)
   - High point: x=1-δε, y=κε (τ=π/2)

4. **Curvature coefficients** (GENERAL formulas, Eq. 3.14e-g):
   ```fortran
   alpha = asin(delta)
   N1 = -(1 + alpha)**2 / (epsilon * kappa**2)
   N2 = (1 - alpha)**2 / (epsilon * kappa**2)
   N3 = -kappa / (epsilon * cos(alpha)**2)
   ```
   These automatically reduce to circular case when κ=1, δ=0!

5. **MATLAB symbolic solve → Fortran LAPACK numeric**:
   - MATLAB uses `solve([eqn1...eqn7])` symbolically
   - We build mat(7,7) and rhs(7) numerically, call dgesv

6. **ITER shaped parameters** (Table 3.1, page 21):
   - R₀=6.2 m, ε=0.32, κ=1.7, δ=0.33, B₀=5.3 T, A=-0.142
   - Expected coefficients (Table 3.2): c₁=0.0673, c₂=-0.1827, c₃=-0.0414,
     c₄=-0.1503, c₅=0.0014, c₆=-0.0106, c₇=-0.0003

7. **Circular limit test**: Same ITER but κ=1, δ=0
   - Coefficients will differ from shaped case
   - Validates general solver reduces correctly

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
