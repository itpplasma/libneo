# DEBUG.md - Systematic Investigation of ntor=0 Field Issues

## Problem Statement

For the single coil test (128-segment nearly-circular coil, tilted), all numerical methods (Fourier, Anvac, Direct) produce:
- **Bφ,n=0 ≈ -0.01 G** at coil center
- This is **50-90% of other components** (BR ≈ 0.001 G, BZ ≈ 0.005 G)

**Expected**: For a circular coil with 128 segments, Bφ,n=0 should be ≈ 0 (< 0.0001 G)

**All methods agree with each other (<1% difference)**, suggesting either:
1. All are computing the same wrong thing (systematic bug)
2. The field is actually correct but expectations are wrong

**User says**: "128 segments is PLENTY" - implying the field should be nearly axisymmetric

## Measured Values

### At Coil Center (R=210.32 cm, Z=78.21 cm)

| Source | BnR (G) | Bnφ (G) | BnZ (G) | \|B\| (G) |
|--------|---------|---------|---------|-----------|
| Fourier (H5) | 1.055e-03 | **-9.509e-03** | 5.028e-03 | 1.081e-02 |
| Anvac (NC) | 1.054e-03 | **-9.459e-03** | 5.031e-03 | 1.076e-02 |
| **Expected** | ~0.001 | **~0** | ~0.005 | ~0.005 |

**Issue**: Bφ is 9x larger than BR, 2x larger than BZ - should be nearly zero!

## Coil Geometry

From `test/magfie/test_data/single_coil.dat`:
- **Segments**: 129 points (128 segments)
- **Shape**: Planar (flat to <4 nm), tilted
- **Plane normal**: (0.350, -0.451, 0.821)
- **Center**: (197.05, 71.84, 78.0) cm
- **Radius**: 35.00 ± 0.19 cm (0.5% variation)
- **Eccentricity**: ~0.02 (nearly circular)

**Conclusion**: Coil is sufficiently circular that Bφ,n=0 should be << 0.001 G

## Test Configuration

From `test/magfie/run_single_coil_case.cmake`:
```
--ntor 0
--axis-origin 197.2682697 72.00853957 78.0  (cm)
--axis-normal 0.35045589 -0.45058614 0.82106808  (unit vector)
--coil-radius 35.0  (cm)
--axis-range 60.0  (cm)
--axis-samples 181
```

**Test compares**: Single n=0 Fourier mode, NOT full 3D field reconstruction

## Hypotheses to Investigate

### Hypothesis 1: Fortran Cartesian → Cylindrical Transformation Bug

**Formula used** (line 469, src/magfie/coil_tools.f90):
```fortran
BR(kphi) = BXYZ(1) * cosphi(kphi) + BXYZ(2) * sinphi(kphi)
Bphi(kphi) = BXYZ(2) * cosphi(kphi) - BXYZ(1) * sinphi(kphi)
BZ(kphi) = BXYZ(3)
```

**Expected**:
```
BR = Bx*cos(φ) + By*sin(φ)  ✓
Bφ = -Bx*sin(φ) + By*cos(φ)  ← CHECK THIS!
BZ = Bz  ✓
```

**Current formula for Bφ**: `By*cos(φ) - Bx*sin(φ)`
**Correct formula**: `-Bx*sin(φ) + By*cos(φ)`

These are algebraically equivalent (just reordered), but **check the sign**!

**Action**: Verify transformation with test case where Bx, By, φ are known

### Hypothesis 2: FFT Normalization Error

**FFT execution** (line 474-476, src/magfie/coil_tools.f90):
```fortran
call fftw_execute_dft_r2c(plan_nphi, BR, BnR)
call fftw_execute_dft_r2c(plan_nphi, Bphi, Bnphi)
call fftw_execute_dft_r2c(plan_nphi, BZ, BnZ)
```

**Normalization** (line 477-479):
```fortran
Bn(0:nmax, 1, kR, kZ, kc) = BnR(1:nmax+1) / dble(nphi)
Bn(0:nmax, 2, kR, kZ, kc) = Bnphi(1:nmax+1) / dble(nphi)
Bn(0:nmax, 3, kR, kZ, kc) = BnZ(1:nmax+1) / dble(nphi)
```

**Question**: Is this the correct normalization for r2c FFT?

FFTW r2c computes: `B̃[k] = Σ B[j] * e^(-2πijk/N)`

For n=0: `B̃[0] = Σ B[j]` (sum of all samples)

To get average: should divide by N (which is `nphi`)

**This looks correct**, but check if there's a factor of 2 issue for n=0 vs n>0 modes.

### Hypothesis 3: Sampling Grid Issue

**φ sampling** (line 422-423, src/magfie/coil_tools.f90):
```fortran
cosphi = [(cos(2.0_dp * pi * real(kphi - 1, dp) / real(nphi, dp)), kphi=1, nphi)]
sinphi = [(sin(2.0_dp * pi * real(kphi - 1, dp) / real(nphi, dp)), kphi=1, nphi)]
```

This samples φ ∈ [0, 2π) at `nphi` points: φ_k = 2π(k-1)/nphi for k=1..nphi

**Question**: Does this include the endpoint 2π? NO - goes from 0 to 2π(nphi-1)/nphi

For FFT, this is correct (periodic, don't duplicate endpoint)

**Check**: Is `nphi` set correctly? From input file: `nphi = 192`

### Hypothesis 4: Curl Formula Error in Python

**Python curl for ntor=0** (line 180-182, python/libneo/biotsavart_fourier.py):
```python
BnR = -dAnphi_dZ
Bnphi = dAnR_dZ - dAnZ_dR
BnZ = dAnphi_dR + Anphi / R[newaxis, :, newaxis]
```

**Expected for axisymmetric** (∂/∂φ = 0):
```
BR = -∂Aφ/∂Z  ✓
Bφ = ∂AR/∂Z - ∂AZ/∂R  ✓
BZ = ∂Aφ/∂R + Aφ/R  ✓
```

**This is correct** - verified with SymPy

### Hypothesis 5: Vector Potential Gauge Issue

For ntor=0, **gauge_Anvac returns AR, AZ unchanged** (line 95-96):
```python
if ntor == 0:
    return AnR, AnZ
```

**Question**: Are AR, AZ themselves correct from Fortran?

The Fortran computes A from Biot-Savart, then transforms to cylindrical, then FFTs.

**Transformation formula** (line 596-598, src/magfie/coil_tools.f90):
```fortran
AR(kphi) = AXYZ(1) * cosphi(kphi) + AXYZ(2) * sinphi(kphi)
Aphi(kphi) = AXYZ(2) * cosphi(kphi) - AXYZ(1) * sinphi(kphi)
AZ(kphi) = AXYZ(3)
```

**Expected**:
```
AR = Ax*cos(φ) + Ay*sin(φ)  ✓
Aφ = -Ax*sin(φ) + Ay*cos(φ)  ← CHECK THIS!
AZ = Az  ✓
```

**Current Aφ formula**: `Ay*cos(φ) - Ax*sin(φ)`

**Verified with SymPy**: This is correct!

### Hypothesis 6: Axis Evaluation Error

**Axis validation evaluates at points along tilted axis** (line 554, plot_biotsavart_fourier.py):
```python
points = origin[newaxis, :] + outer(s_vals, direction)
```

Where `direction = [0.350, -0.451, 0.821]` (normalized axis normal)

For each point `(x, y, z)`:
1. Convert to cylindrical: `R = hypot(x,y)`, `φ = arctan2(y,x)`, `Z = z`
2. Interpolate B̃(R,Z) at that (R,Z)
3. For ntor=0: `B_total = B̃.real` (no φ-dependence)
4. Convert to Cartesian: `Bx = BR*cos(φ) - Bφ*sin(φ)`, etc.
5. Project: `B_parallel = B·direction`

**Question**: Is step 4 correct?

**Check formula** (line 520-521, plot_biotsavart_fourier.py):
```python
Bx = BR_total * cosphi - Bphi_total * sinphi
By = BR_total * sinphi + Bphi_total * cosphi
```

**Expected**:
```
Bx = BR*cos(φ) - Bφ*sin(φ)  ✓
By = BR*sin(φ) + Bφ*cos(φ)  ✓
```

**This is correct!**

### Hypothesis 7: Reference Data is Wrong

**Fourier reference** from `single_reference.h5`:
- Computed by Fortran `GPEC` mode with `Fourier` output
- Directly computes Bnvac Fourier modes

**Question**: Does GPEC use the same (possibly buggy) transformation?

**YES** - Both reference and test use same Fortran code (coil_tools.f90)

So if there's a bug in Fortran, both will be wrong!

### Hypothesis 8: Analytic Formula Misapplied

**Analytic formula** (line 603, plot_biotsavart_fourier.py):
```python
B_analytic = mu0 * current * radius**2 / (2.0 * (radius**2 + s**2)**1.5)
```

This is for a **circular loop** with axis perpendicular to loop, evaluated along axis.

**For tilted coil**: Formula should still apply along the coil's axis

**But**: Formula gives scalar field magnitude, assumes field is purely axial

**Issue**: For slightly elliptical coil, there could be small perpendicular components

**However**: 60-80% error is too large to explain by this!

### Hypothesis 9: Physical vs Covariant Components Confusion

**Cylindrical coordinates** have two conventions:
- **Physical**: (BR, Bφ, BZ) - orthonormal basis, all same units
- **Covariant**: (BR, B_φ, BZ) where B_φ = R·Bφ

**Check**: Are we consistently using physical components?

**In Fortran** (line 469): `Bphi(kphi) = BXYZ(2) * cosphi - BXYZ(1) * sinphi`

This computes Bφ (physical), not B_φ (covariant)

**In Python** (line 182): `Bnphi = dAnR_dZ - dAnZ_dR`

For curl with physical Aφ: Bφ = ∂AR/∂Z - ∂AZ/∂R  ✓

**Conclusion**: All using physical components consistently

## Files to Check Systematically

### Fortran Source (coil_tools.f90)

| Line | Code | Check |
|------|------|-------|
| 422-423 | φ sampling | Correct periodicity? |
| 469 | `Bphi(kphi) = ...` | Sign of transformation? |
| 474-476 | FFT execution | Correct FFTW usage? |
| 477-479 | FFT normalization | Factor of nphi correct? |
| 596-598 | A transformation | Sign of Aφ formula? |
| 609-614 | A FFT and normalize | Same normalization? |

### Python Library (biotsavart_fourier.py)

| Line | Code | Check |
|------|------|-------|
| 95-96 | gauge_Anvac for ntor=0 | Returns unchanged? |
| 164-182 | field_divfree for ntor=0 | Curl formulas correct? |
| 162-163 | Derivative computation | Spline dx/dy correct? |
| 180-182 | Curl formula | Signs correct? |

### Python Script (plot_biotsavart_fourier.py)

| Line | Code | Check |
|------|------|-------|
| 498-499 | R, φ from x,y | arctan2 correct? |
| 508-511 | ntor=0 reconstruction | Just take .real? |
| 520-521 | Cyl → Cart | Signs correct? |
| 554 | Axis points | Direction correct? |
| 603 | Analytic formula | Applied correctly? |

## Specific Tests to Run

### Test 1: Verify Fortran B Field Transformation

**Create test input**:
- Known Cartesian field: Bx=1, By=0, Bz=0
- At φ=0: should give BR=1, Bφ=0, BZ=0
- At φ=π/2: should give BR=0, Bφ=1, BZ=0

**Check**: Does Fortran transformation give these values?

### Test 2: Check FFT Normalization

**Create test input**:
- Constant field: BR=1 everywhere in φ
- After FFT: B̃[0] should equal 1
- Higher modes should be ~0

**Check**: Does normalization give B̃[0]=1?

### Test 3: Manual n=0 Mode Calculation

**From coil geometry**:
- Compute B at one (R,Z) point for φ=0, 2π/192, 4π/192, ...
- Manually compute n=0 mode: B̃_0 = (1/192) Σ B(φ_k)
- Compare to Fortran output

### Test 4: Verify Python Curl

**From Anvac data**:
- Load AR, Aφ, AZ at one (R,Z) point
- Compute derivatives manually (finite difference)
- Apply curl formula manually
- Compare to field_divfree output

### Test 5: Check Cylindrical ↔ Cartesian Consistency

**Round-trip test**:
- Start with known cylindrical: BR=1, Bφ=0, BZ=0
- Convert to Cartesian at φ=π/4
- Convert back to cylindrical
- Should recover original values

## Debug Strategy

1. **Start with Fortran transformation** (Hypothesis 1)
   - Most likely place for sign error
   - Add debug output to coil_tools.f90 at line 469
   - Print BXYZ, BR, Bphi for first few φ points

2. **Check FFT** (Hypothesis 2)
   - Add debug output after line 477
   - Print Bn[0:5] for first (R,Z) point
   - Verify n=0 mode is reasonable

3. **Verify Python independently** (Hypothesis 4)
   - Use test data with known solution
   - Check each step of field_divfree

4. **Compare Fortran vs Python**
   - Load same Anvac data in Python
   - Compute Bφ manually
   - If matches Fortran → bug is in Fortran
   - If doesn't match → bug is in Python

## Expected Outcomes

### If Bug is Found

**Scenario A: Sign error in transformation**
- Fix sign in Fortran line 469 or 597
- Rerun all tests
- Expect Bφ,n=0 → ~0

**Scenario B: FFT normalization error**
- Fix normalization factor
- Expect magnitude change, not sign change

**Scenario C: Coordinate system mismatch**
- Clarify physical vs covariant
- Update formulas consistently

### If No Bug Found

**Then**: The field IS correct, and expectations are wrong

**Possible reasons**:
1. Coil geometry not as circular as thought
2. Gauge choice makes A non-axisymmetric even for axisymmetric B
3. Definition of n=0 mode is different than expected

## Action Items

- [ ] Add debug prints to Fortran coil_tools.f90
- [ ] Create synthetic test cases (constant fields, known transforms)
- [ ] Manually compute n=0 mode from coil geometry
- [ ] Cross-check Python curl calculation
- [ ] Verify coordinate transformations with round-trip tests
- [ ] Document findings in this file

## Notes

- All SymPy validations passed → formulas are mathematically correct
- All methods agree with each other → systematic issue, not random bug
- User says "128 segments is PLENTY" → high confidence that Bφ should be ~0
- Test is intentionally for tilted coil → coordinate transforms are critical
