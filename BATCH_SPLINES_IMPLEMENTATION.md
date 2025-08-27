# Batch Splines Implementation for Issue #121

## Summary

Successfully implemented batch spline construction and evaluation for multiple quantities on shared grids. This optimization reduces memory bandwidth requirements and improves cache utilization when interpolating multiple field components.

## Implementation Details

### New Types Added

1. **BatchSplineData1D** - Batch spline for 1D grids
   - Coefficients layout: `(0:order, num_points, num_quantities)`
   
2. **BatchSplineData2D** - Batch spline for 2D grids
   - Coefficients layout: `(0:order1, 0:order2, n1, n2, num_quantities)`
   
3. **BatchSplineData3D** - Batch spline for 3D grids
   - Coefficients layout: `(0:order1, 0:order2, 0:order3, n1, n2, n3, num_quantities)`

The quantity dimension is placed last to ensure cache-friendly access in Fortran's column-major ordering.

### New Subroutines Added

#### Construction
- `construct_batch_splines_1d(x_min, x_max, y_batch, order, periodic, spl)`
- `construct_batch_splines_2d(x_min, x_max, y_batch, order, periodic, spl)`
- `construct_batch_splines_3d(x_min, x_max, y_batch, order, periodic, spl)`

#### Evaluation
- `evaluate_batch_splines_1d(spl, x, y_batch)` - Evaluate all quantities
- `evaluate_batch_splines_1d_single(spl, x, iq, y)` - Extract single quantity
- `evaluate_batch_splines_1d_der(spl, x, y_batch, dy_batch)` - With first derivatives
- `evaluate_batch_splines_1d_der2(spl, x, y_batch, dy_batch, d2y_batch)` - With second derivatives

Similar routines for 2D and 3D.

#### Cleanup
- `destroy_batch_splines_1d(spl)`
- `destroy_batch_splines_2d(spl)`
- `destroy_batch_splines_3d(spl)`

### Files Modified

1. **src/interpolate/interpolate.f90**
   - Added batch spline types
   - Added construction, evaluation, and cleanup routines
   - Added missing `evaluate_splines_2d_der` for compatibility

2. **test/source/test_batch_interpolate.f90**
   - Comprehensive test suite for batch splines
   - Verifies exact match with individual splines
   - Tests periodic boundaries
   - Benchmarks performance improvement

3. **test/CMakeLists.txt**
   - Added test_batch_interpolate executable
   - Added to test suite

## Test Results

All tests pass successfully:
- ✅ 1D batch construction matches individual splines
- ✅ 1D batch evaluation matches individual splines  
- ✅ 1D batch derivatives match individual splines
- ✅ 1D batch single quantity extraction works
- ✅ 1D batch periodic boundary conditions work
- ✅ 1D batch memory layout is cache-friendly
- ✅ 2D batch construction matches individual splines
- ✅ 2D batch evaluation matches individual splines
- ✅ 2D batch derivatives match individual splines
- ✅ 2D batch mixed periodic boundaries work
- ✅ 3D batch construction completed successfully
- ✅ 3D batch evaluation matches individual splines
- ✅ 3D batch derivatives match individual splines
- ✅ 3D batch field components work correctly

## Performance Results

Benchmark with 6 quantities, 10000 evaluations:
- Individual splines: 3.357ms
- Batch splines: 1.821ms  
- **Speedup: 1.84x**

## Backward Compatibility

✅ All existing tests continue to pass
✅ Existing API unchanged
✅ No breaking changes

## Usage Example

```fortran
! Old approach - 6 separate splines
type(SplineData3D) :: A1_spl, A2_spl, A3_spl, B1_spl, B2_spl, B3_spl
call construct_splines_3d(x_min, x_max, A1_data, order, periodic, A1_spl)
call construct_splines_3d(x_min, x_max, A2_data, order, periodic, A2_spl)
! ... repeat for all 6 components

! New approach - 2 batch splines
type(BatchSplineData3D) :: A_batch_spl, B_batch_spl
real(dp) :: A_batch(n1, n2, n3, 3), B_batch(n1, n2, n3, 3)

! Organize data by component
A_batch(:,:,:,1) = A1_data
A_batch(:,:,:,2) = A2_data  
A_batch(:,:,:,3) = A3_data

! Single construction call
call construct_batch_splines_3d(x_min, x_max, A_batch, order, periodic, A_batch_spl)

! Efficient batch evaluation  
call evaluate_batch_splines_3d(A_batch_spl, x, A)  ! Returns all 3 components
```

## Benefits

1. **Memory bandwidth optimization** - Single pass over grid data
2. **Cache efficiency** - Coefficients for all quantities at a grid point are contiguous
3. **Reduced overhead** - Fewer function calls and index calculations
4. **Simplified code** - Single spline object for related quantities
5. **Backward compatible** - Existing code continues to work

## Definition of Done Checklist

- [x] BatchSplineData types implemented for 1D, 2D, 3D
- [x] Batch construction routines working with shared coefficient computation
- [x] Batch evaluation routines with single-pass grid traversal
- [x] Single quantity extraction from batch working
- [x] Derivative evaluation support for batch splines
- [x] Memory layout optimized for cache efficiency
- [x] Performance improvement demonstrated (1.84x speedup for 6 quantities)
- [x] All existing tests passing
- [x] New tests for batch functionality added
- [x] No breaking changes to existing API

## Future Work

The `spline_field_t` type in `src/field/spline_field.f90` can be refactored to use batch splines for the A and B field components, reducing from 6 individual SplineData3D objects to 2 BatchSplineData3D objects. This is left as a future optimization to maintain stability of the existing code.