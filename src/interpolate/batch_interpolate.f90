module batch_interpolate
    ! Wrapper module that re-exports all batch interpolation functionality
    ! The implementation has been split into separate modules for better maintainability
    
    use batch_interpolate_types
    use batch_interpolate_1d
    use batch_interpolate_2d
    use batch_interpolate_3d
    
    implicit none
    
    ! Re-export types from batch_interpolate_types
    public :: BatchSplineData1D, BatchSplineData2D, BatchSplineData3D
    
    ! Re-export 1D routines from batch_interpolate_1d
    public :: construct_batch_splines_1d, destroy_batch_splines_1d
    public :: evaluate_batch_splines_1d, evaluate_batch_splines_1d_single
    public :: evaluate_batch_splines_1d_der, evaluate_batch_splines_1d_der2
    
    ! Re-export 2D routines from batch_interpolate_2d
    public :: construct_batch_splines_2d, destroy_batch_splines_2d
    public :: evaluate_batch_splines_2d, evaluate_batch_splines_2d_der
    
    ! Re-export 3D routines from batch_interpolate_3d
    public :: construct_batch_splines_3d, destroy_batch_splines_3d
    public :: evaluate_batch_splines_3d, evaluate_batch_splines_3d_der
    public :: evaluate_batch_splines_3d_der2
    
end module batch_interpolate