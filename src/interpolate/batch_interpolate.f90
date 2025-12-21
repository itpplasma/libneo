module batch_interpolate
    ! Wrapper module that re-exports all batch interpolation functionality
    ! The implementation has been split into separate modules for better maintainability

    use batch_interpolate_types, only: BatchSplineData1D, BatchSplineData2D, &
                                       BatchSplineData3D
    use batch_interpolate_1d, only: construct_batch_splines_1d, &
                                    construct_batch_splines_1d_lines, &
                                    construct_batch_splines_1d_resident, &
                                    construct_batch_splines_1d_resident_device, &
                                    destroy_batch_splines_1d, &
                                    evaluate_batch_splines_1d, &
                                    evaluate_batch_splines_1d_single, &
                                    evaluate_batch_splines_1d_many, &
                                    evaluate_batch_splines_1d_many_resident, &
                                    evaluate_batch_splines_1d_der, &
                                    evaluate_batch_splines_1d_der2, &
                                    evaluate_batch_splines_1d_der3
    use batch_interpolate_2d, only: construct_batch_splines_2d, &
                                    construct_batch_splines_2d_lines, &
                                    construct_batch_splines_2d_resident, &
                                    construct_batch_splines_2d_resident_device, &
                                    destroy_batch_splines_2d, &
                                    evaluate_batch_splines_2d, &
                                    evaluate_batch_splines_2d_der, &
                                    evaluate_batch_splines_2d_many, &
                                    evaluate_batch_splines_2d_many_resident
    use batch_interpolate_3d, only: construct_batch_splines_3d, &
                                    construct_batch_splines_3d_lines, &
                                    construct_batch_splines_3d_resident, &
                                    construct_batch_splines_3d_resident_device, &
                                    destroy_batch_splines_3d, &
                                    evaluate_batch_splines_3d, &
                                    evaluate_batch_splines_3d_der, &
                                    evaluate_batch_splines_3d_der2, &
                                    evaluate_batch_splines_3d_many, &
                                    evaluate_batch_splines_3d_many_resident

    implicit none
    
    ! Re-export types from batch_interpolate_types
    public :: BatchSplineData1D, BatchSplineData2D, BatchSplineData3D
    
    ! Re-export 1D routines from batch_interpolate_1d
    public :: construct_batch_splines_1d, destroy_batch_splines_1d
    public :: construct_batch_splines_1d_lines
    public :: construct_batch_splines_1d_resident
    public :: construct_batch_splines_1d_resident_device
    public :: evaluate_batch_splines_1d, evaluate_batch_splines_1d_single
    public :: evaluate_batch_splines_1d_many, evaluate_batch_splines_1d_many_resident
    public :: evaluate_batch_splines_1d_der, evaluate_batch_splines_1d_der2
    public :: evaluate_batch_splines_1d_der3
    
    ! Re-export 2D routines from batch_interpolate_2d
    public :: construct_batch_splines_2d, destroy_batch_splines_2d
    public :: construct_batch_splines_2d_lines
    public :: construct_batch_splines_2d_resident
    public :: construct_batch_splines_2d_resident_device
    public :: evaluate_batch_splines_2d, evaluate_batch_splines_2d_der
    public :: evaluate_batch_splines_2d_many, evaluate_batch_splines_2d_many_resident
    
    ! Re-export 3D routines from batch_interpolate_3d
    public :: construct_batch_splines_3d, destroy_batch_splines_3d
    public :: construct_batch_splines_3d_lines
    public :: construct_batch_splines_3d_resident
    public :: construct_batch_splines_3d_resident_device
    public :: evaluate_batch_splines_3d, evaluate_batch_splines_3d_der
    public :: evaluate_batch_splines_3d_der2
    public :: evaluate_batch_splines_3d_many, evaluate_batch_splines_3d_many_resident
    
end module batch_interpolate
