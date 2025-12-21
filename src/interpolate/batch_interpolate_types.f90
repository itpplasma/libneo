module batch_interpolate_types
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private
    
    ! Export batch spline types
    public :: BatchSplineData1D, BatchSplineData2D, BatchSplineData3D
    ! dp is for internal use only, not exported
    
    ! Batch spline types for multiple quantities on shared grid
    type :: BatchSplineData1D
        ! Shared grid data
        integer :: order
        integer :: num_points  
        logical :: periodic
        real(dp) :: x_min
        real(dp) :: h_step
        
        ! Batch data
        integer :: num_quantities
        ! Memory layout: (n_quantities, 0:order, n_points) for cache efficiency
        real(dp), dimension(:, :, :), allocatable :: coeff
    end type BatchSplineData1D
    
    type :: BatchSplineData2D
        ! Shared grid data
        integer :: order(2)
        integer :: num_points(2)
        logical :: periodic(2)
        real(dp) :: h_step(2)
        real(dp) :: x_min(2)
        
        ! Batch data
        integer :: num_quantities
        ! Memory layout: (n_quantities, 0:order1, 0:order2, n1, n2) for cache efficiency
        real(dp), dimension(:, :, :, :, :), allocatable :: coeff
    end type BatchSplineData2D
    
    type :: BatchSplineData3D
        ! Shared grid data
        integer :: order(3)
        integer :: num_points(3)
        logical :: periodic(3)
        real(dp) :: h_step(3)
        real(dp) :: x_min(3)
        
        ! Batch data
        integer :: num_quantities
        ! Memory layout: (n_quantities, 0:order1, 0:order2, 0:order3, n1, n2, n3) for cache efficiency
        real(dp), dimension(:, :, :, :, :, :, :), allocatable :: coeff
    end type BatchSplineData3D
    
end module batch_interpolate_types
