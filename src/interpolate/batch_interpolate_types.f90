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
      integer :: order = 0
      integer :: num_points = 0
      logical :: periodic = .false.
      real(dp) :: x_min = 0.0_dp
      real(dp) :: h_step = 0.0_dp
      real(dp) :: inv_h_step = 0.0_dp
      real(dp) :: period = 0.0_dp

      ! Batch data
      integer :: num_quantities = 0
      ! Memory layout: (n_quantities, 0:order, n_points) for cache efficiency
      real(dp), dimension(:, :, :), allocatable :: coeff
   end type BatchSplineData1D

   type :: BatchSplineData2D
      ! Shared grid data
      integer :: order(2) = 0
      integer :: num_points(2) = 0
      logical :: periodic(2) = .false.
      real(dp) :: h_step(2) = 0.0_dp
      real(dp) :: x_min(2) = 0.0_dp
      real(dp) :: inv_h_step(2) = 0.0_dp
      real(dp) :: period(2) = 0.0_dp

      ! Batch data
      integer :: num_quantities = 0
      ! Memory layout: (n_quantities, 0:order1, 0:order2, n1, n2) for cache efficiency
      real(dp), dimension(:, :, :, :, :), allocatable :: coeff
   end type BatchSplineData2D

   type :: BatchSplineData3D
      ! Shared grid data
      integer :: order(3) = 0
      integer :: num_points(3) = 0
      logical :: periodic(3) = .false.
      real(dp) :: h_step(3) = 0.0_dp
      real(dp) :: x_min(3) = 0.0_dp
      real(dp) :: inv_h_step(3) = 0.0_dp
      real(dp) :: period(3) = 0.0_dp

      ! Batch data
      integer :: num_quantities = 0
      ! Memory layout: (n_quantities, 0:order1, 0:order2, 0:order3, n1, n2, n3) for
      ! cache efficiency
      real(dp), dimension(:, :, :, :, :, :, :), allocatable :: coeff
   end type BatchSplineData3D

end module batch_interpolate_types
