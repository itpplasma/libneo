module boozer_chartmap_types
    !> Plain data record for an extended Boozer chartmap. No I/O, no field
    !> machinery, so the boozer converter can build splines from it without
    !> pulling in the NetCDF chartmap reader.

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none
    private
    public :: boozer_chartmap_data_t

    type :: boozer_chartmap_data_t
        integer :: n_rho = 0
        integer :: n_s = 0
        integer :: n_theta = 0       !< internal field grid, endpoint-included
        integer :: n_phi = 0         !< internal field grid, endpoint-included
        integer :: nfp = 1
        integer :: signgs = -1       !< sign of the source (s,theta,zeta) Jacobian
        real(dp) :: torflux = 0.0_dp
        real(dp) :: rmajor = 0.0_dp  !< metres, derived from innermost-surface geometry
        real(dp) :: rho_min = 0.0_dp
        real(dp) :: rho_max = 0.0_dp
        real(dp) :: h_rho = 0.0_dp   !< uniform rho step
        real(dp) :: h_s = 0.0_dp
        real(dp) :: h_theta = 0.0_dp !< 2*pi/(n_theta-1)
        real(dp) :: h_phi = 0.0_dp   !< 2*pi/nfp/(n_phi-1)
        real(dp), allocatable :: rho(:)
        real(dp), allocatable :: s(:)
        real(dp), allocatable :: A_phi(:)
        real(dp), allocatable :: B_theta(:)
        real(dp), allocatable :: B_phi(:)
        real(dp), allocatable :: Bmod(:, :, :)  !< (n_rho, n_theta, n_phi)
    end type boozer_chartmap_data_t
end module boozer_chartmap_types
