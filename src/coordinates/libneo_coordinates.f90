module libneo_coordinates
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use interpolate, only: BatchSplineData3D
    implicit none

    type, abstract :: coordinate_system_t
    contains
        procedure(evaluate_point_if), deferred :: evaluate_point
        procedure(covariant_basis_if), deferred :: covariant_basis
        procedure(metric_tensor_if), deferred :: metric_tensor
        procedure(from_cyl_if), deferred :: from_cyl
    end type coordinate_system_t

    abstract interface
        subroutine evaluate_point_if(self, u, x)
            import :: coordinate_system_t, dp
            class(coordinate_system_t), intent(in) :: self
            real(dp), intent(in) :: u(3)
            real(dp), intent(out) :: x(3)
        end subroutine

        subroutine covariant_basis_if(self, u, e_cov)
            import :: coordinate_system_t, dp
            class(coordinate_system_t), intent(in) :: self
            real(dp), intent(in) :: u(3)
            real(dp), intent(out) :: e_cov(3,3)
        end subroutine

        subroutine metric_tensor_if(self, u, g, ginv, sqrtg)
            import :: coordinate_system_t, dp
            class(coordinate_system_t), intent(in) :: self
            real(dp), intent(in) :: u(3)
            real(dp), intent(out) :: g(3,3), ginv(3,3), sqrtg
        end subroutine

        subroutine from_cyl_if(self, xcyl, u, ierr)
            import :: coordinate_system_t, dp
            class(coordinate_system_t), intent(in) :: self
            real(dp), intent(in) :: xcyl(3)
            real(dp), intent(out) :: u(3)
            integer, intent(out) :: ierr
        end subroutine
    end interface

    interface
        module subroutine make_vmec_coordinate_system(cs)
            class(coordinate_system_t), allocatable, intent(out) :: cs
        end subroutine

        module subroutine make_geoflux_coordinate_system(cs)
            class(coordinate_system_t), allocatable, intent(out) :: cs
        end subroutine

        module subroutine make_chartmap_coordinate_system(cs, filename)
            class(coordinate_system_t), allocatable, intent(out) :: cs
            character(len=*), intent(in) :: filename
        end subroutine
    end interface

    type, extends(coordinate_system_t) :: chartmap_coordinate_system_t
        type(BatchSplineData3D) :: spl_xyz
        integer :: nrho = 0
        integer :: ntheta = 0
        integer :: nzeta = 0
        real(dp) :: tol_newton = 1.0e-12_dp
    contains
        procedure :: evaluate_point => chartmap_evaluate_point
        procedure :: covariant_basis => chartmap_covariant_basis
        procedure :: metric_tensor => chartmap_metric_tensor
        procedure :: from_cyl => chartmap_from_cyl
    end type chartmap_coordinate_system_t

    interface
        module subroutine chartmap_evaluate_point(self, u, x)
            class(chartmap_coordinate_system_t), intent(in) :: self
            real(dp), intent(in) :: u(3)
            real(dp), intent(out) :: x(3)
        end subroutine chartmap_evaluate_point

        module subroutine chartmap_covariant_basis(self, u, e_cov)
            class(chartmap_coordinate_system_t), intent(in) :: self
            real(dp), intent(in) :: u(3)
            real(dp), intent(out) :: e_cov(3,3)
        end subroutine chartmap_covariant_basis

        module subroutine chartmap_metric_tensor(self, u, g, ginv, sqrtg)
            class(chartmap_coordinate_system_t), intent(in) :: self
            real(dp), intent(in) :: u(3)
            real(dp), intent(out) :: g(3,3), ginv(3,3), sqrtg
        end subroutine chartmap_metric_tensor

        module subroutine chartmap_from_cyl(self, xcyl, u, ierr)
            class(chartmap_coordinate_system_t), intent(in) :: self
            real(dp), intent(in) :: xcyl(3)
            real(dp), intent(out) :: u(3)
            integer, intent(out) :: ierr
        end subroutine chartmap_from_cyl
    end interface

end module libneo_coordinates
