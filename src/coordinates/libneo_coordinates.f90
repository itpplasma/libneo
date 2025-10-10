module libneo_coordinates
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none

    type, abstract :: coordinate_system_t
    contains
        procedure(evaluate_point_if), deferred :: evaluate_point
        procedure(covariant_basis_if), deferred :: covariant_basis
        procedure(metric_tensor_if), deferred :: metric_tensor
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
    end interface

    interface
        module subroutine make_vmec_coordinate_system(cs)
            class(coordinate_system_t), allocatable, intent(out) :: cs
        end subroutine

        module subroutine make_geoflux_coordinate_system(cs)
            class(coordinate_system_t), allocatable, intent(out) :: cs
        end subroutine
    end interface

end module libneo_coordinates
