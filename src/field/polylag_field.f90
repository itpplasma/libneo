module neo_polylag_field
use, intrinsic :: iso_fortran_env, only: dp => real64
use neo_field_base, only: field_t
implicit none

type, extends(field_t) :: polylag_field_t
    contains
        procedure :: polylag_field_init
        procedure :: compute_abfield
        procedure :: compute_afield
        procedure :: compute_bfield
end type polylag_field_t

contains

subroutine polylag_field_init(self)

    class(polylag_field_t), intent(out) :: self

end subroutine polylag_field_init

subroutine compute_abfield(self, x, A, B)
    class(polylag_field_t), intent(in) :: self
    real(dp), intent(in) :: x(3)
    real(dp), intent(out) :: A(3), B(3)

    call self%compute_afield(x, A)
    call self%compute_bfield(x, B)
end subroutine compute_abfield

subroutine compute_afield(self, x, A)
    class(polylag_field_t), intent(in) :: self
    real(dp), intent(in) :: x(3)
    real(dp), intent(out) :: A(3)

    A = 0.0_dp
end subroutine compute_afield

subroutine compute_bfield(self, x, B)
    class(polylag_field_t), intent(in) :: self
    real(dp), intent(in) :: x(3)
    real(dp), intent(out) :: B(3)

    B = 0.0_dp
end subroutine compute_bfield

end module neo_polylag_field