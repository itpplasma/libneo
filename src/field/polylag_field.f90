module neo_polylag_field
use, intrinsic :: iso_fortran_env, only: dp => real64
use neo_field_base, only: field_t
use neo_field_mesh, only: field_mesh_t
implicit none

type, extends(field_t) :: polylag_field_t
    type(field_mesh_t) :: field_mesh
    contains
        procedure :: polylag_field_init
        procedure :: compute_abfield
        procedure :: compute_afield
        procedure :: compute_bfield
end type polylag_field_t

contains

subroutine polylag_field_init(self, limits, field, n_nodes)
    class(polylag_field_t), intent(out) :: self
    real(dp), dimension(3,2), intent(in):: limits
    class(field_t), intent(in), optional :: field
    integer, dimension(3), intent(in), optional :: n_nodes

    if (.not. present(field)) then
        call self%field_mesh%field_mesh_init_with_field(limits)
    elseif (.not. present(n_nodes)) then
        call self%field_mesh%field_mesh_init_with_field(limits, field)
    else
        call self%field_mesh%field_mesh_init_with_field(limits, field, n_nodes)
    end if
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