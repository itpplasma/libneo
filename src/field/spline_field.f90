module neo_spline_field
use, intrinsic :: iso_fortran_env, only: dp => real64
use neo_field_base, only: field_t
use neo_field_mesh, only: field_mesh_t
use interpolate, only: SplineData3D, construct_splines_3d, evaluate_splines_3d
implicit none

type, extends(field_t) :: spline_field_t
    type(field_mesh_t) :: field_mesh
    type(SplineData3D) :: A1_spline, A2_spline, A3_spline, &
                          B1_spline, B2_spline, B3_spline
    contains
        procedure :: spline_field_init
        procedure :: compute_abfield
        procedure :: compute_afield
        procedure :: compute_bfield
end type spline_field_t

contains

subroutine spline_field_init(self, limits, field, n_points, periodic)
    class(spline_field_t), intent(out) :: self
    real(dp), dimension(3,2), intent(in):: limits
    class(field_t), intent(in), optional :: field
    integer, dimension(3), intent(in), optional :: n_points
    logical, intent(in), optional :: periodic(3)

    call self%field_mesh%field_mesh_init_with_field(limits, field, n_points, periodic)
    call make_spline_from_mesh(self%field_mesh%A1, self%A1_spline)
    call make_spline_from_mesh(self%field_mesh%A2, self%A2_spline)
    call make_spline_from_mesh(self%field_mesh%A3, self%A3_spline)
    call make_spline_from_mesh(self%field_mesh%B1, self%B1_spline)
    call make_spline_from_mesh(self%field_mesh%B2, self%B2_spline)
    call make_spline_from_mesh(self%field_mesh%B3, self%B3_spline)
end subroutine spline_field_init

subroutine make_spline_from_mesh(mesh, spline, order_in)
    use neo_mesh, only: mesh_t

    class(mesh_t), intent(in) :: mesh
    type(SplineData3D), intent(out) :: spline
    integer, optional :: order_in(3)

    real(dp) :: x_min(3), x_max(3)
    integer :: order(3)
        
    if (present(order_in)) then
        order = order_in
    else
        order = 3
    end if
    x_min = [mesh%x1(1), mesh%x2(1), mesh%x3(1)]
    x_max = [mesh%x1(mesh%n1), mesh%x2(mesh%n2), mesh%x3(mesh%n3)]
    call construct_splines_3d(x_min, x_max, mesh%value, order, mesh%periodic, spline)
end subroutine make_spline_from_mesh

subroutine compute_abfield(self, x, A, B)
    class(spline_field_t), intent(in) :: self
    real(dp), intent(in) :: x(3)
    real(dp), intent(out) :: A(3), B(3)

    call self%compute_afield(x, A)
    call self%compute_bfield(x, B)
end subroutine compute_abfield

subroutine compute_afield(self, x, A)
    class(spline_field_t), intent(in) :: self
    real(dp), intent(in) :: x(3)
    real(dp), intent(out) :: A(3)

    call evaluate_splines_3d(self%A1_spline, x, A(1))
    call evaluate_splines_3d(self%A2_spline, x, A(2))
    call evaluate_splines_3d(self%A3_spline, x, A(3))
end subroutine compute_afield

subroutine compute_bfield(self, x, B)
    class(spline_field_t), intent(in) :: self
    real(dp), intent(in) :: x(3)
    real(dp), intent(out) :: B(3)

    call evaluate_splines_3d(self%B1_spline, x, B(1))
    call evaluate_splines_3d(self%B2_spline, x, B(2))
    call evaluate_splines_3d(self%B3_spline, x, B(3))
end subroutine compute_bfield

end module neo_spline_field