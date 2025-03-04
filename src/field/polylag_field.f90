module neo_polylag_field
use neo_field_base, only: field_t
use neo_field_mesh, only: field_mesh_t
implicit none
integer, parameter :: dp = kind(1.0d0)

type, extends(field_t) :: polylag_field_t
    type(field_mesh_t) :: field_mesh
    contains
        procedure :: polylag_field_init
        procedure :: compute_abfield
        procedure :: compute_afield
        procedure :: compute_bfield
end type polylag_field_t

contains

subroutine polylag_field_init(self, limits, field, n_points)
    class(polylag_field_t), intent(out) :: self
    real(dp), dimension(3,2), intent(in):: limits
    class(field_t), intent(in), optional :: field
    integer, dimension(3), intent(in), optional :: n_points

    call self%field_mesh%field_mesh_init_with_field(limits, field, n_points)
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

    A(1) = eval_polylag_on_mesh(x, self%field_mesh%A1)
    A(2) = eval_polylag_on_mesh(x, self%field_mesh%A2)
    A(3) = eval_polylag_on_mesh(x, self%field_mesh%A3)
end subroutine compute_afield

subroutine compute_bfield(self, x, B)
    class(polylag_field_t), intent(in) :: self
    real(dp), intent(in) :: x(3)
    real(dp), intent(out) :: B(3)

    B(1) = eval_polylag_on_mesh(x, self%field_mesh%B1)
    B(2) = eval_polylag_on_mesh(x, self%field_mesh%B2)
    B(3) = eval_polylag_on_mesh(x, self%field_mesh%B3)
end subroutine compute_bfield

function eval_polylag_on_mesh(x, mesh) result(f)
    use neo_polylag_5, only: find_node_index, calc_coefs
    use neo_mesh, only: mesh_t
    integer, parameter :: n_node=6

    real(dp), intent(in) :: x(3)
    type(mesh_t), intent(in) :: mesh
    real(dp) :: f

    integer, dimension(n_node) :: node_idx1, node_idx2, node_idx3
    real(dp), dimension(n_node) :: x1_coefs, x2_coefs, x3_coefs
    integer :: i, j, k
    real(dp) :: fp

    call find_node_index(x(1), mesh%x1(1), mesh%dx1, &
                         mesh%n1, node_idx1)
    call find_node_index(x(2), mesh%x2(1), mesh%dx2, &
                         mesh%n2, node_idx2)
    call find_node_index(x(3), mesh%x3(1), mesh%dx3, &
                         mesh%n3, node_idx3)
    x1_coefs = calc_coefs(x(1), mesh%x1(node_idx1), &
                           mesh%dx1)
    x2_coefs = calc_coefs(x(2), mesh%x2(node_idx2), &
                          mesh%dx2)
    x3_coefs = calc_coefs(x(3), mesh%x3(node_idx3), &
                          mesh%dx3)
    f = 0.0_dp
    do k = 1, n_node
        do j = 1, n_node
            do i = 1, n_node
                fp = mesh%value(node_idx1(i),node_idx2(j),node_idx3(k))
                f = f + fp * x1_coefs(i) * x2_coefs(j) * x3_coefs(k)
            end do
        end do
    end do
end function eval_polylag_on_mesh

end module neo_polylag_field
