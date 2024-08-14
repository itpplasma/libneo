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
    use neo_polylag_5, only: find_node_index, calc_coefs
    
    integer, parameter :: mp=6

    class(polylag_field_t), intent(in) :: self
    real(dp), intent(in) :: x(3)
    real(dp), intent(out) :: A(3)

    integer, dimension(mp) :: index1, index2, index3
    real(dp), dimension(mp) :: x1_coefs, x2_coefs, x3_coefs
    integer :: i, j, k
    real(dp), dimension(3) :: Ap
    
    call find_node_index(x(1), self%field_mesh%x1(1), self%field_mesh%dx1, &
                         self%field_mesh%n1, index1)
    call find_node_index(x(2), self%field_mesh%x2(1), self%field_mesh%dx2, &
                         self%field_mesh%n2, index2)
    call find_node_index(x(3), self%field_mesh%x3(1), self%field_mesh%dx3, &
                         self%field_mesh%n3, index3)
    x1_coefs = calc_coefs(x(1), self%field_mesh%x1(index1), &
                           self%field_mesh%dx1)
    x2_coefs = calc_coefs(x(2), self%field_mesh%x2(index2), &
                          self%field_mesh%dx2)
    x3_coefs = calc_coefs(x(3), self%field_mesh%x3(index3), &
                          self%field_mesh%dx3)
    A = 0.0_dp
    do k = 1, mp
        do j = 1, mp
            do i = 1, mp
                Ap = self%field_mesh%A(:,index1(i),index2(j),index3(k))
                A = A + Ap * x1_coefs(i) * x2_coefs(j) * x3_coefs(k)
            end do
        end do
    end do
end subroutine compute_afield

subroutine compute_bfield(self, x, B)
    use neo_polylag_5, only: find_node_index, calc_coefs
    integer, parameter :: mp=6

    class(polylag_field_t), intent(in) :: self
    real(dp), intent(in) :: x(3)
    real(dp), intent(out) :: B(3)

    integer, dimension(mp) :: index1, index2, index3
    real(dp), dimension(mp) :: x1_coefs, x2_coefs, x3_coefs
    integer :: i, j, k
    real(dp), dimension(3) :: Bp
    
    call find_node_index(x(1), self%field_mesh%x1(1), self%field_mesh%dx1, &
                         self%field_mesh%n1, index1)
    call find_node_index(x(2), self%field_mesh%x2(1), self%field_mesh%dx2, &
                         self%field_mesh%n2, index2)
    call find_node_index(x(3), self%field_mesh%x3(1), self%field_mesh%dx3, &
                         self%field_mesh%n3, index3)
    x1_coefs = calc_coefs(x(1), self%field_mesh%x1(index1), &
                           self%field_mesh%dx1)
    x2_coefs = calc_coefs(x(2), self%field_mesh%x2(index2), &
                          self%field_mesh%dx2)
    x3_coefs = calc_coefs(x(3), self%field_mesh%x3(index3), &
                          self%field_mesh%dx3)
    B = 0.0_dp
    do k = 1, mp
        do j = 1, mp
            do i = 1, mp
                Bp = self%field_mesh%B(:,index1(i),index2(j),index3(k))
                B = B + Bp * x1_coefs(i) * x2_coefs(j) * x3_coefs(k)
            end do
        end do
    end do
end subroutine compute_bfield

end module neo_polylag_field