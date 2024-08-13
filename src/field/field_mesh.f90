module neo_field_mesh
use, intrinsic :: iso_fortran_env, only: dp => real64
implicit none

type :: field_mesh_t
    real(dp), dimension(:), allocatable :: x1, x2, x3
    real(dp), dimension(:,:,:,:), allocatable :: A, B
    integer :: n1, n2, n3
    real(dp) :: dx1, dx2, dx3
    contains
    procedure :: field_mesh_init_with_field
    procedure :: field_mesh_allocate
    procedure :: field_mesh_deallocate
end type field_mesh_t

contains

subroutine field_mesh_init_with_field(self, limits, field, n_nodes)
    use neo_field_base, only: field_t

    class(field_mesh_t), intent(out) :: self
    real(dp), dimension(3,2), intent(in) :: limits
    class(field_t), intent(in), optional :: field
    integer, dimension(3), intent(in), optional :: n_nodes

    integer :: i, j, k
    real(dp), dimension(3) :: x, A, B

    if (present(n_nodes)) then
        call self%field_mesh_allocate(n_nodes(1), n_nodes(2), n_nodes(3))
    else
        call self%field_mesh_allocate(100, 100, 100)
    end if
    self%x1 = linspace(limits(1,1), limits(1,2), size(self%x1))
    self%x2 = linspace(limits(2,1), limits(2,2), size(self%x2))
    self%x3 = linspace(limits(3,1), limits(3,2), size(self%x3))
    self%dx1 = self%x1(2) - self%x1(1)
    self%dx2 = self%x2(2) - self%x2(1)
    self%dx3 = self%x3(2) - self%x3(1)

    if (present(field)) then
        do i = 1, size(self%x1)
            do j = 1, size(self%x2)
                do k = 1, size(self%x3)
                    x = [self%x1(i), self%x2(j), self%x3(k)]
                    call field%compute_abfield(x, A, B)
                    self%A(:,i,j,k) = A
                    self%B(:,i,j,k) = B
                end do
            end do
        end do
    else
        self%A = 0.0_dp
        self%B = 0.0_dp
    end if
end subroutine field_mesh_init_with_field

function linspace(start, stop, n) result(x)
    real(dp), intent(in) :: start, stop
    integer, intent(in) :: n
    real(dp) :: x(n)
    real(dp) :: dx
    integer :: i

    dx = (stop - start) / (n - 1)
    x(1) = start
    do i = 2, n
        x(i) = x(i-1) + dx
    end do
end function linspace

subroutine field_mesh_allocate(self, n1, n2, n3)
    class(field_mesh_t), intent(inout) :: self
    integer, intent(in) :: n1, n2, n3

    allocate(self%x1(n1))
    allocate(self%x2(n2))
    allocate(self%x3(n3))
    self%n1 = n1
    self%n2 = n2
    self%n3 = n3
    allocate(self%A(3,n1,n2,n3))
    allocate(self%B(3,n1,n2,n3))
end subroutine field_mesh_allocate

subroutine field_mesh_deallocate(self)
    class(field_mesh_t), intent(inout) :: self

    deallocate(self%x1)
    deallocate(self%x2)
    deallocate(self%x3)
    deallocate(self%A)
    deallocate(self%B)
    self%n1 = 0
    self%n2 = 0
    self%n3 = 0
    self%dx1 = 0.0_dp
    self%dx2 = 0.0_dp
    self%dx3 = 0.0_dp
end subroutine field_mesh_deallocate

end module neo_field_mesh