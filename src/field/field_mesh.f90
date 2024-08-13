module neo_field_mesh
use, intrinsic :: iso_fortran_env, only: dp => real64
implicit none

type :: field_mesh_t
    real(dp), dimension(:), allocatable :: x1, x2, x3
    real(dp), dimension(:,:,:), allocatable :: A1, A2, A3, B1, B2, B3
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

    if (present(field)) then
        do i = 1, size(self%x1)
            do j = 1, size(self%x2)
                do k = 1, size(self%x3)
                    x = [self%x1(i), self%x2(j), self%x3(k)]
                    call field%compute_abfield(x, A, B)
                    self%A1(i,j,k) = A(1)
                    self%A2(i,j,k) = A(2)
                    self%A3(i,j,k) = A(3)
                    self%B1(i,j,k) = B(1)
                    self%B2(i,j,k) = B(2)
                    self%B3(i,j,k) = B(3)
                end do
            end do
        end do
    else
        self%A1 = 0.0_dp
        self%A2 = 0.0_dp
        self%A3 = 0.0_dp
        self%B1 = 0.0_dp
        self%B2 = 0.0_dp
        self%B3 = 0.0_dp
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

subroutine field_mesh_allocate(self, nx, ny, nz)
    class(field_mesh_t), intent(inout) :: self
    integer, intent(in) :: nx, ny, nz

    allocate(self%x1(nx))
    allocate(self%x2(ny))
    allocate(self%x3(nz))
    allocate(self%A1(nx,ny,nz))
    allocate(self%A2(nx,ny,nz))
    allocate(self%A3(nx,ny,nz))
    allocate(self%B1(nx,ny,nz))
    allocate(self%B2(nx,ny,nz))
    allocate(self%B3(nx,ny,nz))
end subroutine field_mesh_allocate

subroutine field_mesh_deallocate(self)
    class(field_mesh_t), intent(inout) :: self

    deallocate(self%x1)
    deallocate(self%x2)
    deallocate(self%x3)
    deallocate(self%A1)
    deallocate(self%A2)
    deallocate(self%A3)
    deallocate(self%B1)
    deallocate(self%B2)
    deallocate(self%B3)
end subroutine field_mesh_deallocate

end module neo_field_mesh