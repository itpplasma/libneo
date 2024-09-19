module neo_mesh
use, intrinsic :: iso_fortran_env, only: dp => real64
implicit none

type :: mesh_t
    real(dp), dimension(:), allocatable :: x1, x2, x3
    real(dp), dimension(:,:,:), allocatable :: value
    integer :: n1, n2, n3
    real(dp) :: dx1, dx2, dx3
    contains
    procedure :: mesh_init
    procedure :: mesh_deinit
end type mesh_t

contains

subroutine mesh_init(self, x1, x2, x3, value)

    class(mesh_t), intent(out) :: self
    real(dp), dimension(:), intent(in) :: x1, x2, x3
    real(dp), dimension(:,:,:), intent(in), optional :: value

    integer :: n1, n2, n3

    n1 = size(x1)
    n2 = size(x2)
    n3 = size(x3)
    allocate(self%x1(n1), self%x2(n2), self%x3(n3))
    allocate(self%value(n1,n2,n3))
    self%x1 = x1
    self%x2 = x2
    self%x3 = x3
    self%dx1 = self%x1(2) - self%x1(1)
    self%dx2 = self%x2(2) - self%x2(1)
    self%dx3 = self%x3(2) - self%x3(1)
    if (present(value)) then
        self%value = value
    else
        self%value = 0.0_dp
    end if
end subroutine mesh_init

subroutine mesh_deinit(self)
    class(mesh_t), intent(inout) :: self

    deallocate(self%x1, self%x2, self%x3)
    deallocate(self%value)
    self%n1 = 0
    self%n2 = 0
    self%n3 = 0
    self%dx1 = 0.0_dp
    self%dx2 = 0.0_dp
    self%dx3 = 0.0_dp
end subroutine mesh_deinit
    
end module neo_mesh