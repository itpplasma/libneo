module neo_field_mesh
use, intrinsic :: iso_fortran_env, only: dp => real64
use neo_mesh, only: mesh_t, mesh_init
implicit none

type :: field_mesh_t
    type(mesh_t) :: A1, A2, A3, B1, B2, B3
    contains
    procedure :: field_mesh_init_with_field
    procedure :: field_mesh_deinit
end type field_mesh_t

contains

subroutine field_mesh_init_with_field(self, limits, field, n_points)
    use neo_field_base, only: field_t

    class(field_mesh_t), intent(out) :: self
    real(dp), dimension(3,2), intent(in) :: limits
    class(field_t), intent(in), optional :: field
    integer, dimension(3), intent(in), optional :: n_points

    integer :: i, j, k
    real(dp), dimension(:), allocatable :: x1, x2, x3
    real(dp), dimension(3) :: x, A, B

    if (present(n_points)) then
        allocate(x1(n_points(1)), x2(n_points(2)), x3(n_points(3)))
        x1 = linspace(limits(1,1), limits(1,2), n_points(1))
        x2 = linspace(limits(2,1), limits(2,2), n_points(2))
        x3 = linspace(limits(3,1), limits(3,2), n_points(3))
    else
        allocate(x1(10), x2(10), x3(10))
    end if

    call self%A1%mesh_init(x1, x2, x3)
    call self%A2%mesh_init(x1, x2, x3)
    call self%A3%mesh_init(x1, x2, x3)
    call self%B1%mesh_init(x1, x2, x3)
    call self%B2%mesh_init(x1, x2, x3)
    call self%B3%mesh_init(x1, x2, x3)

    if (present(field)) then
        do i = 1, size(x1)
            do j = 1, size(x2)
                do k = 1, size(x3)
                    x = [x1(i), x2(j), x3(k)]
                    call field%compute_abfield(x, A, B)
                    self%A1%value(i,j,k) = A(1)
                    self%A2%value(i,j,k) = A(2)
                    self%A3%value(i,j,k) = A(3)
                    self%B1%value(i,j,k) = B(1)
                    self%B2%value(i,j,k) = B(2)
                    self%B3%value(i,j,k) = B(3)
                end do
            end do
        end do
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

subroutine field_mesh_deinit(self)
    class(field_mesh_t), intent(inout) :: self

    call self%A1%mesh_deinit
    call self%A2%mesh_deinit
    call self%A3%mesh_deinit
    call self%B1%mesh_deinit
    call self%B2%mesh_deinit
    call self%B3%mesh_deinit
end subroutine field_mesh_deinit

end module neo_field_mesh