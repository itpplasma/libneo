program test_mesh
use, intrinsic :: iso_fortran_env, only: dp => real64
use util_for_test, only: print_test, print_ok, print_fail
implicit none

call test_mesh_init
call test_mesh_init_with_default

contains


subroutine test_mesh_init
    use neo_mesh, only: mesh_t

    real(dp), parameter :: tol = 1.0e-12_dp
    integer, dimension(3), parameter :: n = [10, 20, 30]
    real(dp), dimension(3,2), parameter:: limits = transpose( &
                                                    reshape([1.0_dp, 2.0_dp, &
                                                             0.0_dp, 2.0_dp, &
                                                            -1.0_dp, 0.0_dp], [2, 3]))

    type(mesh_t) :: mesh
    real(dp), dimension(:), allocatable :: x1, x2, x3
    real(dp), dimension(:,:,:), allocatable :: values

    call print_test("test_mesh_init")

    allocate(x1(n(1)), x2(n(2)), x3(n(3)))
    x1 = linspace(limits(1,1), limits(1,2), n(1))
    x2 = linspace(limits(2,1), limits(2,2), n(2))
    x3 = linspace(limits(3,1), limits(3,2), n(3))
    allocate(values(n(1), n(2), n(3)))
    values = trial_func(x1, x2, x3)
    call mesh%mesh_init(x1, x2, x3, values)

    if (any(abs(mesh%x1 - x1) > tol)) then
        print *, "mesh%x1 =/= x1"
        call print_fail
        error stop
    end if
    if (any(abs(mesh%x2 - x2) > tol)) then
        print *, "mesh%x2 =/= x2"
        call print_fail
        error stop
    end if
    if (any(abs(mesh%x3 - x3) > tol)) then
        print *, "mesh%x3 =/= x3"
        call print_fail
        error stop
    end if
    if (any(abs(mesh%value - trial_func(x1, x2, x3)) > tol)) then
        print *, "mesh%value =/= trial_func"
        call print_fail
        error stop
    end if

    call print_ok
end subroutine test_mesh_init

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

function trial_func(x1, x2, x3) result(y)
    real(dp), dimension(:), intent(in) :: x1, x2, x3
    real(dp), dimension(size(x1), size(x2), size(x3)) :: y

    integer :: i, j, k

    do i = 1, size(x1)
        do j = 1, size(x2)
            do k = 1, size(x3)
                y(i,j,k) = sin(x1(i)) + cos(x2(j)) + tan(x3(k))
            end do
        end do
    end do
end function trial_func


subroutine test_mesh_init_with_default
    use neo_mesh, only: mesh_t

    real(dp), parameter :: tol = 1.0e-12_dp
    integer, dimension(3), parameter :: n = [10, 20, 30]
    real(dp), dimension(3,2), parameter:: limits = transpose( &
                                                    reshape([1.0_dp, 2.0_dp, &
                                                             0.0_dp, 2.0_dp, &
                                                            -1.0_dp, 0.0_dp], [2, 3]))
    type(mesh_t) :: mesh
    real(dp), dimension(:), allocatable :: x1, x2, x3
    integer, dimension(3) :: mesh_shape, value_shape

    call print_test("test_mesh_init")

    allocate(x1(n(1)), x2(n(2)), x3(n(3)))
    x1 = linspace(limits(1,1), limits(1,2), n(1))
    x2 = linspace(limits(2,1), limits(2,2), n(2))
    x3 = linspace(limits(3,1), limits(3,2), n(3))
    call mesh%mesh_init(x1, x2, x3)

    mesh_shape = [size(x1), size(x2), size(x3)]
    value_shape = shape(mesh%value)

    if (any(mesh_shape /= value_shape)) then
        print *, "shape(mesh%value) =/= shape of mesh"
        call print_fail
        error stop
    end if
    if (any(abs(mesh%value) > tol)) then
        print *, "mesh%value =/= 0"
        call print_fail
        error stop
    end if

    call print_ok
end subroutine test_mesh_init_with_default

end program test_mesh