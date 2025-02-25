program test_field_mesh
use libneo_kinds, only : dp
use util_for_test, only: print_test, print_ok, print_fail
implicit none

call test_field_mesh_init_with_field
call test_field_mesh_init_with_field_default

contains


subroutine test_field_mesh_init_with_field
    use neo_field_mesh, only: field_mesh_t
    use neo_field, only: create_field
    use neo_field_base, only: field_t

    real(dp), parameter :: tol = 1.0e-12_dp

    type(field_mesh_t) :: field_mesh
    class(field_t), allocatable :: field
    real(dp), dimension(3,2) :: limits
    integer, dimension(3) :: n_points
    real(dp), dimension(3) :: x, A, B
    integer, dimension(3) :: node

    call print_test("test_field_mesh_init_with_field")

    call create_field(field, "example")
    n_points = [10, 10, 10]
    limits(1,:) = [1.0_dp, 2.0_dp]
    limits(2,:) = [1.0_dp, 2.0_dp]
    limits(3,:) = [1.0_dp, 2.0_dp]
    call field_mesh%field_mesh_init_with_field(limits, field, n_points)

    node = [5, 5, 5]
    x = [field_mesh%A1%x1(node(1)), &
         field_mesh%A2%x2(node(2)), &
         field_mesh%A3%x3(node(3))]
    call field%compute_abfield(x, A, B)

    if (abs(field_mesh%A1%value(node(1),node(2),node(3)) - A(1)) > tol) then
        print *, "field_mesh%A1 = ", field_mesh%A1%value(node(1),node(2),node(3))
        print *, "A(1) = ", A(1)
        call print_fail
        error stop
    end if
    if (abs(field_mesh%A2%value(node(1),node(2),node(3)) - A(2)) > tol) then
        print *, "field_mesh%A2 = ", field_mesh%A2%value(node(1),node(2),node(3))
        print *, "A(2) = ", A(2)
        call print_fail
        error stop
    end if
    if (abs(field_mesh%A3%value(node(1),node(2),node(3)) - A(3)) > tol) then
        print *, "field_mesh%A3 = ", field_mesh%A3%value(node(1),node(2),node(3))
        print *, "A(3) = ", A(3)
        call print_fail
        error stop
    end if
    if (abs(field_mesh%B1%value(node(1),node(2),node(3)) - B(1)) > tol) then
        print *, "field_mesh%B1 = ", field_mesh%B1%value(node(1),node(2),node(3))
        print *, "B(1) = ", B(1)
        call print_fail
        error stop
    end if
    if (abs(field_mesh%B2%value(node(1),node(2),node(3)) - B(2)) > tol) then
        print *, "field_mesh%B2 = ", field_mesh%B2%value(node(1),node(2),node(3))
        print *, "B(2) = ", B(2)
        call print_fail
        error stop
    end if
    if (abs(field_mesh%B3%value(node(1),node(2),node(3)) - B(3)) > tol) then
        print *, "field_mesh%B3 = ", field_mesh%B3%value(node(1),node(2),node(3))
        print *, "B(3) = ", B(3)
        call print_fail
        error stop
    end if

    call print_ok
end subroutine test_field_mesh_init_with_field


subroutine test_field_mesh_init_with_field_default
    use neo_field_mesh, only: field_mesh_t

    real(dp), parameter :: tol = 1.0e-12_dp

    type(field_mesh_t) :: field_mesh
    real(dp), dimension(3,2) :: limits

    call print_test("test_field_mesh_init_with_field_default")

    limits(1,:) = [1.0_dp, 2.0_dp]
    limits(2,:) = [1.0_dp, 2.0_dp]
    limits(3,:) = [1.0_dp, 2.0_dp]
    call field_mesh%field_mesh_init_with_field(limits)

    if (any(abs(field_mesh%A1%value - 0.0_dp) > tol)) then
        print *, "field_mesh%A1 =/= 0.0_dp"
        call print_fail
        error stop
    end if
    if (any(abs(field_mesh%A2%value - 0.0_dp) > tol)) then
        print *, "field_mesh%A2 =/= 0.0_dp"
        call print_fail
        error stop
    end if
    if (any(abs(field_mesh%A3%value - 0.0_dp) > tol)) then
        print *, "field_mesh%A3 =/= 0.0_dp"
        call print_fail
        error stop
    end if
    if (any(abs(field_mesh%B1%value - 0.0_dp) > tol)) then
        print *, "field_mesh%B1 =/= 0.0_dp"
        call print_fail
        error stop
    end if
    if (any(abs(field_mesh%B2%value - 0.0_dp) > tol)) then
        print *, "field_mesh%B2 =/= 0.0_dp"
        call print_fail
        error stop
    end if
    if (any(abs(field_mesh%B3%value - 0.0_dp) > tol)) then
        print *, "field_mesh%B3 =/= 0.0_dp"
        call print_fail
        error stop
    end if

    call print_ok
end subroutine test_field_mesh_init_with_field_default

end program test_field_mesh
