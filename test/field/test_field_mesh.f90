program test_field_mesh
use, intrinsic :: iso_fortran_env, only: dp => real64
use util_for_test, only: print_test, print_ok, print_fail
implicit none

call test_field_mesh_allocate
call test_field_mesh_init_with_field
call test_field_mesh_init_with_field_default

contains

subroutine test_field_mesh_allocate
    use neo_field_mesh, only: field_mesh_t

    type(field_mesh_t) :: field_mesh

    call print_test("test_field_mesh_allocate")

    call field_mesh%field_mesh_allocate(10, 10, 10)

    call print_ok
end subroutine test_field_mesh_allocate


subroutine test_field_mesh_init_with_field
    use neo_field_mesh, only: field_mesh_t
    use neo_field, only: create_field
    use neo_field_base, only: field_t

    real(dp), parameter :: tol = 1.0e-12_dp

    type(field_mesh_t) :: field_mesh
    class(field_t), allocatable :: field
    real(dp), dimension(3,2) :: limits
    integer, dimension(3) :: n_nodes
    real(dp), dimension(3) :: x, A, B
    integer, dimension(3) :: node

    call print_test("test_field_mesh_init_with_field")

    call create_field(field, "example")
    n_nodes = [10, 10, 10]
    limits(1,:) = [1.0_dp, 2.0_dp]
    limits(2,:) = [1.0_dp, 2.0_dp]
    limits(3,:) = [1.0_dp, 2.0_dp]
    call field_mesh%field_mesh_init_with_field(limits, field, n_nodes)

    node = [5, 5, 5]
    x = [field_mesh%x1(node(1)), &
         field_mesh%x2(node(2)), &
         field_mesh%x3(node(3))]
    call field%compute_abfield(x, A, B)
    if (any(abs(field_mesh%A(:,node(1),node(2),node(3)) - A) > tol)) then
        print *, "field_mesh%A = ", field_mesh%A(:,node(1),node(2),node(3))
        print *, "A = ", A
        call print_fail
        error stop
    end if
    if (any(abs(field_mesh%B(:,node(1),node(2),node(3)) - B) > tol)) then
        print *, "field_mesh%B = ", field_mesh%B(:,node(1),node(2),node(3))
        print *, "B = ", B
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

    if (any(abs(field_mesh%A - 0.0_dp) > tol)) then
        print *, "field_mesh%A =/= 0.0_dp"
        call print_fail
        error stop
    end if
    if (any(abs(field_mesh%B - 0.0_dp) > tol)) then
        print *, "field_mesh%B =/= 0.0_dp"
        call print_fail
        error stop
    end if

    call print_ok
end subroutine test_field_mesh_init_with_field_default

end program test_field_mesh
