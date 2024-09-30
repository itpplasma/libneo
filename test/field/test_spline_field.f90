program test_spline_field
use, intrinsic :: iso_fortran_env, only: dp => real64
use util_for_test, only: print_test, print_ok, print_fail
implicit none

call test_spline_field_init
call test_against_original_field
call test_curla_equal_b

contains

subroutine test_spline_field_init
    use neo_spline_field, only: spline_field_t
    use neo_field, only: create_field
    use neo_field_base, only: field_t
    use neo_field_mesh, only: field_mesh_t

    class(field_t), allocatable :: field
    type(field_mesh_t) :: field_mesh
    type(spline_field_t) :: spline_field
    real(dp), dimension(3,2) :: limits

    call print_test("test_spline_field_init")

    limits(1,:) = [1.0_dp, 2.0_dp]
    limits(2,:) = [1.0_dp, 2.0_dp]
    limits(3,:) = [1.0_dp, 2.0_dp]

    call create_field(field, "example")
    call field_mesh%field_mesh_init_with_field(limits, field)
    call spline_field%spline_field_init(field_mesh)

    call print_ok
end subroutine test_spline_field_init


subroutine test_against_original_field
    use neo_spline_field, only: spline_field_t
    use neo_field, only: create_field
    use neo_field_base, only: field_t
    use neo_field_mesh, only: field_mesh_t

    real(dp), parameter :: tol = 1.0e-7_dp

    class(field_t), allocatable :: field
    type(field_mesh_t) :: field_mesh
    type(spline_field_t) :: spline_field
    real(dp), dimension(3,2) :: limits
    integer, dimension(3) :: n_points
    real(dp) :: x(3), B(3), B_original(3)

    call print_test("test_against_original_field")

    call create_field(field, "example")

    x = [1.5_dp, 1.5_dp, 1.5_dp]

    limits(1,:) = [1.0_dp, 2.0_dp]
    limits(2,:) = [1.0_dp, 2.0_dp]
    limits(3,:) = [1.0_dp, 2.0_dp]
    n_points = int((limits(:,2) - limits(:,1))/tol**(1.0_dp/5.0_dp)) + 1

    call field_mesh%field_mesh_init_with_field(limits, field, n_points)
    call spline_field%spline_field_init(field_mesh)
    call spline_field%compute_bfield(x, B)
    call field%compute_bfield(x, B_original)

    if (maxval(abs(B - B_original)) > tol) then
        print *, "B != B_original"
        print *, "B = ", B
        print *, "B_original = ", B_original
        call print_fail
        error stop
    end if

    call print_ok
end subroutine test_against_original_field

subroutine test_curla_equal_b
    use neo_spline_field, only: spline_field_t
    use neo_field, only: create_field
    use neo_field_base, only: field_t
    use neo_field_mesh, only: field_mesh_t
    use util_for_test_field, only: compute_cartesian_curla

    real(dp), parameter :: tol = 1.0e-7_dp

    class(field_t), allocatable :: field
    type(field_mesh_t) :: field_mesh
    type(spline_field_t) :: spline_field
    real(dp), dimension(3,2) :: limits
    integer, dimension(3) :: n_points
    real(dp) :: x(3), B(3), curl_A(3)

    call print_test("test_curla_equal_b")

    call create_field(field, "example")

    x = [1.5_dp, 1.5_dp, 1.5_dp]

    limits(1,:) = [1.0_dp, 2.0_dp]
    limits(2,:) = [1.0_dp, 2.0_dp]
    limits(3,:) = [1.0_dp, 2.0_dp]
    n_points = int((limits(:,2) - limits(:,1))/tol**(1.0_dp/5.0_dp)) + 1

    call field_mesh%field_mesh_init_with_field(limits, field, n_points)
    call spline_field%spline_field_init(field_mesh)
    call spline_field%compute_bfield(x, B)
    curl_A = compute_cartesian_curla(spline_field, x, tol)

    if (maxval(abs(B - curl_A)) > tol) then
        print *, "curl A != B"
        print *, "B = ", B
        print *, "curl_A = ", curl_A
        call print_fail
        error stop
    end if

    call print_ok
end subroutine test_curla_equal_b

end program test_spline_field