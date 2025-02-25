program test_polylag_field
use libneo_kinds, only : dp
use util_for_test, only: print_test, print_ok, print_fail
implicit none

call test_polylag_field_init
call test_against_original_field
call test_curla_equal_b

contains

subroutine test_polylag_field_init
    use neo_polylag_field, only: polylag_field_t

    type(polylag_field_t) :: polylag_field
    real(dp), dimension(3,2) :: limits

    call print_test("test_polylag_field_init")

    limits(1,:) = [1.0_dp, 2.0_dp]
    limits(2,:) = [1.0_dp, 2.0_dp]
    limits(3,:) = [1.0_dp, 2.0_dp]
    call polylag_field%polylag_field_init(limits)

    call print_ok
end subroutine test_polylag_field_init


subroutine test_against_original_field
    use neo_polylag_field, only: polylag_field_t
    use neo_field, only: create_field
    use neo_field_base, only: field_t

    real(dp), parameter :: tol = 1.0e-9_dp

    class(field_t), allocatable :: field
    type(polylag_field_t) :: polylag_field
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

    call polylag_field%polylag_field_init(limits, field, n_points)
    call polylag_field%compute_bfield(x, B)
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
    use neo_polylag_field, only: polylag_field_t
    use neo_field, only: create_field
    use neo_field_base, only: field_t
    use util_for_test_field, only: compute_cartesian_curla

    real(dp), parameter :: tol = 1.0e-9_dp

    class(field_t), allocatable :: field
    type(polylag_field_t) :: polylag_field
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

    call polylag_field%polylag_field_init(limits, field, n_points)
    call polylag_field%compute_bfield(x, B)
    curl_A = compute_cartesian_curla(polylag_field, x, tol)

    if (maxval(abs(B - curl_A)) > tol) then
        print *, "curl A != B"
        print *, "B = ", B
        print *, "curl_A = ", curl_A
        call print_fail
        error stop
    end if

    call print_ok
end subroutine test_curla_equal_b

end program test_polylag_field
