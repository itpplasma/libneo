program test_polylag_field
use, intrinsic :: iso_fortran_env, only: dp => real64
use util_for_test, only: print_test, print_ok, print_fail
implicit none

call test_polylag_field_init
call test_curla_equal_b

contains

subroutine test_polylag_field_init
    use neo_polylag_field, only: polylag_field_t

    real(dp), parameter :: tol = 1.0e-9_dp

    type(polylag_field_t) :: polylag_field

    call print_test("test_polylag_field_init")

    call polylag_field%polylag_field_init()

    call print_ok
end subroutine test_polylag_field_init


subroutine test_curla_equal_b
    use neo_polylag_field, only: polylag_field_t
    use util_for_test_field, only: compute_cartesian_curla

    real(dp), parameter :: tol = 1.0e-9_dp

    type(polylag_field_t) :: polylag_field
    real(dp) :: x(3), B(3), B_from_A(3)

    call print_test("test_curla_equal_b")

    x = [1.0_dp, 1.0_dp, 1.0_dp]

    call polylag_field%polylag_field_init()
    call polylag_field%compute_bfield(x, B)
    B_from_A = compute_cartesian_curla(polylag_field, x, tol)

    if (maxval(abs(B - B_from_A)) > tol) then
        print *, "curl A != B"
        call print_fail
        error stop
    end if

    call print_ok
end subroutine test_curla_equal_b

end program test_polylag_field
