program test_example_field
use, intrinsic :: iso_fortran_env, only: dp => real64
use test_util, only: print_test, print_ok, print_fail

implicit none


call test_example_field_init
call test_curla_equal_b
call test_divb_0
call test_compute_abfield


contains


subroutine test_example_field_init
    use libneo_example_field, only: example_field_t

    real(dp), parameter :: tol = 1.0e-9_dp

    type(example_field_t) :: example_field

    call print_test("test_example_field_init")

    call example_field%example_field_init(1.0_dp, 2.0_dp)

    if (abs(example_field%ampl - 1.0_dp) > tol) then
        call print_fail
        return
    end if
    if (abs(example_field%ampl2 - 2.0_dp) > tol) then
        call print_fail
        return
    end if

    call print_ok
end subroutine test_example_field_init


subroutine test_curla_equal_b
    use libneo_example_field, only: example_field_t
    use test_libneo_field_util, only: compute_cartesian_curla

    real(dp), parameter :: tol = 1.0e-9_dp

    type(example_field_t) :: example_field
    real(dp) :: x(3), B(3), B_from_A(3)

    call print_test("test_curla_equal_b")

    call example_field%example_field_init(1.0_dp, 2.0_dp)

    x = [1.0_dp, 1.0_dp, 1.0_dp]

    call example_field%compute_bfield(x, B)
    B_from_A = compute_cartesian_curla(example_field, x, tol)

    if (maxval(abs(B - B_from_A)) > tol) then
        print *, "curl A != B"
        call print_fail
        return
    end if

    call print_ok
end subroutine test_curla_equal_b


subroutine test_divb_0
    use libneo_example_field, only: example_field_t
    use test_libneo_field_util, only: compute_cartesian_divb

    real(dp), parameter :: tol = 1.0e-9_dp

    type(example_field_t) :: example_field
    real(dp) :: x(3), divb

    call print_test("test_divb_0")

    call example_field%example_field_init(1.0_dp, 2.0_dp)

    x = [1.0_dp, 1.0_dp, 1.0_dp]

    divb = compute_cartesian_divb(example_field, x, tol)

    if (abs(divb) > tol) then
        print *, "div B = ", divb
        call print_fail
        return
    end if

    call print_ok
end subroutine test_divb_0


subroutine test_compute_abfield
    use libneo_example_field, only: example_field_t

    real(dp), parameter :: tol = 1.0e-9_dp

    type(example_field_t) :: example_field
    real(dp) :: x(3), A(3), B(3), temp(3)

    call print_test("test_compute_abfield")

    call example_field%example_field_init(1.0_dp, 2.0_dp)

    x = [1.0_dp, 1.0_dp, 1.0_dp]

    call example_field%compute_abfield(x, A, B)
    call example_field%compute_afield(x, temp)
    if (maxval(abs(A - temp)) > tol) then
        print *, "A from afield != from abfield"
        call print_fail
        return
    end if
    call example_field%compute_bfield(x, temp)
    if (maxval(abs(B - temp)) > tol) then
        print *, "B from bfield != from abfield"
        call print_fail
        return
    end if

    call print_ok
end subroutine test_compute_abfield
    
end program test_example_field