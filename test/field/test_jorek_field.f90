program test_jorek_field
use, intrinsic :: iso_fortran_env, only: dp => real64
use neo_jorek_field, only: jorek_field_t
use util_for_test, only: print_test, print_ok, print_fail
use util_for_test_jorek_field, only: get_filename, filename_len

implicit none

call test_get_ranges_from_filename
call test_jorek_field_init
call test_jorek_trial_field

contains


subroutine test_get_ranges_from_filename
    use neo_jorek_field, only: get_ranges_from_filename

    real(dp), parameter :: tol = 1.0e-10_dp
    real(dp), parameter :: expected_Rmin = 1.234_dp, expected_Rmax = 4.321_dp
    real(dp), parameter :: expected_Zmin = -1.234_dp, expected_Zmax = 1.234_dp
    real(dp), parameter :: expected_phimin = 1.234_dp, expected_phimax = 4.321_dp

    real(dp) :: Rmin, Rmax, Zmin, Zmax, phimin, phimax
    character(len=200) :: filename

    call print_test("test_get_ranges_from_filename")

    filename = "exprs_Rmin1.234_Rmax4.321_Zmin-1.234_Zmax1.234_phimin1.234_phimax4.321_s40000.h5"
    call get_ranges_from_filename(Rmin, Rmax, Zmin, Zmax, phimin, phimax, filename)
    if (is_not_same(Rmin, expected_Rmin, tol)) then
        print *, "Rmin = ", Rmin, "expected_Rmin = ", expected_Rmin
        call print_fail
        error stop
    end if
    if (is_not_same(Rmax, expected_Rmax, tol)) then
        print *, "Rmax = ", Rmax, "expected_Rmax = ", expected_Rmax
        call print_fail
        error stop
    end if
    if (is_not_same(Zmin, expected_Zmin, tol)) then
        print *, "Zmin = ", Zmin, "expected_Zmin = ", expected_Zmin
        call print_fail
        error stop
    end if
    if (is_not_same(Zmax, expected_Zmax, tol)) then
        print *, "Zmax = ", Zmax, "expected_Zmax = ", expected_Zmax
        call print_fail
        error stop
    end if
    if (is_not_same(phimin, expected_phimin, tol)) then
        print *, "phimin = ", phimin, "expected_phimin = ", expected_phimin
        call print_fail
        error stop
    end if
    if (is_not_same(phimax, expected_phimax, tol)) then
        print *, "phimax = ", phimax, "expected_phimax = ", expected_phimax
        call print_fail
        error stop
    end if

    call print_ok
end subroutine test_get_ranges_from_filename

function is_not_same(a, b, tol)
    real(dp), intent(in) :: a, b, tol
    logical :: is_not_same

    is_not_same = abs(a - b) > tol
end function is_not_same


subroutine test_jorek_field_init
    type(jorek_field_t) :: field
    character(len=filename_len) :: filename

    call print_test("test_jorek_field_init")

    call get_filename(filename)

    call field%jorek_field_init(filename)

    call print_ok
end subroutine test_jorek_field_init


subroutine test_jorek_trial_field
    use util_for_test_jorek_field, only: Rmin, Rmax, Zmin, Zmax, phimin, phimax
    type(jorek_field_t) :: field
    character(len=filename_len) :: trial_filename

    call print_test("test_trial_field")

    call get_filename(trial_filename)

    call field%jorek_field_init(trial_filename)

    call is_trial_field(field, Rmin, Rmax, Zmin, Zmax, phimin, phimax)
    !call is_curla_plus_fluxfunction_equal_b(field, Rmin, Rmax, Zmin, Zmax, phimin, phimax)

    call print_ok
end subroutine test_jorek_trial_field

subroutine is_trial_field(field, Rmin, Rmax, Zmin, Zmax, phimin, phimax)
    type(jorek_field_t), intent(in) :: field
    real(dp), intent(in) :: Rmin, Rmax, Zmin, Zmax, phimin, phimax

    integer, parameter :: n = 1000
    real(dp), parameter :: tol = 1.0e-10_dp
    real(dp) :: A_trial(3), A(3)
    real(dp) :: B_trial(3) = (/0.0_dp, 1.0_dp, 0.0_dp/), B(3)
    real(dp) :: fluxfunction_trial, fluxfunction
    real(dp) :: x(3,n), R
    integer :: idx

    x(1,:) = get_random_numbers(Rmin, Rmax, n)
    x(2,:) = get_random_numbers(phimin, phimax, n)
    x(3,:) = get_random_numbers(Zmin, Zmax, n)

    do idx = 1, n
        call field%compute_abfield(x(:,idx), A, B)
        R = x(1,idx)
        A_trial = (/0.0_dp, 0.0_dp, -0.5_dp/) * R
        call field%compute_fluxfunction(x(:,idx), fluxfunction)
        fluxfunction_trial = 0.5_dp * R
        if (any(abs(A - A_trial) > tol) .or. &
            any(abs(B - B_trial) > tol) .or. &
            abs(fluxfunction - fluxfunction_trial) > tol) then
            print *, "mis-match at x = ", x(:,idx)
            print *, "A = ", A, "B = ", B, "fluxfunction = ", fluxfunction
            print *, "A_trial = ", A_trial, "B_trial = ", B_trial, &
                     "fluxfunction_trial = ", fluxfunction_trial
            call print_fail
            error stop
        end if
    end do
end subroutine is_trial_field

subroutine is_curla_plus_fluxfunction_equal_b(field, &
                                              Rmin, Rmax, &
                                              phimin, phimax, &
                                              Zmin, Zmax)
    use util_for_test_field, only: compute_cylindrical_curla

    type(jorek_field_t), intent(in) :: field
    real(dp), intent(in) :: Rmin, Rmax, phimin, phimax, Zmin, Zmax

    integer, parameter :: n = 1000
    real(dp), parameter :: tol = 1.0e-10_dp
    real(dp) :: A(3), B(3), curla(3), fluxfunction, B_from_a_and_fluxfunction(3)
    real(dp) :: x(3,n), R
    integer :: idx

    x(1,:) = get_random_numbers(Rmin, Rmax, n)
    x(2,:) = get_random_numbers(Zmin, Zmax, n)
    x(3,:) = get_random_numbers(phimin, phimax, n)

    do idx = 1, n
        call field%compute_abfield(x(:,idx), A, B)
        call field%compute_fluxfunction(x(:,idx), fluxfunction)
        curla = compute_cylindrical_curla(field, x(:,idx), tol)
        B_from_a_and_fluxfunction = curla
        R = x(1,idx)
        B_from_a_and_fluxfunction(2) = B_from_a_and_fluxfunction(2) + fluxfunction/R
        if (any(abs(B - B_from_a_and_fluxfunction) > tol)) then
            print *, "mis-match at x = ", x(:,idx)
            print *, "B = ", B
            print *, "B_from_a_and_fluxfunction = ", B_from_a_and_fluxfunction
            print *, "curla = ", curla
            print *, "fluxfunction = ", fluxfunction
            call print_fail
            error stop
        end if
    end do
end subroutine is_curla_plus_fluxfunction_equal_b


function get_random_numbers(xmin, xmax, n, seed) result(x)
    real(dp), intent(in) :: xmin, xmax
    integer, intent(in) :: n
    integer, dimension(:), intent(in), optional :: seed
    real(dp), dimension(:), allocatable :: x

    if (present(seed)) then
        call random_seed(put=seed)
    end if
    allocate(x(n))
    call random_number(x)
    x = xmin + (xmax - xmin) * x
end function get_random_numbers


end program test_jorek_field