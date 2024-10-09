program test_jorek_field
use, intrinsic :: iso_fortran_env, only: dp => real64
use neo_jorek_field, only: jorek_field_t
use util_for_test, only: print_test, print_ok, print_fail
use util_for_test_jorek_field, only: get_filename, filename_len

implicit none


call test_jorek_field_init
call test_jorek_trial_field


contains


subroutine test_jorek_field_init
    type(jorek_field_t) :: field
    character(len=filename_len) :: filename

    call print_test("test_jorek_field_init")

    call get_filename(filename)

    call field%jorek_field_init(filename)

    call print_ok
end subroutine test_jorek_field_init


subroutine test_jorek_trial_field
    type(jorek_field_t) :: field
    character(len=filename_len) :: trial_filename

    call print_test("test_trial_field")

    call get_filename(trial_filename)

    call field%jorek_field_init(trial_filename)

    call is_trial_field(field)

    call print_ok
end subroutine test_jorek_trial_field

subroutine is_trial_field(field)
    use util_for_test_jorek_field, only: Rmin, Rmax, Zmin, Zmax, phimin, phimax
    type(jorek_field_t), intent(in) :: field

    integer, parameter :: n = 1000
    real(dp), parameter :: tol = 1.0e-10_dp

    real(dp) :: A_trial(3), A(3)
    real(dp) :: B_trial(3) = (/0.0_dp, 1.0_dp, 0.0_dp/), B(3)
    real(dp) :: fluxfunction_trial, fluxfunction
    real(dp) :: x(3,n), R
    integer :: idx

    x(1,:) = get_random_numbers(Rmin, Rmax, n)
    x(2,:) = get_random_numbers(Zmin, Zmax, n)
    x(3,:) = get_random_numbers(phimin, phimax, n)

    do idx = 1, n
        call field%compute_abfield(x(:,idx), A, B)
        R = x(1,idx)
        A_trial = (/0.0_dp, 0.0_dp, 1.0_dp/) * R
        if (any(abs(A - A_trial) > tol) .or. any(abs(B - B_trial) > tol)) then
            print *, "mis-match at x = ", x(:,idx)
            print *, "A = ", A, "B = ", B
            print *, "A_trial = ", A_trial, "B_trial = ", B_trial
            call print_fail
            error stop
        end if
        call field%compute_fluxfunction(x(:,idx), fluxfunction)
        fluxfunction_trial = -1.0_dp * R
        if (abs(fluxfunction - fluxfunction_trial) > tol) then
            print *, "mis-match at x = ", x(:,idx)
            print *, "fluxfunction = ", fluxfunction, &
                     "fluxfunction_trial = ", fluxfunction_trial
            call print_fail
            error stop
        end if
    end do
end subroutine is_trial_field

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