module test_gauss_kronrod_callbacks

    use libneo_kinds, only: dp

    implicit none

    integer :: neval

contains

    function f_exp(x) result(fx)
        real(dp), intent(in) :: x
        real(dp) :: fx
        neval = neval + 1
        fx = exp(x)
    end function f_exp

    function f_sin50(x) result(fx)
        real(dp), intent(in) :: x
        real(dp) :: fx
        neval = neval + 1
        fx = sin(50.0d0*x)
    end function f_sin50

    function f_sqrt(x) result(fx)
        real(dp), intent(in) :: x
        real(dp) :: fx
        neval = neval + 1
        fx = sqrt(x)
    end function f_sqrt

    function f_lorentz(x) result(fx)
        real(dp), intent(in) :: x
        real(dp) :: fx
        neval = neval + 1
        fx = 1.0d0/(x*x + 1.0d-8)
    end function f_lorentz

    function f_peak(x) result(fx)
        real(dp), intent(in) :: x
        real(dp) :: fx
        neval = neval + 1
        fx = exp(-400.0d0*(x - 0.3d0)**2)
    end function f_peak

    function f_invsqrt(x) result(fx)
        real(dp), intent(in) :: x
        real(dp) :: fx
        neval = neval + 1
        fx = 1.0d0/sqrt(x)
    end function f_invsqrt

end module test_gauss_kronrod_callbacks

program test_gauss_kronrod

    use libneo_kinds, only: dp
    use neo_gauss_kronrod, only: integrate_gk
    use util_for_test, only: print_test, print_ok, print_fail
    use test_gauss_kronrod_callbacks

    implicit none

    ! references from mpmath.quad at 40 digits
    real(dp), parameter :: ref_exp = 1.71828182845904523536028747135266d0
    real(dp), parameter :: ref_sin50 = 7.0067943015773451862085882197966d-4
    real(dp), parameter :: ref_sqrt = 0.666666666666666666666666666666667d0
    real(dp), parameter :: ref_lorentz = 31413.9265359045990512531004997474d0
    real(dp), parameter :: ref_peak = 0.0886226925452758004113398690035196d0
    real(dp), parameter :: ref_invsqrt = 2.0d0
    real(dp), parameter :: epsrel = 1.0d-10
    integer, parameter :: keys(4) = [15, 21, 31, 61]

    call test_smooth_exp_all_keys
    call test_oscillatory_sin50_all_keys
    call test_sqrt_endpoint_singularity
    call test_near_singular_lorentzian
    call test_gaussian_peak
    call test_inverse_sqrt_singularity
    call test_default_key_and_limit
    call test_invalid_input
    call test_limit_exceeded

contains

    subroutine check_case(name, res, abserr, ierr, ref, epsa, epsr)
        character(*), intent(in) :: name
        real(dp), intent(in) :: res, abserr, ref, epsa, epsr
        integer, intent(in) :: ierr

        real(dp) :: err_true

        err_true = abs(res - ref)
        if (ierr /= 0) then
            call print_fail
            print *, name, ": ierr = ", ierr
            stop "integrate_gk reported failure"
        end if
        if (err_true > max(epsa, epsr*abs(ref))) then
            call print_fail
            print *, name, ": result ", res, " ref ", ref, " err ", err_true
            stop "result outside requested tolerance"
        end if
        if (abserr < err_true) then
            call print_fail
            print *, name, ": abserr ", abserr, " < true error ", err_true
            stop "abserr is not a true error bound"
        end if
        if (abserr > max(epsa, epsr*abs(res))) then
            call print_fail
            print *, name, ": abserr ", abserr, " exceeds requested tolerance"
            stop "claimed error exceeds tolerance despite ierr = 0"
        end if
    end subroutine check_case

    subroutine test_smooth_exp_all_keys
        real(dp) :: res, abserr
        integer :: ierr, k

        call print_test("test_smooth_exp_all_keys")
        do k = 1, size(keys)
            neval = 0
            call integrate_gk(f_exp, 0.0d0, 1.0d0, 0.0d0, epsrel, res, abserr, &
                              ierr, key=keys(k))
            call check_case("exp", res, abserr, ierr, ref_exp, 0.0d0, epsrel)
            print *, "    key ", keys(k), " neval ", neval
        end do
        call print_ok
    end subroutine test_smooth_exp_all_keys

    subroutine test_oscillatory_sin50_all_keys
        real(dp) :: res, abserr
        integer :: ierr, k

        call print_test("test_oscillatory_sin50_all_keys")
        do k = 1, size(keys)
            neval = 0
            call integrate_gk(f_sin50, 0.0d0, 1.0d0, 0.0d0, epsrel, res, abserr, &
                              ierr, key=keys(k))
            call check_case("sin50", res, abserr, ierr, ref_sin50, 0.0d0, epsrel)
            print *, "    key ", keys(k), " neval ", neval
        end do
        call print_ok
    end subroutine test_oscillatory_sin50_all_keys

    subroutine test_sqrt_endpoint_singularity
        real(dp) :: res, abserr
        integer :: ierr, k

        call print_test("test_sqrt_endpoint_singularity")
        do k = 1, size(keys)
            neval = 0
            call integrate_gk(f_sqrt, 0.0d0, 1.0d0, 0.0d0, epsrel, res, abserr, &
                              ierr, key=keys(k))
            call check_case("sqrt", res, abserr, ierr, ref_sqrt, 0.0d0, epsrel)
            print *, "    key ", keys(k), " neval ", neval
        end do
        call print_ok
    end subroutine test_sqrt_endpoint_singularity

    subroutine test_near_singular_lorentzian
        real(dp) :: res, abserr
        integer :: ierr, k

        call print_test("test_near_singular_lorentzian")
        do k = 1, size(keys)
            neval = 0
            call integrate_gk(f_lorentz, -1.0d0, 1.0d0, 0.0d0, epsrel, res, &
                              abserr, ierr, key=keys(k))
            call check_case("lorentz", res, abserr, ierr, ref_lorentz, 0.0d0, &
                            epsrel)
            print *, "    key ", keys(k), " neval ", neval
        end do
        call print_ok
    end subroutine test_near_singular_lorentzian

    subroutine test_gaussian_peak
        real(dp) :: res, abserr
        integer :: ierr, k

        call print_test("test_gaussian_peak")
        do k = 1, size(keys)
            neval = 0
            call integrate_gk(f_peak, 0.0d0, 1.0d0, 0.0d0, epsrel, res, abserr, &
                              ierr, key=keys(k))
            call check_case("peak", res, abserr, ierr, ref_peak, 0.0d0, epsrel)
            print *, "    key ", keys(k), " neval ", neval
        end do
        call print_ok
    end subroutine test_gaussian_peak

    subroutine test_inverse_sqrt_singularity
        real(dp) :: res, abserr
        integer :: ierr

        call print_test("test_inverse_sqrt_singularity")
        neval = 0
        call integrate_gk(f_invsqrt, 0.0d0, 1.0d0, 0.0d0, epsrel, res, abserr, &
                          ierr, key=21, limit=200)
        call check_case("invsqrt", res, abserr, ierr, ref_invsqrt, 0.0d0, epsrel)
        print *, "    key ", 21, " neval ", neval
        call print_ok
    end subroutine test_inverse_sqrt_singularity

    subroutine test_default_key_and_limit
        real(dp) :: res, abserr, res21, abserr21
        integer :: ierr

        call print_test("test_default_key_and_limit")
        call integrate_gk(f_sin50, 0.0d0, 1.0d0, 0.0d0, epsrel, res, abserr, ierr)
        call check_case("defaults", res, abserr, ierr, ref_sin50, 0.0d0, epsrel)
        call integrate_gk(f_sin50, 0.0d0, 1.0d0, 0.0d0, epsrel, res21, abserr21, &
                          ierr, key=21, limit=200)
        if (res /= res21 .or. abserr /= abserr21) then
            call print_fail
            stop "defaults do not match key=21, limit=200"
        end if
        call print_ok
    end subroutine test_default_key_and_limit

    subroutine test_invalid_input
        real(dp) :: res, abserr
        integer :: ierr

        call print_test("test_invalid_input")
        call integrate_gk(f_exp, 0.0d0, 1.0d0, 0.0d0, 1.0d-30, res, abserr, ierr)
        if (ierr /= 6) then
            call print_fail
            stop "expected ierr = 6 for unattainable tolerance"
        end if
        call integrate_gk(f_exp, 0.0d0, 1.0d0, 0.0d0, epsrel, res, abserr, ierr, &
                          key=17)
        if (ierr /= 6) then
            call print_fail
            stop "expected ierr = 6 for invalid key"
        end if
        call integrate_gk(f_exp, 0.0d0, 1.0d0, 0.0d0, epsrel, res, abserr, ierr, &
                          limit=0)
        if (ierr /= 6) then
            call print_fail
            stop "expected ierr = 6 for invalid limit"
        end if
        call print_ok
    end subroutine test_invalid_input

    subroutine test_limit_exceeded
        real(dp) :: res, abserr
        integer :: ierr

        call print_test("test_limit_exceeded")
        call integrate_gk(f_lorentz, -1.0d0, 1.0d0, 0.0d0, epsrel, res, abserr, &
                          ierr, key=15, limit=2)
        if (ierr /= 1) then
            call print_fail
            print *, "ierr = ", ierr
            stop "expected ierr = 1 when subdivision limit is too small"
        end if
        call print_ok
    end subroutine test_limit_exceeded

end program test_gauss_kronrod
