program test_bessel_k_oracle
    ! Compares bessel_kn against gsl_sf_bessel_Kn over the full domain and
    ! benchmarks both. Link with -lgsl.
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use libneo_kinds, only: dp
    use neo_bessel_k, only: bessel_kn
    use util_for_test, only: print_test, print_ok, print_fail

    implicit none

    type, bind(c) :: gsl_sf_result_t
        real(c_double) :: val
        real(c_double) :: err
    end type gsl_sf_result_t

    interface
        function gsl_sf_bessel_kn(n, x) result(k) bind(c, name="gsl_sf_bessel_Kn")
            import :: c_int, c_double
            integer(c_int), value :: n
            real(c_double), value :: x
            real(c_double) :: k
        end function gsl_sf_bessel_kn

        function gsl_sf_bessel_kn_e(n, x, r) result(status) &
            bind(c, name="gsl_sf_bessel_Kn_e")
            import :: c_int, c_double, gsl_sf_result_t
            integer(c_int), value :: n
            real(c_double), value :: x
            type(gsl_sf_result_t), intent(out) :: r
            integer(c_int) :: status
        end function gsl_sf_bessel_kn_e

        function gsl_set_error_handler_off() result(old) &
            bind(c, name="gsl_set_error_handler_off")
            import :: c_ptr
            type(c_ptr) :: old
        end function gsl_set_error_handler_off
    end interface

    type(c_ptr) :: old_handler

    old_handler = gsl_set_error_handler_off()

    call test_against_gsl
    call benchmark

    contains

    subroutine test_against_gsl

        integer, parameter :: n_x = 1200
        integer, parameter :: n_list(15) = [0, 1, 2, 3, 4, 5, 6, 8, 10, 20, 50, &
                                            100, 150, 180, 200]
        real(dp), parameter :: rel_tol = 1.0d-13
        real(dp), parameter :: log_x_lo = -300.0d0
        real(dp), parameter :: log_x_hi = log10(700.0d0)

        type(gsl_sf_result_t) :: r
        real(dp) :: x, k_ours, rel_err, max_rel_err, tol
        integer :: i, j, n_compared, status

        call print_test("test_bessel_kn_against_gsl")

        max_rel_err = 0.0d0
        n_compared = 0
        do i = 1, n_x
            x = 10.0d0**(log_x_lo + (log_x_hi - log_x_lo)*(i - 1)/real(n_x - 1, dp))
            do j = 1, size(n_list)
                status = gsl_sf_bessel_kn_e(n_list(j), x, r)
                if (status /= 0) cycle
                ! GSL returns sub-normal junk where K_n overflows; skip those
                if (.not. ieee_is_finite(r%val)) cycle
                if (r%val < tiny(1.0d0)) cycle
                k_ours = bessel_kn(n_list(j), x)
                rel_err = abs(k_ours/r%val - 1.0d0)
                max_rel_err = max(max_rel_err, rel_err)
                n_compared = n_compared + 1
                ! allow GSL's own error estimate where it exceeds rel_tol
                tol = max(rel_tol, 2.0d0*r%err/r%val)
                if (rel_err > tol) then
                    print *, "n =", n_list(j), " x =", x
                    print *, "ours", k_ours, " gsl", r%val, " rel_err", rel_err
                    call print_fail
                    stop 1
                end if
            end do
        end do
        print *, "    compared", n_compared, "points, max rel err:", max_rel_err

        call print_ok
    end subroutine test_against_gsl


    subroutine benchmark

        integer, parameter :: n_calls = 2000000
        real(dp), parameter :: x_lo = 0.1d0, x_hi = 30.0d0

        real(dp) :: x, acc_ours, acc_gsl, t_ours, t_gsl
        integer :: i, n
        integer(8) :: count0, count1, count_rate

        call print_test("benchmark_bessel_kn_vs_gsl")

        acc_ours = 0.0d0
        call system_clock(count0, count_rate)
        do i = 1, n_calls
            x = x_lo + (x_hi - x_lo)*mod(i*0.6180339887498949d0, 1.0d0)
            n = mod(i, 11)
            acc_ours = acc_ours + bessel_kn(n, x)
        end do
        call system_clock(count1)
        t_ours = real(count1 - count0, dp)/real(count_rate, dp)

        acc_gsl = 0.0d0
        call system_clock(count0)
        do i = 1, n_calls
            x = x_lo + (x_hi - x_lo)*mod(i*0.6180339887498949d0, 1.0d0)
            n = mod(i, 11)
            acc_gsl = acc_gsl + gsl_sf_bessel_kn(n, x)
        end do
        call system_clock(count1)
        t_gsl = real(count1 - count0, dp)/real(count_rate, dp)

        print *, "    ours:", t_ours/n_calls*1.0d9, "ns/call (acc", acc_ours, ")"
        print *, "    gsl: ", t_gsl/n_calls*1.0d9, "ns/call (acc", acc_gsl, ")"

        if (abs(acc_ours/acc_gsl - 1.0d0) > 1.0d-12) then
            call print_fail
            stop 1
        end if

        call print_ok
    end subroutine benchmark

end program test_bessel_k_oracle
