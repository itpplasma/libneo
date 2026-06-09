program test_dawson_oracle

    use, intrinsic :: iso_c_binding, only: c_double
    use libneo_kinds, only: dp
    use neo_dawson, only: dawson
    use util_for_test, only: print_test, print_ok, print_fail

    implicit none

    interface
        function gsl_sf_dawson(x) bind(c, name="gsl_sf_dawson") result(y)
            import :: c_double
            real(c_double), value :: x
            real(c_double) :: y
        end function gsl_sf_dawson
    end interface

    real(dp), parameter :: rel_tol = 1.0d-13

    call test_against_gsl_linear_sweep
    call test_against_gsl_decades
    call benchmark_vs_gsl

contains

    subroutine test_against_gsl_linear_sweep

        integer, parameter :: n = 60001
        real(dp) :: x, f, g, rel_err, worst
        integer :: i

        call print_test("test_against_gsl_linear_sweep")

        worst = 0.0d0
        do i = 1, n
            x = -15.0d0 + 30.0d0*real(i - 1, dp)/real(n - 1, dp)
            f = dawson(x)
            g = gsl_sf_dawson(x)
            if (g == 0.0d0) then
                rel_err = abs(f - g)
            else
                rel_err = abs(f - g)/abs(g)
            end if
            worst = max(worst, rel_err)
            if (rel_err > rel_tol) then
                call print_fail
                print *, "x = ", x, ", ours ", f, ", gsl ", g, ", rel err ", rel_err
                stop "dawson deviates from gsl_sf_dawson on linear sweep"
            end if
        end do
        print *, "    max rel err vs gsl: ", worst
        call print_ok

    end subroutine test_against_gsl_linear_sweep

    subroutine test_against_gsl_decades

        integer, parameter :: n = 601
        real(dp) :: x, f, g, rel_err, worst
        integer :: i

        call print_test("test_against_gsl_decades")

        worst = 0.0d0
        do i = 1, n
            x = 10.0d0**(-300 + i - 1)
            f = dawson(x)
            g = gsl_sf_dawson(x)
            rel_err = abs(f - g)/abs(g)
            worst = max(worst, rel_err)
            if (rel_err > rel_tol) then
                call print_fail
                print *, "x = ", x, ", ours ", f, ", gsl ", g, ", rel err ", rel_err
                stop "dawson deviates from gsl_sf_dawson on decade sweep"
            end if
        end do
        print *, "    max rel err vs gsl: ", worst
        call print_ok

    end subroutine test_against_gsl_decades

    subroutine benchmark_vs_gsl

        integer, parameter :: nx = 4096, nrep = 5000
        real(dp) :: xs(nx)
        real(dp) :: s_ours, s_gsl, t0, t1, ns_ours, ns_gsl
        integer :: i, r
        integer(8) :: count0, count1, count_rate

        call print_test("benchmark_vs_gsl")

        do i = 1, nx
            xs(i) = 12.5d0*real(i, dp)/real(nx, dp)
        end do

        s_ours = 0.0d0
        call system_clock(count0, count_rate)
        do r = 1, nrep
            do i = 1, nx
                s_ours = s_ours + dawson(xs(i))
            end do
            xs(1) = xs(1) + 1.0d-13
        end do
        call system_clock(count1)
        ns_ours = 1.0d9*real(count1 - count0, dp)/real(count_rate, dp) &
                  /real(nx, dp)/real(nrep, dp)

        s_gsl = 0.0d0
        call system_clock(count0, count_rate)
        do r = 1, nrep
            do i = 1, nx
                s_gsl = s_gsl + gsl_sf_dawson(xs(i))
            end do
            xs(1) = xs(1) + 1.0d-13
        end do
        call system_clock(count1)
        ns_gsl = 1.0d9*real(count1 - count0, dp)/real(count_rate, dp) &
                 /real(nx, dp)/real(nrep, dp)

        print *, "    calls each: ", nx*nrep
        print *, "    ours ns/call: ", ns_ours, " (checksum ", s_ours, ")"
        print *, "    gsl  ns/call: ", ns_gsl, " (checksum ", s_gsl, ")"
        call print_ok

    end subroutine benchmark_vs_gsl

end program test_dawson_oracle
