program test_spl_three_to_five
    use spl_three_to_five_sub

    implicit none

    integer, parameter :: dp = kind(1.0d0)

    logical :: all_passed

    all_passed = .true.

    call test_splreg_linear(all_passed)
    call test_splreg_quadratic(all_passed)
    call test_splreg_cubic(all_passed)
    call test_splper_sin(all_passed)
    call test_spl_four_reg_quartic(all_passed)
    call test_spl_five_reg_quintic(all_passed)
    call test_spl_reg_order3(all_passed)
    call test_spl_reg_order4(all_passed)
    call test_spl_reg_order5(all_passed)
    call test_spl_per_order3(all_passed)
    call test_spl_per_order5(all_passed)
    call test_spline_interpolation_accuracy(all_passed)
    call test_spline_derivative_accuracy(all_passed)

    if (.not. all_passed) then
        error stop "One or more spl_three_to_five tests failed"
    end if

contains

    subroutine test_splreg_linear(passed)
        logical, intent(inout) :: passed

        real(dp), parameter :: tol = 1.0d-10
        integer, parameter :: n = 21
        real(dp) :: y(n), bi(n), ci(n), di(n)
        real(dp) :: h, x
        integer :: i

        h = 0.1d0
        do i = 1, n
            x = (i - 1) * h
            y(i) = 2.0d0 * x + 1.0d0
        end do

        call splreg(n, h, y, bi, ci, di)

        do i = 2, n-2
            if (abs(bi(i) - 2.0d0) > tol) then
                write(*,*) "FAIL: splreg_linear bi(", i, ") =", bi(i), &
                    " expected 2.0"
                passed = .false.
                return
            end if
            if (abs(ci(i)) > tol) then
                write(*,*) "FAIL: splreg_linear ci(", i, ") =", ci(i), " expected 0"
                passed = .false.
                return
            end if
            if (abs(di(i)) > tol) then
                write(*,*) "FAIL: splreg_linear di(", i, ") =", di(i), " expected 0"
                passed = .false.
                return
            end if
        end do

        write(*,*) "PASS: test_splreg_linear"
    end subroutine test_splreg_linear

    subroutine test_splreg_quadratic(passed)
        logical, intent(inout) :: passed

        real(dp), parameter :: tol = 1.0d-6
        integer, parameter :: n = 51
        real(dp) :: y(n), bi(n), ci(n), di(n)
        real(dp) :: h, x, x_eval, f_spline, f_exact
        integer :: i, idx

        h = 0.04d0
        do i = 1, n
            x = (i - 1) * h
            y(i) = x**2
        end do

        call splreg(n, h, y, bi, ci, di)

        do i = 1, 5
            x_eval = 0.3d0 + 0.3d0 * (i - 1)
            idx = int(x_eval / h) + 1
            if (idx > n - 1) idx = n - 1

            x = x_eval - (idx - 1) * h
            f_spline = y(idx) + bi(idx) * x + ci(idx) * x**2 + di(idx) * x**3
            f_exact = x_eval**2

            if (abs(f_spline - f_exact) > tol) then
                write(*,*) "FAIL: splreg_quadratic at x=", x_eval, &
                    " f_spline=", f_spline, " f_exact=", f_exact
                passed = .false.
                return
            end if
        end do

        write(*,*) "PASS: test_splreg_quadratic"
    end subroutine test_splreg_quadratic

    subroutine test_splreg_cubic(passed)
        logical, intent(inout) :: passed

        real(dp), parameter :: tol = 1.0d-8
        integer, parameter :: n = 51
        real(dp) :: y(n), bi(n), ci(n), di(n)
        real(dp) :: h, x, x_eval, f_spline, f_exact
        integer :: i, idx

        h = 0.05d0
        do i = 1, n
            x = (i - 1) * h
            y(i) = x**3 - x**2 + x
        end do

        call splreg(n, h, y, bi, ci, di)

        do i = 1, 5
            x_eval = 0.2d0 + 0.3d0 * (i - 1)
            idx = int(x_eval / h) + 1
            if (idx > n - 1) idx = n - 1

            x = x_eval - (idx - 1) * h
            f_spline = y(idx) + bi(idx) * x + ci(idx) * x**2 + di(idx) * x**3
            f_exact = x_eval**3 - x_eval**2 + x_eval

            if (abs(f_spline - f_exact) > tol) then
                write(*,*) "FAIL: splreg_cubic at x=", x_eval, &
                    " f_spline=", f_spline, " f_exact=", f_exact
                passed = .false.
                return
            end if
        end do

        write(*,*) "PASS: test_splreg_cubic"
    end subroutine test_splreg_cubic

    subroutine test_splper_sin(passed)
        logical, intent(inout) :: passed

        real(dp), parameter :: tol = 1.0d-6
        real(dp), parameter :: pi = 3.14159265358979323846d0
        integer, parameter :: n = 65
        real(dp) :: y(n), bi(n), ci(n), di(n)
        real(dp) :: h, x, x_eval, f_spline, f_exact
        integer :: i, idx

        h = 2.0d0 * pi / (n - 1)
        do i = 1, n
            x = (i - 1) * h
            y(i) = sin(x)
        end do

        call splper(n, h, y, bi, ci, di)

        do i = 1, 5
            x_eval = 0.5d0 + 1.0d0 * (i - 1)
            idx = int(x_eval / h) + 1
            if (idx > n - 1) idx = n - 1

            x = x_eval - (idx - 1) * h
            f_spline = y(idx) + bi(idx) * x + ci(idx) * x**2 + di(idx) * x**3
            f_exact = sin(x_eval)

            if (abs(f_spline - f_exact) > tol) then
                write(*,*) "FAIL: splper_sin at x=", x_eval, &
                    " f_spline=", f_spline, " f_exact=", f_exact
                passed = .false.
                return
            end if
        end do

        write(*,*) "PASS: test_splper_sin"
    end subroutine test_splper_sin

    subroutine test_spl_four_reg_quartic(passed)
        logical, intent(inout) :: passed

        real(dp), parameter :: tol = 1.0d-6
        integer, parameter :: n = 51
        real(dp) :: a(n), b(n), c(n), d(n), e(n)
        real(dp) :: h, x, x_eval, f_spline, f_exact
        integer :: i, idx

        h = 0.04d0
        do i = 1, n
            x = (i - 1) * h
            a(i) = x**4 - x**2
        end do

        call spl_four_reg(n, h, a, b, c, d, e)

        do i = 1, 3
            x_eval = 0.2d0 + 0.5d0 * (i - 1)
            idx = int(x_eval / h) + 1
            if (idx > n - 1) idx = n - 1

            x = x_eval - (idx - 1) * h
            f_spline = a(idx) + b(idx) * x + c(idx) * x**2 + d(idx) * x**3 &
                       + e(idx) * x**4
            f_exact = x_eval**4 - x_eval**2

            if (abs(f_spline - f_exact) > tol) then
                write(*,*) "FAIL: spl_four_reg_quartic at x=", x_eval, &
                    " f_spline=", f_spline, " f_exact=", f_exact
                passed = .false.
                return
            end if
        end do

        write(*,*) "PASS: test_spl_four_reg_quartic"
    end subroutine test_spl_four_reg_quartic

    subroutine test_spl_five_reg_quintic(passed)
        logical, intent(inout) :: passed

        real(dp), parameter :: tol = 1.0d-5
        integer, parameter :: n = 51
        real(dp) :: a(n), b(n), c(n), d(n), e(n), f(n)
        real(dp) :: h, x, x_eval, f_spline, f_exact
        integer :: i, idx

        h = 0.04d0
        do i = 1, n
            x = (i - 1) * h
            a(i) = x**5 - x**3 + x
        end do

        call spl_five_reg(n, h, a, b, c, d, e, f)

        do i = 1, 3
            x_eval = 0.2d0 + 0.5d0 * (i - 1)
            idx = int(x_eval / h) + 1
            if (idx > n - 1) idx = n - 1

            x = x_eval - (idx - 1) * h
            f_spline = a(idx) + b(idx) * x + c(idx) * x**2 + d(idx) * x**3 &
                       + e(idx) * x**4 + f(idx) * x**5
            f_exact = x_eval**5 - x_eval**3 + x_eval

            if (abs(f_spline - f_exact) > tol) then
                write(*,*) "FAIL: spl_five_reg_quintic at x=", x_eval, &
                    " f_spline=", f_spline, " f_exact=", f_exact
                passed = .false.
                return
            end if
        end do

        write(*,*) "PASS: test_spl_five_reg_quintic"
    end subroutine test_spl_five_reg_quintic

    subroutine test_spl_reg_order3(passed)
        logical, intent(inout) :: passed

        real(dp), parameter :: tol = 1.0d-6
        integer, parameter :: ns = 3
        integer, parameter :: n = 101
        real(dp) :: splcoe(0:ns, n)
        real(dp) :: h, x, x_eval, f_spline, f_exact
        integer :: i, idx

        h = 0.02d0
        do i = 1, n
            x = (i - 1) * h
            splcoe(0, i) = x**2
        end do

        call spl_reg(ns, n, h, splcoe)

        do i = 1, 5
            x_eval = 0.3d0 + 0.3d0 * (i - 1)
            idx = int(x_eval / h) + 1
            if (idx > n - 1) idx = n - 1

            x = x_eval - (idx - 1) * h
            f_spline = splcoe(0,idx) + splcoe(1,idx) * x + splcoe(2,idx) * x**2 &
                       + splcoe(3,idx) * x**3
            f_exact = x_eval**2

            if (abs(f_spline - f_exact) > tol) then
                write(*,*) "FAIL: spl_reg_order3 at x=", x_eval, &
                    " f_spline=", f_spline, " f_exact=", f_exact
                passed = .false.
                return
            end if
        end do

        write(*,*) "PASS: test_spl_reg_order3"
    end subroutine test_spl_reg_order3

    subroutine test_spl_reg_order4(passed)
        logical, intent(inout) :: passed

        real(dp), parameter :: tol = 1.0d-6
        integer, parameter :: ns = 4
        integer, parameter :: n = 41
        real(dp) :: splcoe(0:ns, n)
        real(dp) :: h, x, x_eval, f_spline, f_exact
        integer :: i, idx

        h = 0.05d0
        do i = 1, n
            x = (i - 1) * h
            splcoe(0, i) = x**3 - x
        end do

        call spl_reg(ns, n, h, splcoe)

        do i = 1, 5
            x_eval = 0.12d0 + 0.3d0 * (i - 1)
            idx = int(x_eval / h) + 1
            if (idx > n - 1) idx = n - 1

            x = x_eval - (idx - 1) * h
            f_spline = splcoe(0,idx) + splcoe(1,idx) * x + splcoe(2,idx) * x**2 &
                       + splcoe(3,idx) * x**3 + splcoe(4,idx) * x**4
            f_exact = x_eval**3 - x_eval

            if (abs(f_spline - f_exact) > tol) then
                write(*,*) "FAIL: spl_reg_order4 at x=", x_eval, &
                    " f_spline=", f_spline, " f_exact=", f_exact
                passed = .false.
                return
            end if
        end do

        write(*,*) "PASS: test_spl_reg_order4"
    end subroutine test_spl_reg_order4

    subroutine test_spl_reg_order5(passed)
        logical, intent(inout) :: passed

        real(dp), parameter :: tol = 1.0d-5
        integer, parameter :: ns = 5
        integer, parameter :: n = 51
        real(dp) :: splcoe(0:ns, n)
        real(dp) :: h, x, x_eval, f_spline, f_exact
        integer :: i, idx

        h = 0.04d0
        do i = 1, n
            x = (i - 1) * h
            splcoe(0, i) = x**4 - x**2 + 1.0d0
        end do

        call spl_reg(ns, n, h, splcoe)

        do i = 1, 5
            x_eval = 0.12d0 + 0.3d0 * (i - 1)
            idx = int(x_eval / h) + 1
            if (idx > n - 1) idx = n - 1

            x = x_eval - (idx - 1) * h
            f_spline = splcoe(0,idx) + splcoe(1,idx) * x + splcoe(2,idx) * x**2 &
                       + splcoe(3,idx) * x**3 + splcoe(4,idx) * x**4 &
                       + splcoe(5,idx) * x**5
            f_exact = x_eval**4 - x_eval**2 + 1.0d0

            if (abs(f_spline - f_exact) > tol) then
                write(*,*) "FAIL: spl_reg_order5 at x=", x_eval, &
                    " f_spline=", f_spline, " f_exact=", f_exact
                passed = .false.
                return
            end if
        end do

        write(*,*) "PASS: test_spl_reg_order5"
    end subroutine test_spl_reg_order5

    subroutine test_spl_per_order3(passed)
        logical, intent(inout) :: passed

        real(dp), parameter :: tol = 1.0d-6
        real(dp), parameter :: pi = 3.14159265358979323846d0
        integer, parameter :: ns = 3
        integer, parameter :: n = 65
        real(dp) :: splcoe(0:ns, n)
        real(dp) :: h, x, x_eval, f_spline, f_exact
        integer :: i, idx

        h = 2.0d0 * pi / (n - 1)
        do i = 1, n
            x = (i - 1) * h
            splcoe(0, i) = cos(x)
        end do

        call spl_per(ns, n, h, splcoe)

        do i = 1, 5
            x_eval = 0.5d0 + 1.0d0 * (i - 1)
            idx = int(x_eval / h) + 1
            if (idx > n - 1) idx = n - 1

            x = x_eval - (idx - 1) * h
            f_spline = splcoe(0,idx) + splcoe(1,idx) * x + splcoe(2,idx) * x**2 &
                       + splcoe(3,idx) * x**3
            f_exact = cos(x_eval)

            if (abs(f_spline - f_exact) > tol) then
                write(*,*) "FAIL: spl_per_order3 at x=", x_eval, &
                    " f_spline=", f_spline, " f_exact=", f_exact
                passed = .false.
                return
            end if
        end do

        write(*,*) "PASS: test_spl_per_order3"
    end subroutine test_spl_per_order3

    subroutine test_spl_per_order5(passed)
        logical, intent(inout) :: passed

        real(dp), parameter :: tol = 1.0d-6
        real(dp), parameter :: pi = 3.14159265358979323846d0
        integer, parameter :: ns = 5
        integer, parameter :: n = 65
        real(dp) :: splcoe(0:ns, n)
        real(dp) :: h, x, x_eval, f_spline, f_exact
        integer :: i, idx

        h = 2.0d0 * pi / (n - 1)
        do i = 1, n
            x = (i - 1) * h
            splcoe(0, i) = sin(2.0d0 * x)
        end do

        call spl_per(ns, n, h, splcoe)

        do i = 1, 5
            x_eval = 0.5d0 + 1.0d0 * (i - 1)
            idx = int(x_eval / h) + 1
            if (idx > n - 1) idx = n - 1

            x = x_eval - (idx - 1) * h
            f_spline = splcoe(0,idx) + splcoe(1,idx) * x + splcoe(2,idx) * x**2 &
                       + splcoe(3,idx) * x**3 + splcoe(4,idx) * x**4 &
                       + splcoe(5,idx) * x**5
            f_exact = sin(2.0d0 * x_eval)

            if (abs(f_spline - f_exact) > tol) then
                write(*,*) "FAIL: spl_per_order5 at x=", x_eval, &
                    " f_spline=", f_spline, " f_exact=", f_exact
                passed = .false.
                return
            end if
        end do

        write(*,*) "PASS: test_spl_per_order5"
    end subroutine test_spl_per_order5

    subroutine test_spline_interpolation_accuracy(passed)
        logical, intent(inout) :: passed

        real(dp), parameter :: tol = 1.0d-10
        integer, parameter :: ns = 3
        integer, parameter :: n = 101
        real(dp) :: splcoe(0:ns, n)
        real(dp) :: h, x
        integer :: i

        h = 0.02d0
        do i = 1, n
            x = (i - 1) * h
            splcoe(0, i) = exp(-x)
        end do

        call spl_reg(ns, n, h, splcoe)

        do i = 1, n
            x = (i - 1) * h
            if (abs(splcoe(0,i) - exp(-x)) > tol) then
                write(*,*) "FAIL: spline_interpolation_accuracy at node i=", i
                passed = .false.
                return
            end if
        end do

        write(*,*) "PASS: test_spline_interpolation_accuracy"
    end subroutine test_spline_interpolation_accuracy

    subroutine test_spline_derivative_accuracy(passed)
        logical, intent(inout) :: passed

        real(dp), parameter :: tol = 1.0d-4
        integer, parameter :: ns = 5
        integer, parameter :: n = 101
        real(dp) :: splcoe(0:ns, n)
        real(dp) :: h, x, df_spline, df_exact
        integer :: i

        h = 0.02d0
        do i = 1, n
            x = (i - 1) * h
            splcoe(0, i) = sin(x)
        end do

        call spl_reg(ns, n, h, splcoe)

        do i = 2, n-1
            x = (i - 1) * h
            df_spline = splcoe(1, i)
            df_exact = cos(x)

            if (abs(df_spline - df_exact) > tol) then
                write(*,*) "FAIL: spline_derivative_accuracy at node i=", i, &
                    " df_spline=", df_spline, " df_exact=", df_exact
                passed = .false.
                return
            end if
        end do

        write(*,*) "PASS: test_spline_derivative_accuracy"
    end subroutine test_spline_derivative_accuracy

end program test_spl_three_to_five
