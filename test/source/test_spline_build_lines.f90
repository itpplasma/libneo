program test_spline_build_lines
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use spline_build_lines, only: spl_three_per_line
    use spl_three_to_five_sub, only: spl_per

    implicit none

    integer, parameter :: sizes(3) = [17, 65, 129]
    integer :: icase, i, n
    real(dp), allocatable :: a(:), b(:), c(:), d(:), reference(:, :)
    real(dp) :: h, x, max_error

    max_error = 0.0_dp
    do icase = 1, size(sizes)
        n = sizes(icase)
        h = 2.0_dp*acos(-1.0_dp)/real(n - 1, dp)
        allocate (a(n), b(n), c(n), d(n), reference(0:3, n))
        do i = 1, n
            x = real(i - 1, dp)*h
            a(i) = cos(x) + 0.25_dp*sin(3.0_dp*x)
        end do
        reference = 0.0_dp
        reference(0, :) = a

        call spl_per(3, n, h, reference)
        call spl_three_per_line(n, h, a, b, c, d)

        max_error = max(max_error, maxval(abs(b - reference(1, :))))
        max_error = max(max_error, maxval(abs(c - reference(2, :))))
        max_error = max(max_error, maxval(abs(d - reference(3, :))))
        deallocate (a, b, c, d, reference)
    end do

    print '(a,es12.4)', 'maximum coefficient error = ', max_error
    if (max_error > 2.0e-12_dp) error stop 'in-place cubic periodic solver differs from spl_per'
end program test_spline_build_lines
