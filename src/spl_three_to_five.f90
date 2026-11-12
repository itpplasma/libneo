module spl_three_to_five_sub
! Equidistant splines of order 3 to 5. The coefficient recurrences live once,
! in spline_build_lines, which the batch splines also use on host and device.

    use spline_build_lines, only: spl_three_reg_line, spl_three_per_line, &
                                  spl_four_reg_line, spl_four_per_line, &
                                  spl_five_reg_line, spl_five_per_line

    implicit none

contains

    recursive subroutine spl_five_reg(n, h, a, b, c, d, e, f)
        integer, parameter :: dp = kind(1.0d0)
        integer, intent(in) :: n
        real(dp), intent(in) :: h
        real(dp), dimension(n), intent(in) :: a
        real(dp), dimension(n), intent(out) :: b, c, d, e, f

        call spl_five_reg_line(n, h, a, b, c, d, e, f)
    end subroutine spl_five_reg

    recursive subroutine spl_five_per(n, h, a, b, c, d, e, f)
        ! First and last values of the function must be the same.
        integer, parameter :: dp = kind(1.0d0)
        integer, intent(in) :: n
        real(dp), intent(in) :: h
        real(dp), dimension(n), intent(in) :: a
        real(dp), dimension(n), intent(out) :: b, c, d, e, f

        call spl_five_per_line(n, h, a, b, c, d, e, f)
    end subroutine spl_five_per

    recursive subroutine spl_four_reg(n, h, a, b, c, d, e)
        integer, parameter :: dp = kind(1.0d0)
        integer, intent(in) :: n
        real(dp), intent(in) :: h
        real(dp), dimension(n), intent(in) :: a
        real(dp), dimension(n), intent(out) :: b, c, d, e

        call spl_four_reg_line(n, h, a, b, c, d, e)
    end subroutine spl_four_reg

    recursive subroutine spl_four_per(n, h, a, b, c, d, e)
        integer, parameter :: dp = kind(1.0d0)
        integer, intent(in) :: n
        real(dp), intent(in) :: h
        real(dp), dimension(n), intent(in) :: a
        real(dp), dimension(n), intent(out) :: b, c, d, e

        call spl_four_per_line(n, h, a, b, c, d, e)
    end subroutine spl_four_per

    recursive subroutine splreg(n, h, y, bi, ci, di)
        integer, parameter :: dp = kind(1.0d0)
        integer, intent(in) :: n
        real(dp), intent(in) :: h
        real(dp), dimension(n), intent(in) :: y
        real(dp), dimension(n), intent(out) :: bi, ci, di

        call spl_three_reg_line(n, h, y, bi, ci, di)
    end subroutine splreg

    recursive subroutine splper(n, h, y, bi, ci, di)
        integer, parameter :: dp = kind(1.0d0)
        integer, intent(in) :: n
        real(dp), intent(in) :: h
        real(dp), dimension(n), intent(in) :: y
        real(dp), dimension(n), intent(out) :: bi, ci, di

        call spl_three_per_line(n, h, y, bi, ci, di)
    end subroutine splper

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    recursive subroutine spl_reg(ns, n, h, splcoe)
        ! Generic routine for equidistant regular spline of order 3,4,5
        !
        ! Input/output parameters:
        ! ns              - spline order (input)
        ! n               - number of data points
        ! h               - abscissa step size
        ! splcoe(0:ns,n)  - spline coefficients; coefficient splcoe(0,:) are the input (function to spline)

        implicit none

        integer, parameter :: dp = kind(1.0d0)

        integer :: ns, n
        real(dp) :: h
        real(dp), dimension(0:ns, n) :: splcoe
        real(dp) :: a(n), b(n), c(n), d(n), e(n), f(n)

        if (ns .eq. 3) then
            a = splcoe(0, :)

            call splreg(n, h, a, b, c, d)

            splcoe(1, :) = b
            splcoe(2, :) = c
            splcoe(3, :) = d
        elseif (ns .eq. 4) then
            a = splcoe(0, :)

            call spl_four_reg(n, h, a, b, c, d, e)

            splcoe(1, :) = b
            splcoe(2, :) = c
            splcoe(3, :) = d
            splcoe(4, :) = e
        elseif (ns .eq. 5) then
            a = splcoe(0, :)

            call spl_five_reg(n, h, a, b, c, d, e, f)

            splcoe(1, :) = b
            splcoe(2, :) = c
            splcoe(3, :) = d
            splcoe(4, :) = e
            splcoe(5, :) = f
        else
            print *, 'wrong spline order'
        end if

    end subroutine spl_reg

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    recursive subroutine spl_per(ns, n, h, splcoe)
        ! Generic routine for equidistant periodic spline of order 3,4,5
        !
        ! Input/output parameters:
        ! ns              - spline order (input)
        ! n               - number of data points
        ! h               - abscissa step size
        ! splcoe(0:ns,n)  - spline coefficients; coefficient splcoe(0,:) are the input (function to spline)

        implicit none

        integer, parameter :: dp = kind(1.0d0)

        integer :: ns, n
        real(dp) :: h
        real(dp), dimension(0:ns, n) :: splcoe
        real(dp) :: a(n), b(n), c(n), d(n), e(n), f(n)

        if (ns .eq. 3) then
            a = splcoe(0, :)

            call splper(n, h, a, b, c, d)

            splcoe(1, :) = b
            splcoe(2, :) = c
            splcoe(3, :) = d
        elseif (ns .eq. 4) then
            a = splcoe(0, :)

            call spl_four_per(n, h, a, b, c, d, e)

            splcoe(1, :) = b
            splcoe(2, :) = c
            splcoe(3, :) = d
            splcoe(4, :) = e
        elseif (ns .eq. 5) then
            a = splcoe(0, :)

            call spl_five_per(n, h, a, b, c, d, e, f)

            splcoe(1, :) = b
            splcoe(2, :) = c
            splcoe(3, :) = d
            splcoe(4, :) = e
            splcoe(5, :) = f
        else
            print *, 'wrong spline order'
        end if

    end subroutine spl_per

end module spl_three_to_five_sub
