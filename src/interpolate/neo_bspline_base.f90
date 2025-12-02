module neo_bspline_base
    !! Core B-spline types and basis function evaluation.
    !! No external dependencies beyond intrinsics.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    type, public :: bspline_1d
        integer :: degree
        integer :: n_ctrl
        real(dp), allocatable :: knots(:)
        real(dp) :: x_min, x_max
    end type bspline_1d

    type, public :: bspline_2d
        type(bspline_1d) :: sx, sy
    end type bspline_2d

    type, public :: bspline_3d
        type(bspline_1d) :: sx, sy, sz
    end type bspline_3d

    public :: bspline_1d_init_uniform, bspline_2d_init_uniform, bspline_3d_init_uniform
    public :: bspline_1d_eval, bspline_2d_eval, bspline_3d_eval
    public :: find_span, basis_funs

contains

    subroutine bspline_1d_init_uniform(spl, degree, n_ctrl, x_min, x_max)
        type(bspline_1d), intent(out) :: spl
        integer, intent(in) :: degree, n_ctrl
        real(dp), intent(in) :: x_min, x_max
        integer :: p, n, m, i
        real(dp) :: h

        p = degree
        n = n_ctrl
        m = n + p + 1
        spl%degree = p
        spl%n_ctrl = n
        spl%x_min = x_min
        spl%x_max = x_max
        allocate(spl%knots(m))
        h = (x_max - x_min)/real(n - p, dp)
        do i = 1, m
            if (i <= p + 1) then
                spl%knots(i) = x_min
            else if (i >= m - p) then
                spl%knots(i) = x_max
            else
                spl%knots(i) = x_min + real(i - (p + 1), dp)*h
            end if
        end do
    end subroutine bspline_1d_init_uniform


    subroutine bspline_2d_init_uniform(spl, degree, n_ctrl, x_min, x_max)
        type(bspline_2d), intent(out) :: spl
        integer, intent(in) :: degree(2), n_ctrl(2)
        real(dp), intent(in) :: x_min(2), x_max(2)

        call bspline_1d_init_uniform(spl%sx, degree(1), n_ctrl(1), x_min(1), &
            x_max(1))
        call bspline_1d_init_uniform(spl%sy, degree(2), n_ctrl(2), x_min(2), &
            x_max(2))
    end subroutine bspline_2d_init_uniform


    subroutine bspline_3d_init_uniform(spl, degree, n_ctrl, x_min, x_max)
        type(bspline_3d), intent(out) :: spl
        integer, intent(in) :: degree(3), n_ctrl(3)
        real(dp), intent(in) :: x_min(3), x_max(3)

        call bspline_1d_init_uniform(spl%sx, degree(1), n_ctrl(1), x_min(1), &
            x_max(1))
        call bspline_1d_init_uniform(spl%sy, degree(2), n_ctrl(2), x_min(2), &
            x_max(2))
        call bspline_1d_init_uniform(spl%sz, degree(3), n_ctrl(3), x_min(3), &
            x_max(3))
    end subroutine bspline_3d_init_uniform


    subroutine find_span(spl, x, span)
        type(bspline_1d), intent(in) :: spl
        real(dp), intent(in) :: x
        integer, intent(out) :: span
        integer :: low, high, mid, p, n, m
        real(dp) :: xx

        p = spl%degree
        n = spl%n_ctrl
        m = n + p + 1
        xx = min(max(x, spl%x_min), spl%x_max)
        if (xx >= spl%knots(m - p)) then
            span = n
            return
        end if
        low = p + 1
        high = n + 1
        do
            mid = (low + high)/2
            if (xx < spl%knots(mid)) then
                high = mid
            else if (xx >= spl%knots(mid + 1)) then
                low = mid
            else
                span = mid
                exit
            end if
        end do
    end subroutine find_span


    subroutine basis_funs(spl, span, x, N)
        type(bspline_1d), intent(in) :: spl
        integer, intent(in) :: span
        real(dp), intent(in) :: x
        real(dp), intent(out) :: N(0:)
        integer :: p, j, r
        real(dp) :: left(0:spl%degree), right(0:spl%degree), saved, temp

        p = spl%degree
        N(0) = 1.0_dp
        do j = 1, p
            left(j) = x - spl%knots(span + 1 - j)
            right(j) = spl%knots(span + j) - x
            saved = 0.0_dp
            do r = 0, j - 1
                temp = N(r)/(right(r + 1) + left(j - r))
                N(r) = saved + right(r + 1)*temp
                saved = left(j - r)*temp
            end do
            N(j) = saved
        end do
    end subroutine basis_funs


    subroutine bspline_1d_eval(spl, coeff, x, y)
        type(bspline_1d), intent(in) :: spl
        real(dp), intent(in) :: coeff(:)
        real(dp), intent(in) :: x
        real(dp), intent(out) :: y
        integer :: span, j, p
        real(dp) :: N(0:spl%degree)

        p = spl%degree
        call find_span(spl, x, span)
        call basis_funs(spl, span, x, N)
        y = 0.0_dp
        !$omp simd reduction(+:y)
        do j = 0, p
            y = y + N(j)*coeff(span - p + j)
        end do
    end subroutine bspline_1d_eval


    subroutine bspline_2d_eval(spl, coeff, x, y)
        type(bspline_2d), intent(in) :: spl
        real(dp), intent(in) :: coeff(:,:)
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: y
        integer :: spanx, spany, px, py, a, b, ix, iy
        real(dp) :: Nx(0:spl%sx%degree), Ny(0:spl%sy%degree)

        px = spl%sx%degree
        py = spl%sy%degree
        call find_span(spl%sx, x(1), spanx)
        call basis_funs(spl%sx, spanx, x(1), Nx)
        call find_span(spl%sy, x(2), spany)
        call basis_funs(spl%sy, spany, x(2), Ny)
        y = 0.0_dp
        do a = 0, px
            ix = spanx - px + a
            do b = 0, py
                iy = spany - py + b
                y = y + Nx(a)*Ny(b)*coeff(ix, iy)
            end do
        end do
    end subroutine bspline_2d_eval


    subroutine bspline_3d_eval(spl, coeff, x, y)
        type(bspline_3d), intent(in) :: spl
        real(dp), intent(in) :: coeff(:,:,:)
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: y
        integer :: spanx, spany, spanz, px, py, pz, a, b, c, ix, iy, iz
        real(dp) :: Nx(0:spl%sx%degree), Ny(0:spl%sy%degree)
        real(dp) :: Nz(0:spl%sz%degree)

        px = spl%sx%degree
        py = spl%sy%degree
        pz = spl%sz%degree
        call find_span(spl%sx, x(1), spanx)
        call basis_funs(spl%sx, spanx, x(1), Nx)
        call find_span(spl%sy, x(2), spany)
        call basis_funs(spl%sy, spany, x(2), Ny)
        call find_span(spl%sz, x(3), spanz)
        call basis_funs(spl%sz, spanz, x(3), Nz)
        y = 0.0_dp
        do a = 0, px
            ix = spanx - px + a
            do b = 0, py
                iy = spany - py + b
                do c = 0, pz
                    iz = spanz - pz + c
                    y = y + Nx(a)*Ny(b)*Nz(c)*coeff(ix, iy, iz)
                end do
            end do
        end do
    end subroutine bspline_3d_eval

end module neo_bspline_base
