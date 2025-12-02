module neo_bspline
    !! Simple 1D B-spline basis and matrix-free LSQ via CGLS.
    !!
    !! Uses textbook Cox–de Boor recursion (Piegl & Tiller) for basis
    !! evaluation and a standard CGLS algorithm for least-squares fitting.
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none
    private

    type :: bspline_1d
        integer :: degree        !! Polynomial degree p
        integer :: n_ctrl        !! Number of basis functions / control points
        real(dp), allocatable :: knots(:)  !! Knot vector, size = n_ctrl+degree+1
        real(dp) :: x_min, x_max
    end type bspline_1d

    public :: bspline_1d
    public :: bspline_1d_init_uniform
    public :: bspline_1d_eval
    public :: bspline_1d_lsq_cgls

contains

    subroutine bspline_1d_init_uniform(spl, degree, n_ctrl, x_min, x_max)
        !! Initialise open-uniform 1D B-spline on [x_min, x_max].
        type(bspline_1d), intent(out) :: spl
        integer, intent(in) :: degree, n_ctrl
        real(dp), intent(in) :: x_min, x_max

        integer :: p, n, m, i
        real(dp) :: h, left, right

        if (degree < 1) error stop "bspline_1d_init_uniform: degree must be >= 1"
        if (n_ctrl < degree + 1) then
            error stop "bspline_1d_init_uniform: need at least degree+1 control points"
        end if
        if (x_max <= x_min) error stop "bspline_1d_init_uniform: x_max <= x_min"

        p = degree
        n = n_ctrl
        m = n + p + 1

        spl%degree = p
        spl%n_ctrl = n
        spl%x_min = x_min
        spl%x_max = x_max

        if (allocated(spl%knots)) deallocate(spl%knots)
        allocate(spl%knots(m))

        ! Open-uniform knot vector:
        ! First p+1 knots at x_min, last p+1 knots at x_max, interior equally spaced.
        h = (x_max - x_min) / real(n - p, dp)
        do i = 1, m
            if (i <= p + 1) then
                spl%knots(i) = x_min
            else if (i >= m - p) then
                spl%knots(i) = x_max
            else
                left = real(i - (p + 1), dp)
                spl%knots(i) = x_min + left*h
            end if
        end do
    end subroutine bspline_1d_init_uniform


    subroutine bspline_1d_eval(spl, coeff, x, y)
        !! Evaluate spline sum_j coeff(j) * N_j(x).
        type(bspline_1d), intent(in) :: spl
        real(dp), intent(in) :: coeff(:)
        real(dp), intent(in) :: x
        real(dp), intent(out) :: y

        integer :: span, j, p
        real(dp), dimension(:), allocatable :: N

        if (size(coeff) /= spl%n_ctrl) then
            error stop "bspline_1d_eval: coeff size mismatch"
        end if

        p = spl%degree
        allocate(N(0:p))

        call find_span(spl, x, span)
        call basis_funs(spl, span, x, N)

        y = 0.0_dp
        do j = 0, p
            y = y + N(j) * coeff(span - p + j)
        end do

        deallocate(N)
    end subroutine bspline_1d_eval


    subroutine bspline_1d_lsq_cgls(spl, x_data, f_data, coeff, max_iter, tol)
        !! Matrix-free CGLS for 1D B-spline LSQ:
        !!   min_c sum_i (S_c(x_i) - f_i)^2
        !!
        !! where S_c(x) is the spline with control points coeff(:).
        type(bspline_1d), intent(in) :: spl
        real(dp), intent(in) :: x_data(:)
        real(dp), intent(in) :: f_data(:)
        real(dp), intent(inout) :: coeff(:)
        integer, intent(in), optional :: max_iter
        real(dp), intent(in), optional :: tol

        integer :: n_data, n_ctrl
        integer :: k, kmax
        real(dp) :: atol, gamma, gamma_new, alpha, beta, denom
        real(dp) :: rhs_norm
        real(dp), allocatable :: r(:), s(:), p(:), q(:)

        n_data = size(x_data)
        if (n_data /= size(f_data)) then
            error stop "bspline_1d_lsq_cgls: x_data and f_data size mismatch"
        end if

        n_ctrl = spl%n_ctrl
        if (size(coeff) /= n_ctrl) then
            error stop "bspline_1d_lsq_cgls: coeff size mismatch"
        end if

        if (present(max_iter)) then
            kmax = max_iter
        else
            kmax = 200
        end if
        if (present(tol)) then
            atol = tol
        else
            atol = 1.0d-10
        end if

        allocate(r(n_data))
        allocate(q(n_data))
        allocate(s(n_ctrl))
        allocate(p(n_ctrl))

        ! Initial guess: coeff = 0
        coeff = 0.0_dp

        ! r = f - A*coeff  (A*coeff = 0 initially)
        r = f_data

        call apply_AT(spl, x_data, r, s)
        p = s
        gamma = dot_product(s, s)
        rhs_norm = sqrt(dot_product(f_data, f_data))

        if (rhs_norm == 0.0_dp) then
            coeff = 0.0_dp
            deallocate(r, q, s, p)
            return
        end if

        do k = 1, kmax
            call apply_A(spl, x_data, p, q)
            denom = dot_product(q, q)
            if (denom <= 0.0_dp) exit

            alpha = gamma/denom
            coeff = coeff + alpha*p
            r = r - alpha*q

            call apply_AT(spl, x_data, r, s)
            gamma_new = dot_product(s, s)

            if (gamma_new <= (atol*rhs_norm)**2) exit

            beta = gamma_new/gamma
            gamma = gamma_new
            p = s + beta*p
        end do

        deallocate(r, q, s, p)
    end subroutine bspline_1d_lsq_cgls


    !-----------------------------------------------------------------
    ! Internal helpers
    !-----------------------------------------------------------------

    subroutine find_span(spl, x, span)
        !! Find knot span index for x (0-based basis index will be span-degree...span).
        type(bspline_1d), intent(in) :: spl
        real(dp), intent(in) :: x
        integer, intent(out) :: span

        integer :: low, high, mid, p, n, m
        real(dp) :: xx

        p = spl%degree
        n = spl%n_ctrl
        m = n + p + 1

        ! Clamp x to [x_min, x_max]
        xx = min(max(x, spl%x_min), spl%x_max)

        ! Special cases at the ends
        if (xx >= spl%knots(m-p)) then
            span = n
            return
        end if

        low = p + 1
        high = n + 1

        do
            mid = (low + high)/2
            if (xx < spl%knots(mid)) then
                high = mid
            else if (xx >= spl%knots(mid+1)) then
                low = mid
            else
                span = mid
                exit
            end if
        end do
    end subroutine find_span


    subroutine basis_funs(spl, span, x, N)
        !! Compute non-zero B-spline basis functions N(0:p) at x.
        type(bspline_1d), intent(in) :: spl
        integer, intent(in) :: span
        real(dp), intent(in) :: x
        real(dp), intent(out) :: N(0:)

        integer :: p, j, r
        real(dp), allocatable :: left(:), right(:)
        real(dp) :: saved, temp

        p = spl%degree
        if (size(N) < p+1) then
            error stop "basis_funs: N has wrong size"
        end if

        allocate(left(0:p))
        allocate(right(0:p))

        N(0) = 1.0_dp
        do j = 1, p
            left(j) = x - spl%knots(span+1-j)
            right(j) = spl%knots(span+j) - x
            saved = 0.0_dp
            do r = 0, j-1
                temp = N(r)/(right(r+1) + left(j-r))
                N(r) = saved + right(r+1)*temp
                saved = left(j-r)*temp
            end do
            N(j) = saved
        end do

        deallocate(left, right)
    end subroutine basis_funs


    subroutine apply_A(spl, x_data, coeff, y)
        !! y = A * coeff  (spline evaluation at all x_data)
        type(bspline_1d), intent(in) :: spl
        real(dp), intent(in) :: x_data(:)
        real(dp), intent(in) :: coeff(:)
        real(dp), intent(out) :: y(:)

        integer :: n_data, i, span, p, j
        real(dp), allocatable :: N(:)

        n_data = size(x_data)
        if (size(y) /= n_data) then
            error stop "apply_A: y size mismatch"
        end if
        if (size(coeff) /= spl%n_ctrl) then
            error stop "apply_A: coeff size mismatch"
        end if

        p = spl%degree
        allocate(N(0:p))

        y = 0.0_dp
        do i = 1, n_data
            call find_span(spl, x_data(i), span)
            call basis_funs(spl, span, x_data(i), N)
            do j = 0, p
                y(i) = y(i) + N(j)*coeff(span - p + j)
            end do
        end do

        deallocate(N)
    end subroutine apply_A


    subroutine apply_AT(spl, x_data, r, g)
        !! g = A^T * r  (adjoint action, accumulating into control points)
        type(bspline_1d), intent(in) :: spl
        real(dp), intent(in) :: x_data(:)
        real(dp), intent(in) :: r(:)
        real(dp), intent(out) :: g(:)

        integer :: n_data, i, span, p, j
        real(dp), allocatable :: N(:)

        n_data = size(x_data)
        if (size(r) /= n_data) then
            error stop "apply_AT: r size mismatch"
        end if
        if (size(g) /= spl%n_ctrl) then
            error stop "apply_AT: g size mismatch"
        end if

        p = spl%degree
        allocate(N(0:p))

        g = 0.0_dp
        do i = 1, n_data
            call find_span(spl, x_data(i), span)
            call basis_funs(spl, span, x_data(i), N)
            do j = 0, p
                g(span - p + j) = g(span - p + j) + N(j)*r(i)
            end do
        end do

        deallocate(N)
    end subroutine apply_AT

end module neo_bspline

