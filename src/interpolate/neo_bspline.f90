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

    type :: bspline_2d
        type(bspline_1d) :: sx
        type(bspline_1d) :: sy
    end type bspline_2d

    public :: bspline_1d
    public :: bspline_1d_init_uniform
    public :: bspline_1d_eval
    public :: bspline_1d_lsq_cgls
    public :: bspline_1d_lsq_cgls_batch

    public :: bspline_2d
    public :: bspline_2d_init_uniform
    public :: bspline_2d_eval
    public :: bspline_2d_lsq_cgls
    public :: bspline_2d_lsq_cgls_batch

    public :: apply_A
    public :: apply_AT
    public :: apply_A2D
    public :: apply_A2D_T
    public :: find_span
    public :: basis_funs

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


    subroutine bspline_1d_lsq_cgls_batch(spl, x_data, f_data, coeff, max_iter, tol)
        !! Batched matrix-free CGLS for 1D B-splines.
        !! Solves independent LSQ problems for multiple right-hand sides:
        !!   min_{c(:,k)} sum_i (S_{c(:,k)}(x_i) - f_i(k))^2
        !! for k = 1..n_rhs.
        type(bspline_1d), intent(in) :: spl
        real(dp), intent(in) :: x_data(:)
        real(dp), intent(in) :: f_data(:,:)
        real(dp), intent(inout) :: coeff(:,:)
        integer, intent(in), optional :: max_iter
        real(dp), intent(in), optional :: tol

        integer :: n_data, n_rhs, n_ctrl, k

        n_data = size(x_data)
        if (n_data /= size(f_data, 1)) then
            error stop "bspline_1d_lsq_cgls_batch: x_data and f_data size mismatch"
        end if

        n_rhs = size(f_data, 2)
        n_ctrl = spl%n_ctrl

        if (size(coeff, 1) /= n_ctrl .or. size(coeff, 2) /= n_rhs) then
            error stop "bspline_1d_lsq_cgls_batch: coeff shape mismatch"
        end if

        do k = 1, n_rhs
            if (present(max_iter) .and. present(tol)) then
                call bspline_1d_lsq_cgls(spl, x_data, f_data(:, k), coeff(:, k), &
                    max_iter=max_iter, tol=tol)
            elseif (present(max_iter)) then
                call bspline_1d_lsq_cgls(spl, x_data, f_data(:, k), coeff(:, k), &
                    max_iter=max_iter)
            elseif (present(tol)) then
                call bspline_1d_lsq_cgls(spl, x_data, f_data(:, k), coeff(:, k), &
                    tol=tol)
            else
                call bspline_1d_lsq_cgls(spl, x_data, f_data(:, k), coeff(:, k))
            end if
        end do
    end subroutine bspline_1d_lsq_cgls_batch


    !-----------------------------------------------------------------
    ! 2D tensor-product B-splines
    !-----------------------------------------------------------------

    subroutine bspline_2d_init_uniform(spl, degree, n_ctrl, x_min, x_max)
        !! Initialise 2D tensor-product B-spline with open-uniform knots
        !! in each dimension on [x_min(j), x_max(j)].
        type(bspline_2d), intent(out) :: spl
        integer, intent(in) :: degree(2), n_ctrl(2)
        real(dp), intent(in) :: x_min(2), x_max(2)

        call bspline_1d_init_uniform(spl%sx, degree(1), n_ctrl(1), x_min(1), x_max(1))
        call bspline_1d_init_uniform(spl%sy, degree(2), n_ctrl(2), x_min(2), x_max(2))
    end subroutine bspline_2d_init_uniform


    subroutine bspline_2d_eval(spl, coeff, x, y)
        !! Evaluate 2D spline S(x,y) = sum_{i,j} coeff(i,j) N_i(x) M_j(y).
        type(bspline_2d), intent(in) :: spl
        real(dp), intent(in) :: coeff(:,:)
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: y

        integer :: nx, ny, spanx, spany, px, py
        integer :: a, b, ix, iy
        real(dp), allocatable :: Nx_b(:), Ny_b(:)

        nx = spl%sx%n_ctrl
        ny = spl%sy%n_ctrl
        if (size(coeff, 1) /= nx .or. size(coeff, 2) /= ny) then
            error stop "bspline_2d_eval: coeff shape mismatch"
        end if

        px = spl%sx%degree
        py = spl%sy%degree
        allocate(Nx_b(0:px))
        allocate(Ny_b(0:py))

        call find_span(spl%sx, x(1), spanx)
        call basis_funs(spl%sx, spanx, x(1), Nx_b)
        call find_span(spl%sy, x(2), spany)
        call basis_funs(spl%sy, spany, x(2), Ny_b)

        y = 0.0_dp
        do a = 0, px
            ix = spanx - px + a
            do b = 0, py
                iy = spany - py + b
                y = y + Nx_b(a)*Ny_b(b)*coeff(ix, iy)
            end do
        end do

        deallocate(Nx_b, Ny_b)
    end subroutine bspline_2d_eval


    subroutine bspline_2d_lsq_cgls(spl, x_data, y_data, f_data, coeff, max_iter, tol)
        !! Matrix-free CGLS for 2D tensor-product B-spline LSQ:
        !!   min_C sum_i (S_C(x_i,y_i) - f_i)^2
        !!
        !! where C is coeff(nx,ny) and S_C is the spline.
        type(bspline_2d), intent(in) :: spl
        real(dp), intent(in) :: x_data(:)
        real(dp), intent(in) :: y_data(:)
        real(dp), intent(in) :: f_data(:)
        real(dp), intent(inout) :: coeff(:,:)
        integer, intent(in), optional :: max_iter
        real(dp), intent(in), optional :: tol

        integer :: n_data, nx, ny
        integer :: k, kmax
        real(dp) :: atol, gamma, gamma_new, alpha, beta, denom
        real(dp) :: rhs_norm
        real(dp), allocatable :: r(:), q(:)
        real(dp), allocatable :: s(:,:), p(:,:)

        n_data = size(x_data)
        if (n_data /= size(y_data) .or. n_data /= size(f_data)) then
            error stop "bspline_2d_lsq_cgls: data size mismatch"
        end if

        nx = spl%sx%n_ctrl
        ny = spl%sy%n_ctrl
        if (size(coeff, 1) /= nx .or. size(coeff, 2) /= ny) then
            error stop "bspline_2d_lsq_cgls: coeff shape mismatch"
        end if

        if (present(max_iter)) then
            kmax = max_iter
        else
            kmax = 400
        end if
        if (present(tol)) then
            atol = tol
        else
            atol = 1.0d-10
        end if

        allocate(r(n_data))
        allocate(q(n_data))
        allocate(s(nx, ny))
        allocate(p(nx, ny))

        coeff = 0.0_dp
        r = f_data

        call apply_A2D_T(spl, x_data, y_data, r, s)
        p = s
        gamma = sum(s*s)
        rhs_norm = sqrt(sum(f_data*f_data))

        if (rhs_norm == 0.0_dp) then
            coeff = 0.0_dp
            deallocate(r, q, s, p)
            return
        end if

        do k = 1, kmax
            call apply_A2D(spl, x_data, y_data, p, q)
            denom = sum(q*q)
            if (denom <= 0.0_dp) exit

            alpha = gamma/denom
            coeff = coeff + alpha*p
            r = r - alpha*q

            call apply_A2D_T(spl, x_data, y_data, r, s)
            gamma_new = sum(s*s)

            if (gamma_new <= (atol*rhs_norm)**2) exit

            beta = gamma_new/gamma
            gamma = gamma_new
            p = s + beta*p
        end do

        deallocate(r, q, s, p)
    end subroutine bspline_2d_lsq_cgls


    subroutine bspline_2d_lsq_cgls_batch(spl, x_data, y_data, f_data, coeff, &
        max_iter, tol)
        !! Batched matrix-free CGLS for 2D tensor-product B-splines.
        !! Solves independent LSQ problems for multiple right-hand sides:
        !!   min_{C(:,:,k)} sum_i (S_{C(:,:,k)}(x_i,y_i) - f_i(k))^2
        type(bspline_2d), intent(in) :: spl
        real(dp), intent(in) :: x_data(:)
        real(dp), intent(in) :: y_data(:)
        real(dp), intent(in) :: f_data(:,:)
        real(dp), intent(inout) :: coeff(:,:,:)
        integer, intent(in), optional :: max_iter
        real(dp), intent(in), optional :: tol

        integer :: n_data, n_rhs, nx, ny, k

        n_data = size(x_data)
        if (n_data /= size(y_data) .or. n_data /= size(f_data, 1)) then
            error stop "bspline_2d_lsq_cgls_batch: data size mismatch"
        end if

        n_rhs = size(f_data, 2)
        nx = spl%sx%n_ctrl
        ny = spl%sy%n_ctrl

        if (size(coeff, 1) /= nx .or. size(coeff, 2) /= ny .or. &
                size(coeff, 3) /= n_rhs) then
            error stop "bspline_2d_lsq_cgls_batch: coeff shape mismatch"
        end if

        do k = 1, n_rhs
            if (present(max_iter) .and. present(tol)) then
                call bspline_2d_lsq_cgls(spl, x_data, y_data, f_data(:, k), &
                    coeff(:, :, k), max_iter=max_iter, tol=tol)
            elseif (present(max_iter)) then
                call bspline_2d_lsq_cgls(spl, x_data, y_data, f_data(:, k), &
                    coeff(:, :, k), max_iter=max_iter)
            elseif (present(tol)) then
                call bspline_2d_lsq_cgls(spl, x_data, y_data, f_data(:, k), &
                    coeff(:, :, k), tol=tol)
            else
                call bspline_2d_lsq_cgls(spl, x_data, y_data, f_data(:, k), &
                    coeff(:, :, k))
            end if
        end do
    end subroutine bspline_2d_lsq_cgls_batch


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
        y = 0.0_dp

!$omp parallel default(shared) private(i, span, j, N) if (n_data > 100)
        allocate(N(0:p))
!$omp do
        do i = 1, n_data
            call find_span(spl, x_data(i), span)
            call basis_funs(spl, span, x_data(i), N)
            do j = 0, p
                y(i) = y(i) + N(j)*coeff(span - p + j)
            end do
        end do
!$omp end do
        deallocate(N)
!$omp end parallel
    end subroutine apply_A


    subroutine apply_AT(spl, x_data, r, g)
        !! g = A^T * r  (adjoint action, accumulating into control points)
        type(bspline_1d), intent(in) :: spl
        real(dp), intent(in) :: x_data(:)
        real(dp), intent(in) :: r(:)
        real(dp), intent(out) :: g(:)

        integer :: n_data, i, span, p, j
        real(dp), allocatable :: N(:)
        real(dp), allocatable :: g_local(:)

        n_data = size(x_data)
        if (size(r) /= n_data) then
            error stop "apply_AT: r size mismatch"
        end if
        if (size(g) /= spl%n_ctrl) then
            error stop "apply_AT: g size mismatch"
        end if

        p = spl%degree
        g = 0.0_dp

!$omp parallel default(shared) private(i, span, j, N, g_local) if (n_data > 100)
        allocate(N(0:p))
        allocate(g_local(size(g)))
        g_local = 0.0_dp
!$omp do
        do i = 1, n_data
            call find_span(spl, x_data(i), span)
            call basis_funs(spl, span, x_data(i), N)
            do j = 0, p
                g_local(span - p + j) = g_local(span - p + j) + N(j)*r(i)
            end do
        end do
!$omp end do
!$omp critical
        g = g + g_local
!$omp end critical
        deallocate(g_local)
        deallocate(N)
!$omp end parallel
    end subroutine apply_AT


    subroutine apply_A2D(spl, x_data, y_data, coeff, f)
        !! f = A * coeff  (2D spline evaluation at all data points)
        type(bspline_2d), intent(in) :: spl
        real(dp), intent(in) :: x_data(:)
        real(dp), intent(in) :: y_data(:)
        real(dp), intent(in) :: coeff(:,:)
        real(dp), intent(out) :: f(:)

        integer :: n_data, nx, ny
        integer :: i, spanx, spany, px, py, a, b, ix, iy
        real(dp), allocatable :: Nx_b(:), Ny_b(:)

        n_data = size(x_data)
        if (n_data /= size(y_data) .or. n_data /= size(f)) then
            error stop "apply_A2D: data size mismatch"
        end if

        nx = spl%sx%n_ctrl
        ny = spl%sy%n_ctrl
        if (size(coeff, 1) /= nx .or. size(coeff, 2) /= ny) then
            error stop "apply_A2D: coeff shape mismatch"
        end if

        px = spl%sx%degree
        py = spl%sy%degree

        f = 0.0_dp

!$omp parallel default(shared) private(i, spanx, spany, a, b, ix, iy, Nx_b, Ny_b) &
!$omp if (n_data > 100)
        allocate(Nx_b(0:px))
        allocate(Ny_b(0:py))
!$omp do
        do i = 1, n_data
            call find_span(spl%sx, x_data(i), spanx)
            call basis_funs(spl%sx, spanx, x_data(i), Nx_b)
            call find_span(spl%sy, y_data(i), spany)
            call basis_funs(spl%sy, spany, y_data(i), Ny_b)
            do a = 0, px
                ix = spanx - px + a
                do b = 0, py
                    iy = spany - py + b
                    f(i) = f(i) + Nx_b(a)*Ny_b(b)*coeff(ix, iy)
                end do
            end do
        end do
!$omp end do
        deallocate(Nx_b, Ny_b)
!$omp end parallel
    end subroutine apply_A2D


    subroutine apply_A2D_T(spl, x_data, y_data, r, g)
        !! g = A^T * r  (adjoint action in 2D coefficient space)
        type(bspline_2d), intent(in) :: spl
        real(dp), intent(in) :: x_data(:)
        real(dp), intent(in) :: y_data(:)
        real(dp), intent(in) :: r(:)
        real(dp), intent(out) :: g(:,:)

        integer :: n_data, nx, ny
        integer :: i, spanx, spany, px, py, a, b, ix, iy
        real(dp), allocatable :: Nx_b(:), Ny_b(:)
        real(dp), allocatable :: g_local(:,:)

        n_data = size(x_data)
        if (n_data /= size(y_data) .or. n_data /= size(r)) then
            error stop "apply_A2D_T: data size mismatch"
        end if

        nx = spl%sx%n_ctrl
        ny = spl%sy%n_ctrl
        if (size(g, 1) /= nx .or. size(g, 2) /= ny) then
            error stop "apply_A2D_T: g shape mismatch"
        end if

        px = spl%sx%degree
        py = spl%sy%degree

        g = 0.0_dp

!$omp parallel default(shared) private(i, spanx, spany, a, b, ix, iy) &
!$omp& private(Nx_b, Ny_b, g_local) if (n_data > 100)
        allocate(Nx_b(0:px))
        allocate(Ny_b(0:py))
        allocate(g_local(nx, ny))
        g_local = 0.0_dp
!$omp do
        do i = 1, n_data
            call find_span(spl%sx, x_data(i), spanx)
            call basis_funs(spl%sx, spanx, x_data(i), Nx_b)
            call find_span(spl%sy, y_data(i), spany)
            call basis_funs(spl%sy, spany, y_data(i), Ny_b)
            do a = 0, px
                ix = spanx - px + a
                do b = 0, py
                    iy = spany - py + b
                    g_local(ix, iy) = g_local(ix, iy) + Nx_b(a)*Ny_b(b)*r(i)
                end do
            end do
        end do
!$omp end do
!$omp critical
        g = g + g_local
!$omp end critical
        deallocate(g_local)
        deallocate(Nx_b, Ny_b)
!$omp end parallel
    end subroutine apply_A2D_T


end module neo_bspline
