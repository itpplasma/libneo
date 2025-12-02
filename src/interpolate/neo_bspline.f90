module neo_bspline
    !! High-performance B-spline module for 1D/2D/3D tensor-product splines.
    !! Optimized for regular grids with OpenMP SIMD and parallel construction.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    type :: bspline_1d
        integer :: degree
        integer :: n_ctrl
        real(dp), allocatable :: knots(:)
        real(dp) :: x_min, x_max
    end type bspline_1d

    type :: bspline_2d
        type(bspline_1d) :: sx, sy
    end type bspline_2d

    type :: bspline_3d
        type(bspline_1d) :: sx, sy, sz
    end type bspline_3d

    type :: grid_cache_1d
        integer :: ng1, px
        integer, allocatable :: span1(:)
        real(dp), allocatable :: B1(:,:)
    end type grid_cache_1d

    type :: grid_cache_2d
        integer :: ng1, ng2, px, py
        integer, allocatable :: span1(:), span2(:)
        real(dp), allocatable :: B1(:,:), B2(:,:)
    end type grid_cache_2d

    type :: grid_cache_3d
        integer :: ng1, ng2, ng3, px, py, pz
        integer, allocatable :: span1(:), span2(:), span3(:)
        real(dp), allocatable :: B1(:,:), B2(:,:), B3(:,:)
    end type grid_cache_3d

    public :: bspline_1d, bspline_2d, bspline_3d
    public :: bspline_1d_init_uniform, bspline_2d_init_uniform, bspline_3d_init_uniform
    public :: bspline_1d_eval, bspline_2d_eval, bspline_3d_eval
    public :: bspline_1d_lsq_cgls, bspline_2d_lsq_cgls, bspline_3d_lsq_cgls
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

        call bspline_1d_init_uniform(spl%sx, degree(1), n_ctrl(1), x_min(1), x_max(1))
        call bspline_1d_init_uniform(spl%sy, degree(2), n_ctrl(2), x_min(2), x_max(2))
    end subroutine bspline_2d_init_uniform


    subroutine bspline_3d_init_uniform(spl, degree, n_ctrl, x_min, x_max)
        type(bspline_3d), intent(out) :: spl
        integer, intent(in) :: degree(3), n_ctrl(3)
        real(dp), intent(in) :: x_min(3), x_max(3)

        call bspline_1d_init_uniform(spl%sx, degree(1), n_ctrl(1), x_min(1), x_max(1))
        call bspline_1d_init_uniform(spl%sy, degree(2), n_ctrl(2), x_min(2), x_max(2))
        call bspline_1d_init_uniform(spl%sz, degree(3), n_ctrl(3), x_min(3), x_max(3))
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
        real(dp) :: Nx(0:spl%sx%degree), Ny(0:spl%sy%degree), Nz(0:spl%sz%degree)

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


    subroutine precompute_grid_1d(spl, x1, cache)
        type(bspline_1d), intent(in) :: spl
        real(dp), intent(in) :: x1(:)
        type(grid_cache_1d), intent(out) :: cache
        integer :: i, n1, px

        n1 = size(x1)
        px = spl%degree
        cache%ng1 = n1
        cache%px = px
        allocate(cache%span1(n1), cache%B1(0:px, n1))

        !$omp parallel do simd private(i) schedule(simd:static)
        do i = 1, n1
            call find_span(spl, x1(i), cache%span1(i))
        end do

        !$omp parallel do private(i) schedule(static)
        do i = 1, n1
            call basis_funs(spl, cache%span1(i), x1(i), cache%B1(:, i))
        end do
    end subroutine precompute_grid_1d


    subroutine precompute_grid_2d(spl, x1, x2, cache)
        type(bspline_2d), intent(in) :: spl
        real(dp), intent(in) :: x1(:), x2(:)
        type(grid_cache_2d), intent(out) :: cache
        integer :: i, n1, n2, px, py

        n1 = size(x1)
        n2 = size(x2)
        px = spl%sx%degree
        py = spl%sy%degree
        cache%ng1 = n1
        cache%ng2 = n2
        cache%px = px
        cache%py = py
        allocate(cache%span1(n1), cache%span2(n2))
        allocate(cache%B1(0:px, n1), cache%B2(0:py, n2))

        !$omp parallel do simd private(i) schedule(simd:static)
        do i = 1, n1
            call find_span(spl%sx, x1(i), cache%span1(i))
        end do
        !$omp parallel do private(i) schedule(static)
        do i = 1, n1
            call basis_funs(spl%sx, cache%span1(i), x1(i), cache%B1(:, i))
        end do

        !$omp parallel do simd private(i) schedule(simd:static)
        do i = 1, n2
            call find_span(spl%sy, x2(i), cache%span2(i))
        end do
        !$omp parallel do private(i) schedule(static)
        do i = 1, n2
            call basis_funs(spl%sy, cache%span2(i), x2(i), cache%B2(:, i))
        end do
    end subroutine precompute_grid_2d


    subroutine precompute_grid_3d(spl, x1, x2, x3, cache)
        type(bspline_3d), intent(in) :: spl
        real(dp), intent(in) :: x1(:), x2(:), x3(:)
        type(grid_cache_3d), intent(out) :: cache
        integer :: i, n1, n2, n3, px, py, pz

        n1 = size(x1)
        n2 = size(x2)
        n3 = size(x3)
        px = spl%sx%degree
        py = spl%sy%degree
        pz = spl%sz%degree
        cache%ng1 = n1
        cache%ng2 = n2
        cache%ng3 = n3
        cache%px = px
        cache%py = py
        cache%pz = pz
        allocate(cache%span1(n1), cache%span2(n2), cache%span3(n3))
        allocate(cache%B1(0:px, n1), cache%B2(0:py, n2), cache%B3(0:pz, n3))

        !$omp parallel do simd private(i) schedule(simd:static)
        do i = 1, n1
            call find_span(spl%sx, x1(i), cache%span1(i))
        end do
        !$omp parallel do private(i) schedule(static)
        do i = 1, n1
            call basis_funs(spl%sx, cache%span1(i), x1(i), cache%B1(:, i))
        end do

        !$omp parallel do simd private(i) schedule(simd:static)
        do i = 1, n2
            call find_span(spl%sy, x2(i), cache%span2(i))
        end do
        !$omp parallel do private(i) schedule(static)
        do i = 1, n2
            call basis_funs(spl%sy, cache%span2(i), x2(i), cache%B2(:, i))
        end do

        !$omp parallel do simd private(i) schedule(simd:static)
        do i = 1, n3
            call find_span(spl%sz, x3(i), cache%span3(i))
        end do
        !$omp parallel do private(i) schedule(static)
        do i = 1, n3
            call basis_funs(spl%sz, cache%span3(i), x3(i), cache%B3(:, i))
        end do
    end subroutine precompute_grid_3d


    subroutine apply_A_1d_grid(cache, coeff, f)
        type(grid_cache_1d), intent(in) :: cache
        real(dp), intent(in) :: coeff(:)
        real(dp), intent(out) :: f(:)
        integer :: i, a, ix, n1, px, s1

        n1 = cache%ng1
        px = cache%px
        !$omp parallel do simd private(i, a, ix, s1) schedule(simd:static)
        do i = 1, n1
            s1 = cache%span1(i)
            f(i) = 0.0_dp
            do a = 0, px
                ix = s1 - px + a
                f(i) = f(i) + cache%B1(a, i)*coeff(ix)
            end do
        end do
    end subroutine apply_A_1d_grid


    subroutine apply_AT_1d_grid(cache, r, g, n_ctrl)
        type(grid_cache_1d), intent(in) :: cache
        real(dp), intent(in) :: r(:)
        real(dp), intent(out) :: g(:)
        integer, intent(in) :: n_ctrl
        integer :: i, a, ix, n1, px, s1
        real(dp), allocatable :: g_local(:)

        n1 = cache%ng1
        px = cache%px
        g = 0.0_dp

        !$omp parallel private(g_local, i, a, ix, s1) default(shared)
        allocate(g_local(n_ctrl))
        g_local = 0.0_dp
        !$omp do schedule(static)
        do i = 1, n1
            s1 = cache%span1(i)
            do a = 0, px
                ix = s1 - px + a
                g_local(ix) = g_local(ix) + cache%B1(a, i)*r(i)
            end do
        end do
        !$omp end do
        !$omp critical
        g = g + g_local
        !$omp end critical
        deallocate(g_local)
        !$omp end parallel
    end subroutine apply_AT_1d_grid


    subroutine apply_A_2d_grid(cache, coeff, f)
        type(grid_cache_2d), intent(in) :: cache
        real(dp), intent(in) :: coeff(:,:)
        real(dp), intent(out) :: f(:,:)
        integer :: i1, i2, a, b, ix, iy, n1, n2, px, py, s1, s2
        real(dp) :: val

        n1 = cache%ng1
        n2 = cache%ng2
        px = cache%px
        py = cache%py

        !$omp parallel do collapse(2) private(i1, i2, a, b, ix, iy, s1, s2, val) &
        !$omp& schedule(static)
        do i2 = 1, n2
            do i1 = 1, n1
                s1 = cache%span1(i1)
                s2 = cache%span2(i2)
                val = 0.0_dp
                do a = 0, px
                    ix = s1 - px + a
                    do b = 0, py
                        iy = s2 - py + b
                        val = val + cache%B1(a, i1)*cache%B2(b, i2)*coeff(ix, iy)
                    end do
                end do
                f(i1, i2) = val
            end do
        end do
    end subroutine apply_A_2d_grid


    subroutine apply_AT_2d_grid(cache, r, g, nx, ny)
        type(grid_cache_2d), intent(in) :: cache
        real(dp), intent(in) :: r(:,:)
        real(dp), intent(out) :: g(:,:)
        integer, intent(in) :: nx, ny
        integer :: i1, i2, a, b, ix, iy, n1, n2, px, py, s1, s2
        real(dp), allocatable :: g_local(:,:)

        n1 = cache%ng1
        n2 = cache%ng2
        px = cache%px
        py = cache%py
        g = 0.0_dp

        !$omp parallel private(g_local, i1, i2, a, b, ix, iy, s1, s2) default(shared)
        allocate(g_local(nx, ny))
        g_local = 0.0_dp
        !$omp do collapse(2) schedule(static)
        do i2 = 1, n2
            do i1 = 1, n1
                s1 = cache%span1(i1)
                s2 = cache%span2(i2)
                do a = 0, px
                    ix = s1 - px + a
                    do b = 0, py
                        iy = s2 - py + b
                        g_local(ix, iy) = g_local(ix, iy) + &
                            cache%B1(a, i1)*cache%B2(b, i2)*r(i1, i2)
                    end do
                end do
            end do
        end do
        !$omp end do
        !$omp critical
        g = g + g_local
        !$omp end critical
        deallocate(g_local)
        !$omp end parallel
    end subroutine apply_AT_2d_grid


    subroutine apply_A_3d_grid(cache, coeff, f)
        type(grid_cache_3d), intent(in) :: cache
        real(dp), intent(in) :: coeff(:,:,:)
        real(dp), intent(out) :: f(:,:,:)
        integer :: i1, i2, i3, a, b, c, ix, iy, iz
        integer :: n1, n2, n3, px, py, pz, s1, s2, s3
        real(dp) :: val, Nab

        n1 = cache%ng1
        n2 = cache%ng2
        n3 = cache%ng3
        px = cache%px
        py = cache%py
        pz = cache%pz

        !$omp parallel do collapse(3) &
        !$omp& private(i1, i2, i3, a, b, c, ix, iy, iz, s1, s2, s3, val, Nab) &
        !$omp& schedule(static)
        do i3 = 1, n3
            do i2 = 1, n2
                do i1 = 1, n1
                    s1 = cache%span1(i1)
                    s2 = cache%span2(i2)
                    s3 = cache%span3(i3)
                    val = 0.0_dp
                    do a = 0, px
                        ix = s1 - px + a
                        do b = 0, py
                            iy = s2 - py + b
                            Nab = cache%B1(a, i1)*cache%B2(b, i2)
                            !$omp simd reduction(+:val)
                            do c = 0, pz
                                iz = s3 - pz + c
                                val = val + Nab*cache%B3(c, i3)*coeff(ix, iy, iz)
                            end do
                        end do
                    end do
                    f(i1, i2, i3) = val
                end do
            end do
        end do
    end subroutine apply_A_3d_grid


    subroutine apply_AT_3d_grid(cache, r, g, nx, ny, nz)
        type(grid_cache_3d), intent(in) :: cache
        real(dp), intent(in) :: r(:,:,:)
        real(dp), intent(out) :: g(:,:,:)
        integer, intent(in) :: nx, ny, nz
        integer :: i1, i2, i3, a, b, c, ix, iy, iz
        integer :: n1, n2, n3, px, py, pz, s1, s2, s3
        real(dp) :: ri, Nab
        real(dp), allocatable :: g_local(:,:,:)

        n1 = cache%ng1
        n2 = cache%ng2
        n3 = cache%ng3
        px = cache%px
        py = cache%py
        pz = cache%pz
        g = 0.0_dp

        !$omp parallel private(g_local, i1, i2, i3, a, b, c, ix, iy, iz, s1, s2, s3, ri, Nab) &
        !$omp& default(shared)
        allocate(g_local(nx, ny, nz))
        g_local = 0.0_dp
        !$omp do collapse(3) schedule(static)
        do i3 = 1, n3
            do i2 = 1, n2
                do i1 = 1, n1
                    s1 = cache%span1(i1)
                    s2 = cache%span2(i2)
                    s3 = cache%span3(i3)
                    ri = r(i1, i2, i3)
                    do a = 0, px
                        ix = s1 - px + a
                        do b = 0, py
                            iy = s2 - py + b
                            Nab = cache%B1(a, i1)*cache%B2(b, i2)*ri
                            do c = 0, pz
                                iz = s3 - pz + c
                                g_local(ix, iy, iz) = g_local(ix, iy, iz) + &
                                    Nab*cache%B3(c, i3)
                            end do
                        end do
                    end do
                end do
            end do
        end do
        !$omp end do
        !$omp critical
        g = g + g_local
        !$omp end critical
        deallocate(g_local)
        !$omp end parallel
    end subroutine apply_AT_3d_grid


    subroutine bspline_1d_lsq_cgls(spl, x_data, f_data, coeff, max_iter, tol)
        type(bspline_1d), intent(in) :: spl
        real(dp), intent(in) :: x_data(:)
        real(dp), intent(in) :: f_data(:)
        real(dp), intent(inout) :: coeff(:)
        integer, intent(in), optional :: max_iter
        real(dp), intent(in), optional :: tol

        type(grid_cache_1d) :: cache
        integer :: n_data, n_ctrl, k, kmax
        real(dp) :: atol, gamma, gamma_new, alpha, beta, denom, rhs_norm
        real(dp), allocatable :: r(:), q(:), s(:), p(:)

        n_data = size(x_data)
        n_ctrl = spl%n_ctrl
        kmax = 200
        atol = 1.0d-10
        if (present(max_iter)) kmax = max_iter
        if (present(tol)) atol = tol

        call precompute_grid_1d(spl, x_data, cache)

        allocate(r(n_data), q(n_data), s(n_ctrl), p(n_ctrl))
        coeff = 0.0_dp
        r = f_data
        call apply_AT_1d_grid(cache, r, s, n_ctrl)
        p = s
        gamma = dot_product(s, s)
        rhs_norm = sqrt(dot_product(f_data, f_data))
        if (rhs_norm == 0.0_dp) return

        do k = 1, kmax
            call apply_A_1d_grid(cache, p, q)
            denom = dot_product(q, q)
            if (denom <= 0.0_dp) exit
            alpha = gamma/denom
            coeff = coeff + alpha*p
            r = r - alpha*q
            call apply_AT_1d_grid(cache, r, s, n_ctrl)
            gamma_new = dot_product(s, s)
            if (gamma_new <= (atol*rhs_norm)**2) exit
            beta = gamma_new/gamma
            gamma = gamma_new
            p = s + beta*p
        end do
    end subroutine bspline_1d_lsq_cgls


    subroutine bspline_2d_lsq_cgls(spl, x1, x2, f_grid, coeff, max_iter, tol)
        type(bspline_2d), intent(in) :: spl
        real(dp), intent(in) :: x1(:), x2(:)
        real(dp), intent(in) :: f_grid(:,:)
        real(dp), intent(inout) :: coeff(:,:)
        integer, intent(in), optional :: max_iter
        real(dp), intent(in), optional :: tol

        type(grid_cache_2d) :: cache
        integer :: n1, n2, nx, ny, k, kmax
        real(dp) :: atol, gamma, gamma_new, alpha, beta, denom, rhs_norm
        real(dp), allocatable :: r(:,:), q(:,:), s(:,:), p(:,:)

        n1 = size(x1)
        n2 = size(x2)
        nx = spl%sx%n_ctrl
        ny = spl%sy%n_ctrl
        kmax = 400
        atol = 1.0d-10
        if (present(max_iter)) kmax = max_iter
        if (present(tol)) atol = tol

        call precompute_grid_2d(spl, x1, x2, cache)

        allocate(r(n1, n2), q(n1, n2), s(nx, ny), p(nx, ny))
        coeff = 0.0_dp
        r = f_grid
        call apply_AT_2d_grid(cache, r, s, nx, ny)
        p = s
        gamma = sum(s*s)
        rhs_norm = sqrt(sum(f_grid*f_grid))
        if (rhs_norm == 0.0_dp) return

        do k = 1, kmax
            call apply_A_2d_grid(cache, p, q)
            denom = sum(q*q)
            if (denom <= 0.0_dp) exit
            alpha = gamma/denom
            coeff = coeff + alpha*p
            r = r - alpha*q
            call apply_AT_2d_grid(cache, r, s, nx, ny)
            gamma_new = sum(s*s)
            if (gamma_new <= (atol*rhs_norm)**2) exit
            beta = gamma_new/gamma
            gamma = gamma_new
            p = s + beta*p
        end do
    end subroutine bspline_2d_lsq_cgls


    subroutine bspline_3d_lsq_cgls(spl, x1, x2, x3, f_grid, coeff, max_iter, tol)
        type(bspline_3d), intent(in) :: spl
        real(dp), intent(in) :: x1(:), x2(:), x3(:)
        real(dp), intent(in) :: f_grid(:,:,:)
        real(dp), intent(inout) :: coeff(:,:,:)
        integer, intent(in), optional :: max_iter
        real(dp), intent(in), optional :: tol

        type(grid_cache_3d) :: cache
        integer :: n1, n2, n3, nx, ny, nz, k, kmax
        real(dp) :: atol, gamma, gamma_new, alpha, beta, denom, rhs_norm
        real(dp), allocatable :: r(:,:,:), q(:,:,:), s(:,:,:), p(:,:,:)

        n1 = size(x1)
        n2 = size(x2)
        n3 = size(x3)
        nx = spl%sx%n_ctrl
        ny = spl%sy%n_ctrl
        nz = spl%sz%n_ctrl
        kmax = 400
        atol = 1.0d-10
        if (present(max_iter)) kmax = max_iter
        if (present(tol)) atol = tol

        call precompute_grid_3d(spl, x1, x2, x3, cache)

        allocate(r(n1, n2, n3), q(n1, n2, n3), s(nx, ny, nz), p(nx, ny, nz))
        coeff = 0.0_dp
        r = f_grid
        call apply_AT_3d_grid(cache, r, s, nx, ny, nz)
        p = s
        gamma = sum(s*s)
        rhs_norm = sqrt(sum(f_grid*f_grid))
        if (rhs_norm == 0.0_dp) return

        do k = 1, kmax
            call apply_A_3d_grid(cache, p, q)
            denom = sum(q*q)
            if (denom <= 0.0_dp) exit
            alpha = gamma/denom
            coeff = coeff + alpha*p
            r = r - alpha*q
            call apply_AT_3d_grid(cache, r, s, nx, ny, nz)
            gamma_new = sum(s*s)
            if (gamma_new <= (atol*rhs_norm)**2) exit
            beta = gamma_new/gamma
            gamma = gamma_new
            p = s + beta*p
        end do
    end subroutine bspline_3d_lsq_cgls


end module neo_bspline
