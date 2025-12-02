module neo_bspline_3d
    !! 3D tensor-product B-splines and matrix-free LSQ via CGLS.
    use, intrinsic :: iso_fortran_env, only : dp => real64
    use neo_bspline, only : bspline_1d, bspline_1d_init_uniform, find_span, basis_funs
    implicit none
    private

    type :: bspline_3d
        type(bspline_1d) :: sx
        type(bspline_1d) :: sy
        type(bspline_1d) :: sz
    end type bspline_3d

    public :: bspline_3d
    public :: bspline_3d_init_uniform
    public :: bspline_3d_eval
    public :: bspline_3d_lsq_cgls
    public :: bspline_3d_lsq_cgls_batch
    public :: apply_A3D
    public :: apply_A3D_T

contains

    subroutine bspline_3d_init_uniform(spl, degree, n_ctrl, x_min, x_max)
        !! Initialise 3D tensor-product B-spline with open-uniform knots
        !! in each dimension on [x_min(j), x_max(j)].
        type(bspline_3d), intent(out) :: spl
        integer, intent(in) :: degree(3), n_ctrl(3)
        real(dp), intent(in) :: x_min(3), x_max(3)

        call bspline_1d_init_uniform(spl%sx, degree(1), n_ctrl(1), x_min(1), x_max(1))
        call bspline_1d_init_uniform(spl%sy, degree(2), n_ctrl(2), x_min(2), x_max(2))
        call bspline_1d_init_uniform(spl%sz, degree(3), n_ctrl(3), x_min(3), x_max(3))
    end subroutine bspline_3d_init_uniform


    subroutine bspline_3d_eval(spl, coeff, x, y)
        !! Evaluate 3D spline
        !!   S(x,y,z) = sum_{i,j,k} coeff(i,j,k) N_i(x) M_j(y) L_k(z).
        type(bspline_3d), intent(in) :: spl
        real(dp), intent(in) :: coeff(:,:,:)
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: y

        integer :: nx, ny, nz
        integer :: spanx, spany, spanz
        integer :: px, py, pz
        integer :: a, b, c
        integer :: ix, iy, iz
        real(dp), allocatable :: Nx_b(:), Ny_b(:), Nz_b(:)

        nx = spl%sx%n_ctrl
        ny = spl%sy%n_ctrl
        nz = spl%sz%n_ctrl
        if (size(coeff, 1) /= nx .or. size(coeff, 2) /= ny .or. &
                size(coeff, 3) /= nz) then
            error stop "bspline_3d_eval: coeff shape mismatch"
        end if

        px = spl%sx%degree
        py = spl%sy%degree
        pz = spl%sz%degree

        allocate(Nx_b(0:px))
        allocate(Ny_b(0:py))
        allocate(Nz_b(0:pz))

        call find_span(spl%sx, x(1), spanx)
        call basis_funs(spl%sx, spanx, x(1), Nx_b)
        call find_span(spl%sy, x(2), spany)
        call basis_funs(spl%sy, spany, x(2), Ny_b)
        call find_span(spl%sz, x(3), spanz)
        call basis_funs(spl%sz, spanz, x(3), Nz_b)

        y = 0.0_dp
        do a = 0, px
            ix = spanx - px + a
            do b = 0, py
                iy = spany - py + b
                do c = 0, pz
                    iz = spanz - pz + c
                    y = y + Nx_b(a)*Ny_b(b)*Nz_b(c)*coeff(ix, iy, iz)
                end do
            end do
        end do

        deallocate(Nx_b, Ny_b, Nz_b)
    end subroutine bspline_3d_eval


    subroutine bspline_3d_lsq_cgls(spl, x_data, y_data, z_data, f_data, coeff, &
        max_iter, tol)
        !! Matrix-free CGLS for 3D tensor-product B-spline LSQ:
        !!   min_C sum_i (S_C(x_i,y_i,z_i) - f_i)^2
        !!
        !! where C is coeff(nx,ny,nz) and S_C is the spline.
        type(bspline_3d), intent(in) :: spl
        real(dp), intent(in) :: x_data(:)
        real(dp), intent(in) :: y_data(:)
        real(dp), intent(in) :: z_data(:)
        real(dp), intent(in) :: f_data(:)
        real(dp), intent(inout) :: coeff(:,:,:)
        integer, intent(in), optional :: max_iter
        real(dp), intent(in), optional :: tol

        integer :: n_data, nx, ny, nz
        integer :: k, kmax
        real(dp) :: atol, gamma, gamma_new, alpha, beta, denom
        real(dp) :: rhs_norm
        real(dp), allocatable :: r(:), q(:)
        real(dp), allocatable :: s(:,:,:), p(:,:,:)

        n_data = size(x_data)
        if (n_data /= size(y_data) .or. n_data /= size(z_data) .or. &
                n_data /= size(f_data)) then
            error stop "bspline_3d_lsq_cgls: data size mismatch"
        end if

        nx = spl%sx%n_ctrl
        ny = spl%sy%n_ctrl
        nz = spl%sz%n_ctrl
        if (size(coeff, 1) /= nx .or. size(coeff, 2) /= ny .or. &
                size(coeff, 3) /= nz) then
            error stop "bspline_3d_lsq_cgls: coeff shape mismatch"
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
        allocate(s(nx, ny, nz))
        allocate(p(nx, ny, nz))

        coeff = 0.0_dp
        r = f_data

        call apply_A3D_T(spl, x_data, y_data, z_data, r, s)
        p = s
        gamma = sum(s*s)
        rhs_norm = sqrt(sum(f_data*f_data))

        if (rhs_norm == 0.0_dp) then
            coeff = 0.0_dp
            deallocate(r, q, s, p)
            return
        end if

        do k = 1, kmax
            call apply_A3D(spl, x_data, y_data, z_data, p, q)
            denom = sum(q*q)
            if (denom <= 0.0_dp) exit

            alpha = gamma/denom
            coeff = coeff + alpha*p
            r = r - alpha*q

            call apply_A3D_T(spl, x_data, y_data, z_data, r, s)
            gamma_new = sum(s*s)

            if (gamma_new <= (atol*rhs_norm)**2) exit

            beta = gamma_new/gamma
            gamma = gamma_new
            p = s + beta*p
        end do

        deallocate(r, q, s, p)
    end subroutine bspline_3d_lsq_cgls


    subroutine bspline_3d_lsq_cgls_batch(spl, x_data, y_data, z_data, f_data, &
        coeff, max_iter, tol)
        !! Batched matrix-free CGLS for 3D tensor-product B-splines.
        !! Solves independent LSQ problems for multiple right-hand sides.
        type(bspline_3d), intent(in) :: spl
        real(dp), intent(in) :: x_data(:)
        real(dp), intent(in) :: y_data(:)
        real(dp), intent(in) :: z_data(:)
        real(dp), intent(in) :: f_data(:,:)
        real(dp), intent(inout) :: coeff(:,:,:,:)
        integer, intent(in), optional :: max_iter
        real(dp), intent(in), optional :: tol

        integer :: n_data, n_rhs, nx, ny, nz, k

        n_data = size(x_data)
        if (n_data /= size(y_data) .or. n_data /= size(z_data) .or. &
                n_data /= size(f_data, 1)) then
            error stop "bspline_3d_lsq_cgls_batch: data size mismatch"
        end if

        n_rhs = size(f_data, 2)
        nx = spl%sx%n_ctrl
        ny = spl%sy%n_ctrl
        nz = spl%sz%n_ctrl

        if (size(coeff, 1) /= nx .or. size(coeff, 2) /= ny .or. &
                size(coeff, 3) /= nz .or. size(coeff, 4) /= n_rhs) then
            error stop "bspline_3d_lsq_cgls_batch: coeff shape mismatch"
        end if

        do k = 1, n_rhs
            if (present(max_iter) .and. present(tol)) then
                call bspline_3d_lsq_cgls(spl, x_data, y_data, z_data, f_data(:, k), &
                    coeff(:, :, :, k), max_iter=max_iter, tol=tol)
            elseif (present(max_iter)) then
                call bspline_3d_lsq_cgls(spl, x_data, y_data, z_data, f_data(:, k), &
                    coeff(:, :, :, k), max_iter=max_iter)
            elseif (present(tol)) then
                call bspline_3d_lsq_cgls(spl, x_data, y_data, z_data, f_data(:, k), &
                    coeff(:, :, :, k), tol=tol)
            else
                call bspline_3d_lsq_cgls(spl, x_data, y_data, z_data, f_data(:, k), &
                    coeff(:, :, :, k))
            end if
        end do
    end subroutine bspline_3d_lsq_cgls_batch


    subroutine apply_A3D(spl, x_data, y_data, z_data, coeff, f)
        !! f = A * coeff  (3D spline evaluation at all data points)
        type(bspline_3d), intent(in) :: spl
        real(dp), intent(in) :: x_data(:)
        real(dp), intent(in) :: y_data(:)
        real(dp), intent(in) :: z_data(:)
        real(dp), intent(in) :: coeff(:,:,:)
        real(dp), intent(out) :: f(:)

        integer :: n_data, nx, ny, nz
        integer :: i, spanx, spany, spanz
        integer :: px, py, pz
        integer :: a, b, c, ix, iy, iz
        real(dp), allocatable :: Nx_b(:), Ny_b(:), Nz_b(:)

        n_data = size(x_data)
        if (n_data /= size(y_data) .or. n_data /= size(z_data) .or. &
                n_data /= size(f)) then
            error stop "apply_A3D: data size mismatch"
        end if

        nx = spl%sx%n_ctrl
        ny = spl%sy%n_ctrl
        nz = spl%sz%n_ctrl
        if (size(coeff, 1) /= nx .or. size(coeff, 2) /= ny .or. &
                size(coeff, 3) /= nz) then
            error stop "apply_A3D: coeff shape mismatch"
        end if

        px = spl%sx%degree
        py = spl%sy%degree
        pz = spl%sz%degree

        f = 0.0_dp

!$omp parallel default(shared) private(i, spanx, spany, spanz, a, b, c, ix, iy, iz) &
!$omp& private(Nx_b, Ny_b, Nz_b) if (n_data > 100)
        allocate(Nx_b(0:px))
        allocate(Ny_b(0:py))
        allocate(Nz_b(0:pz))
!$omp do
        do i = 1, n_data
            call find_span(spl%sx, x_data(i), spanx)
            call basis_funs(spl%sx, spanx, x_data(i), Nx_b)
            call find_span(spl%sy, y_data(i), spany)
            call basis_funs(spl%sy, spany, y_data(i), Ny_b)
            call find_span(spl%sz, z_data(i), spanz)
            call basis_funs(spl%sz, spanz, z_data(i), Nz_b)

            do a = 0, px
                ix = spanx - px + a
                do b = 0, py
                    iy = spany - py + b
                    do c = 0, pz
                        iz = spanz - pz + c
                        f(i) = f(i) + Nx_b(a)*Ny_b(b)*Nz_b(c)*coeff(ix, iy, iz)
                    end do
                end do
            end do
        end do
!$omp end do
        deallocate(Nx_b, Ny_b, Nz_b)
!$omp end parallel
    end subroutine apply_A3D


    subroutine apply_A3D_T(spl, x_data, y_data, z_data, r, g)
        !! g = A^T * r  (adjoint action in 3D coefficient space)
        type(bspline_3d), intent(in) :: spl
        real(dp), intent(in) :: x_data(:)
        real(dp), intent(in) :: y_data(:)
        real(dp), intent(in) :: z_data(:)
        real(dp), intent(in) :: r(:)
        real(dp), intent(out) :: g(:,:,:)

        integer :: n_data, nx, ny, nz
        integer :: i, spanx, spany, spanz
        integer :: px, py, pz
        integer :: a, b, c, ix, iy, iz
        real(dp), allocatable :: Nx_b(:), Ny_b(:), Nz_b(:)
        real(dp), allocatable :: g_local(:,:,:)

        n_data = size(x_data)
        if (n_data /= size(y_data) .or. n_data /= size(z_data) .or. &
                n_data /= size(r)) then
            error stop "apply_A3D_T: data size mismatch"
        end if

        nx = spl%sx%n_ctrl
        ny = spl%sy%n_ctrl
        nz = spl%sz%n_ctrl
        if (size(g, 1) /= nx .or. size(g, 2) /= ny .or. size(g, 3) /= nz) then
            error stop "apply_A3D_T: g shape mismatch"
        end if

        px = spl%sx%degree
        py = spl%sy%degree
        pz = spl%sz%degree

        g = 0.0_dp

!$omp parallel default(shared) private(i, spanx, spany, spanz, a, b, c, ix, iy, iz) &
!$omp& private(Nx_b, Ny_b, Nz_b, g_local) if (n_data > 100)
        allocate(Nx_b(0:px))
        allocate(Ny_b(0:py))
        allocate(Nz_b(0:pz))
        allocate(g_local(nx, ny, nz))
        g_local = 0.0_dp
!$omp do
        do i = 1, n_data
            call find_span(spl%sx, x_data(i), spanx)
            call basis_funs(spl%sx, spanx, x_data(i), Nx_b)
            call find_span(spl%sy, y_data(i), spany)
            call basis_funs(spl%sy, spany, y_data(i), Ny_b)
            call find_span(spl%sz, z_data(i), spanz)
            call basis_funs(spl%sz, spanz, z_data(i), Nz_b)

            do a = 0, px
                ix = spanx - px + a
                do b = 0, py
                    iy = spany - py + b
                    do c = 0, pz
                        iz = spanz - pz + c
                        g_local(ix, iy, iz) = g_local(ix, iy, iz) + &
                            Nx_b(a)*Ny_b(b)*Nz_b(c)*r(i)
                    end do
                end do
            end do
        end do
!$omp end do
!$omp critical
        g = g + g_local
!$omp end critical
        deallocate(g_local)
        deallocate(Nx_b, Ny_b, Nz_b)
!$omp end parallel
    end subroutine apply_A3D_T

end module neo_bspline_3d
