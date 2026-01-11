module cgls_dense
    !! Conjugate gradient for normal equations using an explicit
    !! design matrix stored in column-major form.  OpenMP and SIMD
    !! are used inside the matrix-vector products.
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none
    private

    public :: cgls_dense_solve

contains

    subroutine matvec_phi_core(phi, x, y)
        !! y = A*x  where A == phi (n_rows x n_cols)
        real(dp), intent(in) :: phi(:,:)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)

        integer :: i, j, n_rows, n_cols
        real(dp) :: acc

        n_rows = size(phi, 1)
        n_cols = size(phi, 2)

        if (size(x) /= n_cols) then
            error stop "matvec_phi_core: x size mismatch"
        end if
        if (size(y) /= n_rows) then
            error stop "matvec_phi_core: y size mismatch"
        end if

        y = 0.0_dp

        !$omp parallel do default(shared) private(i, j, acc)
        do i = 1, n_rows
            acc = 0.0_dp
            !$omp simd reduction(+:acc)
            do j = 1, n_cols
                acc = acc + phi(i, j)*x(j)
            end do
            y(i) = acc
        end do
        !$omp end parallel do
    end subroutine matvec_phi_core


    subroutine matvec_phi_t_core(phi, r, y)
        !! y = A^T*r  where A == phi (n_rows x n_cols)
        real(dp), intent(in) :: phi(:,:)
        real(dp), intent(in) :: r(:)
        real(dp), intent(out) :: y(:)

        integer :: i, j, n_rows, n_cols
        real(dp) :: rloc

        n_rows = size(phi, 1)
        n_cols = size(phi, 2)

        if (size(r) /= n_rows) then
            error stop "matvec_phi_t_core: r size mismatch"
        end if
        if (size(y) /= n_cols) then
            error stop "matvec_phi_t_core: y size mismatch"
        end if

        y = 0.0_dp

        !$omp parallel do default(shared) private(j, i, rloc) schedule(static)
        do j = 1, n_cols
            rloc = 0.0_dp
            !$omp simd reduction(+:rloc)
            do i = 1, n_rows
                rloc = rloc + phi(i, j)*r(i)
            end do
            y(j) = rloc
        end do
        !$omp end parallel do
    end subroutine matvec_phi_t_core


    subroutine cgls_dense_core(phi, rhs, x, max_iter, tol)
        !! Conjugate Gradient for normal equations using
        !! explicit matrix-vector products (unweighted case).
        real(dp), intent(in) :: phi(:,:)
        real(dp), intent(in) :: rhs(:)
        real(dp), intent(inout) :: x(:)
        integer, intent(in) :: max_iter
        real(dp), intent(in) :: tol

        integer :: iter, n_rows, n_cols
        real(dp) :: gamma, gamma_new, alpha, beta, denom
        real(dp) :: rhs_norm

        real(dp), allocatable :: r(:), s(:), p(:), q(:)

        n_rows = size(phi, 1)
        n_cols = size(phi, 2)

        if (size(rhs) /= n_rows) then
            error stop "cgls_dense_core: rhs size mismatch"
        end if
        if (size(x) /= n_cols) then
            error stop "cgls_dense_core: x size mismatch"
        end if

        allocate(r(n_rows))
        allocate(s(n_cols))
        allocate(p(n_cols))
        allocate(q(n_rows))

        r = rhs

        call matvec_phi_t_core(phi, r, s)
        p = s
        gamma = sum(s*s)
        rhs_norm = sqrt(sum(r*r))

        if (rhs_norm == 0.0_dp) then
            x = 0.0_dp
            deallocate(r, s, p, q)
            return
        end if

        do iter = 1, max_iter
            call matvec_phi_core(phi, p, q)
            denom = sum(q*q)
            if (denom <= 0.0_dp) exit

            alpha = gamma/denom
            x = x + alpha*p
            r = r - alpha*q

            call matvec_phi_t_core(phi, r, s)
            gamma_new = sum(s*s)

            if (gamma_new <= (tol*rhs_norm)**2) exit

            beta = gamma_new/gamma
            gamma = gamma_new
            p = s + beta*p
        end do

        deallocate(r, s, p, q)
    end subroutine cgls_dense_core


    subroutine cgls_dense_solve(phi, rhs, x, w, max_iter, tol)
        !! Conjugate Gradient for normal equations using
        !! explicit matrix-vector products.
        !! If w is present, each row i is scaled by w(i).
        real(dp), intent(in) :: phi(:,:)
        real(dp), intent(in) :: rhs(:)
        real(dp), intent(out) :: x(:)
        real(dp), intent(in), optional :: w(:)
        integer, intent(in), optional :: max_iter
        real(dp), intent(in), optional :: tol

        integer :: n_rows, n_cols, kmax, i, j
        real(dp) :: atol
        real(dp), allocatable :: phi_eff(:,:), rhs_eff(:)

        n_rows = size(phi, 1)
        n_cols = size(phi, 2)

        if (size(rhs) /= n_rows) then
            error stop "cgls_dense_solve: rhs size mismatch"
        end if
        if (size(x) /= n_cols) then
            error stop "cgls_dense_solve: x size mismatch"
        end if
        if (present(w)) then
            if (size(w) /= n_rows) then
                error stop "cgls_dense_solve: w size mismatch"
            end if
        end if

        if (present(max_iter)) then
            kmax = max_iter
        else
            kmax = 200
        end if
        if (present(tol)) then
            atol = tol
        else
            atol = 1.0d-12
        end if

        x = 0.0_dp

        if (present(w)) then
            allocate(phi_eff(n_rows, n_cols))
            allocate(rhs_eff(n_rows))

            !$omp parallel do default(shared) private(i, j)
            do i = 1, n_rows
                rhs_eff(i) = w(i)*rhs(i)
                !$omp simd
                do j = 1, n_cols
                    phi_eff(i, j) = w(i)*phi(i, j)
                end do
            end do
            !$omp end parallel do

            call cgls_dense_core(phi_eff, rhs_eff, x, kmax, atol)

            deallocate(phi_eff, rhs_eff)
        else
            call cgls_dense_core(phi, rhs, x, kmax, atol)
        end if
    end subroutine cgls_dense_solve

end module cgls_dense
