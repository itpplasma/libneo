module neo_bspline_interp
    !! Direct B-spline interpolation via collocation.
    !! Solves the linear system A*c = f where A_ij = B_j(x_i).
    !! Uses separable tensor-product structure for efficiency in 2D/3D.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neo_bspline_base, only: bspline_1d, bspline_2d, bspline_3d, &
        find_span, basis_funs
    implicit none
    private

    public :: bspline_1d_interp, bspline_2d_interp, bspline_3d_interp

contains

    subroutine bspline_1d_interp(spl, x_data, f_data, coeff)
        !! Direct 1D B-spline interpolation: solve collocation system A*c = f.
        type(bspline_1d), intent(in) :: spl
        real(dp), intent(in) :: x_data(:)
        real(dp), intent(in) :: f_data(:)
        real(dp), intent(out) :: coeff(:)

        integer :: n, p, i, j, span, info
        real(dp), allocatable :: A(:,:), work_f(:), Nvals(:)
        integer, allocatable :: ipiv(:)

        n = size(x_data)
        p = spl%degree
        allocate(A(n, n), work_f(n), ipiv(n), Nvals(0:p))
        A = 0.0_dp

        do i = 1, n
            call find_span(spl, x_data(i), span)
            call basis_funs(spl, span, x_data(i), Nvals)
            do j = 0, p
                A(i, span - p + j) = Nvals(j)
            end do
        end do

        work_f = f_data
        call dgesv(n, 1, A, n, ipiv, work_f, n, info)
        if (info /= 0) then
            coeff = 0.0_dp
            return
        end if
        coeff = work_f
    end subroutine bspline_1d_interp


    subroutine bspline_2d_interp(spl, x1, x2, f_grid, coeff)
        !! Direct 2D B-spline interpolation using separable tensor product.
        !! Solves dimension-by-dimension: first along x1, then along x2.
        type(bspline_2d), intent(in) :: spl
        real(dp), intent(in) :: x1(:), x2(:)
        real(dp), intent(in) :: f_grid(:,:)
        real(dp), intent(out) :: coeff(:,:)

        integer :: n1, n2, i, j, span, p1, p2, info, pmax
        real(dp), allocatable :: A1(:,:), A2(:,:), temp(:,:), work(:), Nvals(:)
        integer, allocatable :: ipiv(:)

        n1 = size(x1)
        n2 = size(x2)
        p1 = spl%sx%degree
        p2 = spl%sy%degree
        pmax = max(p1, p2)
        allocate(A1(n1, n1), A2(n2, n2), temp(n1, n2), ipiv(max(n1, n2)))
        allocate(work(max(n1, n2)), Nvals(0:pmax))

        A1 = 0.0_dp
        do i = 1, n1
            call find_span(spl%sx, x1(i), span)
            call basis_funs(spl%sx, span, x1(i), Nvals(0:p1))
            do j = 0, p1
                A1(i, span - p1 + j) = Nvals(j)
            end do
        end do

        A2 = 0.0_dp
        do i = 1, n2
            call find_span(spl%sy, x2(i), span)
            call basis_funs(spl%sy, span, x2(i), Nvals(0:p2))
            do j = 0, p2
                A2(i, span - p2 + j) = Nvals(j)
            end do
        end do

        temp = f_grid
        call dgetrf(n1, n1, A1, n1, ipiv, info)
        if (info /= 0) then
            coeff = 0.0_dp
            return
        end if
        do j = 1, n2
            work(1:n1) = temp(:, j)
            call dgetrs('N', n1, 1, A1, n1, ipiv, work, n1, info)
            temp(:, j) = work(1:n1)
        end do

        call dgetrf(n2, n2, A2, n2, ipiv, info)
        if (info /= 0) then
            coeff = 0.0_dp
            return
        end if
        do i = 1, n1
            work(1:n2) = temp(i, :)
            call dgetrs('N', n2, 1, A2, n2, ipiv, work, n2, info)
            coeff(i, :) = work(1:n2)
        end do
    end subroutine bspline_2d_interp


    subroutine bspline_3d_interp(spl, x1, x2, x3, f_grid, coeff)
        !! Direct 3D B-spline interpolation using separable tensor product.
        !! Solves dimension-by-dimension: x1 -> x2 -> x3.
        type(bspline_3d), intent(in) :: spl
        real(dp), intent(in) :: x1(:), x2(:), x3(:)
        real(dp), intent(in) :: f_grid(:,:,:)
        real(dp), intent(out) :: coeff(:,:,:)

        integer :: n1, n2, n3, i, j, k, span, p1, p2, p3, info, pmax, nmax
        real(dp), allocatable :: A1(:,:), A2(:,:), A3(:,:), temp(:,:,:)
        real(dp), allocatable :: Nvals(:)
        integer, allocatable :: ipiv1(:), ipiv2(:), ipiv3(:)

        n1 = size(x1)
        n2 = size(x2)
        n3 = size(x3)
        p1 = spl%sx%degree
        p2 = spl%sy%degree
        p3 = spl%sz%degree
        pmax = max(p1, p2, p3)
        nmax = max(n1, n2, n3)
        allocate(A1(n1, n1), A2(n2, n2), A3(n3, n3))
        allocate(temp(n1, n2, n3), ipiv1(n1), ipiv2(n2), ipiv3(n3))
        allocate(Nvals(0:pmax))

        A1 = 0.0_dp
        do i = 1, n1
            call find_span(spl%sx, x1(i), span)
            call basis_funs(spl%sx, span, x1(i), Nvals(0:p1))
            do j = 0, p1
                A1(i, span - p1 + j) = Nvals(j)
            end do
        end do

        A2 = 0.0_dp
        do i = 1, n2
            call find_span(spl%sy, x2(i), span)
            call basis_funs(spl%sy, span, x2(i), Nvals(0:p2))
            do j = 0, p2
                A2(i, span - p2 + j) = Nvals(j)
            end do
        end do

        A3 = 0.0_dp
        do i = 1, n3
            call find_span(spl%sz, x3(i), span)
            call basis_funs(spl%sz, span, x3(i), Nvals(0:p3))
            do j = 0, p3
                A3(i, span - p3 + j) = Nvals(j)
            end do
        end do

        temp = f_grid

        call dgetrf(n1, n1, A1, n1, ipiv1, info)
        if (info /= 0) then
            coeff = 0.0_dp
            return
        end if
        call dgetrf(n2, n2, A2, n2, ipiv2, info)
        if (info /= 0) then
            coeff = 0.0_dp
            return
        end if
        call dgetrf(n3, n3, A3, n3, ipiv3, info)
        if (info /= 0) then
            coeff = 0.0_dp
            return
        end if

        !$omp parallel do collapse(2) private(j, k) schedule(static)
        do k = 1, n3
            do j = 1, n2
                call dgetrs('N', n1, 1, A1, n1, ipiv1, temp(:, j, k), n1, info)
            end do
        end do

        !$omp parallel do collapse(2) private(i, k) schedule(static)
        do k = 1, n3
            do i = 1, n1
                call dgetrs('N', n2, 1, A2, n2, ipiv2, temp(i, :, k), n2, info)
            end do
        end do

        !$omp parallel do collapse(2) private(i, j) schedule(static)
        do j = 1, n2
            do i = 1, n1
                call dgetrs('N', n3, 1, A3, n3, ipiv3, temp(i, j, :), n3, info)
            end do
        end do

        coeff = temp
    end subroutine bspline_3d_interp

end module neo_bspline_interp
