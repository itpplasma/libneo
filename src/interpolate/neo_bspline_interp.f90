module neo_bspline_interp
    !! Direct B-spline interpolation via collocation.
    !! Solves the linear system A*c = f where A_ij = B_j(x_i).
    !! Uses separable tensor-product structure for efficiency in 2D/3D.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neo_bspline_base, only: bspline_1d, bspline_2d, bspline_3d, &
        find_span, basis_funs
    implicit none
    private

    interface
        subroutine dgbsv(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info)
            import :: dp
            integer :: n, kl, ku, nrhs, ldab, ldb, info
            integer :: ipiv(*)
            real(dp) :: ab(ldab, *), b(ldb, *)
        end subroutine dgbsv
    end interface

    public :: bspline_1d_interp, bspline_2d_interp, bspline_3d_interp

contains

    subroutine bspline_1d_interp(spl, x_data, f_data, coeff)
        !! Direct 1D B-spline interpolation: solve collocation system A*c = f.
        type(bspline_1d), intent(in) :: spl
        real(dp), intent(in) :: x_data(:)
        real(dp), intent(in) :: f_data(:)
        real(dp), intent(out) :: coeff(:)

        integer :: n, p, i, j, col, span, info, kl, ku, ldab
        real(dp), allocatable :: ab(:,:), work_f(:), Nvals(:)
        integer, allocatable :: ipiv(:)

        n = size(x_data)
        p = spl%degree
        kl = p
        ku = p
        ldab = 2*kl + ku + 1
        allocate(ab(ldab, n), work_f(n), ipiv(n), Nvals(0:p))
        ab = 0.0_dp

        do i = 1, n
            call find_span(spl, x_data(i), span)
            call basis_funs(spl, span, x_data(i), Nvals)
            do j = 0, p
                col = span - p + j
                ab(kl + ku + 1 + i - col, col) = Nvals(j)
            end do
        end do

        work_f = f_data
        call dgbsv(n, kl, ku, 1, ab, ldab, ipiv, work_f, n, info)
        if (info /= 0) then
            coeff = 0.0_dp
            return
        end if
        coeff = work_f
    end subroutine bspline_1d_interp


    subroutine bspline_2d_interp(spl, x1, x2, f_grid, coeff)
        !! Direct 2D B-spline interpolation using separable tensor product.
        !! Solves dimension-by-dimension: first along x1, then along x2.
        !! Each collocation system is banded (bandwidth = degree), so dgbsv
        !! solves it in O(n) instead of the O(n^3) dense dgetrf/dgetrs.
        type(bspline_2d), intent(in) :: spl
        real(dp), intent(in) :: x1(:), x2(:)
        real(dp), intent(in) :: f_grid(:,:)
        real(dp), intent(out) :: coeff(:,:)

        integer :: n1, n2, i, j, col, span, info
        integer :: p1, p2, kl1, ku1, ldab1, kl2, ku2, ldab2
        real(dp), allocatable :: ab1(:,:), ab2(:,:), temp(:,:), temp_t(:,:)
        real(dp), allocatable :: Nvals(:)
        integer, allocatable :: ipiv(:)

        n1 = size(x1)
        n2 = size(x2)
        p1 = spl%sx%degree
        p2 = spl%sy%degree
        kl1 = p1
        ku1 = p1
        ldab1 = 2*kl1 + ku1 + 1
        kl2 = p2
        ku2 = p2
        ldab2 = 2*kl2 + ku2 + 1
        allocate(ab1(ldab1, n1), ab2(ldab2, n2))
        allocate(temp(n1, n2), temp_t(n2, n1))
        allocate(ipiv(max(n1, n2)), Nvals(0:max(p1, p2)))

        ab1 = 0.0_dp
        do i = 1, n1
            call find_span(spl%sx, x1(i), span)
            call basis_funs(spl%sx, span, x1(i), Nvals(0:p1))
            do j = 0, p1
                col = span - p1 + j
                ab1(kl1 + ku1 + 1 + i - col, col) = Nvals(j)
            end do
        end do

        ab2 = 0.0_dp
        do i = 1, n2
            call find_span(spl%sy, x2(i), span)
            call basis_funs(spl%sy, span, x2(i), Nvals(0:p2))
            do j = 0, p2
                col = span - p2 + j
                ab2(kl2 + ku2 + 1 + i - col, col) = Nvals(j)
            end do
        end do

        temp = f_grid
        call dgbsv(n1, kl1, ku1, n2, ab1, ldab1, ipiv(1:n1), temp, n1, info)
        if (info /= 0) then
            coeff = 0.0_dp
            return
        end if

        temp_t = transpose(temp)
        call dgbsv(n2, kl2, ku2, n1, ab2, ldab2, ipiv(1:n2), temp_t, n2, info)
        if (info /= 0) then
            coeff = 0.0_dp
            return
        end if

        coeff = transpose(temp_t)
    end subroutine bspline_2d_interp


    subroutine bspline_3d_interp(spl, x1, x2, x3, f_grid, coeff)
        !! Direct 3D B-spline interpolation using separable tensor product.
        !! Solves dimension-by-dimension: x1 -> x2 -> x3.
        !! Each collocation system is banded (bandwidth = degree), so dgbsv
        !! solves it in O(n) instead of the O(n^3) dense dgetrf/dgetrs.
        type(bspline_3d), intent(in) :: spl
        real(dp), intent(in) :: x1(:), x2(:), x3(:)
        real(dp), intent(in) :: f_grid(:,:,:)
        real(dp), intent(out) :: coeff(:,:,:)

        integer :: n1, n2, n3, i, j, k, col, span, info
        integer :: p1, p2, p3, kl1, ku1, ldab1, kl2, ku2, ldab2, kl3, ku3, ldab3
        real(dp), allocatable :: ab1(:,:), ab2(:,:), ab3(:,:)
        real(dp), allocatable :: temp(:,:,:), perm(:,:,:), Nvals(:)
        integer, allocatable :: ipiv(:)

        n1 = size(x1)
        n2 = size(x2)
        n3 = size(x3)
        p1 = spl%sx%degree
        p2 = spl%sy%degree
        p3 = spl%sz%degree
        kl1 = p1; ku1 = p1; ldab1 = 2*kl1 + ku1 + 1
        kl2 = p2; ku2 = p2; ldab2 = 2*kl2 + ku2 + 1
        kl3 = p3; ku3 = p3; ldab3 = 2*kl3 + ku3 + 1
        allocate(ab1(ldab1, n1), ab2(ldab2, n2), ab3(ldab3, n3))
        allocate(temp(n1, n2, n3))
        allocate(ipiv(max(n1, n2, n3)), Nvals(0:max(p1, p2, p3)))

        ab1 = 0.0_dp
        do i = 1, n1
            call find_span(spl%sx, x1(i), span)
            call basis_funs(spl%sx, span, x1(i), Nvals(0:p1))
            do j = 0, p1
                col = span - p1 + j
                ab1(kl1 + ku1 + 1 + i - col, col) = Nvals(j)
            end do
        end do

        ab2 = 0.0_dp
        do i = 1, n2
            call find_span(spl%sy, x2(i), span)
            call basis_funs(spl%sy, span, x2(i), Nvals(0:p2))
            do j = 0, p2
                col = span - p2 + j
                ab2(kl2 + ku2 + 1 + i - col, col) = Nvals(j)
            end do
        end do

        ab3 = 0.0_dp
        do i = 1, n3
            call find_span(spl%sz, x3(i), span)
            call basis_funs(spl%sz, span, x3(i), Nvals(0:p3))
            do j = 0, p3
                col = span - p3 + j
                ab3(kl3 + ku3 + 1 + i - col, col) = Nvals(j)
            end do
        end do

        temp = f_grid

        call dgbsv(n1, kl1, ku1, n2*n3, ab1, ldab1, ipiv(1:n1), temp, n1, info)
        if (info /= 0) then
            coeff = 0.0_dp
            return
        end if

        allocate(perm(n2, n1, n3))
        do k = 1, n3
            perm(:, :, k) = transpose(temp(:, :, k))
        end do
        call dgbsv(n2, kl2, ku2, n1*n3, ab2, ldab2, ipiv(1:n2), perm, n2, info)
        if (info /= 0) then
            coeff = 0.0_dp
            return
        end if
        do k = 1, n3
            temp(:, :, k) = transpose(perm(:, :, k))
        end do
        deallocate(perm)

        allocate(perm(n3, n1, n2))
        do j = 1, n2
            do i = 1, n1
                perm(:, i, j) = temp(i, j, :)
            end do
        end do
        call dgbsv(n3, kl3, ku3, n1*n2, ab3, ldab3, ipiv(1:n3), perm, n3, info)
        if (info /= 0) then
            coeff = 0.0_dp
            return
        end if
        do j = 1, n2
            do i = 1, n1
                coeff(i, j, :) = perm(:, i, j)
            end do
        end do
    end subroutine bspline_3d_interp

end module neo_bspline_interp
