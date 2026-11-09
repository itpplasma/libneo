!> The legacy uniform-mesh Lagrange interfaces (polylag_3, neo_polylag_5 and
!> the *_bdf externals) all evaluate through plag_coeff. They are checked
!> against the closed-form uniform-mesh weights they used to carry and against
!> polynomials of their own degree, which they reproduce exactly.
program test_lagrange_legacy
    use libneo_kinds, only: dp
    implicit none

    real(dp), parameter :: tol = 1.0e-12_dp
    logical :: ok

    ok = .true.
    call check_cubic_weights(ok)
    call check_quintic_weights(ok)
    call check_plag3d_cubic(ok)
    call check_plag3d_quintic(ok)
    call check_bdf(ok)
    call check_stencils(ok)
    if (.not. ok) error stop "legacy Lagrange interfaces disagree"

contains

    subroutine expect(ok, name, got, ref)
        logical, intent(inout) :: ok
        character(*), intent(in) :: name
        real(dp), intent(in) :: got, ref

        if (abs(got - ref) > tol*max(1.0_dp, abs(ref))) then
            write (*, '(a, 2es24.16)') "FAIL "//name, got, ref
            ok = .false.
        end if
    end subroutine expect

    pure subroutine closed_form(u, up, du, c, dc)
        real(dp), intent(in) :: u, up(:), du
        real(dp), intent(out) :: c(size(up)), dc(size(up))

        integer :: i, j, k, n
        real(dp) :: denom, term

        n = size(up)
        do i = 1, n
            denom = du**(n - 1)
            do j = 1, n
                if (j /= i) denom = denom*real(i - j, dp)
            end do
            c(i) = 1.0_dp
            dc(i) = 0.0_dp
            do j = 1, n
                if (j == i) cycle
                c(i) = c(i)*(u - up(j))
                term = 1.0_dp
                do k = 1, n
                    if (k /= i .and. k /= j) term = term*(u - up(k))
                end do
                dc(i) = dc(i) + term
            end do
            c(i) = c(i)/denom
            dc(i) = dc(i)/denom
        end do
    end subroutine closed_form

    subroutine check_cubic_weights(ok)
        use polylag_3, only: mp, coefs, coefs1, plag1d
        logical, intent(inout) :: ok

        real(dp), parameter :: h = 0.37_dp, u0 = -1.3_dp
        real(dp) :: up(mp), c(mp), dc(mp), cref(mp), dcref(mp), fp(mp), f, df
        integer :: i, s

        up = [(u0 + h*(i - 1), i=1, mp)]
        do s = 0, 30
            call closed_form(u0 - 0.5_dp*h + s*0.13_dp, up, h, cref, dcref)
            call coefs(u0 - 0.5_dp*h + s*0.13_dp, up, 1.0_dp/h, c)
            call coefs1(u0 - 0.5_dp*h + s*0.13_dp, up, 1.0_dp/h, dc)
            do i = 1, mp
                call expect(ok, "coefs", c(i), cref(i))
                call expect(ok, "coefs1", dc(i), dcref(i))
            end do
            fp = cubic1(up)
            call plag1d(u0 + s*0.05_dp, fp, 1.0_dp/h, up, f, df)
            call expect(ok, "plag1d f", f, cubic1(u0 + s*0.05_dp))
            call expect(ok, "plag1d df", df, dcubic1(u0 + s*0.05_dp))
        end do
    end subroutine check_cubic_weights

    subroutine check_quintic_weights(ok)
        use neo_polylag_5, only: mp, calc_coefs, calc_derivative_coefs, plag1d
        logical, intent(inout) :: ok

        real(dp), parameter :: h = 0.21_dp, u0 = 0.4_dp
        real(dp) :: up(mp), cref(mp), dcref(mp), c(mp), dc(mp), fp(mp), f, df, u
        integer :: i, s

        up = [(u0 + h*(i - 1), i=1, mp)]
        do s = 0, 30
            u = u0 - 0.3_dp*h + s*0.04_dp
            call closed_form(u, up, h, cref, dcref)
            c = calc_coefs(u, up, h)
            dc = calc_derivative_coefs(u, up, h)
            do i = 1, mp
                call expect(ok, "calc_coefs", c(i), cref(i))
                call expect(ok, "calc_derivative_coefs", dc(i), dcref(i))
            end do
            fp = quintic1(up)
            call plag1d(u, up, fp, h, f, df)
            call expect(ok, "plag1d_5 f", f, quintic1(u))
            call expect(ok, "plag1d_5 df", df, dquintic1(u))
        end do
    end subroutine check_quintic_weights

    subroutine check_plag3d_cubic(ok)
        use polylag_3, only: mp, plag3d
        logical, intent(inout) :: ok

        real(dp) :: xp(mp), yp(mp), zp(mp), fp(mp, mp, mp)
        real(dp) :: x, y, z, f, fx, fy, fz
        integer :: i, j, k

        xp = [(0.1_dp*i, i=1, mp)]
        yp = [(-0.3_dp + 0.2_dp*i, i=1, mp)]
        zp = [(1.0_dp + 0.15_dp*i, i=1, mp)]
        do k = 1, mp
            do j = 1, mp
                do i = 1, mp
                    fp(i, j, k) = cubic1(xp(i))*cubic1(yp(j))*cubic1(zp(k))
                end do
            end do
        end do
        x = 0.27_dp
        y = 0.18_dp
        z = 1.41_dp
        call plag3d(x, y, z, fp, 10.0_dp, 5.0_dp, 1.0_dp/0.15_dp, xp, yp, zp, &
                    f, fx, fy, fz)
        call expect(ok, "plag3d f", f, cubic1(x)*cubic1(y)*cubic1(z))
        call expect(ok, "plag3d fx", fx, dcubic1(x)*cubic1(y)*cubic1(z))
        call expect(ok, "plag3d fy", fy, cubic1(x)*dcubic1(y)*cubic1(z))
        call expect(ok, "plag3d fz", fz, cubic1(x)*cubic1(y)*dcubic1(z))
    end subroutine check_plag3d_cubic

    subroutine check_plag3d_quintic(ok)
        use neo_polylag_5, only: mp, get_plag3d
        logical, intent(inout) :: ok

        real(dp) :: xp(mp), yp(mp), zp(mp), fp(mp, mp, mp)
        real(dp) :: x, y, z, f, fx, fy, fz
        integer :: i, j, k

        xp = [(0.1_dp*i, i=1, mp)]
        yp = [(-0.3_dp + 0.2_dp*i, i=1, mp)]
        zp = [(0.05_dp*i, i=1, mp)]
        do k = 1, mp
            do j = 1, mp
                do i = 1, mp
                    fp(i, j, k) = quintic1(xp(i))*cubic1(yp(j))*quintic1(zp(k))
                end do
            end do
        end do
        x = 0.33_dp
        y = 0.41_dp
        z = 0.12_dp
        call get_plag3d(x, y, z, xp, yp, zp, fp, 0.1_dp, 0.2_dp, 0.05_dp, &
                        f, fx, fy, fz)
        call expect(ok, "get_plag3d f", f, quintic1(x)*cubic1(y)*quintic1(z))
        call expect(ok, "get_plag3d fx", fx, dquintic1(x)*cubic1(y)*quintic1(z))
        call expect(ok, "get_plag3d fy", fy, quintic1(x)*dcubic1(y)*quintic1(z))
        call expect(ok, "get_plag3d fz", fz, quintic1(x)*cubic1(y)*dquintic1(z))
    end subroutine check_plag3d_quintic

    subroutine check_bdf(ok)
        logical, intent(inout) :: ok

        interface
            subroutine coefs_bdf(u, up, dum1, cu)
                import :: dp
                real(dp), intent(in) :: u, up(4), dum1
                real(dp), intent(out) :: cu(4)
            end subroutine coefs_bdf
            subroutine plag2d_bdf(x, y, fp, dxm1, dym1, xp, yp, polyl2d)
                import :: dp
                real(dp), intent(in) :: x, y, fp(4, 4), dxm1, dym1, xp(4), yp(4)
                real(dp), intent(out) :: polyl2d
            end subroutine plag2d_bdf
        end interface

        real(dp) :: xp(4), yp(4), fp(4, 4), c(4), cref(4), dcref(4), p
        integer :: i, j

        xp = [(0.5_dp + 0.25_dp*i, i=1, 4)]
        yp = [(-1.0_dp + 0.4_dp*i, i=1, 4)]
        call closed_form(0.93_dp, xp, 0.25_dp, cref, dcref)
        call coefs_bdf(0.93_dp, xp, 4.0_dp, c)
        do i = 1, 4
            call expect(ok, "coefs_bdf", c(i), cref(i))
        end do
        do j = 1, 4
            do i = 1, 4
                fp(i, j) = cubic1(xp(i))*cubic1(yp(j))
            end do
        end do
        call plag2d_bdf(0.93_dp, -0.13_dp, fp, 4.0_dp, 2.5_dp, xp, yp, p)
        call expect(ok, "plag2d_bdf", p, cubic1(0.93_dp)*cubic1(-0.13_dp))
    end subroutine check_bdf

    subroutine check_stencils(ok)
        use polylag_3, only: indef, indsmp
        use neo_polylag_5, only: find_node_index
        logical, intent(inout) :: ok

        interface
            subroutine indef_bdf(u, umin, dum1, nup, indu)
                import :: dp
                real(dp), intent(in) :: u, umin, dum1
                integer, intent(in) :: nup
                integer, intent(out) :: indu(4)
            end subroutine indef_bdf
        end interface

        integer, parameter :: nup = 11
        real(dp), parameter :: umin = 2.0_dp, h = 0.5_dp
        integer :: i4(4), j4(4), i6(6), s
        real(dp) :: u

        do s = -4, 60
            u = umin + s*0.1_dp
            call indef(u, umin, 1.0_dp/h, nup, i4)
            if (any(i4 /= reference(int((u - umin)/h), 4, nup))) call flag(ok, "indef")
            call indef_bdf(u, umin, 1.0_dp/h, nup, j4)
            if (any(j4 /= i4)) call flag(ok, "indef_bdf")
            call find_node_index(u, umin, h, nup, i6)
            if (any(i6 /= reference(int((u - umin)/h) + 1, 6, nup))) &
                call flag(ok, "find_node_index")
        end do
        do s = -1, nup + 2
            call indsmp(s, nup, i4)
            if (any(i4 /= reference(s - 1, 4, nup))) call flag(ok, "indsmp")
        end do

    end subroutine check_stencils

    subroutine flag(ok, name)
        logical, intent(inout) :: ok
        character(*), intent(in) :: name

        write (*, '(a)') "FAIL "//name
        ok = .false.
    end subroutine flag

    pure function reference(first, m, nup) result(ind)
        integer, intent(in) :: first, m, nup
        integer :: ind(m), i

        ind(1) = first
        if (ind(1) <= 0) ind(1) = 1
        ind(m) = ind(1) + m - 1
        if (ind(m) > nup) then
            ind(m) = nup
            ind(1) = ind(m) - m + 1
        end if
        do i = 2, m - 1
            ind(i) = ind(i - 1) + 1
        end do
    end function reference

    elemental real(dp) function cubic1(x)
        real(dp), intent(in) :: x

        cubic1 = 2.0_dp*x**3 - x**2 + 0.5_dp*x - 3.0_dp
    end function cubic1

    elemental real(dp) function dcubic1(x)
        real(dp), intent(in) :: x

        dcubic1 = 6.0_dp*x**2 - 2.0_dp*x + 0.5_dp
    end function dcubic1

    elemental real(dp) function quintic1(x)
        real(dp), intent(in) :: x

        quintic1 = x**5 - 2.0_dp*x**4 + 0.3_dp*x**2 - x + 1.0_dp
    end function quintic1

    elemental real(dp) function dquintic1(x)
        real(dp), intent(in) :: x

        dquintic1 = 5.0_dp*x**4 - 8.0_dp*x**3 + 0.6_dp*x - 1.0_dp
    end function dquintic1

end program test_lagrange_legacy
