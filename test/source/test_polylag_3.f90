program test_polylag_3
    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    logical :: all_passed

    all_passed = .true.

    call test_indef_middle(all_passed)
    call test_indef_boundary_left(all_passed)
    call test_indef_boundary_right(all_passed)
    call test_plag1d_polynomial(all_passed)
    call test_plag1d_derivative(all_passed)
    call test_coefs_partition_of_unity(all_passed)
    call test_plag3d_polynomial(all_passed)
    call test_plag3d_derivatives(all_passed)

    if (.not. all_passed) then
        error stop "One or more polylag_3 tests failed"
    end if

contains

    subroutine test_indef_middle(passed)
        use polylag_3, only: indef, mp
        logical, intent(inout) :: passed

        integer :: indu(mp)
        real(dp) :: u, umin, dum1
        integer :: nup

        umin = 0.0d0
        dum1 = 10.0d0
        nup = 20
        u = 0.55d0

        call indef(u, umin, dum1, nup, indu)

        if (indu(1) /= 5) then
            write(*,*) "FAIL: indef middle. indu(1)=", indu(1), " expected 5"
            passed = .false.
            return
        end if

        if (indu(mp) /= 8) then
            write(*,*) "FAIL: indef middle. indu(mp)=", indu(mp), " expected 8"
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_indef_middle"
    end subroutine test_indef_middle

    subroutine test_indef_boundary_left(passed)
        use polylag_3, only: indef, mp
        logical, intent(inout) :: passed

        integer :: indu(mp)
        real(dp) :: u, umin, dum1
        integer :: nup

        umin = 0.0d0
        dum1 = 10.0d0
        nup = 20
        u = 0.05d0

        call indef(u, umin, dum1, nup, indu)

        if (indu(1) /= 1) then
            write(*,*) "FAIL: indef left boundary. indu(1)=", indu(1), " expected 1"
            passed = .false.
            return
        end if

        if (indu(mp) /= 4) then
            write(*,*) "FAIL: indef left boundary. indu(mp)=", indu(mp), " expected 4"
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_indef_boundary_left"
    end subroutine test_indef_boundary_left

    subroutine test_indef_boundary_right(passed)
        use polylag_3, only: indef, mp
        logical, intent(inout) :: passed

        integer :: indu(mp)
        real(dp) :: u, umin, dum1
        integer :: nup

        umin = 0.0d0
        dum1 = 10.0d0
        nup = 20
        u = 1.95d0

        call indef(u, umin, dum1, nup, indu)

        if (indu(mp) /= 20) then
            write(*,*) "FAIL: indef right boundary. indu(mp)=", indu(mp), " expected 20"
            passed = .false.
            return
        end if

        if (indu(1) /= 17) then
            write(*,*) "FAIL: indef right boundary. indu(1)=", indu(1), " expected 17"
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_indef_boundary_right"
    end subroutine test_indef_boundary_right

    subroutine test_plag1d_polynomial(passed)
        use polylag_3, only: plag1d, mp
        logical, intent(inout) :: passed

        real(dp) :: x, fp(mp), xp(mp), polyl1d, p1x1d
        real(dp) :: dxm1, dx
        real(dp) :: expected
        integer :: i

        dx = 0.1d0
        dxm1 = 1.0d0 / dx

        do i = 1, mp
            xp(i) = 0.5d0 + (i - 1) * dx
            fp(i) = xp(i)**3
        end do

        x = 0.65d0
        expected = x**3

        call plag1d(x, fp, dxm1, xp, polyl1d, p1x1d)

        if (abs(polyl1d - expected) > 1.0d-10) then
            write(*,*) "FAIL: plag1d cubic polynomial. Got:", polyl1d, " expected:", expected
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_plag1d_polynomial"
    end subroutine test_plag1d_polynomial

    subroutine test_plag1d_derivative(passed)
        use polylag_3, only: plag1d, mp
        logical, intent(inout) :: passed

        real(dp) :: x, fp(mp), xp(mp), polyl1d, p1x1d
        real(dp) :: dxm1, dx
        real(dp) :: expected_deriv
        integer :: i

        dx = 0.1d0
        dxm1 = 1.0d0 / dx

        do i = 1, mp
            xp(i) = 0.5d0 + (i - 1) * dx
            fp(i) = xp(i)**3
        end do

        x = 0.65d0
        expected_deriv = 3.0d0 * x**2

        call plag1d(x, fp, dxm1, xp, polyl1d, p1x1d)

        if (abs(p1x1d - expected_deriv) > 1.0d-8) then
            write(*,*) "FAIL: plag1d derivative. Got:", p1x1d, " expected:", expected_deriv
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_plag1d_derivative"
    end subroutine test_plag1d_derivative

    subroutine test_coefs_partition_of_unity(passed)
        use polylag_3, only: coefs, mp
        logical, intent(inout) :: passed

        real(dp) :: u, up(mp), cu(mp), dum1, dx
        real(dp) :: coef_sum
        integer :: i

        dx = 0.1d0
        dum1 = 1.0d0 / dx

        do i = 1, mp
            up(i) = 0.5d0 + (i - 1) * dx
        end do

        u = 0.65d0
        call coefs(u, up, dum1, cu)

        coef_sum = sum(cu)

        if (abs(coef_sum - 1.0d0) > 1.0d-12) then
            write(*,*) "FAIL: coefs partition of unity. Sum:", coef_sum
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_coefs_partition_of_unity"
    end subroutine test_coefs_partition_of_unity

    subroutine test_plag3d_polynomial(passed)
        use polylag_3, only: plag3d, mp
        logical, intent(inout) :: passed

        real(dp) :: x, y, z
        real(dp) :: fp(mp, mp, mp), xp(mp), yp(mp), zp(mp)
        real(dp) :: polyl3d, poly1x, poly1y, poly1z
        real(dp) :: dxm1, dym1, dzm1, dx
        real(dp) :: expected
        integer :: i, j, k

        dx = 0.1d0
        dxm1 = 1.0d0 / dx
        dym1 = 1.0d0 / dx
        dzm1 = 1.0d0 / dx

        do i = 1, mp
            xp(i) = 0.5d0 + (i - 1) * dx
            yp(i) = 0.5d0 + (i - 1) * dx
            zp(i) = 0.5d0 + (i - 1) * dx
        end do

        do i = 1, mp
            do j = 1, mp
                do k = 1, mp
                    fp(i, j, k) = xp(i) * yp(j) * zp(k)
                end do
            end do
        end do

        x = 0.65d0
        y = 0.55d0
        z = 0.72d0
        expected = x * y * z

        call plag3d(x, y, z, fp, dxm1, dym1, dzm1, xp, yp, zp, &
                    polyl3d, poly1x, poly1y, poly1z)

        if (abs(polyl3d - expected) > 1.0d-10) then
            write(*,*) "FAIL: plag3d trilinear. Got:", polyl3d, " expected:", expected
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_plag3d_polynomial"
    end subroutine test_plag3d_polynomial

    subroutine test_plag3d_derivatives(passed)
        use polylag_3, only: plag3d, mp
        logical, intent(inout) :: passed

        real(dp) :: x, y, z
        real(dp) :: fp(mp, mp, mp), xp(mp), yp(mp), zp(mp)
        real(dp) :: polyl3d, poly1x, poly1y, poly1z
        real(dp) :: dxm1, dym1, dzm1, dx
        real(dp) :: expected_dx, expected_dy, expected_dz
        integer :: i, j, k

        dx = 0.1d0
        dxm1 = 1.0d0 / dx
        dym1 = 1.0d0 / dx
        dzm1 = 1.0d0 / dx

        do i = 1, mp
            xp(i) = 0.5d0 + (i - 1) * dx
            yp(i) = 0.5d0 + (i - 1) * dx
            zp(i) = 0.5d0 + (i - 1) * dx
        end do

        do i = 1, mp
            do j = 1, mp
                do k = 1, mp
                    fp(i, j, k) = xp(i) * yp(j) * zp(k)
                end do
            end do
        end do

        x = 0.65d0
        y = 0.55d0
        z = 0.72d0

        expected_dx = y * z
        expected_dy = x * z
        expected_dz = x * y

        call plag3d(x, y, z, fp, dxm1, dym1, dzm1, xp, yp, zp, &
                    polyl3d, poly1x, poly1y, poly1z)

        if (abs(poly1x - expected_dx) > 1.0d-8) then
            write(*,*) "FAIL: plag3d dx. Got:", poly1x, " expected:", expected_dx
            passed = .false.
            return
        end if

        if (abs(poly1y - expected_dy) > 1.0d-8) then
            write(*,*) "FAIL: plag3d dy. Got:", poly1y, " expected:", expected_dy
            passed = .false.
            return
        end if

        if (abs(poly1z - expected_dz) > 1.0d-8) then
            write(*,*) "FAIL: plag3d dz. Got:", poly1z, " expected:", expected_dz
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_plag3d_derivatives"
    end subroutine test_plag3d_derivatives

end program test_polylag_3
