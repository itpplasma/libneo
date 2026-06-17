program test_flux_pseudocartesian
    !> Pseudo-Cartesian near-axis chart (X,Y) = (sqrt(s)cos th, sqrt(s)sin th):
    !> round-trip identity, analytic Jacobians vs finite differences, and the
    !> defining property that the inverse map (X,Y)->(s,theta) stays smooth at
    !> the axis where the forward flux chart is singular.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use flux_pseudocartesian, only: flux_to_pseudocart, pseudocart_to_flux
    implicit none

    call test_round_trip()
    call test_inverse_jacobian_fd()
    call test_forward_jacobian_fd()
    call test_inverse_smooth_at_axis()

    print *, 'flux pseudo-Cartesian chart: all checks passed'

contains

    subroutine test_round_trip()
        real(dp) :: u(3), xy(3), u2(3)
        real(dp) :: svals(5)
        integer :: it
        svals = [1.0e-8_dp, 1.0e-3_dp, 1.0e-2_dp, 0.1_dp, 0.9_dp]
        do it = 1, 5
            u = [svals(it), 0.7_dp, 0.3_dp]
            call flux_to_pseudocart(u, xy)
            call pseudocart_to_flux(xy, u2)
            if (abs(u2(1) - u(1)) + abs(u2(2) - u(2)) + abs(u2(3) - u(3)) > 1.0e-12_dp) then
                print *, 's =', svals(it), ' round-trip residual too large'
                error stop 'pseudo-Cartesian round trip failed'
            end if
        end do
    end subroutine test_round_trip

    subroutine test_inverse_jacobian_fd()
        real(dp) :: xy(3), u(3), J(3, 3), Jfd(3, 3)
        real(dp) :: xp(3), xm(3), up(3), um(3), h
        integer :: it, k
        real(dp) :: xys(2, 4)
        xys = reshape([0.2_dp, 0.1_dp, 0.4_dp, -0.3_dp, 0.05_dp, 0.0_dp, &
                       -0.1_dp, 0.25_dp], [2, 4])
        h = 1.0e-6_dp
        do it = 1, 4
            xy = [xys(1, it), xys(2, it), 0.3_dp]
            call pseudocart_to_flux(xy, u, J)
            do k = 1, 3
                xp = xy; xm = xy
                xp(k) = xy(k) + h; xm(k) = xy(k) - h
                call pseudocart_to_flux(xp, up)
                call pseudocart_to_flux(xm, um)
                Jfd(:, k) = (up - um)/(2.0_dp*h)
            end do
            if (maxval(abs(J - Jfd)) > 1.0e-4_dp) then
                print *, 'inverse Jacobian FD mismatch =', maxval(abs(J - Jfd))
                error stop 'pseudocart_to_flux Jacobian wrong'
            end if
        end do
    end subroutine test_inverse_jacobian_fd

    subroutine test_forward_jacobian_fd()
        real(dp) :: u(3), xy(3), J(3, 3), Jfd(3, 3)
        real(dp) :: up(3), um(3), xp(3), xm(3), h
        integer :: k
        ! away from the axis the forward Jacobian is also finite and correct
        u = [0.25_dp, 0.6_dp, 0.2_dp]
        h = 1.0e-6_dp
        call flux_to_pseudocart(u, xy, J)
        do k = 1, 3
            up = u; um = u
            up(k) = u(k) + h; um(k) = u(k) - h
            call flux_to_pseudocart(up, xp)
            call flux_to_pseudocart(um, xm)
            Jfd(:, k) = (xp - xm)/(2.0_dp*h)
        end do
        if (maxval(abs(J - Jfd)) > 1.0e-4_dp) then
            print *, 'forward Jacobian FD mismatch =', maxval(abs(J - Jfd))
            error stop 'flux_to_pseudocart Jacobian wrong'
        end if
    end subroutine test_forward_jacobian_fd

    subroutine test_inverse_smooth_at_axis()
        real(dp) :: xy(3), u(3), J(3, 3)
        real(dp) :: rr
        integer :: i
        ! ds/dX, ds/dY stay bounded as (X,Y) -> 0; that is the whole point.
        do i = 1, 6
            rr = 10.0_dp**(-i)
            xy = [rr, 0.5_dp*rr, 0.0_dp]
            call pseudocart_to_flux(xy, u, J)
            if (abs(J(1, 1)) > 1.0_dp .or. abs(J(1, 2)) > 1.0_dp) then
                print *, 'radius =', rr, ' ds/dX =', J(1, 1), ' ds/dY =', J(1, 2)
                error stop 'inverse radial Jacobian not bounded near axis'
            end if
        end do
    end subroutine test_inverse_smooth_at_axis

end program test_flux_pseudocartesian
