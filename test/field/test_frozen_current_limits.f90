!> Limit cases of the frozen-current field with an artificial boundary in a
!> three-dimensional coil field.  Ten rectangular toroidal-field coils give a
!> field with a strong n = 10 ripple; a torus inside them is taken as the
!> "equilibrium" volume V with the coil field as its field.
!>   Vacuum: the source carries no current, so A_p must vanish to the
!>   quadrature error (it falls with the source resolution).
!>   Known current: ten azimuthal ball currents, one per period, added to
!>   the coil field must give their dipole field outside the balls.
program test_frozen_current_limits
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neo_biotsavart, only: coils_t, coils_init, clight
    use neo_equilibrium_source, only: source_points_t, source_resolution_t, &
                                      flux_quadrature
    use neo_free_boundary_coils, only: coil_field_cyl
    use neo_frozen_current_field, only: frozen_current_patch_t, &
                                        frozen_current_plasma_t, solve_frozen_current, &
                                        frozen_current_plasma_cyl
    use neo_frozen_current_leaves, only: cart_to_cyl
    implicit none

    real(dp), parameter :: pi = 3.141592653589793238462643383279502884_dp
    integer, parameter :: nfp = 10
    real(dp), parameter :: r_axis = 3.0_dp, b_minor = 1.0_dp
    real(dp), parameter :: a_ball = 0.5_dp, m_dip = 1.0_dp

    type(coils_t) :: coils
    type(source_resolution_t) :: res
    type(source_points_t) :: vol, surf
    type(frozen_current_patch_t) :: patch
    type(frozen_current_plasma_t) :: plasma
    real(dp) :: probes(3, 4), a(3), b(3), bref(3), err
    integer :: status, i

    call tf_coils(coils)
    res%ns = 32
    res%ntheta = 96
    res%nzeta = 32
    call torus_points(res, vol, surf)
    patch%r0 = 1.5_dp
    patch%r1 = 4.5_dp
    patch%z0 = -1.5_dp
    patch%z1 = 1.5_dp
    patch%nr = 24
    patch%nz = 24
    patch%degree = 3
    patch%kmax = 4

    ! two points inside V, between the balls, and two in the vacuum around V
    probes(:, 1) = [3.0_dp, pi/nfp, 0.6_dp]
    probes(:, 2) = [3.7_dp, 0.25_dp, -0.3_dp]
    probes(:, 3) = [4.3_dp, 0.05_dp, 0.4_dp]
    probes(:, 4) = [2.4_dp, 0.4_dp, 1.3_dp]

    call set_source(vol, .false.)
    call set_source(surf, .false.)
    call solve_frozen_current(patch, vol, surf, plasma, status)
    call check(status == 0, "vacuum solve")
    do i = 1, size(probes, 2)
        call frozen_current_plasma_cyl(plasma, probes(:, i), a, b, status)
        err = norm2(b)/norm2(scaled_coil_field(probes(:, i)))
        print '(a,3f7.3,a,es10.2)', " vacuum probe", probes(:, i), &
            "  |curl A_p|/|B_coils|", err
        call check(err < merge(5.0e-3_dp, 1.0e-4_dp, i <= 2), &
                   "zero current gives no field")
    end do

    call set_source(vol, .true.)
    call set_source(surf, .true.)
    patch%kmax = 6
    call solve_frozen_current(patch, vol, surf, plasma, status)
    call check(status == 0, "ball solve")
    print *, "max gauge fraction", maxval(plasma%gauge)
    call check(maxval(plasma%gauge) < 1.0e-3_dp, "small multiplier")
    do i = 1, size(probes, 2)
        call frozen_current_plasma_cyl(plasma, probes(:, i), a, b, status)
        bref = balls_cyl(probes(:, i))
        err = norm2(b - bref)/norm2(bref)
        print '(a,3f7.3,a,es10.2)', " ball probe", probes(:, i), "  |dB|/|B|", err
        call check(err < 5.0e-3_dp, "ball currents give their dipole fields")
    end do
    print *, "PASS"

contains

    !> nfp coils with legs at R = 1.2 m and 4.8 m, |Z| <= 2.5 m, 1 MA each,
    !> in cm and statA.
    subroutine tf_coils(coils)
        type(coils_t), intent(out) :: coils
        real(dp) :: x(5*nfp), y(5*nfp), z(5*nfp), c(5*nfp), rz(2, 5), phi
        integer :: k, j, m

        rz = reshape([120.0_dp, -250.0_dp, 480.0_dp, -250.0_dp, 480.0_dp, 250.0_dp, &
                      120.0_dp, 250.0_dp, 120.0_dp, -250.0_dp], [2, 5])
        m = 0
        do k = 1, nfp
            phi = 2.0_dp*pi*(k - 0.5_dp)/nfp
            do j = 1, 5
                m = m + 1
                x(m) = rz(1, j)*cos(phi)
                y(m) = rz(1, j)*sin(phi)
                z(m) = rz(2, j)
                c(m) = merge(0.0_dp, 1.0e6_dp*clight/10.0_dp, j == 5)
            end do
        end do
        call coils_init(x, y, z, c, coils)
    end subroutine tf_coils

    !> Torus R = R0 + rho cos t, Z = rho sin t, rho = b sqrt(s), one period;
    !> lengths in m (the coil field is evaluated at 100 x).
    subroutine torus_points(res, vol, surf)
        type(source_resolution_t), intent(in) :: res
        type(source_points_t), intent(out) :: vol, surf
        real(dp), allocatable :: s(:), ws(:), th(:), ze(:)
        real(dp) :: dth, dze, rho, r
        integer :: i, j, k, m, n

        call flux_quadrature(res, nfp, s, ws, th, ze, dth, dze)
        n = res%ns*res%ntheta*res%nzeta
        vol%nfp = nfp
        surf%nfp = nfp
        allocate (vol%x(3, n), vol%w(n), vol%b(3, n))
        allocate (surf%x(3, res%ntheta*res%nzeta), surf%b(3, res%ntheta*res%nzeta))
        allocate (surf%area(3, res%ntheta*res%nzeta))
        do k = 1, res%nzeta
            do j = 1, res%ntheta
                do i = 1, res%ns
                    m = i + res%ns*((j - 1) + res%ntheta*(k - 1))
                    rho = b_minor*sqrt(s(i))
                    r = r_axis + rho*cos(th(j))
                    vol%x(:, m) = [r, ze(k), rho*sin(th(j))]
                    vol%w(m) = 0.5_dp*r*b_minor**2*ws(i)*dth*dze
                end do
                m = j + res%ntheta*(k - 1)
                r = r_axis + b_minor*cos(th(j))
                surf%x(:, m) = [r, ze(k), b_minor*sin(th(j))]
                surf%area(:, m) = [cos(th(j)), 0.0_dp, sin(th(j))]*r*b_minor*dth*dze
            end do
        end do
    end subroutine torus_points

    !> Source field: the scaled coil field, plus the balls when asked.
    subroutine set_source(pts, with_balls)
        type(source_points_t), intent(inout) :: pts
        logical, intent(in) :: with_balls
        integer :: m

        do m = 1, size(pts%x, 2)
            pts%b(:, m) = scaled_coil_field(pts%x(:, m))
            if (with_balls) pts%b(:, m) = pts%b(:, m) + balls_cyl(pts%x(:, m))
        end do
    end subroutine set_source

    !> Coil field at x (m), scaled to 0.1 at R = 3 m on the midplane.
    function scaled_coil_field(x) result(b)
        real(dp), intent(in) :: x(3)
        real(dp) :: b(3)

        b = coil_field_cyl(coils, x*[100.0_dp, 1.0_dp, 100.0_dp])* &
            (0.1_dp/norm2(coil_field_cyl(coils, [300.0_dp, 0.0_dp, 0.0_dp])))
    end function scaled_coil_field

    !> Field of the ball currents centred at R0, phi = 2 pi k/nfp, Z = 0.
    function balls_cyl(x) result(b)
        real(dp), intent(in) :: x(3)
        real(dp) :: b(3), d(3), ap(3), bp(3), bc(3), phi
        integer :: k

        bc = 0.0_dp
        do k = 0, nfp - 1
            phi = 2.0_dp*pi*k/nfp
            d = [x(1)*cos(x(2)) - r_axis*cos(phi), x(1)*sin(x(2)) - r_axis*sin(phi), x(3)]
            if (norm2(d) < a_ball) then
                call mms_ball_inner(d(1), d(2), d(3), m_dip, a_ball, ap(1), ap(2), &
                                    ap(3), bp(1), bp(2), bp(3))
            else
                call mms_ball_outer(d(1), d(2), d(3), m_dip, a_ball, ap(1), ap(2), &
                                    ap(3), bp(1), bp(2), bp(3))
            end if
            bc = bc + bp
        end do
        call cart_to_cyl(cos(x(2)), sin(x(2)), bc(1), bc(2), b(1), b(2))
        b(3) = bc(3)
    end function balls_cyl

    subroutine check(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what

        if (.not. ok) then
            print *, "FAIL: ", what
            error stop 1
        end if
    end subroutine check

    include '../../src/generated/frozen_current_mms_leaves.inc'

end program test_frozen_current_limits
