!> Frozen-current field of a manufactured plasma current: an azimuthal
!> current in a ball of radius a centred off the torus axis, inside a torus
!> volume V that supplies the source points (nfp = 1, so every toroidal mode
!> enters).  Outside the ball its field is that of a dipole, so curl A_p must
!> match the dipole inside V outside the ball and in the vacuum around V, and
!> the field of a coil, which is curl-free in V, must give no A_p.  The
!> multiplier is small because the current is divergence-free and vanishes
!> on the boundary.
program test_frozen_current_field
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neo_biotsavart, only: coils_t, coils_init
    use neo_equilibrium_source, only: source_points_t, source_resolution_t, &
                                      flux_quadrature
    use neo_free_boundary_coils, only: coil_field_cyl
    use neo_frozen_current_field, only: frozen_current_patch_t, &
                                        frozen_current_plasma_t, solve_frozen_current, &
                                        frozen_current_plasma_cyl, &
                                        frozen_current_field_t
    use neo_frozen_current_leaves, only: cart_to_cyl, cyl_to_cart
    use neo_frozen_current_io, only: save_frozen_current_field, load_frozen_current_field
    implicit none

    real(dp), parameter :: pi = 3.141592653589793238462643383279502884_dp
    real(dp), parameter :: r_axis = 3.0_dp, b_minor = 1.0_dp
    real(dp), parameter :: c_ball(3) = [3.0_dp, 0.0_dp, 0.05_dp], a_ball = 0.9_dp
    real(dp), parameter :: m_dip = 1.0_dp

    type(source_resolution_t) :: res
    type(source_points_t) :: vol, surf
    type(frozen_current_patch_t) :: patch
    type(frozen_current_plasma_t) :: plasma
    type(frozen_current_field_t) :: field, loaded
    real(dp) :: bl(3), dbl(3, 3)
    type(coils_t) :: loop
    real(dp) :: probes(3, 4), err, emax, bref(3), a(3), b(3), ac(3), bc(3)
    real(dp) :: xc(3), dbc(3, 3), fd(3, 3), bp(3), bm(3), e(3), scale
    real(dp), parameter :: h = 1.0e-6_dp
    integer :: status, i, j

    res%ns = 40
    res%ntheta = 128
    res%nzeta = 96
    call torus_points(res, vol, surf)
    call ball_field(vol)
    call ball_field(surf)

    patch%r0 = 1.5_dp
    patch%r1 = 4.5_dp
    patch%z0 = -1.5_dp
    patch%z1 = 1.5_dp
    patch%nr = 24
    patch%nz = 24
    patch%degree = 3
    patch%kmax = 24
    call solve_frozen_current(patch, vol, surf, plasma, status)
    call check(status == 0, "solve")
    print *, "max gauge fraction", maxval(plasma%gauge)
    call check(maxval(plasma%gauge) < 1.0e-3_dp, "small multiplier")

    ! two points inside V outside the ball, two in the vacuum around V
    probes(:, 1) = [3.0_dp, 0.0_dp, -0.97_dp]
    probes(:, 2) = [3.8_dp, 0.4_dp, -0.2_dp]
    probes(:, 3) = [4.3_dp, 0.05_dp, 0.4_dp]
    probes(:, 4) = [2.4_dp, -0.6_dp, 1.3_dp]
    do i = 1, size(probes, 2)
        call frozen_current_plasma_cyl(plasma, probes(:, i), a, b, status)
        call check(status == 0, "evaluation in the patch")
        bref = dipole_cyl(probes(:, i))
        err = norm2(b - bref)/norm2(bref)
        print '(a,3f7.3,a,es10.2)', " probe", probes(:, i), "  |dB|/|B|", err
        emax = merge(5.0e-2_dp, 1.0e-2_dp, i <= 2)
        call check(err < emax, "curl A_p matches the dipole field")
    end do

    ! field_t without coils: Cartesian plasma part
    field%plasma = plasma
    call field%compute_abfield(cart(probes(:, 2)), ac, bc)
    call frozen_current_plasma_cyl(plasma, probes(:, 2), a, b, status)
    call check(norm2(bc - cart_vec(probes(:, 2), b)) < 1.0e-12_dp*norm2(b), &
               "field_t matches the cylindrical evaluation")

    ! exact Jacobian of the plasma part plus a coil loop
    call circular_loop(loop, 3.0_dp, 2.0_dp, 1.0e3_dp)
    field%coils = loop
    do i = 1, size(probes, 2)
        xc = cart(probes(:, i))
        call field%compute_field_jacobian(xc, bc, dbc)
        call field%compute_abfield(xc, ac, bref)
        call check(norm2(bc - bref) < 1.0e-13_dp*norm2(bref), &
                   "the Jacobian route gives the same B")
        do j = 1, 3
            e = 0.0_dp
            e(j) = h
            call field%compute_abfield(xc + e, ac, bp)
            call field%compute_abfield(xc - e, ac, bm)
            fd(:, j) = (bp - bm)/(2*h)
        end do
        scale = maxval(abs(dbc))
        print '(a,i2,a,es10.2,a,es10.2)', " probe", i, "  |dB - fd|", &
            maxval(abs(dbc - fd))/scale, "  |div B|", &
            abs(dbc(1, 1) + dbc(2, 2) + dbc(3, 3))/scale
        call check(maxval(abs(dbc - fd)) < 1.0e-6_dp*scale, &
                   "dB matches central differences")
        call check(abs(dbc(1, 1) + dbc(2, 2) + dbc(3, 3)) < 1.0e-12_dp*scale, &
                   "div B vanishes")
    end do

    ! a saved and loaded field evaluates identically
    call save_frozen_current_field('frozen_current_test.h5', patch, field)
    call load_frozen_current_field('frozen_current_test.h5', loaded, status)
    call check(status == 0, "load")
    do i = 1, size(probes, 2)
        xc = cart(probes(:, i))
        call field%compute_field_jacobian(xc, bc, dbc)
        call loaded%compute_field_jacobian(xc, bl, dbl)
        call check(all(bl == bc) .and. all(dbl == dbc), &
                   "the loaded field equals the saved one")
    end do

    ! a curl-free source: the field of a coil loop around V
    call circular_loop(loop, 3.0_dp, 2.0_dp, 1.0e3_dp)
    do i = 1, size(vol%x, 2)
        vol%b(:, i) = coil_field_cyl(loop, vol%x(:, i))
    end do
    do i = 1, size(surf%x, 2)
        surf%b(:, i) = coil_field_cyl(loop, surf%x(:, i))
    end do
    patch%kmax = 2
    call solve_frozen_current(patch, vol, surf, plasma, status)
    call check(status == 0, "solve")
    call frozen_current_plasma_cyl(plasma, probes(:, 2), a, b, status)
    bref = coil_field_cyl(loop, probes(:, 2))
    print *, "coil source: |curl A_p|/|B_coil|", norm2(b)/norm2(bref)
    call check(norm2(b) < 1.0e-3_dp*norm2(bref), "a curl-free source gives no field")
    print *, "PASS"

contains

    !> Torus R = R0 + rho cos t, Z = rho sin t, rho = b sqrt(s), phi = zeta.
    subroutine torus_points(res, vol, surf)
        type(source_resolution_t), intent(in) :: res
        type(source_points_t), intent(out) :: vol, surf
        real(dp), allocatable :: s(:), ws(:), th(:), ze(:)
        real(dp) :: dth, dze, rho, r
        integer :: i, j, k, m

        call flux_quadrature(res, 1, s, ws, th, ze, dth, dze)
        allocate (vol%x(3, res%ns*res%ntheta*res%nzeta), &
                  vol%w(res%ns*res%ntheta*res%nzeta))
        allocate (surf%x(3, res%ntheta*res%nzeta), surf%area(3, res%ntheta*res%nzeta))
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

    subroutine ball_field(pts)
        type(source_points_t), intent(inout) :: pts
        integer :: m

        allocate (pts%b(3, size(pts%x, 2)))
        do m = 1, size(pts%x, 2)
            pts%b(:, m) = ball_cyl(pts%x(:, m))
        end do
    end subroutine ball_field

    function ball_cyl(x) result(b)
        real(dp), intent(in) :: x(3)
        real(dp) :: b(3), d(3), ap(3), bp(3)

        d = cart(x) - cart(c_ball)
        if (norm2(d) < a_ball) then
            call mms_ball_inner(d(1), d(2), d(3), m_dip, a_ball, ap(1), ap(2), ap(3), &
                                bp(1), bp(2), bp(3))
        else
            call mms_ball_outer(d(1), d(2), d(3), m_dip, a_ball, ap(1), ap(2), ap(3), &
                                bp(1), bp(2), bp(3))
        end if
        b = cyl_vec(x, bp)
    end function ball_cyl

    function dipole_cyl(x) result(b)
        real(dp), intent(in) :: x(3)
        real(dp) :: b(3), d(3), ap(3), bp(3)

        d = cart(x) - cart(c_ball)
        call mms_ball_outer(d(1), d(2), d(3), m_dip, a_ball, ap(1), ap(2), ap(3), &
                            bp(1), bp(2), bp(3))
        b = cyl_vec(x, bp)
    end function dipole_cyl

    function cart(x) result(y)
        real(dp), intent(in) :: x(3)
        real(dp) :: y(3)

        y = [x(1)*cos(x(2)), x(1)*sin(x(2)), x(3)]
    end function cart

    function cyl_vec(x, v) result(w)
        real(dp), intent(in) :: x(3), v(3)
        real(dp) :: w(3)

        call cart_to_cyl(cos(x(2)), sin(x(2)), v(1), v(2), w(1), w(2))
        w(3) = v(3)
    end function cyl_vec

    function cart_vec(x, v) result(w)
        real(dp), intent(in) :: x(3), v(3)
        real(dp) :: w(3)

        call cyl_to_cart(cos(x(2)), sin(x(2)), v(1), v(2), w(1), w(2))
        w(3) = v(3)
    end function cart_vec

    subroutine circular_loop(coils, radius, z, current)
        type(coils_t), intent(out) :: coils
        real(dp), intent(in) :: radius, z, current
        integer, parameter :: n = 720
        real(dp) :: x(n + 1), y(n + 1), zz(n + 1), c(n + 1)
        integer :: i

        do i = 1, n + 1
            x(i) = radius*cos(2.0_dp*pi*(i - 1)/n)
            y(i) = radius*sin(2.0_dp*pi*(i - 1)/n)
        end do
        zz = z
        c = current
        c(n + 1) = 0.0_dp
        call coils_init(x, y, zz, c, coils)
    end subroutine circular_loop

    subroutine check(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what

        if (.not. ok) then
            print *, "FAIL: ", what
            error stop 1
        end if
    end subroutine check

    include '../../src/generated/frozen_current_mms_leaves.inc'

end program test_frozen_current_field
