program test_geoflux
    use, intrinsic :: iso_fortran_env, only : dp => real64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use geoflux_coordinates, only : geoflux_to_cyl, cyl_to_geoflux, &
                                    geoflux_get_axis, geoflux_get_flux_profiles
    use geoflux_field, only : spline_geoflux_data, splint_geoflux_field
    use field_sub, only : field_eq, psif
    use field_eq_mod, only : psi_axis, psi_sep
    use math_constants, only : pi

    implicit none

    character(len=*), parameter :: fallback_geqdsk = '../../python/tests/test.geqdsk'
    character(len=*), parameter :: env_geqdsk = 'LIBNEO_TEST_GEQDSK'
    real(dp), parameter :: tol_roundtrip = 2.0d-2

    character(len=512) :: geqdsk_file
    character(len=512) :: arg_buffer
    integer :: arg_status
    integer :: arg_len
    real(dp) :: x_geo(3), x_geo_back(3)
    real(dp) :: x_cyl(3)
    real(dp) :: max_diff
    real(dp) :: Acov(3), hcov(3), Bmod
    real(dp) :: jac(3,3)
    real(dp) :: Br, Bphi, Bz
    real(dp) :: dBrdR, dBrdp, dBrdZ
    real(dp) :: dBpdR, dBpdp, dBpdZ
    real(dp) :: dBzdR, dBzdp, dBzdZ
    real(dp) :: Bmod_expected
    real(dp) :: hcov_expected(3)
    real(dp) :: psi_expected
    real(dp) :: tol_field
    real(dp) :: q_profile, dq_ds, psi_pol, dpsi_pol_ds, psi_tor_edge
    real(dp) :: R_axis, Z_axis, radius_inner, radius_outer, radius_ratio
    real(dp) :: jac_left(3,3), jac_right(3,3), derivative_jump
    real(dp) :: volume_inner, volume_outer, volume_ratio
    integer, parameter :: ns_cache = 64, ntheta_cache = 128
    real(dp), parameter :: cache_rho_step = 1.0_dp/real(ns_cache - 1, dp)
    real(dp), parameter :: s_inner = (0.025_dp*cache_rho_step)**2
    real(dp), parameter :: s_outer = (0.05_dp*cache_rho_step)**2
    real(dp), parameter :: tol_axis_scaling = 1.0e-3_dp

    geqdsk_file = fallback_geqdsk
    arg_buffer = ''
    call get_command_argument(1, value=arg_buffer, length=arg_len, status=arg_status)
    if (arg_status == 0 .and. arg_len > 0 .and. len_trim(arg_buffer) > 0) then
        geqdsk_file = trim(arg_buffer)
    else
        call get_environment_variable(env_geqdsk, value=arg_buffer, &
                                      length=arg_len, status=arg_status)
        if (arg_status == 0 .and. arg_len > 0 .and. len_trim(arg_buffer) > 0) then
            geqdsk_file = trim(arg_buffer)
        end if
    end if

    call spline_geoflux_data(trim(geqdsk_file), ns_cache, ntheta_cache)

    call geoflux_get_axis(R_axis, Z_axis)
    x_geo = [s_inner, 0.4_dp, 0.0_dp]
    call geoflux_to_cyl(x_geo, x_cyl)
    radius_inner = hypot(x_cyl(1) - R_axis, x_cyl(3) - Z_axis)
    x_geo(1) = s_outer
    call geoflux_to_cyl(x_geo, x_cyl)
    radius_outer = hypot(x_cyl(1) - R_axis, x_cyl(3) - Z_axis)
    radius_ratio = radius_outer/radius_inner
    if (abs(radius_ratio - sqrt(s_outer/s_inner)) > tol_axis_scaling) then
        write(*,*) 'Near-axis radius must scale as sqrt(s): ', radius_ratio
        error stop
    end if

    ! The orbit integrator crosses cache-cell boundaries.  A piecewise-linear
    ! surface has a discontinuous metric there even when its values converge;
    ! the spline representation must provide a continuous first derivative.
    x_geo = [(20.0_dp*cache_rho_step - 1.0e-8_dp)**2, 0.37_dp, 0.0_dp]
    call geoflux_to_cyl(x_geo, x_cyl, jac_left)
    x_geo(1) = (20.0_dp*cache_rho_step + 1.0e-8_dp)**2
    call geoflux_to_cyl(x_geo, x_cyl, jac_right)
    derivative_jump = maxval(abs(jac_right(:,1:2) - jac_left(:,1:2))) / &
        max(1.0_dp, maxval(abs(jac_left(:,1:2))))
    if (derivative_jump > 1.0e-5_dp) then
        write(*,*) 'Cached surface metric is discontinuous: ', derivative_jump
        error stop
    end if

    volume_inner = toroidal_volume(s_inner)
    volume_outer = toroidal_volume(s_outer)
    volume_ratio = volume_outer/volume_inner
    if (abs(volume_ratio - s_outer/s_inner) > 1.0e-2_dp) then
        write(*,*) 'Near-axis enclosed volume must scale as s: ', volume_ratio
        error stop
    end if

    call geoflux_get_flux_profiles(0.3_dp, q_profile, dq_ds, psi_pol, &
                                   dpsi_pol_ds, psi_tor_edge)
    if (.not. all(ieee_is_finite([q_profile, dq_ds, psi_pol, &
                                  dpsi_pol_ds, psi_tor_edge]))) then
        write(*,*) 'Non-finite geoflux profile metadata.'
        error stop
    end if
    if (abs(q_profile) <= tiny(q_profile) .or. &
        abs(dpsi_pol_ds) <= tiny(dpsi_pol_ds) .or. &
        abs(psi_tor_edge) <= tiny(psi_tor_edge)) then
        write(*,*) 'Degenerate geoflux profile metadata.'
        error stop
    end if

    x_geo = [0.3_dp, 0.5_dp, 0.0_dp]
    call geoflux_to_cyl(x_geo, x_cyl)
    call cyl_to_geoflux(x_cyl, x_geo_back)
    max_diff = maxval(abs(x_geo - x_geo_back))
    if (max_diff > tol_roundtrip) then
        write(*,*) 'Round-trip mismatch for mid-surface point:', max_diff
        error stop
    end if

    x_geo = [1.0d-4, -1.2_dp, 1.7_dp]
    call geoflux_to_cyl(x_geo, x_cyl)
    call cyl_to_geoflux(x_cyl, x_geo_back)
    max_diff = maxval(abs(x_geo - x_geo_back))
    if (max_diff > tol_roundtrip) then
        write(*,*) 'Round-trip mismatch near axis:', max_diff
        error stop
    end if

    x_geo = [0.8_dp, 0.9_dp, 2.3_dp]
    call geoflux_to_cyl(x_geo, x_cyl, jac)
    if (x_cyl(1) <= 0.0_dp) then
        write(*,*) 'Non-positive R coordinate on outer surface.'
        error stop
    end if
    if (abs(x_cyl(2) - x_geo(3)) > 1.0d-12) then
        write(*,*) 'Phi was not preserved by geoflux_to_cyl.'
        error stop
    end if

    call field_eq(x_cyl(1), x_cyl(2), x_cyl(3), Br, Bphi, Bz, &
                  dBrdR, dBrdp, dBrdZ, dBpdR, dBpdp, dBpdZ, &
                  dBzdR, dBzdp, dBzdZ)

    if (psi_axis /= 0.0_dp) then
        write(*,*) 'Shifted field_eq flux gauge must have psi_axis=0: ', psi_axis
        error stop
    end if
    if (.not. ieee_is_finite(psi_sep) .or. psi_sep == 0.0_dp) then
        write(*,*) 'Shifted field_eq psi_sep must be finite and nonzero: ', psi_sep
        error stop
    end if

    Bmod_expected = sqrt(Br*Br + Bphi*Bphi + Bz*Bz)
    if (Bmod_expected <= 0.0_dp) then
        write(*,*) 'Reference Bmod not positive: ', Bmod_expected
        error stop
    end if

    hcov_expected(1) = (Br*jac(1,1) + Bz*jac(3,1)) / Bmod_expected
    hcov_expected(2) = (Br*jac(1,2) + Bz*jac(3,2)) / Bmod_expected
    hcov_expected(3) = (Bphi * x_cyl(1)) / Bmod_expected

    psi_expected = psif

    tol_field = 5.0d-10

    call splint_geoflux_field(x_geo(1), x_geo(2), x_geo(3), Acov, hcov, Bmod)

    if (.not. ieee_is_finite(Bmod) .or. Bmod <= 0.0_dp) then
        write(*,*) 'Bmod not positive/finite: ', Bmod
        error stop
    end if

    if (abs(Bmod - Bmod_expected) > tol_field) then
        write(*,*) 'Bmod mismatch. expected=', Bmod_expected, ' got=', Bmod
        error stop
    end if

    if (maxval(abs(hcov - hcov_expected)) > tol_field) then
        write(*,*) 'Covariant unit vector mismatch.'
        write(*,*) 'expected=', hcov_expected
        write(*,*) 'got=', hcov
        error stop
    end if

    if (.not. ieee_is_finite(Acov(3))) then
        write(*,*) 'A_phi is not finite: ', Acov(3)
        error stop
    end if

    if (abs(Acov(3) - psi_expected) > tol_field) then
        write(*,*) 'A_phi does not match psi difference.'
        write(*,*) 'expected=', psi_expected, ' got=', Acov(3)
        error stop
    end if

contains

    function toroidal_volume(s_val) result(volume)
        real(dp), intent(in) :: s_val
        real(dp) :: volume
        integer, parameter :: ntheta = 256
        real(dp) :: points(2, ntheta), x_geo_local(3), x_cyl_local(3)
        real(dp) :: cross, moment
        integer :: i, next_i

        do i = 1, ntheta
            x_geo_local = [s_val, 2.0_dp*pi*real(i - 1, dp)/real(ntheta, dp), &
                           0.0_dp]
            call geoflux_to_cyl(x_geo_local, x_cyl_local)
            points(:, i) = [x_cyl_local(1), x_cyl_local(3)]
        end do

        moment = 0.0_dp
        do i = 1, ntheta
            next_i = modulo(i, ntheta) + 1
            cross = points(1, i)*points(2, next_i) &
                    - points(1, next_i)*points(2, i)
            moment = moment + (points(1, i) + points(1, next_i))*cross
        end do
        volume = pi*abs(moment)/3.0_dp
    end function toroidal_volume

end program test_geoflux
