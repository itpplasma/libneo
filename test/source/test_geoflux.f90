program test_geoflux
    use, intrinsic :: iso_fortran_env, only : dp => real64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use geoflux_coordinates, only : geoflux_to_cyl, cyl_to_geoflux
    use geoflux_field, only : spline_geoflux_data, splint_geoflux_field
    use field_sub, only : field_eq, psif
    use field_eq_mod, only : psi_axis

    implicit none

    character(len=*), parameter :: geqdsk_file = '../../python/tests/test.geqdsk'
    real(dp), parameter :: tol_roundtrip = 2.0d-2

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

    call spline_geoflux_data(geqdsk_file, 64, 128)

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

    Bmod_expected = sqrt(Br*Br + Bphi*Bphi + Bz*Bz)
    if (Bmod_expected <= 0.0_dp) then
        write(*,*) 'Reference Bmod not positive: ', Bmod_expected
        error stop
    end if

    hcov_expected(1) = (Br*jac(1,1) + Bz*jac(3,1)) / Bmod_expected
    hcov_expected(2) = (Br*jac(1,2) + Bz*jac(3,2)) / Bmod_expected
    hcov_expected(3) = (Bphi * x_cyl(1)) / Bmod_expected

    psi_expected = psif - psi_axis

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

end program test_geoflux
