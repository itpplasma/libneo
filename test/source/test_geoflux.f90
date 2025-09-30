program test_geoflux
    use, intrinsic :: iso_fortran_env, only : dp => real64
    use geoflux_coordinates, only : geoflux_to_cyl, cyl_to_geoflux
    use geoflux_field, only : spline_geoflux_data, splint_geoflux_field
    use field_sub, only : psif
    use field_eq_mod, only : psi_axis

    implicit none

    character(len=*), parameter :: geqdsk_file = '../../python/tests/test.geqdsk'
    real(dp), parameter :: tol_roundtrip = 2.0d-2

    real(dp) :: x_geo(3), x_geo_back(3)
    real(dp) :: x_cyl(3)
    real(dp) :: max_diff
    real(dp) :: Acov(3), hcov(3), Bmod
    real(dp) :: unit_length

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
    call geoflux_to_cyl(x_geo, x_cyl)
    if (x_cyl(1) <= 0.0_dp) then
        write(*,*) 'Non-positive R coordinate on outer surface.'
        error stop
    end if
    if (abs(x_cyl(2) - x_geo(3)) > 1.0d-12) then
        write(*,*) 'Phi was not preserved by geoflux_to_cyl.'
        error stop
    end if

    call splint_geoflux_field(x_geo(1), x_geo(2), x_geo(3), Acov, hcov, Bmod)

    if (Bmod <= 0.0_dp) then
        write(*,*) 'Bmod should be positive but is ', Bmod
        error stop
    end if

    unit_length = maxval(abs(hcov))
    if (.not. (unit_length < 1.0d6 .and. unit_length > 0.0d0)) then
        write(*,*) 'Unexpected hcov magnitude: ', unit_length
        error stop
    end if

    if (abs(Acov(3) - (psif - psi_axis)) > 1.0d-8) then
        write(*,*) 'Unexpected A_phi covariant component.'
        error stop
    end if

end program test_geoflux
