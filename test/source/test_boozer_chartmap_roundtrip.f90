!> Pin field equivalence across the Boozer chartmap export/import roundtrip.
!>
!> Builds the Boozer field from the LandremanPaul2021 QA reference wout via
!> get_boozer_coordinates, records |B| at interior test points through
!> splint_boozer_coord, exports a chartmap NetCDF file, reimports it with
!> load_boozer_from_chartmap, and re-evaluates |B| at the same points. The
!> maximum relative |B| error over the interior points must stay below 8e-5.
program test_boozer_chartmap_roundtrip
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use boozer_coordinates_mod, only: use_B_r
    use new_vmec_stuff_mod, only: nper
    use boozer_sub, only: get_boozer_coordinates, splint_boozer_coord
    use boozer_chartmap, only: export_boozer_chartmap, load_boozer_from_chartmap
    implicit none

    real(dp), parameter :: PI = 3.14159265358979_dp
    real(dp), parameter :: TWOPI = 2.0_dp*PI
    real(dp), parameter :: field_tol = 8.0e-5_dp
    character(len=*), parameter :: wout_file = "wout.nc"
    character(len=*), parameter :: chartmap_file = "roundtrip_test.nc"

    integer, parameter :: n_test = 50
    real(dp) :: s_test(n_test), th_test(n_test), ph_test(n_test)
    real(dp) :: bmod_ref(n_test), bmod_new(n_test)

    real(dp) :: A_theta, A_phi_val, dA_theta_dr, dA_phi_dr
    real(dp) :: d2A_phi_dr2, d3A_phi_dr3
    real(dp) :: Bth, dBth, d2Bth, Bph, dBph, d2Bph
    real(dp) :: Bmod, dBmod(3), d2Bmod(6), Br, dBr(3), d2Br(6)
    real(dp) :: phi_period, rel_err, max_err_bmod
    integer :: i

    use_B_r = .false.
    call get_boozer_coordinates(wout_file, radial_spline_order=5, &
                                angular_spline_order=5, grid_refinment=3)

    phi_period = TWOPI/real(nper, dp)
    call init_test_points(phi_period)

    do i = 1, n_test
        call splint_boozer_coord(s_test(i), th_test(i), ph_test(i), 0, &
                                 A_theta, A_phi_val, dA_theta_dr, dA_phi_dr, &
                                 d2A_phi_dr2, d3A_phi_dr3, Bth, dBth, d2Bth, &
                                 Bph, dBph, d2Bph, Bmod, dBmod, d2Bmod, &
                                 Br, dBr, d2Br)
        bmod_ref(i) = Bmod
    end do

    call export_boozer_chartmap(chartmap_file)
    call load_boozer_from_chartmap(chartmap_file)

    do i = 1, n_test
        call splint_boozer_coord(s_test(i), th_test(i), ph_test(i), 0, &
                                 A_theta, A_phi_val, dA_theta_dr, dA_phi_dr, &
                                 d2A_phi_dr2, d3A_phi_dr3, Bth, dBth, d2Bth, &
                                 Bph, dBph, d2Bph, Bmod, dBmod, d2Bmod, &
                                 Br, dBr, d2Br)
        bmod_new(i) = Bmod
    end do

    max_err_bmod = 0.0_dp
    do i = 1, n_test
        if (abs(bmod_ref(i)) > 0.0_dp) then
            rel_err = abs(bmod_new(i) - bmod_ref(i))/abs(bmod_ref(i))
        else
            rel_err = abs(bmod_new(i))
        end if
        max_err_bmod = max(max_err_bmod, rel_err)
    end do

    print *, "max relative |B| error over chartmap roundtrip:", max_err_bmod
    if (max_err_bmod > field_tol) then
        error stop "Bmod chartmap roundtrip error too large"
    end if
    print *, "Boozer chartmap export/import preserves |B| within tolerance"

contains

    !> Interior test points: s in (0.1, 0.8), full theta period, scattered phi.
    subroutine init_test_points(phi_per)
        real(dp), intent(in) :: phi_per
        integer :: j
        real(dp) :: frac

        do j = 1, n_test
            frac = real(j, dp)/real(n_test + 1, dp)
            s_test(j) = 0.1_dp + 0.7_dp*frac
            th_test(j) = TWOPI*frac
            ph_test(j) = phi_per*mod(real(2*j + 1, dp), real(n_test + 1, dp)) &
                         /real(n_test + 1, dp)
        end do
    end subroutine init_test_points

end program test_boozer_chartmap_roundtrip
