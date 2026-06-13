!> Pin the field-object entry get_boozer_coordinates_with_field against the
!> global-VMEC path. Splines the LandremanPaul2021 QA reference wout via
!> get_boozer_coordinates and records the full splint_boozer_coord output at
!> interior points. Then builds a vmec_field_t through create_vmec_field, runs
!> get_boozer_coordinates_with_field on it, and re-evaluates the same points.
!> A VMEC field object must reproduce the global-VMEC path to ~1e-10.
program test_boozer_with_field
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use boozer_coordinates_mod, only: use_B_r
    use new_vmec_stuff_mod, only: nper
    use field_vmec, only: vmec_field_t, create_vmec_field
    use boozer_sub, only: get_boozer_coordinates, &
                          get_boozer_coordinates_with_field, splint_boozer_coord
    implicit none

    real(dp), parameter :: PI = 3.14159265358979_dp
    real(dp), parameter :: TWOPI = 2.0_dp*PI
    real(dp), parameter :: tol = 1.0e-10_dp
    character(len=*), parameter :: wout_file = "wout.nc"

    integer, parameter :: n_test = 20
    real(dp) :: s_test(n_test), th_test(n_test), ph_test(n_test)
    real(dp) :: bmod_ref(n_test), bth_ref(n_test), bph_ref(n_test)
    real(dp) :: aphi_ref(n_test)
    real(dp) :: bmod_new(n_test), bth_new(n_test), bph_new(n_test)
    real(dp) :: aphi_new(n_test)

    type(vmec_field_t) :: field
    real(dp) :: phi_period, max_err
    integer :: i

    use_B_r = .false.
    call get_boozer_coordinates(wout_file, radial_spline_order=5, &
                                angular_spline_order=5, grid_refinment=3)

    phi_period = TWOPI/real(nper, dp)
    call init_test_points(phi_period)

    do i = 1, n_test
        call eval_point(s_test(i), th_test(i), ph_test(i), &
                        bmod_ref(i), bth_ref(i), bph_ref(i), aphi_ref(i))
    end do

    ! Rebuild through the field-object entry. The VMEC splines are still loaded
    ! from get_boozer_coordinates above, which the vmec_field path reuses.
    call create_vmec_field(field)
    call get_boozer_coordinates_with_field(field)

    do i = 1, n_test
        call eval_point(s_test(i), th_test(i), ph_test(i), &
                        bmod_new(i), bth_new(i), bph_new(i), aphi_new(i))
    end do

    max_err = 0.0_dp
    do i = 1, n_test
        max_err = max(max_err, rel_diff(bmod_new(i), bmod_ref(i)))
        max_err = max(max_err, rel_diff(bth_new(i), bth_ref(i)))
        max_err = max(max_err, rel_diff(bph_new(i), bph_ref(i)))
        max_err = max(max_err, rel_diff(aphi_new(i), aphi_ref(i)))
    end do

    print *, "max relative error field-object vs global VMEC path:", max_err
    if (max_err > tol) then
        error stop "get_boozer_coordinates_with_field deviates from global path"
    end if
    print *, "get_boozer_coordinates_with_field reproduces global VMEC path"

contains

    subroutine eval_point(s, th, ph, bmod, bth, bph, aphi)
        real(dp), intent(in) :: s, th, ph
        real(dp), intent(out) :: bmod, bth, bph, aphi

        real(dp) :: A_theta, dA_theta_dr, dA_phi_dr
        real(dp) :: d2A_phi_dr2, d3A_phi_dr3
        real(dp) :: dBth, d2Bth, dBph, d2Bph
        real(dp) :: dBmod(3), d2Bmod(6), sqrt_g_ss, Br, dBr(3), d2Br(6)

        call splint_boozer_coord(s, th, ph, 0, &
                                 A_theta, aphi, dA_theta_dr, dA_phi_dr, &
                                 d2A_phi_dr2, d3A_phi_dr3, bth, dBth, d2Bth, &
                                 bph, dBph, d2Bph, bmod, dBmod, d2Bmod, &
                                 Br, dBr, d2Br, sqrt_g_ss)
    end subroutine eval_point

    pure function rel_diff(a, b) result(d)
        real(dp), intent(in) :: a, b
        real(dp) :: d
        if (abs(b) > 0.0_dp) then
            d = abs(a - b)/abs(b)
        else
            d = abs(a)
        end if
    end function rel_diff

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

end program test_boozer_with_field
