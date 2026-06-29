!> Pin the invertibility of the VMEC <-> Boozer angle transforms.
!>
!> Builds the Boozer field from the LandremanPaul2021 QA reference wout via
!> get_boozer_coordinates, then for several flux surfaces and angle pairs
!> checks that vmec_to_boozer followed by boozer_to_vmec recovers the VMEC
!> angles, and that the reverse composition recovers the Boozer angles, both
!> to high precision modulo the angular periods.
program test_boozer_angle_roundtrip
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use boozer_coordinates_mod, only: use_B_r
    use new_vmec_stuff_mod, only: nper
    use boozer_sub, only: get_boozer_coordinates, vmec_to_boozer, boozer_to_vmec
    implicit none

    real(dp), parameter :: PI = 3.14159265358979_dp
    real(dp), parameter :: TWOPI = 2.0_dp*PI
    real(dp), parameter :: tol = 1.0e-10_dp
    character(len=*), parameter :: wout_file = "wout.nc"

    integer, parameter :: n_s = 3, n_ang = 4
    real(dp) :: stor(n_s), vt(n_ang), vp(n_ang)
    real(dp) :: phi_period
    integer :: is, ia
    logical :: test_failed

    use_B_r = .true.
    call get_boozer_coordinates(wout_file, radial_spline_order=5, &
                                angular_spline_order=5, grid_refinment=3)

    phi_period = TWOPI/real(nper, dp)

    stor = [0.2_dp, 0.5_dp, 0.8_dp]
    vt = [0.3_dp, 1.5_dp, 3.0_dp, 5.0_dp]
    vp = [0.1_dp, 0.7_dp, 1.4_dp, 2.5_dp]

    test_failed = .false.
    do is = 1, n_s
        do ia = 1, n_ang
            call check_vmec_roundtrip(stor(is), vt(ia), vp(ia), test_failed)
            call check_boozer_roundtrip(stor(is), vt(ia), vp(ia), test_failed)
        end do
    end do

    if (test_failed) error stop "boozer angle transform is not invertible"
    print *, "vmec<->boozer angle transform round-trips at all points"

contains

    !> Forward then inverse: VMEC -> Boozer -> VMEC must recover the input.
    subroutine check_vmec_roundtrip(s, theta, varphi, test_failed)
        real(dp), intent(in) :: s, theta, varphi
        logical, intent(inout) :: test_failed

        real(dp) :: vartheta_B, varphi_B, theta_rec, varphi_rec

        call vmec_to_boozer(s, theta, varphi, vartheta_B, varphi_B)
        call boozer_to_vmec(s, vartheta_B, varphi_B, theta_rec, varphi_rec)

        call check_angle("vmec theta", s, theta, theta_rec, TWOPI, test_failed)
        call check_angle("vmec varphi", s, varphi, varphi_rec, phi_period, test_failed)
    end subroutine check_vmec_roundtrip

    !> Inverse then forward: Boozer -> VMEC -> Boozer must recover the input.
    subroutine check_boozer_roundtrip(s, vartheta_B, varphi_B, test_failed)
        real(dp), intent(in) :: s, vartheta_B, varphi_B
        logical, intent(inout) :: test_failed

        real(dp) :: theta, varphi, vt_rec, vp_rec

        call boozer_to_vmec(s, vartheta_B, varphi_B, theta, varphi)
        call vmec_to_boozer(s, theta, varphi, vt_rec, vp_rec)

        call check_angle("boozer vartheta", s, vartheta_B, vt_rec, TWOPI, test_failed)
        call check_angle("boozer varphi", s, varphi_B, vp_rec, phi_period, test_failed)
    end subroutine check_boozer_roundtrip

    !> Compare two angles modulo a period, taking the shortest wrapped distance.
    subroutine check_angle(name, s, expected, got, period, test_failed)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: s, expected, got, period
        logical, intent(inout) :: test_failed

        real(dp) :: diff

        diff = modulo(got - expected, period)
        if (diff > 0.5_dp*period) diff = diff - period

        if (abs(diff) > tol) then
            print *, "---------------------------------------------------"
            print *, "round-trip mismatch in ", name, " at s = ", s
            print *, "expected: ", expected
            print *, "got     : ", got
            print *, "wrapped err: ", abs(diff)
            test_failed = .true.
        end if
    end subroutine check_angle

end program test_boozer_angle_roundtrip
