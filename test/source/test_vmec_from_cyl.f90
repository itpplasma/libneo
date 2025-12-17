program test_vmec_from_cyl
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_coordinates, only: vmec_coordinate_system_t
    use math_constants, only: TWOPI
    use new_vmec_stuff_mod, only: netcdffile, nper
    use spline_vmec_sub, only: spline_vmec_data
    implicit none

    type(vmec_coordinate_system_t) :: vmec
    integer :: ierr, nerrors
    integer :: i
    real(dp) :: u(3), u_back(3), xcyl(3), xcyl_back(3)
    real(dp) :: zeta_period
    real(dp) :: max_dr, max_dz

    nerrors = 0
    max_dr = 0.0_dp
    max_dz = 0.0_dp

    netcdffile = "wout.nc"
    call spline_vmec_data

    zeta_period = TWOPI/real(max(1, nper), dp)

    do i = 1, 20
        u(1) = real(i, dp)/real(21, dp)
        u(2) = modulo(real(37*i, dp)*0.017_dp*TWOPI, TWOPI)
        u(3) = modulo(real(29*i, dp)*0.013_dp*zeta_period, zeta_period)

        call vmec%evaluate_cyl(u, xcyl)
        call vmec%from_cyl(xcyl, u_back, ierr)
        if (ierr /= 0) then
            print *, "  FAIL: vmec_from_cyl ierr=", ierr
            print *, "    u     =", u
            print *, "    xcyl  =", xcyl
            nerrors = nerrors + 1
            cycle
        end if

        call vmec%evaluate_cyl(u_back, xcyl_back)
        max_dr = max(max_dr, abs(xcyl_back(1) - xcyl(1)))
        max_dz = max(max_dz, abs(xcyl_back(3) - xcyl(3)))
    end do

    if (max_dr > 1.0e-6_dp .or. max_dz > 1.0e-6_dp) then
        print *, "  FAIL: vmec from_cyl roundtrip mismatch"
        print *, "    max|dR|=", max_dr
        print *, "    max|dZ|=", max_dz
        nerrors = nerrors + 1
    else
        print *, "  PASS: vmec from_cyl roundtrip max|dR|=", max_dr, " max|dZ|=", max_dz
    end if

    if (nerrors > 0) then
        print *, "FAILED: ", nerrors, " error(s) in VMEC from_cyl test"
        error stop 1
    end if
end program test_vmec_from_cyl

