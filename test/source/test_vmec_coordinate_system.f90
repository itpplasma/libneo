program test_vmec_coordinate_system
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_coordinates, only: vmec_coordinate_system_t
    use math_constants, only: TWOPI
    use new_vmec_stuff_mod, only: netcdffile, nper
    use spline_vmec_sub, only: spline_vmec_data
    implicit none

    type(vmec_coordinate_system_t) :: vmec
    real(dp) :: u(3), xcyl(3)
    real(dp) :: e_cov(3, 3)
    real(dp) :: g(3, 3), ginv(3, 3), sqrtg
    real(dp) :: zeta_period
    integer :: nerrors

    nerrors = 0

    netcdffile = "wout.nc"
    call spline_vmec_data
    if (nper < 1) then
        print *, "  FAIL: VMEC nper invalid"
        nerrors = nerrors + 1
    end if

    zeta_period = TWOPI/real(max(1, nper), dp)
    u = [0.25_dp, 0.13_dp*TWOPI, 0.2_dp*zeta_period]
    call vmec%evaluate_point(u, xcyl)

    if (.not. (xcyl(1) > 100.0_dp .and. xcyl(1) < 2.0e5_dp)) then
        print *, "  FAIL: VMEC R not in cm range, R=", xcyl(1)
        nerrors = nerrors + 1
    end if
    if (abs(xcyl(2) - u(3)) > 1.0e-14_dp) then
        print *, "  FAIL: VMEC phi mismatch"
        nerrors = nerrors + 1
    end if
    if (.not. (abs(xcyl(3)) < 2.0e5_dp)) then
        print *, "  FAIL: VMEC Z not finite/plausible, Z=", xcyl(3)
        nerrors = nerrors + 1
    end if

    call vmec%covariant_basis(u, e_cov)
    if (any(.not. (abs(e_cov) < huge(1.0_dp)))) then
        print *, "  FAIL: VMEC covariant basis contains invalid values"
        nerrors = nerrors + 1
    end if

    call vmec%metric_tensor(u, g, ginv, sqrtg)
    if (.not. (sqrtg > 0.0_dp .and. sqrtg < huge(1.0_dp))) then
        print *, "  FAIL: VMEC sqrtg invalid, sqrtg=", sqrtg
        nerrors = nerrors + 1
    end if
    if (maxval(abs(g - transpose(g))) > 1.0e-10_dp) then
        print *, "  FAIL: VMEC metric tensor not symmetric"
        nerrors = nerrors + 1
    end if

    if (nerrors > 0) then
        print *, "FAILED: ", nerrors, " error(s) in VMEC coordinate system test"
        error stop 1
    end if

    print *, "  PASS: VMEC coordinate system basic checks"
end program test_vmec_coordinate_system

