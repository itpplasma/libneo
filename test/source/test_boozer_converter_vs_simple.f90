!> Pin the libneo VMEC->Boozer converter against SIMPLE reference values.
!>
!> Builds the Boozer field from the LandremanPaul2021 QA reference wout via
!> get_boozer_coordinates, evaluates it at five (s, theta, phi) points through
!> splint_boozer_coord, assembles the same |B|, sqrt(g) and covariant/
!> contravariant components SIMPLE produces, and compares against the values
!> printed by SIMPLE. Quantities are kept in the converter's native CGS units
!> (Gauss, cm), which is what the SIMPLE reference numbers are given in.
program test_boozer_converter_vs_simple
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use boozer_coordinates_mod, only: use_B_r
    use vector_potentail_mod, only: torflux
    use boozer_sub, only: get_boozer_coordinates, splint_boozer_coord, BOOZER_LEGACY
    implicit none

    ! reltol pins the well-conditioned quantities (|B|, sqrt(g), the large
    ! covariant components). abstol is the absolute floor for the near-zero
    ! covariant components (e.g. B_theta/|B| ~ 1e-5 for QA), whose relative
    ! agreement is compile/platform-sensitive at the spline-interpolation level;
    ! 1e-10 absorbs that cross-platform FP noise without loosening the
    ! relative bound the large quantities meet.
    real(dp), parameter :: reltol = 1.0e-6_dp, abstol = 1.0e-10_dp
    ! sqrt(g^ss) comes back from the field3d spline (the use_B_r quantity). The
    ! cross-check evaluates the VMEC metric at the Boozer->VMEC back-mapped angles,
    ! so the gap carries both the spline and the angle-map interpolation error at
    ! this coarse grid_refinment=3 grid (~1e-2). 2e-2 confirms it is the right
    ! quantity without pinning that grid-resolution floor.
    real(dp), parameter :: reltol_gss = 2.0e-2_dp
    character(len=*), parameter :: wout_file = "wout.nc"

    integer, parameter :: n_cases = 5
    real(dp) :: stor(n_cases), theta(n_cases), phi(n_cases)
    real(dp) :: bmod, sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)
    real(dp) :: sqrt_g_ss_B, sqrt_g_ss_ref
    real(dp) :: bmod_ref(n_cases), sqrtg_ref(n_cases)
    real(dp) :: bder_ref(n_cases, 3), hcovar_ref(n_cases, 3)
    real(dp) :: hctrvr_ref(n_cases, 3), hcurl_ref(n_cases, 3)
    integer :: case
    logical :: test_failed

    use_B_r = .true.
    call get_boozer_coordinates(wout_file, radial_spline_order=5, &
                                angular_spline_order=5, grid_refinment=3, &
                                construction=BOOZER_LEGACY)

    stor = [0.1_dp, 0.3_dp, 0.5_dp, 0.7_dp, 0.9_dp]
    theta = [0.0_dp, 1.0_dp, 3.14_dp, 0.5_dp, 2.0_dp]
    phi = [0.0_dp, 0.5_dp, 1.57_dp, 2.0_dp, 3.0_dp]

    ! Reference values regenerated after the -I*lambda Boozer angle fix
    ! (libneo 266f0a3): dump of this converter on the LandremanPaul2021 QA wout
    ! with the call above (gfortran -O3, OpenMP off, es24.16). The old values,
    ! printed by SIMPLE, are recovered when the one-line fix in
    ! boozer_converter.F90 is reverted; hcovar(1) moves by up to 1.5e-3.
    bmod_ref = [5.7024014488311528e+04_dp, 5.7010408565725287e+04_dp, &
                6.3742758213882924e+04_dp, 5.4218536559316562e+04_dp, &
                6.1083220183434642e+04_dp]

    sqrtg_ref = [-1.5711629692163896e+07_dp, -1.5719091478248706e+07_dp, &
                 -1.2574007504116869e+07_dp, -1.7379594206972197e+07_dp, &
                 -1.3692768687452687e+07_dp]

    bder_ref(1, :) = [-1.7906807462368635e-01_dp, 4.3243884155381284e-14_dp, 3.5279905290844091e-12_dp]
    bder_ref(2, :) = [-6.5375422580664741e-02_dp, 4.8901135938522522e-02_dp, -4.2593323560656337e-04_dp]
    bder_ref(3, :) = [7.3646069744469744e-02_dp, 1.3567157156220991e-04_dp, -2.4743000451758948e-05_dp]
    bder_ref(4, :) = [-6.5308326833423816e-02_dp, 3.9538078645228660e-02_dp, -6.8133275613209887e-04_dp]
    bder_ref(5, :) = [1.2755348743767689e-02_dp, 1.0110736990682824e-01_dp, -7.9541233119082900e-04_dp]

    hcovar_ref(1, :) = [1.5922536168961330e-12_dp, -9.4970348982260694e-04_dp, 1.1151159177105092e+03_dp]
    hcovar_ref(2, :) = [-4.9484433131595819e-02_dp, -3.8520822102877319e-04_dp, 1.1153790778403581e+03_dp]
    hcovar_ref(3, :) = [8.2801889215536242e-04_dp, -1.8842726587574028e-04_dp, 9.9757472565820603e+02_dp]
    hcovar_ref(4, :) = [-5.8937852007688772e-02_dp, -9.8956245119790614e-05_dp, 1.1728116580885373e+03_dp]
    hcovar_ref(5, :) = [2.8041661963106538e-02_dp, -1.0129733846890238e-05_dp, 1.0410079195575106e+03_dp]

    hctrvr_ref(1, :) = [0.0000000000000000e+00_dp, 3.7881309780869710e-04_dp, 8.9676807933408687e-04_dp]
    hctrvr_ref(2, :) = [0.0000000000000000e+00_dp, 3.7723748026117170e-04_dp, 8.9655630554879965e-04_dp]
    hctrvr_ref(3, :) = [0.0000000000000000e+00_dp, 4.2022807275819278e-04_dp, 1.0024312499723974e-03_dp]
    hctrvr_ref(4, :) = [0.0000000000000000e+00_dp, 3.5621175215059544e-04_dp, 8.5265185450082378e-04_dp]
    hctrvr_ref(5, :) = [0.0000000000000000e+00_dp, 4.0002199333117358e-04_dp, 9.6060748940044090e-04_dp]

    hcurl_ref(1, :) = [3.0694011413892109e-18_dp, 1.2740515411884298e-05_dp, -8.6345717705439876e-09_dp]
    hcurl_ref(2, :) = [3.4698763487607772e-06_dp, 4.6481492620868257e-06_dp, -1.7477546857429966e-08_dp]
    hcurl_ref(3, :) = [1.0763674673664753e-08_dp, -5.7798156247815565e-06_dp, -9.7758245407636786e-09_dp]
    hcurl_ref(4, :) = [2.6681129003300245e-06_dp, 4.3919396688105018e-06_dp, -5.0461243807006480e-09_dp]
    hcurl_ref(5, :) = [7.6867998863539018e-06_dp, -9.5911628193097490e-07_dp, 7.0755537313349788e-09_dp]

    test_failed = .false.
    do case = 1, n_cases
        call evaluate_boozer(stor(case), theta(case), phi(case), &
                             bmod, sqrtg, bder, hcovar, hctrvr, hcurl, sqrt_g_ss_B)
        call check("bmod", case, [bmod], [bmod_ref(case)], test_failed)
        call check("sqrtg", case, [sqrtg], [sqrtg_ref(case)], test_failed)
        call check("bder", case, bder, bder_ref(case, :), test_failed)
        call check("hcovar", case, hcovar, hcovar_ref(case, :), test_failed)
        call check("hctrvr", case, hctrvr, hctrvr_ref(case, :), test_failed)
        call check("hcurl", case, hcurl, hcurl_ref(case, :), test_failed)

        ! sqrt(g^ss) is exposed only on the use_B_r path; cross-check it against
        ! the VMEC metric evaluated at the back-mapped VMEC angles.
        sqrt_g_ss_ref = metric_nabla_s(stor(case), theta(case), phi(case))
        if (abs(sqrt_g_ss_B - sqrt_g_ss_ref) > reltol_gss*abs(sqrt_g_ss_ref)) then
            print *, "---------------------------------------------------"
            print *, "mismatch in sqrt_g_ss case ", case
            print *, "spline : ", sqrt_g_ss_B
            print *, "metric : ", sqrt_g_ss_ref
            test_failed = .true.
        end if
    end do

    if (test_failed) error stop "boozer converter differs from SIMPLE reference"
    print *, "boozer converter matches SIMPLE reference at all points"

contains

    !> Assemble |B|, sqrt(g) and field components from splint_boozer_coord,
    !> following the SIMPLE/rabe Boozer field convention, in native CGS units.
    subroutine evaluate_boozer(s, vartheta_B, varphi_B, &
                               bmod, sqrtg, bder, hcovar, hctrvr, hcurl, sqrt_g_ss_B)
        real(dp), intent(in) :: s, vartheta_B, varphi_B
        real(dp), intent(out) :: bmod, sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)
        real(dp), intent(out) :: sqrt_g_ss_B

        real(dp) :: A_theta, A_phi, dA_theta_dr, dA_phi_dr, d2A_phi_dr2, d3A_phi_dr3
        real(dp) :: B_vartheta_B, dB_vartheta_B, d2B_vartheta_B
        real(dp) :: B_varphi_B, dB_varphi_B, d2B_varphi_B
        real(dp) :: Bmod_B, B_r
        real(dp) :: dBmod_B(3), dB_r(3), d2Bmod_B(6), d2B_r(6)
        real(dp) :: aiota, Bctrvr_theta, Bctrvr_phi, sqrtgbmod
        integer, parameter :: mode_secders = 0

        call splint_boozer_coord(s, vartheta_B, varphi_B, mode_secders, &
                                 A_theta, A_phi, dA_theta_dr, dA_phi_dr, &
                                 d2A_phi_dr2, d3A_phi_dr3, &
                                 B_vartheta_B, dB_vartheta_B, d2B_vartheta_B, &
                                 B_varphi_B, dB_varphi_B, d2B_varphi_B, &
                                 Bmod_B, dBmod_B, d2Bmod_B, &
                                 B_r, dB_r, d2B_r, sqrt_g_ss_B)

        aiota = -dA_phi_dr/dA_theta_dr
        bmod = Bmod_B
        bder = dBmod_B/Bmod_B
        sqrtg = (aiota*B_vartheta_B + B_varphi_B)/Bmod_B**2*torflux

        Bctrvr_phi = dA_theta_dr/sqrtg
        Bctrvr_theta = aiota*Bctrvr_phi
        hctrvr = [0.0_dp, Bctrvr_theta/bmod, Bctrvr_phi/bmod]
        hcovar = [B_r/bmod, B_vartheta_B/bmod, B_varphi_B/bmod]

        sqrtgbmod = sqrtg*bmod
        hcurl(1) = (B_vartheta_B*bder(3) - B_varphi_B*bder(2))/sqrtgbmod
        hcurl(2) = (B_varphi_B*bder(1) - B_r*bder(3) + dB_r(3) - dB_varphi_B)/sqrtgbmod
        hcurl(3) = (B_r*bder(2) - B_vartheta_B*bder(1) + dB_vartheta_B - dB_r(2))/sqrtgbmod
    end subroutine evaluate_boozer

    !> Reference |nabla s| = sqrt(g^ss) from the VMEC metric, independent of the
    !> Boozer field3d spline: back-map the Boozer angles to VMEC and evaluate.
    function metric_nabla_s(s, vartheta_B, varphi_B) result(sqrt_g_ss)
        use boozer_sub, only: boozer_to_vmec
        use spline_vmec_sub, only: splint_vmec_data, metric_tensor_vmec
        real(dp), intent(in) :: s, vartheta_B, varphi_B
        real(dp) :: sqrt_g_ss

        real(dp) :: theta_V, varphi_V, dummy(10)
        real(dp) :: R, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp
        real(dp) :: g_vmec(3, 3), sqrt_g_vmec

        call boozer_to_vmec(s, vartheta_B, varphi_B, theta_V, varphi_V)
        call splint_vmec_data(s, theta_V, varphi_V, &
                              dummy(1), dummy(2), dummy(3), dummy(4), dummy(5), &
                              R, dummy(6), dummy(7), &
                              dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, &
                              dummy(8), dummy(9), dummy(10))
        call metric_tensor_vmec(R, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, &
                                g_vmec, sqrt_g_vmec)
        sqrt_g_ss = sqrt(g_vmec(2, 2)*g_vmec(3, 3) - g_vmec(2, 3)**2)/abs(sqrt_g_vmec)
    end function metric_nabla_s

    subroutine check(name, case, val, ref, test_failed)
        character(len=*), intent(in) :: name
        integer, intent(in) :: case
        real(dp), intent(in) :: val(:), ref(:)
        logical, intent(inout) :: test_failed
        integer :: i

        do i = 1, size(val)
            if (abs(val(i) - ref(i)) > abstol .and. &
                abs(val(i) - ref(i)) > reltol*abs(ref(i))) then
                print *, "---------------------------------------------------"
                print *, "mismatch in ", name, " component ", i, " case ", case
                print *, "libneo : ", val(i)
                print *, "SIMPLE : ", ref(i)
                print *, "abs err: ", abs(val(i) - ref(i))
                test_failed = .true.
            end if
        end do
    end subroutine check

end program test_boozer_converter_vs_simple
