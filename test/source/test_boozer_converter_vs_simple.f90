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
    use boozer_sub, only: get_boozer_coordinates, splint_boozer_coord
    implicit none

    ! reltol pins the well-conditioned quantities (|B|, sqrt(g), the large
    ! covariant components). abstol is the absolute floor for the near-zero
    ! covariant components (e.g. B_theta/|B| ~ 1e-5 for QA), whose relative
    ! agreement is compile/platform-sensitive at the spline-interpolation level;
    ! 1e-10 absorbs that cross-platform FP noise without loosening the
    ! relative bound the large quantities meet.
    real(dp), parameter :: reltol = 1.0e-6_dp, abstol = 1.0e-10_dp
    character(len=*), parameter :: wout_file = "wout.nc"

    integer, parameter :: n_cases = 5
    real(dp) :: stor(n_cases), theta(n_cases), phi(n_cases)
    real(dp) :: bmod, sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)
    real(dp) :: bmod_ref(n_cases), sqrtg_ref(n_cases)
    real(dp) :: bder_ref(n_cases, 3), hcovar_ref(n_cases, 3)
    real(dp) :: hctrvr_ref(n_cases, 3), hcurl_ref(n_cases, 3)
    integer :: case
    logical :: test_failed

    use_B_r = .true.
    call get_boozer_coordinates(wout_file, radial_spline_order=5, &
                                angular_spline_order=5, grid_refinment=3)

    stor = [0.1_dp, 0.3_dp, 0.5_dp, 0.7_dp, 0.9_dp]
    theta = [0.0_dp, 1.0_dp, 3.14_dp, 0.5_dp, 2.0_dp]
    phi = [0.0_dp, 0.5_dp, 1.57_dp, 2.0_dp, 3.0_dp]

    bmod_ref = [5.7024014488290690e+04_dp, 5.7010408592594038e+04_dp, &
                6.3742758213874906e+04_dp, 5.4218536572575875e+04_dp, &
                6.1083220187449493e+04_dp]

    sqrtg_ref = [-1.5711629692175472e+07_dp, -1.5719091463432141e+07_dp, &
                 -1.2574007504120229e+07_dp, -1.7379594198471960e+07_dp, &
                 -1.3692768685651677e+07_dp]

    bder_ref(1, :) = [-1.7906807464749691e-01_dp, -7.4166439213088501e-14_dp, -1.0958457393698770e-13_dp]
    bder_ref(2, :) = [-6.5375422685029619e-02_dp, 4.8901137158983657e-02_dp, -4.2593419771991015e-04_dp]
    bder_ref(3, :) = [7.3646069743459427e-02_dp, 1.3567156911686458e-04_dp, -2.4743000679278641e-05_dp]
    bder_ref(4, :) = [-6.5308328098055424e-02_dp, 3.9538079263501176e-02_dp, -6.8133243699882112e-04_dp]
    bder_ref(5, :) = [1.2755344859923703e-02_dp, 1.0110736983359410e-01_dp, -7.9541221413346563e-04_dp]

    hcovar_ref(1, :) = [-9.4330617146665213e-13_dp, -9.4970349047515505e-04_dp, 1.1151159177109234e+03_dp]
    hcovar_ref(2, :) = [-4.9526222447404823e-02_dp, -3.8520822192834354e-04_dp, 1.1153790773146934e+03_dp]
    hcovar_ref(3, :) = [8.2784211500102169e-04_dp, -1.8842726486691493e-04_dp, 9.9757472565834667e+02_dp]
    hcovar_ref(4, :) = [-5.9028395703390649e-02_dp, -9.8956245690572466e-05_dp, 1.1728116578017384e+03_dp]
    hcovar_ref(5, :) = [2.8027264710120063e-02_dp, -1.0129726085459463e-05_dp, 1.0410079194890063e+03_dp]

    hctrvr_ref(1, :) = [0.0000000000000000e+00_dp, 3.7881309780855648e-04_dp, 8.9676807933375402e-04_dp]
    hctrvr_ref(2, :) = [0.0000000000000000e+00_dp, 3.7723748043895933e-04_dp, 8.9655630597133623e-04_dp]
    hctrvr_ref(3, :) = [0.0000000000000000e+00_dp, 4.2022807275813337e-04_dp, 1.0024312499722556e-03_dp]
    hctrvr_ref(4, :) = [0.0000000000000000e+00_dp, 3.5621175223770338e-04_dp, 8.5265185470933104e-04_dp]
    hctrvr_ref(5, :) = [0.0000000000000000e+00_dp, 4.0002199335749606e-04_dp, 9.6060748946365118e-04_dp]

    hcurl_ref(1, :) = [-5.2638894003780156e-18_dp, 1.2740530997096499e-05_dp, -8.6197541867643010e-09_dp]
    hcurl_ref(2, :) = [3.4698764369963150e-06_dp, 4.6481422398506964e-06_dp, -1.7483693891917036e-08_dp]
    hcurl_ref(3, :) = [1.0763674479658840e-08_dp, -5.7798304929806501e-06_dp, -9.7743856721982719e-09_dp]
    hcurl_ref(4, :) = [2.6681129427048551e-06_dp, 4.3919431197038171e-06_dp, -5.0485824545366451e-09_dp]
    hcurl_ref(5, :) = [7.6867998812914078e-06_dp, -9.5911521268646427e-07_dp, 7.0764342860889441e-09_dp]

    test_failed = .false.
    do case = 1, n_cases
        call evaluate_boozer(stor(case), theta(case), phi(case), &
                             bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        call check("bmod", case, [bmod], [bmod_ref(case)], test_failed)
        call check("sqrtg", case, [sqrtg], [sqrtg_ref(case)], test_failed)
        call check("bder", case, bder, bder_ref(case, :), test_failed)
        call check("hcovar", case, hcovar, hcovar_ref(case, :), test_failed)
        call check("hctrvr", case, hctrvr, hctrvr_ref(case, :), test_failed)
        call check("hcurl", case, hcurl, hcurl_ref(case, :), test_failed)
    end do

    if (test_failed) error stop "boozer converter differs from SIMPLE reference"
    print *, "boozer converter matches SIMPLE reference at all points"

contains

    !> Assemble |B|, sqrt(g) and field components from splint_boozer_coord,
    !> following the SIMPLE/rabe Boozer field convention, in native CGS units.
    subroutine evaluate_boozer(s, vartheta_B, varphi_B, &
                               bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        real(dp), intent(in) :: s, vartheta_B, varphi_B
        real(dp), intent(out) :: bmod, sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)

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
                                 B_r, dB_r, d2B_r)

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
