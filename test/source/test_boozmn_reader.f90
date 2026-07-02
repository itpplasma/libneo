!> Pin the boozmn NetCDF reader and the boozmn -> Boozer-spline loader.
!>
!> Writes synthetic booz_xform boozmn files (symmetric and asymmetric),
!> checks the raw fields read_boozmn returns, then loads the symmetric file
!> via load_boozer_from_boozmn and compares |B| from splint_boozer_coord
!> against the analytic Fourier sum of the fixture modes.
program test_boozmn_reader
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use netcdf
    use boozmn_file, only: boozmn_data_t, read_boozmn
    use boozmn_reader, only: load_boozer_from_boozmn
    use boozer_sub, only: splint_boozer_coord

    implicit none

    character(len=*), parameter :: sym_file = 'boozmn_reader_sym.nc'
    character(len=*), parameter :: asym_file = 'boozmn_reader_asym.nc'
    integer, parameter :: ns_full = 5
    integer, parameter :: nsurf = 3
    integer, parameter :: nmode = 4
    integer, dimension(nmode), parameter :: ixm = [0, 1, 1, 2]
    integer, dimension(nmode), parameter :: ixn = [0, 1, -1, 0]
    real(dp), dimension(nmode), parameter :: bmn = &
        [2.0_dp, 0.1_dp, -0.07_dp, 0.02_dp]
    integer, dimension(nsurf), parameter :: jlist = [2, 3, 4]
    real(dp), dimension(ns_full), parameter :: iota_full = &
        [0.5_dp, 0.55_dp, 0.6_dp, 0.65_dp, 0.7_dp]
    real(dp), dimension(ns_full), parameter :: phi_full = &
        [0.0_dp, 0.2_dp, 0.5_dp, 0.9_dp, 1.3_dp]

    real(dp), parameter :: TWOPI = 8.0_dp*atan(1.0_dp)
    real(dp), parameter :: GAUSS_PER_T = 1.0e4_dp
    real(dp), parameter :: stor = 0.5_dp
    real(dp), parameter :: theta = TWOPI*7.0_dp/48.0_dp
    real(dp), parameter :: zeta = TWOPI*11.0_dp/96.0_dp
    real(dp), parameter :: reltol = 1.0e-6_dp

    type(boozmn_data_t) :: d
    real(dp) :: bmod, expected

    call write_boozmn_fixture(sym_file, asym=.false.)
    call write_boozmn_fixture(asym_file, asym=.true.)

    call read_boozmn(sym_file, d)
    call check_raw_symmetric(d)

    call read_boozmn(asym_file, d)
    if (.not. d%asym) error stop 'asym flag not read'
    if (.not. allocated(d%bmns)) error stop 'bmns not allocated for asym file'
    call assert_close(d%bmns(2, 1), 0.5_dp*bmn(2), 'bmns spot value')
    call assert_close(d%pmnc(3, 2), 0.25_dp*bmn(3), 'pmnc spot value')
    print *, 'OK raw asymmetric read'

    call load_boozer_from_boozmn(sym_file)
    call eval_bmod(stor, theta, zeta, bmod)
    expected = expected_bmod(theta, zeta)*GAUSS_PER_T
    print *, 'boozmn Bmod:', bmod, ' expected:', expected
    if (abs(bmod - expected) > reltol*abs(expected)) then
        error stop 'boozmn Bmod mismatch'
    end if
    print *, 'OK boozmn -> Boozer splines'

contains

    subroutine check_raw_symmetric(b)
        type(boozmn_data_t), intent(in) :: b

        if (b%ns /= ns_full) error stop 'ns mismatch'
        if (b%nfp /= 1) error stop 'nfp mismatch'
        if (b%nmodes /= nmode) error stop 'nmodes mismatch'
        if (b%nsurf /= nsurf) error stop 'nsurf mismatch'
        if (b%asym) error stop 'asym flag set for symmetric file'
        if (any(b%jlist /= jlist)) error stop 'jlist mismatch'
        if (any(b%ixm /= ixm)) error stop 'ixm mismatch'
        if (any(b%ixn /= ixn)) error stop 'ixn mismatch'
        call assert_close(b%iota(3), iota_full(3), 'iota spot value')
        call assert_close(b%phi(ns_full), phi_full(ns_full), 'phi spot value')
        call assert_close(b%bmnc(2, 1), bmn(2), 'bmnc spot value')
        call assert_close(b%zmns(4, 2), 0.1_dp*bmn(4), 'zmns spot value')
        call assert_close(b%pmns(1, 3), 0.2_dp*bmn(1), 'pmns spot value')
        print *, 'OK raw symmetric read'
    end subroutine check_raw_symmetric

    subroutine assert_close(actual, ref, what)
        real(dp), intent(in) :: actual, ref
        character(len=*), intent(in) :: what

        if (abs(actual - ref) > 1.0e-12_dp*max(1.0_dp, abs(ref))) then
            print *, 'FAIL: ', what, ' actual:', actual, ' ref:', ref
            error stop 'raw boozmn value mismatch'
        end if
    end subroutine assert_close

    function expected_bmod(theta_in, zeta_in) result(value)
        real(dp), intent(in) :: theta_in, zeta_in
        real(dp) :: value
        integer :: imode

        value = 0.0_dp
        do imode = 1, nmode
            value = value + bmn(imode)*cos(real(ixm(imode), dp)*theta_in &
                                           - real(ixn(imode), dp)*zeta_in)
        end do
    end function expected_bmod

    subroutine eval_bmod(s_in, theta_in, zeta_in, bmod_out)
        real(dp), intent(in) :: s_in, theta_in, zeta_in
        real(dp), intent(out) :: bmod_out

        real(dp) :: A_theta, A_phi, dA_theta_dr, dA_phi_dr
        real(dp) :: d2A_phi_dr2, d3A_phi_dr3
        real(dp) :: Bth, dBth, d2Bth, Bph, dBph, d2Bph
        real(dp) :: dBmod(3), d2Bmod(6), Br, dBr(3), d2Br(6)

        call splint_boozer_coord(s_in, theta_in, zeta_in, 0, &
                                 A_theta, A_phi, dA_theta_dr, dA_phi_dr, &
                                 d2A_phi_dr2, d3A_phi_dr3, Bth, dBth, d2Bth, &
                                 Bph, dBph, d2Bph, bmod_out, dBmod, d2Bmod, &
                                 Br, dBr, d2Br)
    end subroutine eval_bmod

    subroutine write_boozmn_fixture(path, asym)
        character(len=*), intent(in) :: path
        logical, intent(in) :: asym

        integer :: ncid, dim_radius, dim_mode, dim_surf
        integer :: var_ns, var_nfp, var_lasym
        integer :: var_jlist, var_ixm, var_ixn
        integer :: var_iota, var_buco, var_bvco, var_phi
        integer :: var_bmnc, var_rmnc, var_zmns, var_pmns
        integer :: var_bmns, var_rmns, var_zmnc, var_pmnc
        integer :: lasym_int
        real(dp), dimension(ns_full) :: buco, bvco
        real(dp), dimension(nmode, nsurf) :: coeff

        call check_nc(nf90_create(path, nf90_clobber, ncid), 'create boozmn')
        call check_nc(nf90_def_dim(ncid, 'radius', ns_full, dim_radius), &
                      'define radius')
        call check_nc(nf90_def_dim(ncid, 'mn_mode', nmode, dim_mode), &
                      'define mn_mode')
        call check_nc(nf90_def_dim(ncid, 'comput_surfs', nsurf, dim_surf), &
                      'define comput_surfs')

        call check_nc(nf90_def_var(ncid, 'ns_b', nf90_int, varid=var_ns), &
                      'define ns_b')
        call check_nc(nf90_def_var(ncid, 'nfp_b', nf90_int, varid=var_nfp), &
                      'define nfp_b')
        call check_nc(nf90_def_var(ncid, 'lasym__logical__', nf90_int, &
                                   varid=var_lasym), 'define lasym')
        call check_nc(nf90_def_var(ncid, 'jlist', nf90_int, [dim_surf], &
                                   var_jlist), 'define jlist')
        call check_nc(nf90_def_var(ncid, 'ixm_b', nf90_int, [dim_mode], &
                                   var_ixm), 'define ixm_b')
        call check_nc(nf90_def_var(ncid, 'ixn_b', nf90_int, [dim_mode], &
                                   var_ixn), 'define ixn_b')
        call check_nc(nf90_def_var(ncid, 'iota_b', nf90_double, [dim_radius], &
                                   var_iota), 'define iota_b')
        call check_nc(nf90_def_var(ncid, 'buco_b', nf90_double, [dim_radius], &
                                   var_buco), 'define buco_b')
        call check_nc(nf90_def_var(ncid, 'bvco_b', nf90_double, [dim_radius], &
                                   var_bvco), 'define bvco_b')
        call check_nc(nf90_def_var(ncid, 'phi_b', nf90_double, [dim_radius], &
                                   var_phi), 'define phi_b')
        call check_nc(nf90_def_var(ncid, 'bmnc_b', nf90_double, &
                                   [dim_mode, dim_surf], var_bmnc), 'define bmnc_b')
        call check_nc(nf90_def_var(ncid, 'rmnc_b', nf90_double, &
                                   [dim_mode, dim_surf], var_rmnc), 'define rmnc_b')
        call check_nc(nf90_def_var(ncid, 'zmns_b', nf90_double, &
                                   [dim_mode, dim_surf], var_zmns), 'define zmns_b')
        call check_nc(nf90_def_var(ncid, 'pmns_b', nf90_double, &
                                   [dim_mode, dim_surf], var_pmns), 'define pmns_b')
        if (asym) then
            call check_nc(nf90_def_var(ncid, 'bmns_b', nf90_double, &
                                       [dim_mode, dim_surf], var_bmns), &
                          'define bmns_b')
            call check_nc(nf90_def_var(ncid, 'rmns_b', nf90_double, &
                                       [dim_mode, dim_surf], var_rmns), &
                          'define rmns_b')
            call check_nc(nf90_def_var(ncid, 'zmnc_b', nf90_double, &
                                       [dim_mode, dim_surf], var_zmnc), &
                          'define zmnc_b')
            call check_nc(nf90_def_var(ncid, 'pmnc_b', nf90_double, &
                                       [dim_mode, dim_surf], var_pmnc), &
                          'define pmnc_b')
        end if
        call check_nc(nf90_enddef(ncid), 'end definitions')

        buco = [0.01_dp, 0.02_dp, 0.03_dp, 0.04_dp, 0.05_dp]
        bvco = [0.2_dp, 0.21_dp, 0.22_dp, 0.23_dp, 0.24_dp]
        lasym_int = 0
        if (asym) lasym_int = 1

        call check_nc(nf90_put_var(ncid, var_ns, ns_full), 'write ns_b')
        call check_nc(nf90_put_var(ncid, var_nfp, 1), 'write nfp_b')
        call check_nc(nf90_put_var(ncid, var_lasym, lasym_int), 'write lasym')
        call check_nc(nf90_put_var(ncid, var_jlist, jlist), 'write jlist')
        call check_nc(nf90_put_var(ncid, var_ixm, ixm), 'write ixm_b')
        call check_nc(nf90_put_var(ncid, var_ixn, ixn), 'write ixn_b')
        call check_nc(nf90_put_var(ncid, var_iota, iota_full), 'write iota_b')
        call check_nc(nf90_put_var(ncid, var_buco, buco), 'write buco_b')
        call check_nc(nf90_put_var(ncid, var_bvco, bvco), 'write bvco_b')
        call check_nc(nf90_put_var(ncid, var_phi, phi_full), 'write phi_b')

        coeff = spread(bmn, dim=2, ncopies=nsurf)
        call check_nc(nf90_put_var(ncid, var_bmnc, coeff), 'write bmnc_b')
        coeff = 0.0_dp
        coeff(1, :) = [1.8_dp, 1.9_dp, 2.0_dp]
        call check_nc(nf90_put_var(ncid, var_rmnc, coeff), 'write rmnc_b')
        coeff = 0.1_dp*spread(bmn, dim=2, ncopies=nsurf)
        call check_nc(nf90_put_var(ncid, var_zmns, coeff), 'write zmns_b')
        coeff = 0.2_dp*spread(bmn, dim=2, ncopies=nsurf)
        call check_nc(nf90_put_var(ncid, var_pmns, coeff), 'write pmns_b')
        if (asym) then
            coeff = 0.5_dp*spread(bmn, dim=2, ncopies=nsurf)
            call check_nc(nf90_put_var(ncid, var_bmns, coeff), 'write bmns_b')
            coeff = 0.3_dp*spread(bmn, dim=2, ncopies=nsurf)
            call check_nc(nf90_put_var(ncid, var_rmns, coeff), 'write rmns_b')
            coeff = 0.4_dp*spread(bmn, dim=2, ncopies=nsurf)
            call check_nc(nf90_put_var(ncid, var_zmnc, coeff), 'write zmnc_b')
            coeff = 0.25_dp*spread(bmn, dim=2, ncopies=nsurf)
            call check_nc(nf90_put_var(ncid, var_pmnc, coeff), 'write pmnc_b')
        end if
        call check_nc(nf90_close(ncid), 'close boozmn')
    end subroutine write_boozmn_fixture

    subroutine check_nc(status, context)
        integer, intent(in) :: status
        character(len=*), intent(in) :: context

        if (status /= nf90_noerr) then
            print *, trim(context), ': ', trim(nf90_strerror(status))
            error stop 'boozmn fixture NetCDF error'
        end if
    end subroutine check_nc

end program test_boozmn_reader
