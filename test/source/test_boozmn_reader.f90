!> Pin the boozmn NetCDF writer/reader and the boozmn -> Boozer-spline loader.
!>
!> Writes synthetic booz_xform boozmn files (symmetric and asymmetric) via
!> write_boozmn, checks the raw fields read_boozmn returns (round trip), then
!> loads the symmetric file via load_boozer_from_boozmn and compares |B| from
!> splint_boozer_coord against the analytic Fourier sum of the fixture modes.
program test_boozmn_reader
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use boozmn_file, only: boozmn_data_t, read_boozmn, write_boozmn
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

        type(boozmn_data_t) :: b

        b%ns = ns_full
        b%nfp = 1
        b%nmodes = nmode
        b%nsurf = nsurf
        b%asym = asym
        allocate (b%jlist(nsurf), b%ixm(nmode), b%ixn(nmode))
        allocate (b%iota(ns_full), b%buco(ns_full), b%bvco(ns_full), &
            b%phi(ns_full))
        allocate (b%bmnc(nmode, nsurf), b%zmns(nmode, nsurf), &
            b%pmns(nmode, nsurf))
        allocate (b%rmnc(nmode, nsurf), source=0.0_dp)
        b%jlist = jlist
        b%ixm = ixm
        b%ixn = ixn
        b%iota = iota_full
        b%buco = [0.01_dp, 0.02_dp, 0.03_dp, 0.04_dp, 0.05_dp]
        b%bvco = [0.2_dp, 0.21_dp, 0.22_dp, 0.23_dp, 0.24_dp]
        b%phi = phi_full
        b%bmnc = spread(bmn, dim=2, ncopies=nsurf)
        b%rmnc(1, :) = [1.8_dp, 1.9_dp, 2.0_dp]
        b%zmns = 0.1_dp*b%bmnc
        b%pmns = 0.2_dp*b%bmnc
        if (asym) then
            allocate (b%bmns(nmode, nsurf), b%rmns(nmode, nsurf), &
                b%zmnc(nmode, nsurf), b%pmnc(nmode, nsurf))
            b%bmns = 0.5_dp*b%bmnc
            b%rmns = 0.3_dp*b%bmnc
            b%zmnc = 0.4_dp*b%bmnc
            b%pmnc = 0.25_dp*b%bmnc
        end if

        call write_boozmn(path, b)
    end subroutine write_boozmn_fixture

end program test_boozmn_reader
