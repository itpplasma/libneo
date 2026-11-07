!> Round trip small synthetic Boozer .bc files, in the six-column
!> (stellarator-symmetric) and the ten-column layout, through read_bc_file and
!> check every field of bc_data_t against the values used to write them.
program test_bc_file
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use bc_file, only: bc_data_t, read_bc_file

    implicit none

    character(len=*), parameter :: bc_filename = "test_bc_file_synthetic.bc"
    character(len=*), parameter :: bc_asym_filename = "test_bc_file_asymmetric.bc"
    real(dp), parameter :: tol = 1.0e-10_dp

    integer, parameter :: nsurf = 3, nmode = 4
    integer, parameter :: m0b = 2, n0b = 1, nper = 5
    real(dp), parameter :: flux = 1.5_dp, a = 0.5_dp, R = 6.0_dp

    integer, parameter :: m(nmode) = [0, 1, 0, 2]
    integer, parameter :: n(nmode) = [0, 0, 1, 1]

    real(dp) :: s(nsurf), iota(nsurf), Jpol_over_nper(nsurf), Itor(nsurf)
    real(dp) :: rmnc(nsurf, nmode), zmns(nsurf, nmode)
    real(dp) :: vmns(nsurf, nmode), bmnc(nsurf, nmode)
    real(dp) :: rmns(nsurf, nmode), zmnc(nsurf, nmode)
    real(dp) :: vmnc(nsurf, nmode), bmns(nsurf, nmode)

    type(bc_data_t) :: d
    logical :: failed
    integer :: isurf, imode

    failed = .false.

    do isurf = 1, nsurf
        s(isurf) = 0.25_dp*real(isurf, dp)
        iota(isurf) = 0.25_dp + 0.0625_dp*real(isurf, dp)
        Jpol_over_nper(isurf) = 99.5_dp + 1.0_dp*real(isurf, dp)
        Itor(isurf) = 0.25_dp + 10.0_dp*real(isurf, dp)
        do imode = 1, nmode
            rmnc(isurf, imode) = 1.0_dp + 0.25_dp*real(isurf, dp) &
                                  + 0.0625_dp*real(imode, dp)
            zmns(isurf, imode) = -0.5_dp + 0.125_dp*real(isurf, dp) &
                                  + 0.03125_dp*real(imode, dp)
            vmns(isurf, imode) = 0.0625_dp*real(isurf, dp) &
                                  - 0.015625_dp*real(imode, dp)
            bmnc(isurf, imode) = 2.0_dp + 0.5_dp*real(isurf, dp) &
                                  + 0.125_dp*real(imode, dp)
            rmns(isurf, imode) = 0.01_dp*real(isurf + 2*imode, dp)
            zmnc(isurf, imode) = -0.02_dp*real(isurf + imode, dp)
            vmnc(isurf, imode) = 0.003_dp*real(2*isurf - imode, dp)
            bmns(isurf, imode) = 0.004_dp*real(isurf*imode, dp)
        end do
    end do

    call write_synthetic_bc_file(bc_filename)
    call read_bc_file(bc_filename, d)

    call check_int(d%m0b, m0b, "m0b", failed)
    call check_int(d%n0b, n0b, "n0b", failed)
    call check_int(d%nsurf, nsurf, "nsurf", failed)
    call check_int(d%nper, nper, "nper", failed)
    call check_int(d%nmode, nmode, "nmode", failed)
    call check_real(d%flux, flux, "flux", failed)
    call check_real(d%a, a, "a", failed)
    call check_real(d%R, R, "R", failed)

    if (.not. allocated(d%m) .or. .not. allocated(d%n)) then
        print *, "FAIL: mode arrays not allocated"
        failed = .true.
    else if (size(d%m) /= nmode .or. size(d%n) /= nmode) then
        print *, "FAIL: mode array size mismatch"
        failed = .true.
    else
        do imode = 1, nmode
            call check_int(d%m(imode), m(imode), "m(imode)", failed)
            call check_int(d%n(imode), n(imode), "n(imode)", failed)
        end do
    end if

    do isurf = 1, nsurf
        call check_real(d%s(isurf), s(isurf), "s(isurf)", failed)
        call check_real(d%iota(isurf), iota(isurf), "iota(isurf)", failed)
        call check_real(d%Jpol_over_nper(isurf), Jpol_over_nper(isurf), &
                        "Jpol_over_nper(isurf)", failed)
        call check_real(d%Itor(isurf), Itor(isurf), "Itor(isurf)", failed)
        do imode = 1, nmode
            call check_real(d%rmnc(isurf, imode), rmnc(isurf, imode), &
                            "rmnc(isurf,imode)", failed)
            call check_real(d%zmns(isurf, imode), zmns(isurf, imode), &
                            "zmns(isurf,imode)", failed)
            call check_real(d%vmns(isurf, imode), vmns(isurf, imode), &
                            "vmns(isurf,imode)", failed)
            call check_real(d%bmnc(isurf, imode), bmnc(isurf, imode), &
                            "bmnc(isurf,imode)", failed)
        end do
    end do

    if (.not. d%symmetric) then
        print *, "FAIL: six-column file read as asymmetric"
        failed = .true.
    end if
    if (any(d%rmns /= 0) .or. any(d%zmnc /= 0) .or. any(d%vmnc /= 0) &
        .or. any(d%bmns /= 0)) then
        print *, "FAIL: asymmetric terms of a six-column file are not zero"
        failed = .true.
    end if

    call write_synthetic_bc_file(bc_asym_filename, asymmetric=.true.)
    call read_bc_file(bc_asym_filename, d)
    if (d%symmetric) then
        print *, "FAIL: ten-column file read as symmetric"
        failed = .true.
    end if
    call check_int(d%nmode, nmode, "nmode (asymmetric)", failed)
    do isurf = 1, nsurf
        call check_real(d%iota(isurf), iota(isurf), "iota (asymmetric)", failed)
        do imode = 1, nmode
            call check_real(d%rmnc(isurf, imode), rmnc(isurf, imode), &
                            "rmnc (asymmetric)", failed)
            call check_real(d%rmns(isurf, imode), rmns(isurf, imode), "rmns", failed)
            call check_real(d%zmnc(isurf, imode), zmnc(isurf, imode), "zmnc", failed)
            call check_real(d%zmns(isurf, imode), zmns(isurf, imode), &
                            "zmns (asymmetric)", failed)
            call check_real(d%vmnc(isurf, imode), vmnc(isurf, imode), "vmnc", failed)
            call check_real(d%vmns(isurf, imode), vmns(isurf, imode), &
                            "vmns (asymmetric)", failed)
            call check_real(d%bmnc(isurf, imode), bmnc(isurf, imode), &
                            "bmnc (asymmetric)", failed)
            call check_real(d%bmns(isurf, imode), bmns(isurf, imode), "bmns", failed)
        end do
    end do

    if (failed) then
        error stop "test_bc_file failed"
    end if
    print *, "test_bc_file passed"

contains

    subroutine write_synthetic_bc_file(filename, asymmetric)
        character(len=*), intent(in) :: filename
        logical, intent(in), optional :: asymmetric

        integer :: iunit, is, k
        logical :: asym

        asym = .false.
        if (present(asymmetric)) asym = asymmetric

        open (newunit=iunit, file=trim(filename), status='replace', &
              action='write')

        write (iunit, '(A)') 'CC Boozer-coordinate data file'
        write (iunit, '(A)') 'CC Version:'
        write (iunit, '(A)') 'CC Author:'
        write (iunit, '(A)') 'CC shot:    0'
        write (iunit, '(A)') ' m0b   n0b  nsurf  nper    flux [Tm^2]'// &
            '        a [m]          R [m]'
        write (iunit, '(2I6, 2I6, 3ES16.8)') m0b, n0b, nsurf, nper, flux, a, R

        do is = 1, nsurf
            write (iunit, '(A)') '        s               iota'// &
                '           Jpol/nper          Itor            pprime'// &
                '         sqrt g(0,0)'
            write (iunit, '(A)') '                                          [A]'// &
                '           [A]             [Pa]         (dV/ds)/nper'
            write (iunit, '(6ES16.8)') s(is), iota(is), Jpol_over_nper(is), &
                Itor(is), 0.0_dp, 0.0_dp

            if (asym) then
                write (iunit, '(A)') '    m    n      rmnc [m]         rmns [m]'// &
                    '         zmnc [m]         zmns [m]         vmnc [ ]'// &
                    '         vmns [ ]         bmnc [T]         bmns [T]'
                do k = 1, nmode
                    write (iunit, '(2I6, 8ES16.8)') m(k), n(k), rmnc(is, k), &
                        rmns(is, k), zmnc(is, k), zmns(is, k), vmnc(is, k), &
                        vmns(is, k), bmnc(is, k), bmns(is, k)
                end do
            else
                write (iunit, '(A)') '    m    n      rmnc [m]'// &
                    '         zmns [m]         vmns [ ]         bmnc [T]'
                do k = 1, nmode
                    write (iunit, '(2I6, 4ES16.8)') m(k), n(k), rmnc(is, k), &
                        zmns(is, k), vmns(is, k), bmnc(is, k)
                end do
            end if
        end do

        close (iunit)
    end subroutine write_synthetic_bc_file

    subroutine check_int(actual, expected, label, failed)
        integer, intent(in) :: actual, expected
        character(len=*), intent(in) :: label
        logical, intent(inout) :: failed

        if (actual /= expected) then
            print *, "FAIL: ", label, " = ", actual, " expected ", expected
            failed = .true.
        end if
    end subroutine check_int

    subroutine check_real(actual, expected, label, failed)
        real(dp), intent(in) :: actual, expected
        character(len=*), intent(in) :: label
        logical, intent(inout) :: failed

        if (abs(actual - expected) > tol*max(1.0_dp, abs(expected))) then
            print *, "FAIL: ", label, " = ", actual, " expected ", expected
            failed = .true.
        end if
    end subroutine check_real

end program test_bc_file
