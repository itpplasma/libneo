module boozmn_file
!> Raw reader for booz_xform boozmn NetCDF output.
!>
!> Layout follows the booz_xform convention: Fourier harmonics live on the
!> half-grid surfaces listed in jlist, the profiles iota, buco, bvco, phi on
!> the full grid of ns surfaces, and ixn already contains the nfp factor.
!> Angles enter as cos(m*theta - n*zeta) for the *mnc and
!> sin(m*theta - n*zeta) for the *mns coefficients. All values keep the file
!> units (SI: T, m, Wb).

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none
    private

    public :: boozmn_data_t, read_boozmn

    type :: boozmn_data_t
        integer :: ns = 0
        integer :: nfp = 0
        integer :: nmodes = 0
        integer :: nsurf = 0
        logical :: asym = .false.
        integer, allocatable :: jlist(:)
        integer, allocatable :: ixm(:), ixn(:)
        real(dp), allocatable :: iota(:), buco(:), bvco(:), phi(:)
        real(dp), allocatable :: bmnc(:, :), rmnc(:, :), zmns(:, :), pmns(:, :)
        real(dp), allocatable :: bmns(:, :), rmns(:, :), zmnc(:, :), pmnc(:, :)
    end type boozmn_data_t

contains

    subroutine read_boozmn(filename, d)
        use nctools_module, only: nc_open, nc_close, nc_inq_dim, nc_get
        use netcdf, only: nf90_inq_varid, nf90_noerr

        character(len=*), intent(in) :: filename
        type(boozmn_data_t), intent(out) :: d

        integer :: ncid, varid, lasym

        call nc_open(trim(filename), ncid)
        call nc_get(ncid, 'ns_b', d%ns)
        call nc_get(ncid, 'nfp_b', d%nfp)
        call nc_inq_dim(ncid, 'ixm_b', d%nmodes)
        call nc_inq_dim(ncid, 'jlist', d%nsurf)

        lasym = 0
        if (nf90_inq_varid(ncid, 'lasym__logical__', varid) == nf90_noerr) then
            call nc_get(ncid, 'lasym__logical__', lasym)
        end if
        d%asym = lasym /= 0

        allocate (d%jlist(d%nsurf), d%ixm(d%nmodes), d%ixn(d%nmodes))
        allocate (d%iota(d%ns), d%buco(d%ns), d%bvco(d%ns), d%phi(d%ns))
        allocate (d%bmnc(d%nmodes, d%nsurf), d%rmnc(d%nmodes, d%nsurf), &
                  d%zmns(d%nmodes, d%nsurf), d%pmns(d%nmodes, d%nsurf))

        call nc_get(ncid, 'jlist', d%jlist)
        call nc_get(ncid, 'ixm_b', d%ixm)
        call nc_get(ncid, 'ixn_b', d%ixn)
        call nc_get(ncid, 'iota_b', d%iota)
        call nc_get(ncid, 'buco_b', d%buco)
        call nc_get(ncid, 'bvco_b', d%bvco)
        call nc_get(ncid, 'phi_b', d%phi)
        call nc_get(ncid, 'bmnc_b', d%bmnc)
        call nc_get(ncid, 'rmnc_b', d%rmnc)
        call nc_get(ncid, 'zmns_b', d%zmns)
        call nc_get(ncid, 'pmns_b', d%pmns)

        if (d%asym) then
            allocate (d%bmns(d%nmodes, d%nsurf), d%rmns(d%nmodes, d%nsurf), &
                      d%zmnc(d%nmodes, d%nsurf), d%pmnc(d%nmodes, d%nsurf))
            call nc_get(ncid, 'bmns_b', d%bmns)
            call nc_get(ncid, 'rmns_b', d%rmns)
            call nc_get(ncid, 'zmnc_b', d%zmnc)
            call nc_get(ncid, 'pmnc_b', d%pmnc)
        end if

        call nc_close(ncid)
    end subroutine read_boozmn

end module boozmn_file
