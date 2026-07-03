module boozmn_file
!> Raw reader and writer for booz_xform boozmn NetCDF output.
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

    public :: boozmn_data_t, read_boozmn, write_boozmn

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

    subroutine write_boozmn(filename, d)
        use netcdf, only: nf90_create, nf90_def_dim, nf90_enddef, nf90_close, &
                          nf90_clobber, nf90_int, nf90_double

        character(len=*), intent(in) :: filename
        type(boozmn_data_t), intent(in) :: d

        integer :: ncid, dim_radius, dim_mode, dim_surf, lasym

        call require_components(d)

        call chk(nf90_create(trim(filename), nf90_clobber, ncid), &
                 'create '//trim(filename))
        call chk(nf90_def_dim(ncid, 'radius', d%ns, dim_radius), &
                 'define radius')
        call chk(nf90_def_dim(ncid, 'mn_mode', d%nmodes, dim_mode), &
                 'define mn_mode')
        call chk(nf90_def_dim(ncid, 'comput_surfs', d%nsurf, dim_surf), &
                 'define comput_surfs')

        call def_var(ncid, 'ns_b', nf90_int)
        call def_var(ncid, 'nfp_b', nf90_int)
        call def_var(ncid, 'mnboz_b', nf90_int)
        call def_var(ncid, 'lasym__logical__', nf90_int)
        call def_var(ncid, 'jlist', nf90_int, [dim_surf])
        call def_var(ncid, 'ixm_b', nf90_int, [dim_mode])
        call def_var(ncid, 'ixn_b', nf90_int, [dim_mode])
        call def_var(ncid, 'iota_b', nf90_double, [dim_radius])
        call def_var(ncid, 'buco_b', nf90_double, [dim_radius])
        call def_var(ncid, 'bvco_b', nf90_double, [dim_radius])
        call def_var(ncid, 'phi_b', nf90_double, [dim_radius])
        call def_var(ncid, 'bmnc_b', nf90_double, [dim_mode, dim_surf])
        call def_var(ncid, 'rmnc_b', nf90_double, [dim_mode, dim_surf])
        call def_var(ncid, 'zmns_b', nf90_double, [dim_mode, dim_surf])
        call def_var(ncid, 'pmns_b', nf90_double, [dim_mode, dim_surf])
        if (d%asym) then
            call def_var(ncid, 'bmns_b', nf90_double, [dim_mode, dim_surf])
            call def_var(ncid, 'rmns_b', nf90_double, [dim_mode, dim_surf])
            call def_var(ncid, 'zmnc_b', nf90_double, [dim_mode, dim_surf])
            call def_var(ncid, 'pmnc_b', nf90_double, [dim_mode, dim_surf])
        end if
        call chk(nf90_enddef(ncid), 'end definitions')

        lasym = 0
        if (d%asym) lasym = 1
        call put_int_0(ncid, 'ns_b', d%ns)
        call put_int_0(ncid, 'nfp_b', d%nfp)
        call put_int_0(ncid, 'mnboz_b', d%nmodes)
        call put_int_0(ncid, 'lasym__logical__', lasym)
        call put_int_1(ncid, 'jlist', d%jlist)
        call put_int_1(ncid, 'ixm_b', d%ixm)
        call put_int_1(ncid, 'ixn_b', d%ixn)
        call put_double_1(ncid, 'iota_b', d%iota)
        call put_double_1(ncid, 'buco_b', d%buco)
        call put_double_1(ncid, 'bvco_b', d%bvco)
        call put_double_1(ncid, 'phi_b', d%phi)
        call put_double_2(ncid, 'bmnc_b', d%bmnc)
        call put_double_2(ncid, 'rmnc_b', d%rmnc)
        call put_double_2(ncid, 'zmns_b', d%zmns)
        call put_double_2(ncid, 'pmns_b', d%pmns)
        if (d%asym) then
            call put_double_2(ncid, 'bmns_b', d%bmns)
            call put_double_2(ncid, 'rmns_b', d%rmns)
            call put_double_2(ncid, 'zmnc_b', d%zmnc)
            call put_double_2(ncid, 'pmnc_b', d%pmnc)
        end if
        call chk(nf90_close(ncid), 'close '//trim(filename))
    end subroutine write_boozmn

    subroutine require_components(d)
        type(boozmn_data_t), intent(in) :: d

        call require_allocated(allocated(d%jlist), 'jlist')
        call require_allocated(allocated(d%ixm), 'ixm')
        call require_allocated(allocated(d%ixn), 'ixn')
        call require_allocated(allocated(d%iota), 'iota')
        call require_allocated(allocated(d%buco), 'buco')
        call require_allocated(allocated(d%bvco), 'bvco')
        call require_allocated(allocated(d%phi), 'phi')
        call require_allocated(allocated(d%bmnc), 'bmnc')
        call require_allocated(allocated(d%rmnc), 'rmnc')
        call require_allocated(allocated(d%zmns), 'zmns')
        call require_allocated(allocated(d%pmns), 'pmns')
        if (d%asym) then
            call require_allocated(allocated(d%bmns), 'bmns')
            call require_allocated(allocated(d%rmns), 'rmns')
            call require_allocated(allocated(d%zmnc), 'zmnc')
            call require_allocated(allocated(d%pmnc), 'pmnc')
        end if
    end subroutine require_components

    subroutine require_allocated(is_allocated, name)
        logical, intent(in) :: is_allocated
        character(len=*), intent(in) :: name

        if (.not. is_allocated) then
            print *, 'write_boozmn: component not allocated: ', name
            error stop 'write_boozmn: missing boozmn_data_t component'
        end if
    end subroutine require_allocated

    subroutine def_var(ncid, name, xtype, dimids)
        use netcdf, only: nf90_def_var

        integer, intent(in) :: ncid, xtype
        character(len=*), intent(in) :: name
        integer, intent(in), optional :: dimids(:)

        integer :: varid

        if (present(dimids)) then
            call chk(nf90_def_var(ncid, name, xtype, dimids, varid), &
                     'define '//name)
        else
            call chk(nf90_def_var(ncid, name, xtype, varid=varid), &
                     'define '//name)
        end if
    end subroutine def_var

    subroutine put_int_0(ncid, name, value)
        use netcdf, only: nf90_inq_varid, nf90_put_var

        integer, intent(in) :: ncid, value
        character(len=*), intent(in) :: name

        integer :: varid

        call chk(nf90_inq_varid(ncid, name, varid), 'find '//name)
        call chk(nf90_put_var(ncid, varid, value), 'write '//name)
    end subroutine put_int_0

    subroutine put_int_1(ncid, name, values)
        use netcdf, only: nf90_inq_varid, nf90_put_var

        integer, intent(in) :: ncid
        character(len=*), intent(in) :: name
        integer, intent(in) :: values(:)

        integer :: varid

        call chk(nf90_inq_varid(ncid, name, varid), 'find '//name)
        call chk(nf90_put_var(ncid, varid, values), 'write '//name)
    end subroutine put_int_1

    subroutine put_double_1(ncid, name, values)
        use netcdf, only: nf90_inq_varid, nf90_put_var

        integer, intent(in) :: ncid
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: values(:)

        integer :: varid

        call chk(nf90_inq_varid(ncid, name, varid), 'find '//name)
        call chk(nf90_put_var(ncid, varid, values), 'write '//name)
    end subroutine put_double_1

    subroutine put_double_2(ncid, name, values)
        use netcdf, only: nf90_inq_varid, nf90_put_var

        integer, intent(in) :: ncid
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: values(:, :)

        integer :: varid

        call chk(nf90_inq_varid(ncid, name, varid), 'find '//name)
        call chk(nf90_put_var(ncid, varid, values), 'write '//name)
    end subroutine put_double_2

    subroutine chk(status, context)
        use netcdf, only: nf90_noerr, nf90_strerror

        integer, intent(in) :: status
        character(len=*), intent(in) :: context

        if (status /= nf90_noerr) then
            print *, trim(context), ': ', trim(nf90_strerror(status))
            error stop 'write_boozmn NetCDF error'
        end if
    end subroutine chk

end module boozmn_file
