module neo_perturbation_field_netcdf
    ! NetCDF classic file for neo_perturbation_field (Fortio nf90 interface).
    !
    ! Dimensions: R(nR), Z(nZ), mode(n_modes).
    ! Variables (Fortran index order; C/ncdump order is reversed):
    !   R(R), Z(Z)                       grid, equidistant, R(1) > 0
    !   ntor(mode)                       toroidal mode numbers n (int)
    !   dA_R_re, dA_R_im (R, Z, mode)    real/imaginary parts of the amplitudes,
    !   dA_phi_re, dA_phi_im, dA_Z_re, dA_Z_im   physical cylindrical components
    !   dPhi_re, dPhi_im (R, Z, mode)    optional scalar potential
    ! Global attributes:
    !   format = "libneo_perturbation_field", format_version = 1
    !   units = "SI" | "Gaussian" (no conversion is applied on read)
    !   length_unit, vector_potential_unit, scalar_potential_unit, field_unit
    !   coordinates, convention (human-readable statement of the sign convention)
    !   has_potential = 0 | 1
    use netcdf, only: nf90_create, nf90_open, nf90_close, nf90_def_dim, nf90_def_var, &
        nf90_enddef, nf90_put_var, nf90_get_var, nf90_put_att, &
        nf90_get_att, nf90_inq_varid, nf90_inq_dimid, &
        nf90_inquire_dimension, NF90_NOERR, NF90_CLOBBER, &
        NF90_NOWRITE, NF90_GLOBAL, NF90_DOUBLE, NF90_INT
    use neo_perturbation_field, only: perturbation_field_t, UNITS_SI, UNITS_GAUSSIAN, &
        PERTFIELD_OK, PERTFIELD_BAD_INPUT
    implicit none
    private

    integer, parameter :: dp = kind(1.0d0)

    character(len=*), parameter, public :: &
        PERTFIELD_FORMAT = 'libneo_perturbation_field'
    integer, parameter, public :: PERTFIELD_FORMAT_VERSION = 1
    integer, parameter, public :: PERTFIELD_IO_ERROR = 10

    character(len=*), parameter :: COORDINATES = 'right-handed (R, phi, Z); phi is '// &
        'the geometric toroidal angle, counterclockwise '// &
        'seen from +Z; physical cylindrical components'
    character(len=*), parameter :: CONVENTION = 'real field = Re(sum_n amplitude_n'// &
        '(R, Z) * exp(+i n phi))'
    character(len=6), parameter :: COMP(3) = ['dA_R  ', 'dA_phi', 'dA_Z  ']

    public :: write_perturbation_field_netcdf, read_perturbation_field_netcdf

contains

    subroutine write_perturbation_field_netcdf(field, filename, ierr)
        type(perturbation_field_t), intent(in) :: field
        character(len=*), intent(in) :: filename
        integer, intent(out) :: ierr

        integer :: ncid, dims(3), c, vid, st

        ierr = PERTFIELD_BAD_INPUT
        if (field%n_modes < 1) return
        ierr = PERTFIELD_IO_ERROR
        if (nf90_create(filename, NF90_CLOBBER, ncid) /= NF90_NOERR) return
        st = define_header(ncid, field, dims)
        do c = 1, 3
            if (st == NF90_NOERR) st = def_complex(ncid, trim(COMP(c)), dims)
        end do
        if (field%has_potential) then
            if (st == NF90_NOERR) st = def_complex(ncid, 'dPhi', dims)
        end if
        if (st == NF90_NOERR) st = nf90_enddef(ncid)
        if (st == NF90_NOERR) st = put_r1(ncid, 'R', field%R)
        if (st == NF90_NOERR) st = put_r1(ncid, 'Z', field%Z)
        if (st == NF90_NOERR) st = nf90_inq_varid(ncid, 'ntor', vid)
        if (st == NF90_NOERR) st = nf90_put_var(ncid, vid, field%ntor)
        do c = 1, 3
            if (st == NF90_NOERR) st = put_complex(ncid, trim(COMP(c)), &
                field%dA(c, :, :, :))
        end do
        if (field%has_potential) then
            if (st == NF90_NOERR) st = put_complex(ncid, 'dPhi', field%dPhi)
        end if
        if (nf90_close(ncid) /= NF90_NOERR) return
        if (st /= NF90_NOERR) return
        ierr = PERTFIELD_OK
    end subroutine write_perturbation_field_netcdf

    integer function define_header(ncid, field, dims) result(st)
        integer, intent(in) :: ncid
        type(perturbation_field_t), intent(in) :: field
        integer, intent(out) :: dims(3)

        integer :: vid, has_phi
        character(len=8) :: u_len, u_a, u_phi, u_b

        if (field%units == UNITS_SI) then
            u_len = 'm'; u_a = 'T m'; u_phi = 'V'; u_b = 'T'
        else
            u_len = 'cm'; u_a = 'G cm'; u_phi = 'statV'; u_b = 'G'
        end if
        has_phi = merge(1, 0, field%has_potential)
        st = nf90_def_dim(ncid, 'R', size(field%R), dims(1))
        if (st == NF90_NOERR) st = nf90_def_dim(ncid, 'Z', size(field%Z), dims(2))
        if (st == NF90_NOERR) st = nf90_def_dim(ncid, 'mode', field%n_modes, dims(3))
        if (st == NF90_NOERR) st = nf90_def_var(ncid, 'R', NF90_DOUBLE, dims(1), vid)
        if (st == NF90_NOERR) st = nf90_put_att(ncid, vid, 'units', trim(u_len))
        if (st == NF90_NOERR) st = nf90_def_var(ncid, 'Z', NF90_DOUBLE, dims(2), vid)
        if (st == NF90_NOERR) st = nf90_put_att(ncid, vid, 'units', trim(u_len))
        if (st == NF90_NOERR) st = nf90_def_var(ncid, 'ntor', NF90_INT, dims(3), vid)
        if (st == NF90_NOERR) st = gatt(ncid, 'format', PERTFIELD_FORMAT)
        if (st == NF90_NOERR) st = nf90_put_att(ncid, NF90_GLOBAL, 'format_version', &
            PERTFIELD_FORMAT_VERSION)
        if (st == NF90_NOERR) st = gatt(ncid, 'units', trim(field%units))
        if (st == NF90_NOERR) st = gatt(ncid, 'length_unit', trim(u_len))
        if (st == NF90_NOERR) st = gatt(ncid, 'vector_potential_unit', trim(u_a))
        if (st == NF90_NOERR) st = gatt(ncid, 'scalar_potential_unit', trim(u_phi))
        if (st == NF90_NOERR) st = gatt(ncid, 'field_unit', trim(u_b))
        if (st == NF90_NOERR) st = gatt(ncid, 'coordinates', COORDINATES)
        if (st == NF90_NOERR) st = gatt(ncid, 'convention', CONVENTION)
        if (st == NF90_NOERR) st = nf90_put_att(ncid, NF90_GLOBAL, 'has_potential', &
            has_phi)
    end function define_header

    integer function gatt(ncid, name, value) result(st)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: name, value

        st = nf90_put_att(ncid, NF90_GLOBAL, name, value)
    end function gatt

    integer function def_complex(ncid, name, dims) result(st)
        integer, intent(in) :: ncid, dims(3)
        character(len=*), intent(in) :: name

        integer :: vid

        st = nf90_def_var(ncid, name//'_re', NF90_DOUBLE, dims, vid)
        if (st == NF90_NOERR) st = nf90_def_var(ncid, name//'_im', NF90_DOUBLE, dims, &
            vid)
    end function def_complex

    integer function put_r1(ncid, name, x) result(st)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: x(:)

        integer :: vid

        st = nf90_inq_varid(ncid, name, vid)
        if (st == NF90_NOERR) st = nf90_put_var(ncid, vid, x)
    end function put_r1

    integer function put_complex(ncid, name, a) result(st)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: name
        complex(dp), intent(in) :: a(:, :, :)

        real(dp) :: part(size(a, 1), size(a, 2), size(a, 3))
        integer :: vid

        st = nf90_inq_varid(ncid, name//'_re', vid)
        part = a%re
        if (st == NF90_NOERR) st = nf90_put_var(ncid, vid, part)
        if (st == NF90_NOERR) st = nf90_inq_varid(ncid, name//'_im', vid)
        part = a%im
        if (st == NF90_NOERR) st = nf90_put_var(ncid, vid, part)
    end function put_complex

    subroutine read_perturbation_field_netcdf(filename, field, ierr)
        character(len=*), intent(in) :: filename
        type(perturbation_field_t), intent(inout) :: field
        integer, intent(out) :: ierr

        integer :: ncid

        ierr = PERTFIELD_IO_ERROR
        if (nf90_open(filename, NF90_NOWRITE, ncid) /= NF90_NOERR) return
        call read_open_file(ncid, field, ierr)
        if (nf90_close(ncid) /= NF90_NOERR) ierr = PERTFIELD_IO_ERROR
    end subroutine read_perturbation_field_netcdf

    subroutine read_open_file(ncid, field, ierr)
        integer, intent(in) :: ncid
        type(perturbation_field_t), intent(inout) :: field
        integer, intent(out) :: ierr

        character(len=64) :: fmt, units
        integer :: version, has_phi, nR, nZ, nm, vid, c, st
        real(dp), allocatable :: R(:), Z(:)
        integer, allocatable :: ntor(:)
        complex(dp), allocatable :: dA(:, :, :, :), dPhi(:, :, :)

        ierr = PERTFIELD_IO_ERROR
        if (nf90_get_att(ncid, NF90_GLOBAL, 'format', fmt) /= NF90_NOERR) return
        if (fmt /= PERTFIELD_FORMAT) return
        st = nf90_get_att(ncid, NF90_GLOBAL, 'format_version', version)
        if (st /= NF90_NOERR) return
        if (version /= PERTFIELD_FORMAT_VERSION) return
        if (nf90_get_att(ncid, NF90_GLOBAL, 'units', units) /= NF90_NOERR) return
        if (units /= UNITS_SI .and. units /= UNITS_GAUSSIAN) return
        st = nf90_get_att(ncid, NF90_GLOBAL, 'has_potential', has_phi)
        if (st /= NF90_NOERR) return
        if (dim_len(ncid, 'R', nR) /= NF90_NOERR) return
        if (dim_len(ncid, 'Z', nZ) /= NF90_NOERR) return
        if (dim_len(ncid, 'mode', nm) /= NF90_NOERR) return

        allocate (R(nR), Z(nZ), ntor(nm), dA(3, nR, nZ, nm), dPhi(nR, nZ, nm))
        st = get_r1(ncid, 'R', R)
        if (st == NF90_NOERR) st = get_r1(ncid, 'Z', Z)
        if (st == NF90_NOERR) st = nf90_inq_varid(ncid, 'ntor', vid)
        if (st == NF90_NOERR) st = nf90_get_var(ncid, vid, ntor)
        do c = 1, 3
            if (st == NF90_NOERR) st = get_complex(ncid, trim(COMP(c)), dA(c, :, :, :))
        end do
        if (has_phi == 1) then
            if (st == NF90_NOERR) st = get_complex(ncid, 'dPhi', dPhi)
        end if
        if (st /= NF90_NOERR) return

        if (has_phi == 1) then
            call field%init(R, Z, ntor, dA, trim(units), ierr, dPhi=dPhi)
        else
            call field%init(R, Z, ntor, dA, trim(units), ierr)
        end if
    end subroutine read_open_file

    integer function dim_len(ncid, name, n) result(st)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: name
        integer, intent(out) :: n

        integer :: dimid

        n = 0
        st = nf90_inq_dimid(ncid, name, dimid)
        if (st == NF90_NOERR) st = nf90_inquire_dimension(ncid, dimid, len=n)
    end function dim_len

    integer function get_r1(ncid, name, x) result(st)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: name
        real(dp), intent(out) :: x(:)

        integer :: vid

        st = nf90_inq_varid(ncid, name, vid)
        if (st == NF90_NOERR) st = nf90_get_var(ncid, vid, x)
    end function get_r1

    integer function get_complex(ncid, name, a) result(st)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: name
        complex(dp), intent(out) :: a(:, :, :)

        real(dp) :: re(size(a, 1), size(a, 2), size(a, 3))
        real(dp) :: im(size(a, 1), size(a, 2), size(a, 3))
        integer :: vid

        st = nf90_inq_varid(ncid, name//'_re', vid)
        if (st == NF90_NOERR) st = nf90_get_var(ncid, vid, re)
        if (st == NF90_NOERR) st = nf90_inq_varid(ncid, name//'_im', vid)
        if (st == NF90_NOERR) st = nf90_get_var(ncid, vid, im)
        a = cmplx(re, im, dp)
    end function get_complex

end module neo_perturbation_field_netcdf
