module neo_field_netcdf
! On-disk netCDF format for a field mesh. Stores the vector potential A on a
! regular grid; B is optional and written only for interchange with codes that
! expect it. The grid is cylindrical (R, phi, Z) or Cartesian (x, y, z),
! selected by the coordinate_system global attribute. The splined field
! reconstructs B = curl A from the A-spline, so div B = 0 holds independent of
! grid resolution; stored B is for external compatibility only.
use netcdf
use neo_field_mesh, only: field_mesh_t
use neo_spline_field, only: COORD_CARTESIAN, COORD_CYLINDRICAL
implicit none
integer, parameter :: dp = kind(1.0d0)

private
public :: write_field_mesh_netcdf, read_field_mesh_netcdf

contains

subroutine write_field_mesh_netcdf(filename, field_mesh, coordinate_system, nfp, write_b)
    character(*), intent(in) :: filename
    type(field_mesh_t), intent(in) :: field_mesh
    integer, intent(in) :: coordinate_system
    integer, intent(in), optional :: nfp
    logical, intent(in), optional :: write_b

    character(len=8) :: axis(3), acomp(3), bcomp(3)
    integer :: ncid, d1, d2, d3, ddim, n1, n2, n3, nfp_loc, i
    integer :: vx(3), vper, va(3), vb(3)
    logical :: do_b

    nfp_loc = 1
    if (present(nfp)) nfp_loc = nfp
    do_b = .false.
    if (present(write_b)) do_b = write_b
    call component_names(coordinate_system, axis, acomp, bcomp)

    n1 = field_mesh%A1%n1
    n2 = field_mesh%A1%n2
    n3 = field_mesh%A1%n3

    call nc_check('create', nf90_create(filename, NF90_NETCDF4, ncid))
    if (coordinate_system == COORD_CYLINDRICAL) then
        call nc_check('att', nf90_put_att(ncid, NF90_GLOBAL, 'coordinate_system', 'cylindrical'))
    else
        call nc_check('att', nf90_put_att(ncid, NF90_GLOBAL, 'coordinate_system', 'cartesian'))
    end if
    call nc_check('att', nf90_put_att(ncid, NF90_GLOBAL, 'nfp', nfp_loc))

    call nc_check('dim', nf90_def_dim(ncid, 'n1', n1, d1))
    call nc_check('dim', nf90_def_dim(ncid, 'n2', n2, d2))
    call nc_check('dim', nf90_def_dim(ncid, 'n3', n3, d3))
    call nc_check('dim', nf90_def_dim(ncid, 'three', 3, ddim))

    call nc_check('var', nf90_def_var(ncid, trim(axis(1)), NF90_DOUBLE, [d1], vx(1)))
    call nc_check('var', nf90_def_var(ncid, trim(axis(2)), NF90_DOUBLE, [d2], vx(2)))
    call nc_check('var', nf90_def_var(ncid, trim(axis(3)), NF90_DOUBLE, [d3], vx(3)))
    call nc_check('var', nf90_def_var(ncid, 'is_periodic', NF90_INT, [ddim], vper))
    do i = 1, 3
        call nc_check('var', nf90_def_var(ncid, trim(acomp(i)), NF90_DOUBLE, [d1, d2, d3], va(i)))
    end do
    if (do_b) then
        do i = 1, 3
            call nc_check('var', nf90_def_var(ncid, trim(bcomp(i)), NF90_DOUBLE, [d1, d2, d3], vb(i)))
        end do
    end if
    call nc_check('enddef', nf90_enddef(ncid))

    call nc_check('put', nf90_put_var(ncid, vx(1), field_mesh%A1%x1))
    call nc_check('put', nf90_put_var(ncid, vx(2), field_mesh%A1%x2))
    call nc_check('put', nf90_put_var(ncid, vx(3), field_mesh%A1%x3))
    call nc_check('put', nf90_put_var(ncid, vper, logical_to_int(field_mesh%A1%is_periodic)))
    call nc_check('put', nf90_put_var(ncid, va(1), field_mesh%A1%value))
    call nc_check('put', nf90_put_var(ncid, va(2), field_mesh%A2%value))
    call nc_check('put', nf90_put_var(ncid, va(3), field_mesh%A3%value))
    if (do_b) then
        call nc_check('put', nf90_put_var(ncid, vb(1), field_mesh%B1%value))
        call nc_check('put', nf90_put_var(ncid, vb(2), field_mesh%B2%value))
        call nc_check('put', nf90_put_var(ncid, vb(3), field_mesh%B3%value))
    end if
    call nc_check('close', nf90_close(ncid))
end subroutine write_field_mesh_netcdf

subroutine read_field_mesh_netcdf(filename, field_mesh, coordinate_system, has_b)
    character(*), intent(in) :: filename
    type(field_mesh_t), intent(out) :: field_mesh
    integer, intent(out) :: coordinate_system
    logical, intent(out) :: has_b

    character(len=8) :: axis(3), acomp(3), bcomp(3)
    character(len=32) :: cs_name
    integer :: ncid, d1, d2, d3, n1, n2, n3, vid, iper(3), i
    real(dp), allocatable :: x1(:), x2(:), x3(:), v1(:,:,:), v2(:,:,:), v3(:,:,:)
    logical :: is_periodic(3)

    call nc_check('open', nf90_open(filename, NF90_NOWRITE, ncid))
    cs_name = ''
    call nc_check('att', nf90_get_att(ncid, NF90_GLOBAL, 'coordinate_system', cs_name))
    coordinate_system = coordinate_system_id(cs_name)
    call component_names(coordinate_system, axis, acomp, bcomp)

    call nc_check('dim', nf90_inq_dimid(ncid, 'n1', d1))
    call nc_check('dim', nf90_inq_dimid(ncid, 'n2', d2))
    call nc_check('dim', nf90_inq_dimid(ncid, 'n3', d3))
    call nc_check('dim', nf90_inquire_dimension(ncid, d1, len=n1))
    call nc_check('dim', nf90_inquire_dimension(ncid, d2, len=n2))
    call nc_check('dim', nf90_inquire_dimension(ncid, d3, len=n3))

    allocate(x1(n1), x2(n2), x3(n3), v1(n1,n2,n3), v2(n1,n2,n3), v3(n1,n2,n3))
    call nc_check('var', nf90_inq_varid(ncid, trim(axis(1)), vid))
    call nc_check('get', nf90_get_var(ncid, vid, x1))
    call nc_check('var', nf90_inq_varid(ncid, trim(axis(2)), vid))
    call nc_check('get', nf90_get_var(ncid, vid, x2))
    call nc_check('var', nf90_inq_varid(ncid, trim(axis(3)), vid))
    call nc_check('get', nf90_get_var(ncid, vid, x3))
    call nc_check('var', nf90_inq_varid(ncid, 'is_periodic', vid))
    call nc_check('get', nf90_get_var(ncid, vid, iper))
    is_periodic = iper /= 0

    call nc_check('var', nf90_inq_varid(ncid, trim(acomp(1)), vid))
    call nc_check('get', nf90_get_var(ncid, vid, v1))
    call nc_check('var', nf90_inq_varid(ncid, trim(acomp(2)), vid))
    call nc_check('get', nf90_get_var(ncid, vid, v2))
    call nc_check('var', nf90_inq_varid(ncid, trim(acomp(3)), vid))
    call nc_check('get', nf90_get_var(ncid, vid, v3))

    call field_mesh%A1%mesh_init(x1, x2, x3, v1, is_periodic)
    call field_mesh%A2%mesh_init(x1, x2, x3, v2, is_periodic)
    call field_mesh%A3%mesh_init(x1, x2, x3, v3, is_periodic)

    has_b = nf90_inq_varid(ncid, trim(bcomp(1)), vid) == NF90_NOERR
    if (has_b) then
        do i = 1, 3
            call nc_check('var', nf90_inq_varid(ncid, trim(bcomp(i)), vid))
            if (i == 1) call nc_check('get', nf90_get_var(ncid, vid, v1))
            if (i == 2) call nc_check('get', nf90_get_var(ncid, vid, v2))
            if (i == 3) call nc_check('get', nf90_get_var(ncid, vid, v3))
        end do
    else
        v1 = 0.0_dp; v2 = 0.0_dp; v3 = 0.0_dp
    end if
    call field_mesh%B1%mesh_init(x1, x2, x3, v1, is_periodic)
    call field_mesh%B2%mesh_init(x1, x2, x3, v2, is_periodic)
    call field_mesh%B3%mesh_init(x1, x2, x3, v3, is_periodic)

    call nc_check('close', nf90_close(ncid))
end subroutine read_field_mesh_netcdf

subroutine component_names(coordinate_system, axis, acomp, bcomp)
    integer, intent(in) :: coordinate_system
    character(len=8), intent(out) :: axis(3), acomp(3), bcomp(3)

    if (coordinate_system == COORD_CYLINDRICAL) then
        axis = ['R   ', 'phi ', 'Z   ']
        acomp = ['A_R  ', 'A_phi', 'A_Z  ']
        bcomp = ['B_R  ', 'B_phi', 'B_Z  ']
    else
        axis = ['x   ', 'y   ', 'z   ']
        acomp = ['A_x ', 'A_y ', 'A_z ']
        bcomp = ['B_x ', 'B_y ', 'B_z ']
    end if
end subroutine component_names

function coordinate_system_id(name) result(id)
    character(*), intent(in) :: name
    integer :: id

    if (trim(name) == 'cylindrical') then
        id = COORD_CYLINDRICAL
    else if (trim(name) == 'cartesian') then
        id = COORD_CARTESIAN
    else
        print *, 'unknown coordinate_system: ', trim(name)
        error stop
    end if
end function coordinate_system_id

function logical_to_int(flags) result(ints)
    logical, intent(in) :: flags(3)
    integer :: ints(3), i

    do i = 1, 3
        if (flags(i)) then
            ints(i) = 1
        else
            ints(i) = 0
        end if
    end do
end function logical_to_int

subroutine nc_check(context, status)
    character(*), intent(in) :: context
    integer, intent(in) :: status

    if (status /= NF90_NOERR) then
        print *, 'netCDF error (', context, '): ', trim(nf90_strerror(status))
        error stop
    end if
end subroutine nc_check

end module neo_field_netcdf
