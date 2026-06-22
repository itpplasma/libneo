module neo_mgrid
! makegrid (STELLOPT/VMEC) mgrid netCDF format and its generalization.
!
! Standard mgrid stores the vacuum field as physical cylindrical components
! br, bp, bz per coil group on a regular grid covering ONE field period in phi
! (kp planes at phi = (k-1)*2*pi/(nfp*kp); the endpoint 2*pi/nfp is NOT stored).
! read_mgrid sums the coil groups (scaling by raw_coil_cur in 'R'aw mode,
! summing directly in 'S'caled mode), appends the wrap plane at 2*pi/nfp in
! memory so the phi spline is periodic, and returns a field_mesh on (R, phi, Z).
!
! The same schema is reused for the generalized form: optional vector-potential
! groups ar, ap, az. These are makegrid's own potential variables and use the
! same convention as the field: physical components in the orthonormal (R, phi,
! Z) basis (a_p = A . phi_hat, NOT the covariant R*A_phi). Storing only A is the
! reduced form; the splined field then reconstructs B = curl A with the
! physical-component cylindrical curl, divergence-free by construction (it is a
! sin to store a H(div) field independently). A cartesian variant
! (coordinate_system = "cartesian", components ax/ay/az, bx/by/bz on a linear
! x/y/z grid with an extra ymin/ymax) shares all the I/O below.
use, intrinsic :: iso_fortran_env, only: error_unit
use netcdf
use neo_mesh, only: mesh_t
use neo_field_mesh, only: field_mesh_t
use neo_spline_field, only: COORD_CARTESIAN, COORD_CYLINDRICAL
implicit none
integer, parameter :: dp = kind(1.0d0)
real(dp), parameter :: twopi = 8.0_dp*atan(1.0_dp)

private
public :: read_mgrid, write_mgrid, is_mgrid_file

contains

subroutine read_mgrid(filename, field_mesh, nfp, coordinate_system, has_a, has_b, verbose)
    character(*), intent(in) :: filename
    type(field_mesh_t), intent(out) :: field_mesh
    integer, intent(out) :: nfp
    integer, intent(out), optional :: coordinate_system
    logical, intent(out), optional :: has_a, has_b
    logical, intent(in), optional :: verbose

    integer :: ncid, vid, ir, jz, kp, nextcur, cs, np2, k
    real(dp) :: rmin, rmax, zmin, zmax, ymin, ymax
    character(len=1) :: mode_c
    character(len=2) :: aname(3), bname(3)
    logical :: per(3), got_a, got_b, is_raw, do_report
    real(dp), allocatable :: cur(:)
    real(dp), allocatable :: ax(:,:,:), ay(:,:,:), az(:,:,:)
    real(dp), allocatable :: bx(:,:,:), by(:,:,:), bz(:,:,:)
    real(dp), allocatable :: x1(:), x2(:), x3(:)
    real(dp), allocatable :: A1(:,:,:), A2(:,:,:), A3(:,:,:)
    real(dp), allocatable :: B1(:,:,:), B2(:,:,:), B3(:,:,:)

    call nc_check('open', nf90_open(filename, NF90_NOWRITE, ncid))
    cs = read_coordinate_system(ncid)
    call comp_names(cs, aname, bname)
    call dim_len(ncid, 'rad', ir)
    call dim_len(ncid, 'zee', jz)
    call dim_len(ncid, 'phi', kp)
    call get_int(ncid, 'nfp', nfp)
    call get_int(ncid, 'nextcur', nextcur)
    call get_dbl(ncid, 'rmin', rmin); call get_dbl(ncid, 'rmax', rmax)
    call get_dbl(ncid, 'zmin', zmin); call get_dbl(ncid, 'zmax', zmax)
    call nc_check('mode', nf90_inq_varid(ncid, 'mgrid_mode', vid))
    call nc_check('mode', nf90_get_var(ncid, vid, mode_c))
    call warn_unknown_mode(mode_c)
    is_raw = mode_is_raw(mode_c)
    allocate(cur(nextcur))
    call nc_check('cur', nf90_inq_varid(ncid, 'raw_coil_cur', vid))
    call nc_check('cur', nf90_get_var(ncid, vid, cur))

    got_a = group_present(ncid, aname(1))
    got_b = group_present(ncid, bname(1))

    allocate(ax(ir,jz,kp), ay(ir,jz,kp), az(ir,jz,kp))
    allocate(bx(ir,jz,kp), by(ir,jz,kp), bz(ir,jz,kp))
    ax = 0; ay = 0; az = 0; bx = 0; by = 0; bz = 0
    if (got_a) call sum_groups(ncid, aname, nextcur, is_raw, cur, ax, ay, az)
    if (got_b) call sum_groups(ncid, bname, nextcur, is_raw, cur, bx, by, bz)

    do_report = .false.
    if (present(verbose)) do_report = verbose
    if (do_report .and. got_b) call report_scaling(ncid, bname, nextcur, mode_c, cur, bx, by, bz)

    if (cs == COORD_CYLINDRICAL) then
        np2 = kp + 1                ! append wrap plane for the periodic phi spline
        allocate(x2(np2))
        do k = 1, np2
            x2(k) = real(k-1, dp)*twopi/(real(nfp, dp)*real(kp, dp))
        end do
        per = [.false., .true., .false.]
    else
        np2 = kp
        call get_dbl(ncid, 'ymin', ymin); call get_dbl(ncid, 'ymax', ymax)
        allocate(x2(np2))
        x2 = linspace(ymin, ymax, np2)
        per = [.false., .false., .false.]
    end if
    allocate(x1(ir), x3(jz))
    x1 = linspace(rmin, rmax, ir)
    x3 = linspace(zmin, zmax, jz)
    call nc_check('close', nf90_close(ncid))

    allocate(A1(ir,np2,jz), A2(ir,np2,jz), A3(ir,np2,jz))
    allocate(B1(ir,np2,jz), B2(ir,np2,jz), B3(ir,np2,jz))
    call to_mesh(ax, A1, kp, cs); call to_mesh(ay, A2, kp, cs); call to_mesh(az, A3, kp, cs)
    call to_mesh(bx, B1, kp, cs); call to_mesh(by, B2, kp, cs); call to_mesh(bz, B3, kp, cs)

    call field_mesh%A1%mesh_init(x1, x2, x3, A1, per)
    call field_mesh%A2%mesh_init(x1, x2, x3, A2, per)
    call field_mesh%A3%mesh_init(x1, x2, x3, A3, per)
    call field_mesh%B1%mesh_init(x1, x2, x3, B1, per)
    call field_mesh%B2%mesh_init(x1, x2, x3, B2, per)
    call field_mesh%B3%mesh_init(x1, x2, x3, B3, per)

    if (present(coordinate_system)) coordinate_system = cs
    if (present(has_a)) has_a = got_a
    if (present(has_b)) has_b = got_b
end subroutine read_mgrid

subroutine write_mgrid(filename, field_mesh, nfp, coordinate_system, write_a, write_b, group_name)
    character(*), intent(in) :: filename
    type(field_mesh_t), intent(in) :: field_mesh
    integer, intent(in) :: nfp
    integer, intent(in), optional :: coordinate_system
    logical, intent(in), optional :: write_a, write_b
    character(*), intent(in), optional :: group_name

    integer, parameter :: ss = 30
    integer :: ncid, ir, jz, kp, n2, cs, v
    integer :: d_ss, d_grp, d_one, d_cur, d_r, d_z, d_p
    logical :: do_a, do_b
    real(dp) :: rmin, rmax, zmin, zmax, period, span
    character(len=2) :: aname(3), bname(3)
    character(len=ss) :: gnames(1)
    character(len=1) :: mode_c

    cs = COORD_CYLINDRICAL
    if (present(coordinate_system)) cs = coordinate_system
    do_b = .true.; if (present(write_b)) do_b = write_b
    do_a = .false.; if (present(write_a)) do_a = write_a
    call comp_names(cs, aname, bname)

    ir = field_mesh%B1%n1
    n2 = field_mesh%B1%n2
    jz = field_mesh%B1%n3
    if (cs == COORD_CYLINDRICAL) then
        period = twopi/real(nfp, dp)
        span = field_mesh%B1%x2(n2) - field_mesh%B1%x2(1)
        if (abs(span - period) < 1.0e-6_dp*period) then
            kp = n2 - 1            ! drop the in-memory wrap plane; mgrid stores one period
        else
            kp = n2
        end if
    else
        kp = n2
    end if
    rmin = field_mesh%B1%x1(1); rmax = field_mesh%B1%x1(ir)
    zmin = field_mesh%B1%x3(1); zmax = field_mesh%B1%x3(jz)

    gnames(1) = 'field'
    if (present(group_name)) gnames(1) = group_name

    call nc_check('create', nf90_create(filename, NF90_CLOBBER, ncid))
    if (cs == COORD_CARTESIAN) call nc_check('att', &
        nf90_put_att(ncid, NF90_GLOBAL, 'coordinate_system', 'cartesian'))
    call nc_check('dim', nf90_def_dim(ncid, 'stringsize', ss, d_ss))
    call nc_check('dim', nf90_def_dim(ncid, 'external_coil_groups', 1, d_grp))
    call nc_check('dim', nf90_def_dim(ncid, 'dim_00001', 1, d_one))
    call nc_check('dim', nf90_def_dim(ncid, 'external_coils', 1, d_cur))
    call nc_check('dim', nf90_def_dim(ncid, 'rad', ir, d_r))
    call nc_check('dim', nf90_def_dim(ncid, 'zee', jz, d_z))
    call nc_check('dim', nf90_def_dim(ncid, 'phi', kp, d_p))

    call def_int(ncid, 'ir'); call def_int(ncid, 'jz'); call def_int(ncid, 'kp')
    call def_int(ncid, 'nfp'); call def_int(ncid, 'nextcur')
    call def_dbl(ncid, 'rmin'); call def_dbl(ncid, 'rmax')
    call def_dbl(ncid, 'zmin'); call def_dbl(ncid, 'zmax')
    if (cs == COORD_CARTESIAN) then
        call def_dbl(ncid, 'ymin'); call def_dbl(ncid, 'ymax')
    end if
    call nc_check('def', nf90_def_var(ncid, 'coil_group', NF90_CHAR, [d_ss, d_grp], v))
    call nc_check('def', nf90_def_var(ncid, 'mgrid_mode', NF90_CHAR, [d_one], v))
    call nc_check('def', nf90_def_var(ncid, 'raw_coil_cur', NF90_DOUBLE, [d_cur], v))
    if (do_b) call def_group(ncid, bname, d_r, d_z, d_p)
    if (do_a) call def_group(ncid, aname, d_r, d_z, d_p)
    call nc_check('enddef', nf90_enddef(ncid))

    call put_int(ncid, 'ir', ir); call put_int(ncid, 'jz', jz); call put_int(ncid, 'kp', kp)
    call put_int(ncid, 'nfp', nfp); call put_int(ncid, 'nextcur', 1)
    call put_dbl(ncid, 'rmin', rmin); call put_dbl(ncid, 'rmax', rmax)
    call put_dbl(ncid, 'zmin', zmin); call put_dbl(ncid, 'zmax', zmax)
    if (cs == COORD_CARTESIAN) then
        call put_dbl(ncid, 'ymin', field_mesh%B1%x2(1))
        call put_dbl(ncid, 'ymax', field_mesh%B1%x2(n2))
    end if
    mode_c = 'S'
    call nc_check('v', nf90_inq_varid(ncid, 'coil_group', v))
    call nc_check('p', nf90_put_var(ncid, v, gnames))
    call nc_check('v', nf90_inq_varid(ncid, 'mgrid_mode', v))
    call nc_check('p', nf90_put_var(ncid, v, mode_c))
    call nc_check('v', nf90_inq_varid(ncid, 'raw_coil_cur', v))
    call nc_check('p', nf90_put_var(ncid, v, [1.0_dp]))
    if (do_b) call put_group(ncid, bname, field_mesh%B1, field_mesh%B2, field_mesh%B3, kp)
    if (do_a) call put_group(ncid, aname, field_mesh%A1, field_mesh%A2, field_mesh%A3, kp)
    call nc_check('close', nf90_close(ncid))
end subroutine write_mgrid

function is_mgrid_file(filename) result(is_mgrid)
    character(*), intent(in) :: filename
    logical :: is_mgrid
    integer :: ncid, vid, stat

    is_mgrid = .false.
    if (nf90_open(filename, NF90_NOWRITE, ncid) /= NF90_NOERR) return
    is_mgrid = nf90_inq_varid(ncid, 'coil_group', vid) == NF90_NOERR
    stat = nf90_close(ncid)
end function is_mgrid_file

! Permute a stored group (ir, jz, kp) into mesh order (ir, np2, jz), appending
! the wrap plane in the cylindrical case so np2 = kp+1.
subroutine to_mesh(g, v, kp, cs)
    real(dp), intent(in) :: g(:,:,:)
    real(dp), intent(out) :: v(:,:,:)
    integer, intent(in) :: kp, cs
    integer :: i, j, k, ir, jz

    ir = size(g, 1); jz = size(g, 2)
    do j = 1, jz
        do k = 1, kp
            do i = 1, ir
                v(i,k,j) = g(i,j,k)
            end do
        end do
    end do
    if (cs == COORD_CYLINDRICAL) then
        do j = 1, jz
            v(:,kp+1,j) = v(:,1,j)
        end do
    end if
end subroutine to_mesh

subroutine sum_groups(ncid, base, nextcur, apply_current, cur, c1, c2, c3)
    integer, intent(in) :: ncid, nextcur
    character(len=2), intent(in) :: base(3)
    logical, intent(in) :: apply_current
    real(dp), intent(in) :: cur(:)
    real(dp), intent(inout) :: c1(:,:,:), c2(:,:,:), c3(:,:,:)

    integer :: g, vid
    real(dp) :: scale
    real(dp), allocatable :: gb(:,:,:)

    allocate(gb(size(c1,1), size(c1,2), size(c1,3)))
    do g = 1, nextcur
        scale = 1.0_dp
        if (apply_current) scale = cur(g)
        call nc_check('g', nf90_inq_varid(ncid, group_var_name(base(1), g), vid))
        call nc_check('g', nf90_get_var(ncid, vid, gb)); c1 = c1 + scale*gb
        call nc_check('g', nf90_inq_varid(ncid, group_var_name(base(2), g), vid))
        call nc_check('g', nf90_get_var(ncid, vid, gb)); c2 = c2 + scale*gb
        call nc_check('g', nf90_inq_varid(ncid, group_var_name(base(3), g), vid))
        call nc_check('g', nf90_get_var(ncid, vid, gb)); c3 = c3 + scale*gb
    end do
end subroutine sum_groups

subroutine def_group(ncid, base, d_r, d_z, d_p)
    integer, intent(in) :: ncid, d_r, d_z, d_p
    character(len=2), intent(in) :: base(3)
    integer :: c, vid

    do c = 1, 3
        call nc_check('def', nf90_def_var(ncid, group_var_name(base(c), 1), &
            NF90_DOUBLE, [d_r, d_z, d_p], vid))
    end do
end subroutine def_group

! Permute the first kp phi planes of a mesh component (ir, n2, jz) into stored
! group order (ir, jz, kp) and write the single coil group.
subroutine put_group(ncid, base, m1, m2, m3, kp)
    integer, intent(in) :: ncid, kp
    character(len=2), intent(in) :: base(3)
    type(mesh_t), intent(in) :: m1, m2, m3

    integer :: ir, jz, i, j, k, vid
    real(dp), allocatable :: g1(:,:,:), g2(:,:,:), g3(:,:,:)

    ir = m1%n1; jz = m1%n3
    allocate(g1(ir,jz,kp), g2(ir,jz,kp), g3(ir,jz,kp))
    do k = 1, kp
        do j = 1, jz
            do i = 1, ir
                g1(i,j,k) = m1%value(i,k,j)
                g2(i,j,k) = m2%value(i,k,j)
                g3(i,j,k) = m3%value(i,k,j)
            end do
        end do
    end do
    call nc_check('v', nf90_inq_varid(ncid, group_var_name(base(1), 1), vid))
    call nc_check('p', nf90_put_var(ncid, vid, g1))
    call nc_check('v', nf90_inq_varid(ncid, group_var_name(base(2), 1), vid))
    call nc_check('p', nf90_put_var(ncid, vid, g2))
    call nc_check('v', nf90_inq_varid(ncid, group_var_name(base(3), 1), vid))
    call nc_check('p', nf90_put_var(ncid, vid, g3))
end subroutine put_group

function mode_is_raw(mode_c) result(raw)
    character(len=1), intent(in) :: mode_c
    logical :: raw

    raw = mode_c == 'R' .or. mode_c == 'r'
end function mode_is_raw

! Standard makegrid tags: 'R' raw (scale by raw_coil_cur), 'S' scaled (sum
! directly). 'N' (no scaling) is read as direct-sum. Anything else is treated
! as direct-sum and warned about: the convention cannot be inferred from the
! field (it is invariant under a global current scale), so the tag is trusted.
subroutine warn_unknown_mode(mode_c)
    character(len=1), intent(in) :: mode_c

    select case (mode_c)
    case ('R', 'r', 'S', 's', 'N', 'n')
    case default
        write(error_unit, '(3a)') 'neo_mgrid WARNING: unrecognized mgrid_mode "', &
            mode_c, '"; treating as already scaled (direct group sum).'
    end select
end subroutine warn_unknown_mode

! Diagnostic only: report peak |B| as read versus under the opposite scaling, so
! a mislabelled mode is visible. Never changes the field that is returned.
subroutine report_scaling(ncid, base, nextcur, mode_c, cur, c1, c2, c3)
    integer, intent(in) :: ncid, nextcur
    character(len=2), intent(in) :: base(3)
    character(len=1), intent(in) :: mode_c
    real(dp), intent(in) :: cur(:)
    real(dp), intent(in) :: c1(:,:,:), c2(:,:,:), c3(:,:,:)

    real(dp), allocatable :: a1(:,:,:), a2(:,:,:), a3(:,:,:)

    allocate(a1(size(c1,1), size(c1,2), size(c1,3)), &
             a2(size(c1,1), size(c1,2), size(c1,3)), &
             a3(size(c1,1), size(c1,2), size(c1,3)))
    a1 = 0; a2 = 0; a3 = 0
    call sum_groups(ncid, base, nextcur, .not. mode_is_raw(mode_c), cur, a1, a2, a3)

    write(error_unit, '(3a)') &
        'neo_mgrid scaling diagnostic (mgrid_mode = "', mode_c, '"):'
    write(error_unit, '(a,es12.4)') &
        '  max|B| as read                = ', field_norm_max(c1, c2, c3)
    write(error_unit, '(a,es12.4)') &
        '  max|B| under opposite scaling = ', field_norm_max(a1, a2, a3)
    write(error_unit, '(a)') &
        '  verify the intended mode (raw scales by raw_coil_cur, scaled does not).'
end subroutine report_scaling

function field_norm_max(c1, c2, c3) result(m)
    real(dp), intent(in) :: c1(:,:,:), c2(:,:,:), c3(:,:,:)
    real(dp) :: m

    m = sqrt(maxval(c1**2 + c2**2 + c3**2))
end function field_norm_max

subroutine comp_names(cs, aname, bname)
    integer, intent(in) :: cs
    character(len=2), intent(out) :: aname(3), bname(3)

    if (cs == COORD_CARTESIAN) then
        aname = ['ax', 'ay', 'az']; bname = ['bx', 'by', 'bz']
    else
        aname = ['ar', 'ap', 'az']; bname = ['br', 'bp', 'bz']
    end if
end subroutine comp_names

function group_var_name(base, g) result(name)
    character(len=2), intent(in) :: base
    integer, intent(in) :: g
    character(len=6) :: name

    write(name, '(a2,a1,i3.3)') base, '_', g
end function group_var_name

function group_present(ncid, base) result(found)
    integer, intent(in) :: ncid
    character(len=2), intent(in) :: base
    logical :: found
    integer :: vid

    found = nf90_inq_varid(ncid, group_var_name(base, 1), vid) == NF90_NOERR
end function group_present

function read_coordinate_system(ncid) result(cs)
    integer, intent(in) :: ncid
    integer :: cs
    character(len=32) :: name

    name = ''
    if (nf90_get_att(ncid, NF90_GLOBAL, 'coordinate_system', name) == NF90_NOERR) then
        if (trim(name) == 'cartesian') then
            cs = COORD_CARTESIAN
            return
        end if
    end if
    cs = COORD_CYLINDRICAL
end function read_coordinate_system

subroutine dim_len(ncid, name, n)
    integer, intent(in) :: ncid
    character(*), intent(in) :: name
    integer, intent(out) :: n
    integer :: did

    call nc_check('dimid', nf90_inq_dimid(ncid, name, did))
    call nc_check('dimlen', nf90_inquire_dimension(ncid, did, len=n))
end subroutine dim_len

subroutine get_int(ncid, name, val)
    integer, intent(in) :: ncid
    character(*), intent(in) :: name
    integer, intent(out) :: val
    integer :: vid

    call nc_check('iget', nf90_inq_varid(ncid, name, vid))
    call nc_check('iget', nf90_get_var(ncid, vid, val))
end subroutine get_int

subroutine get_dbl(ncid, name, val)
    integer, intent(in) :: ncid
    character(*), intent(in) :: name
    real(dp), intent(out) :: val
    integer :: vid

    call nc_check('dget', nf90_inq_varid(ncid, name, vid))
    call nc_check('dget', nf90_get_var(ncid, vid, val))
end subroutine get_dbl

subroutine def_int(ncid, name)
    integer, intent(in) :: ncid
    character(*), intent(in) :: name
    integer :: vid

    call nc_check('idef', nf90_def_var(ncid, name, NF90_INT, vid))
end subroutine def_int

subroutine def_dbl(ncid, name)
    integer, intent(in) :: ncid
    character(*), intent(in) :: name
    integer :: vid

    call nc_check('ddef', nf90_def_var(ncid, name, NF90_DOUBLE, vid))
end subroutine def_dbl

subroutine put_int(ncid, name, val)
    integer, intent(in) :: ncid
    character(*), intent(in) :: name
    integer, intent(in) :: val
    integer :: vid

    call nc_check('iput', nf90_inq_varid(ncid, name, vid))
    call nc_check('iput', nf90_put_var(ncid, vid, val))
end subroutine put_int

subroutine put_dbl(ncid, name, val)
    integer, intent(in) :: ncid
    character(*), intent(in) :: name
    real(dp), intent(in) :: val
    integer :: vid

    call nc_check('dput', nf90_inq_varid(ncid, name, vid))
    call nc_check('dput', nf90_put_var(ncid, vid, val))
end subroutine put_dbl

function linspace(start, stop, n) result(x)
    real(dp), intent(in) :: start, stop
    integer, intent(in) :: n
    real(dp) :: x(n)
    integer :: i

    do i = 1, n
        x(i) = start + (stop - start)*real(i-1, dp)/real(n-1, dp)
    end do
end function linspace

subroutine nc_check(context, status)
    character(*), intent(in) :: context
    integer, intent(in) :: status

    if (status /= NF90_NOERR) then
        print *, 'mgrid netCDF error (', context, '): ', trim(nf90_strerror(status))
        error stop
    end if
end subroutine nc_check

end module neo_mgrid
