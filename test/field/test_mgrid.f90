program test_mgrid
use util_for_test, only: print_test, print_ok, print_fail
use netcdf
use neo_field_mesh, only: field_mesh_t
use neo_spline_field, only: COORD_CARTESIAN, COORD_CYLINDRICAL
use neo_mgrid, only: read_mgrid, write_mgrid
implicit none

integer, parameter :: dp = kind(1.0d0)
real(dp), parameter :: twopi = 8.0_dp*atan(1.0_dp)

call test_roundtrip_b_cyl
call test_raw_multigroup
call test_curl_cyl_factory
call test_curl_cart_factory
call test_curl_order5

contains

! Periodic-in-phi vacuum field for the B round-trip (period 2*pi/nfp).
subroutine analytic_b_cyl(x, nfp, B)
    real(dp), intent(in) :: x(3)
    integer, intent(in) :: nfp
    real(dp), intent(out) :: B(3)

    B(1) = cos(real(nfp, dp)*x(2))
    B(2) = x(1)
    B(3) = sin(real(nfp, dp)*x(2)) + x(3)
end subroutine analytic_b_cyl

! Cylindrical potential A = (0, R Z, 0), physical components, x = (R, phi, Z),
! so curl A = (-R, 0, 2 Z), divergence-free.
subroutine analytic_a_cyl(x, A, B)
    real(dp), intent(in) :: x(3)
    real(dp), intent(out) :: A(3), B(3)

    A = [0.0_dp, x(1)*x(3), 0.0_dp]
    B = [-x(1), 0.0_dp, 2.0_dp*x(3)]
end subroutine analytic_a_cyl

! Cartesian potential A = (y^2, z^2, x^2), curl A = (-2z, -2x, -2y).
subroutine analytic_a_cart(x, A, B)
    real(dp), intent(in) :: x(3)
    real(dp), intent(out) :: A(3), B(3)

    A = [x(2)**2, x(3)**2, x(1)**2]
    B = [-2.0_dp*x(3), -2.0_dp*x(1), -2.0_dp*x(2)]
end subroutine analytic_a_cart

! field_mesh on (R, phi, Z), phi inclusive of the wrap endpoint at 2*pi/nfp
! (n2 = kp+1), as the periodic spline and write_mgrid expect.
subroutine build_cyl_mesh(field_mesh, nfp, ir, kp, jz, store_a)
    type(field_mesh_t), intent(out) :: field_mesh
    integer, intent(in) :: nfp, ir, kp, jz
    logical, intent(in) :: store_a

    real(dp), allocatable :: x1(:), x2(:), x3(:)
    real(dp), allocatable :: A1(:,:,:), A2(:,:,:), A3(:,:,:)
    real(dp), allocatable :: B1(:,:,:), B2(:,:,:), B3(:,:,:)
    real(dp) :: x(3), A(3), B(3)
    logical :: per(3)
    integer :: i, j, k

    allocate(x1(ir), x2(kp+1), x3(jz))
    x1 = linspace(0.8_dp, 1.6_dp, ir)
    x2 = linspace(0.0_dp, twopi/real(nfp, dp), kp+1)
    x3 = linspace(-0.5_dp, 0.5_dp, jz)
    per = [.false., .true., .false.]

    allocate(A1(ir,kp+1,jz), A2(ir,kp+1,jz), A3(ir,kp+1,jz))
    allocate(B1(ir,kp+1,jz), B2(ir,kp+1,jz), B3(ir,kp+1,jz))
    A1 = 0; A2 = 0; A3 = 0; B1 = 0; B2 = 0; B3 = 0
    do i = 1, ir
        do k = 1, kp+1
            do j = 1, jz
                x = [x1(i), x2(k), x3(j)]
                if (store_a) then
                    call analytic_a_cyl(x, A, B)
                    A1(i,k,j) = A(1); A2(i,k,j) = A(2); A3(i,k,j) = A(3)
                else
                    call analytic_b_cyl(x, nfp, B)
                    B1(i,k,j) = B(1); B2(i,k,j) = B(2); B3(i,k,j) = B(3)
                end if
            end do
        end do
    end do

    call field_mesh%A1%mesh_init(x1, x2, x3, A1, per)
    call field_mesh%A2%mesh_init(x1, x2, x3, A2, per)
    call field_mesh%A3%mesh_init(x1, x2, x3, A3, per)
    call field_mesh%B1%mesh_init(x1, x2, x3, B1, per)
    call field_mesh%B2%mesh_init(x1, x2, x3, B2, per)
    call field_mesh%B3%mesh_init(x1, x2, x3, B3, per)
end subroutine build_cyl_mesh

subroutine build_cart_a_mesh(field_mesh, n)
    type(field_mesh_t), intent(out) :: field_mesh
    integer, intent(in) :: n

    real(dp), allocatable :: x1(:), x2(:), x3(:)
    real(dp), allocatable :: A1(:,:,:), A2(:,:,:), A3(:,:,:), zero(:,:,:)
    real(dp) :: x(3), A(3), B(3)
    logical :: per(3)
    integer :: i, j, k

    allocate(x1(n), x2(n), x3(n))
    x1 = linspace(1.0_dp, 2.0_dp, n)
    x2 = linspace(1.0_dp, 2.0_dp, n)
    x3 = linspace(1.0_dp, 2.0_dp, n)
    per = [.false., .false., .false.]

    allocate(A1(n,n,n), A2(n,n,n), A3(n,n,n), zero(n,n,n))
    zero = 0.0_dp
    do i = 1, n
        do k = 1, n
            do j = 1, n
                x = [x1(i), x2(k), x3(j)]
                call analytic_a_cart(x, A, B)
                A1(i,k,j) = A(1); A2(i,k,j) = A(2); A3(i,k,j) = A(3)
            end do
        end do
    end do

    call field_mesh%A1%mesh_init(x1, x2, x3, A1, per)
    call field_mesh%A2%mesh_init(x1, x2, x3, A2, per)
    call field_mesh%A3%mesh_init(x1, x2, x3, A3, per)
    call field_mesh%B1%mesh_init(x1, x2, x3, zero, per)
    call field_mesh%B2%mesh_init(x1, x2, x3, zero, per)
    call field_mesh%B3%mesh_init(x1, x2, x3, zero, per)
end subroutine build_cart_a_mesh

function linspace(start, stop, n) result(x)
    real(dp), intent(in) :: start, stop
    integer, intent(in) :: n
    real(dp) :: x(n)
    integer :: i

    do i = 1, n
        x(i) = start + (stop - start)*real(i-1, dp)/real(n-1, dp)
    end do
end function linspace

subroutine test_roundtrip_b_cyl
    character(*), parameter :: fname = "test_mgrid_b.nc"
    integer, parameter :: nfp = 5, ir = 8, kp = 12, jz = 6
    type(field_mesh_t) :: mesh_out, mesh_in
    integer :: nfp_in, kp_disk
    real(dp) :: dev

    call print_test("test_roundtrip_b_cyl")
    call build_cyl_mesh(mesh_out, nfp, ir, kp, jz, store_a=.false.)
    call write_mgrid(fname, mesh_out, nfp)

    ! one field period is stored without the wrap plane
    kp_disk = disk_dim(fname, 'phi')
    call read_mgrid(fname, mesh_in, nfp_in)
    call unlink(fname)

    dev = maxval(abs(mesh_in%B1%value - mesh_out%B1%value)) &
        + maxval(abs(mesh_in%B2%value - mesh_out%B2%value)) &
        + maxval(abs(mesh_in%B3%value - mesh_out%B3%value)) &
        + maxval(abs(mesh_in%B1%x2 - mesh_out%B1%x2))

    if (kp_disk /= kp .or. nfp_in /= nfp .or. mesh_in%B1%n2 /= kp+1 &
        .or. dev > 1.0e-12_dp) then
        print *, "kp_disk=", kp_disk, " nfp_in=", nfp_in, &
            " n2=", mesh_in%B1%n2, " dev=", dev
        call print_fail
        error stop
    end if
    call print_ok
end subroutine test_roundtrip_b_cyl

! Hand-built raw-mode two-group mgrid: total B = sum_g cur(g)*field_g.
subroutine test_raw_multigroup
    character(*), parameter :: fname = "test_mgrid_raw.nc"
    integer, parameter :: nfp = 3, ir = 5, jz = 4, kp = 6
    real(dp), parameter :: c1 = 2.0_dp, c2 = -3.0_dp
    type(field_mesh_t) :: mesh
    integer :: nfp_in
    real(dp) :: dev

    call print_test("test_raw_multigroup")
    call build_raw_mgrid(fname, nfp, ir, jz, kp, c1, c2)
    call read_mgrid(fname, mesh, nfp_in)
    call unlink(fname)

    dev = maxval(abs(mesh%B1%value - (c1*1.0_dp + c2*2.0_dp))) &
        + maxval(abs(mesh%B2%value - (c1*10.0_dp + c2*20.0_dp))) &
        + maxval(abs(mesh%B3%value - (c1*100.0_dp + c2*200.0_dp)))

    if (nfp_in /= nfp .or. mesh%B1%n2 /= kp+1 .or. dev > 1.0e-12_dp &
        .or. abs(mesh%B1%x2(kp+1) - twopi/real(nfp, dp)) > 1.0e-12_dp) then
        print *, "nfp_in=", nfp_in, " n2=", mesh%B1%n2, " dev=", dev
        call print_fail
        error stop
    end if
    call print_ok
end subroutine test_raw_multigroup

! A-only cylindrical mgrid through the factory: B = curl A.
subroutine test_curl_cyl_factory
    use neo_field, only: create_field
    use neo_field_base, only: field_t

    character(*), parameter :: fname = "test_mgrid_acyl.nc"
    integer, parameter :: nfp = 4, ir = 40, kp = 12, jz = 40
    real(dp), parameter :: tol = 1.0e-7_dp
    type(field_mesh_t) :: mesh
    class(field_t), allocatable :: field
    real(dp) :: x(3), B(3), A_ref(3), B_ref(3)

    call print_test("test_curl_cyl_factory")
    call build_cyl_mesh(mesh, nfp, ir, kp, jz, store_a=.true.)
    call write_mgrid(fname, mesh, nfp, coordinate_system=COORD_CYLINDRICAL, &
        write_a=.true., write_b=.false.)
    call create_field(field, "spline", filename=fname)
    call unlink(fname)

    x = [1.13_dp, 0.37_dp*twopi/real(nfp, dp), 0.21_dp]
    call field%compute_bfield(x, B)
    call analytic_a_cyl(x, A_ref, B_ref)
    if (maxval(abs(B - B_ref)) > tol) then
        print *, "B=", B, " ref=", B_ref
        call print_fail
        error stop
    end if
    call print_ok
end subroutine test_curl_cyl_factory

! A-only cartesian mgrid through the factory: B = curl A.
subroutine test_curl_cart_factory
    use neo_field, only: create_field
    use neo_field_base, only: field_t

    character(*), parameter :: fname = "test_mgrid_acart.nc"
    integer, parameter :: n = 40
    real(dp), parameter :: tol = 1.0e-8_dp
    type(field_mesh_t) :: mesh
    class(field_t), allocatable :: field
    real(dp) :: x(3), B(3), A_ref(3), B_ref(3)

    call print_test("test_curl_cart_factory")
    call build_cart_a_mesh(mesh, n)
    call write_mgrid(fname, mesh, 1, coordinate_system=COORD_CARTESIAN, &
        write_a=.true., write_b=.false.)
    call create_field(field, "spline", filename=fname)
    call unlink(fname)

    x = [1.37_dp, 1.51_dp, 1.62_dp]
    call field%compute_bfield(x, B)
    call analytic_a_cart(x, A_ref, B_ref)
    if (maxval(abs(B - B_ref)) > tol) then
        print *, "B=", B, " ref=", B_ref
        call print_fail
        error stop
    end if
    call print_ok
end subroutine test_curl_cart_factory

! Both coordinate variants through the factory at quintic spline order: the
! cylindrical case exercises the periodic phi spline (spl_per), the cartesian
! case the regular spline (spl_reg), both with order-5 derivatives in the curl.
subroutine test_curl_order5
    use neo_field, only: create_field
    use neo_field_base, only: field_t

    character(*), parameter :: fcyl = "test_mgrid_o5_cyl.nc"
    character(*), parameter :: fcart = "test_mgrid_o5_cart.nc"
    integer, parameter :: nfp = 4
    real(dp), parameter :: tol = 1.0e-7_dp
    type(field_mesh_t) :: mesh
    class(field_t), allocatable :: field
    real(dp) :: x(3), B(3), A_ref(3), B_ref(3)

    call print_test("test_curl_order5")

    call build_cyl_mesh(mesh, nfp, 40, 12, 40, store_a=.true.)
    call write_mgrid(fcyl, mesh, nfp, coordinate_system=COORD_CYLINDRICAL, &
        write_a=.true., write_b=.false.)
    call create_field(field, "spline", filename=fcyl, order=5)
    call unlink(fcyl)
    x = [1.13_dp, 0.37_dp*twopi/real(nfp, dp), 0.21_dp]
    call field%compute_bfield(x, B)
    call analytic_a_cyl(x, A_ref, B_ref)
    if (maxval(abs(B - B_ref)) > tol) then
        print *, "cyl order5 B=", B, " ref=", B_ref
        call print_fail
        error stop
    end if
    deallocate(field)

    call build_cart_a_mesh(mesh, 40)
    call write_mgrid(fcart, mesh, 1, coordinate_system=COORD_CARTESIAN, &
        write_a=.true., write_b=.false.)
    call create_field(field, "spline", filename=fcart, order=5)
    call unlink(fcart)
    x = [1.37_dp, 1.51_dp, 1.62_dp]
    call field%compute_bfield(x, B)
    call analytic_a_cart(x, A_ref, B_ref)
    if (maxval(abs(B - B_ref)) > tol) then
        print *, "cart order5 B=", B, " ref=", B_ref
        call print_fail
        error stop
    end if
    call print_ok
end subroutine test_curl_order5

function disk_dim(fname, name) result(n)
    character(*), intent(in) :: fname, name
    integer :: n, ncid, did, stat

    call nc(nf90_open(fname, NF90_NOWRITE, ncid))
    call nc(nf90_inq_dimid(ncid, name, did))
    call nc(nf90_inquire_dimension(ncid, did, len=n))
    stat = nf90_close(ncid)
end function disk_dim

subroutine build_raw_mgrid(fname, nfp, ir, jz, kp, c1, c2)
    character(*), intent(in) :: fname
    integer, intent(in) :: nfp, ir, jz, kp
    real(dp), intent(in) :: c1, c2

    integer :: ncid, d_r, d_z, d_p, d1, d2, v
    real(dp) :: f(ir,jz,kp)
    character(len=1) :: mode_c

    call nc(nf90_create(fname, NF90_CLOBBER, ncid))
    call nc(nf90_def_dim(ncid, 'rad', ir, d_r))
    call nc(nf90_def_dim(ncid, 'zee', jz, d_z))
    call nc(nf90_def_dim(ncid, 'phi', kp, d_p))
    call nc(nf90_def_dim(ncid, 'dim_00001', 1, d1))
    call nc(nf90_def_dim(ncid, 'external_coils', 2, d2))
    call nc(nf90_def_var(ncid, 'nfp', NF90_INT, v))
    call nc(nf90_def_var(ncid, 'nextcur', NF90_INT, v))
    call nc(nf90_def_var(ncid, 'rmin', NF90_DOUBLE, v))
    call nc(nf90_def_var(ncid, 'rmax', NF90_DOUBLE, v))
    call nc(nf90_def_var(ncid, 'zmin', NF90_DOUBLE, v))
    call nc(nf90_def_var(ncid, 'zmax', NF90_DOUBLE, v))
    call nc(nf90_def_var(ncid, 'coil_group', NF90_CHAR, [d1], v))
    call nc(nf90_def_var(ncid, 'mgrid_mode', NF90_CHAR, [d1], v))
    call nc(nf90_def_var(ncid, 'raw_coil_cur', NF90_DOUBLE, [d2], v))
    call def_grp(ncid, 'br_001', d_r, d_z, d_p)
    call def_grp(ncid, 'bp_001', d_r, d_z, d_p)
    call def_grp(ncid, 'bz_001', d_r, d_z, d_p)
    call def_grp(ncid, 'br_002', d_r, d_z, d_p)
    call def_grp(ncid, 'bp_002', d_r, d_z, d_p)
    call def_grp(ncid, 'bz_002', d_r, d_z, d_p)
    call nc(nf90_enddef(ncid))

    call put_i(ncid, 'nfp', nfp); call put_i(ncid, 'nextcur', 2)
    call put_d(ncid, 'rmin', 0.8_dp); call put_d(ncid, 'rmax', 1.6_dp)
    call put_d(ncid, 'zmin', -0.5_dp); call put_d(ncid, 'zmax', 0.5_dp)
    mode_c = 'R'
    call nc(nf90_inq_varid(ncid, 'mgrid_mode', v)); call nc(nf90_put_var(ncid, v, mode_c))
    call nc(nf90_inq_varid(ncid, 'raw_coil_cur', v)); call nc(nf90_put_var(ncid, v, [c1, c2]))
    f = 1.0_dp;   call put_grp(ncid, 'br_001', f)
    f = 10.0_dp;  call put_grp(ncid, 'bp_001', f)
    f = 100.0_dp; call put_grp(ncid, 'bz_001', f)
    f = 2.0_dp;   call put_grp(ncid, 'br_002', f)
    f = 20.0_dp;  call put_grp(ncid, 'bp_002', f)
    f = 200.0_dp; call put_grp(ncid, 'bz_002', f)
    call nc(nf90_close(ncid))
end subroutine build_raw_mgrid

subroutine def_grp(ncid, name, d_r, d_z, d_p)
    integer, intent(in) :: ncid, d_r, d_z, d_p
    character(*), intent(in) :: name
    integer :: v
    call nc(nf90_def_var(ncid, name, NF90_DOUBLE, [d_r, d_z, d_p], v))
end subroutine def_grp

subroutine put_i(ncid, name, val)
    integer, intent(in) :: ncid, val
    character(*), intent(in) :: name
    integer :: v
    call nc(nf90_inq_varid(ncid, name, v)); call nc(nf90_put_var(ncid, v, val))
end subroutine put_i

subroutine put_d(ncid, name, val)
    integer, intent(in) :: ncid
    character(*), intent(in) :: name
    real(dp), intent(in) :: val
    integer :: v
    call nc(nf90_inq_varid(ncid, name, v)); call nc(nf90_put_var(ncid, v, val))
end subroutine put_d

subroutine put_grp(ncid, name, f)
    integer, intent(in) :: ncid
    character(*), intent(in) :: name
    real(dp), intent(in) :: f(:,:,:)
    integer :: v
    call nc(nf90_inq_varid(ncid, name, v)); call nc(nf90_put_var(ncid, v, f))
end subroutine put_grp

subroutine nc(status)
    integer, intent(in) :: status
    if (status /= NF90_NOERR) then
        print *, 'netCDF error: ', trim(nf90_strerror(status))
        call print_fail
        error stop
    end if
end subroutine nc

end program test_mgrid
