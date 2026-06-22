program test_field_netcdf
use util_for_test, only: print_test, print_ok, print_fail
use neo_spline_field, only: spline_field_t, COORD_CARTESIAN, COORD_CYLINDRICAL
use neo_field_mesh, only: field_mesh_t
use neo_field_netcdf, only: write_field_mesh_netcdf, read_field_mesh_netcdf
implicit none

integer, parameter :: dp = kind(1.0d0)
real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)

call test_roundtrip_cartesian
call test_roundtrip_cylindrical_with_b
call test_curl_cartesian
call test_curl_cylindrical
call test_factory_netcdf

contains

! Analytic Cartesian potential A and its curl B. A = (y^2, z^2, x^2),
! so B = curl A = (-2z, -2x, -2y). Quadratic, reproduced exactly by cubics.
subroutine analytic_cart(x, A, B)
    real(dp), intent(in) :: x(3)
    real(dp), intent(out) :: A(3), B(3)

    A = [x(2)**2, x(3)**2, x(1)**2]
    B = [-2.0_dp*x(3), -2.0_dp*x(1), -2.0_dp*x(2)]
end subroutine analytic_cart

! Analytic cylindrical potential A = (0, R Z, 0), physical components,
! x = (R, phi, Z). Curl gives B = (-R, 0, 2 Z), which is divergence-free.
subroutine analytic_cyl(x, A, B)
    real(dp), intent(in) :: x(3)
    real(dp), intent(out) :: A(3), B(3)

    A = [0.0_dp, x(1)*x(3), 0.0_dp]
    B = [-x(1), 0.0_dp, 2.0_dp*x(3)]
end subroutine analytic_cyl

subroutine build_mesh(field_mesh, coordinate_system, n, with_b)
    type(field_mesh_t), intent(out) :: field_mesh
    integer, intent(in) :: coordinate_system, n(3)
    logical, intent(in) :: with_b

    real(dp), allocatable :: x1(:), x2(:), x3(:)
    real(dp), allocatable :: a1(:,:,:), a2(:,:,:), a3(:,:,:)
    real(dp), allocatable :: b1(:,:,:), b2(:,:,:), b3(:,:,:)
    logical :: per(3)
    real(dp) :: x(3), A(3), B(3)
    integer :: i, j, k

    allocate(x1(n(1)), x2(n(2)), x3(n(3)))
    allocate(a1(n(1),n(2),n(3)), a2(n(1),n(2),n(3)), a3(n(1),n(2),n(3)))
    allocate(b1(n(1),n(2),n(3)), b2(n(1),n(2),n(3)), b3(n(1),n(2),n(3)))

    if (coordinate_system == COORD_CYLINDRICAL) then
        x1 = linspace(1.0_dp, 2.0_dp, n(1))
        x2 = linspace(0.0_dp, 2.0_dp*pi, n(2))
        x3 = linspace(-1.0_dp, 1.0_dp, n(3))
        per = [.false., .true., .false.]
    else
        x1 = linspace(1.0_dp, 2.0_dp, n(1))
        x2 = linspace(1.0_dp, 2.0_dp, n(2))
        x3 = linspace(1.0_dp, 2.0_dp, n(3))
        per = [.false., .false., .false.]
    end if

    do i = 1, n(1)
        do j = 1, n(2)
            do k = 1, n(3)
                x = [x1(i), x2(j), x3(k)]
                if (coordinate_system == COORD_CYLINDRICAL) then
                    call analytic_cyl(x, A, B)
                else
                    call analytic_cart(x, A, B)
                end if
                a1(i,j,k) = A(1); a2(i,j,k) = A(2); a3(i,j,k) = A(3)
                b1(i,j,k) = B(1); b2(i,j,k) = B(2); b3(i,j,k) = B(3)
            end do
        end do
    end do

    call field_mesh%A1%mesh_init(x1, x2, x3, a1, per)
    call field_mesh%A2%mesh_init(x1, x2, x3, a2, per)
    call field_mesh%A3%mesh_init(x1, x2, x3, a3, per)
    if (with_b) then
        call field_mesh%B1%mesh_init(x1, x2, x3, b1, per)
        call field_mesh%B2%mesh_init(x1, x2, x3, b2, per)
        call field_mesh%B3%mesh_init(x1, x2, x3, b3, per)
    else
        call field_mesh%B1%mesh_init(x1, x2, x3, is_periodic=per)
        call field_mesh%B2%mesh_init(x1, x2, x3, is_periodic=per)
        call field_mesh%B3%mesh_init(x1, x2, x3, is_periodic=per)
    end if
end subroutine build_mesh

function linspace(start, stop, n) result(x)
    real(dp), intent(in) :: start, stop
    integer, intent(in) :: n
    real(dp) :: x(n)
    integer :: i

    do i = 1, n
        x(i) = start + (stop - start)*real(i-1, dp)/real(n-1, dp)
    end do
end function linspace

subroutine test_roundtrip_cartesian
    character(*), parameter :: fname = "test_field_cart.nc"
    type(field_mesh_t) :: mesh_out, mesh_in
    integer :: n(3), cs
    logical :: has_b
    real(dp) :: dev

    call print_test("test_roundtrip_cartesian")
    n = [12, 12, 12]
    call build_mesh(mesh_out, COORD_CARTESIAN, n, with_b=.false.)
    call write_field_mesh_netcdf(fname, mesh_out, COORD_CARTESIAN)
    call read_field_mesh_netcdf(fname, mesh_in, cs, has_b)
    call unlink(fname)

    dev = maxval(abs(mesh_in%A1%value - mesh_out%A1%value)) &
        + maxval(abs(mesh_in%A2%value - mesh_out%A2%value)) &
        + maxval(abs(mesh_in%A3%value - mesh_out%A3%value)) &
        + maxval(abs(mesh_in%A1%x1 - mesh_out%A1%x1))

    if (cs /= COORD_CARTESIAN .or. has_b .or. dev > 1.0e-12_dp) then
        print *, "cs=", cs, " has_b=", has_b, " dev=", dev
        call print_fail
        error stop
    end if
    call print_ok
end subroutine test_roundtrip_cartesian

subroutine test_roundtrip_cylindrical_with_b
    character(*), parameter :: fname = "test_field_cyl.nc"
    type(field_mesh_t) :: mesh_out, mesh_in
    integer :: n(3), cs
    logical :: has_b
    real(dp) :: dev

    call print_test("test_roundtrip_cylindrical_with_b")
    n = [12, 16, 12]
    call build_mesh(mesh_out, COORD_CYLINDRICAL, n, with_b=.true.)
    call write_field_mesh_netcdf(fname, mesh_out, COORD_CYLINDRICAL, nfp=5, write_b=.true.)
    call read_field_mesh_netcdf(fname, mesh_in, cs, has_b)
    call unlink(fname)

    dev = maxval(abs(mesh_in%B1%value - mesh_out%B1%value)) &
        + maxval(abs(mesh_in%B3%value - mesh_out%B3%value)) &
        + maxval(abs(mesh_in%A2%value - mesh_out%A2%value))

    if (cs /= COORD_CYLINDRICAL .or. .not. has_b .or. dev > 1.0e-12_dp) then
        print *, "cs=", cs, " has_b=", has_b, " dev=", dev
        call print_fail
        error stop
    end if
    call print_ok
end subroutine test_roundtrip_cylindrical_with_b

subroutine test_curl_cartesian
    real(dp), parameter :: tol = 1.0e-8_dp
    type(field_mesh_t) :: field_mesh
    type(spline_field_t) :: spline_field
    real(dp) :: x(3), B(3), A_ref(3), B_ref(3)

    call print_test("test_curl_cartesian")
    call build_mesh(field_mesh, COORD_CARTESIAN, [40, 40, 40], with_b=.false.)
    call spline_field%spline_field_init_acurl(field_mesh, COORD_CARTESIAN)

    x = [1.37_dp, 1.51_dp, 1.62_dp]
    call spline_field%compute_bfield(x, B)
    call analytic_cart(x, A_ref, B_ref)

    if (maxval(abs(B - B_ref)) > tol) then
        print *, "B = ", B, " B_ref = ", B_ref
        call print_fail
        error stop
    end if
    call print_ok
end subroutine test_curl_cartesian

subroutine test_curl_cylindrical
    real(dp), parameter :: tol = 1.0e-7_dp
    type(field_mesh_t) :: field_mesh
    type(spline_field_t) :: spline_field
    real(dp) :: x(3), B(3), A_ref(3), B_ref(3)

    call print_test("test_curl_cylindrical")
    call build_mesh(field_mesh, COORD_CYLINDRICAL, [48, 48, 48], with_b=.false.)
    call spline_field%spline_field_init_acurl(field_mesh, COORD_CYLINDRICAL)

    x = [1.33_dp, 0.7_dp, 0.21_dp]
    call spline_field%compute_bfield(x, B)
    call analytic_cyl(x, A_ref, B_ref)

    if (maxval(abs(B - B_ref)) > tol) then
        print *, "B = ", B, " B_ref = ", B_ref
        call print_fail
        error stop
    end if
    call print_ok
end subroutine test_curl_cylindrical

subroutine test_factory_netcdf
    use neo_field, only: create_field
    use neo_field_base, only: field_t

    real(dp), parameter :: tol = 1.0e-8_dp
    character(*), parameter :: fname = "test_field_factory.nc"
    type(field_mesh_t) :: field_mesh
    class(field_t), allocatable :: field
    real(dp) :: x(3), B(3), A_ref(3), B_ref(3)

    call print_test("test_factory_netcdf")
    call build_mesh(field_mesh, COORD_CARTESIAN, [40, 40, 40], with_b=.false.)
    call write_field_mesh_netcdf(fname, field_mesh, COORD_CARTESIAN)

    call create_field(field, "spline", filename=fname)
    call unlink(fname)

    x = [1.42_dp, 1.33_dp, 1.58_dp]
    call field%compute_bfield(x, B)
    call analytic_cart(x, A_ref, B_ref)

    if (maxval(abs(B - B_ref)) > tol) then
        print *, "B = ", B, " B_ref = ", B_ref
        call print_fail
        error stop
    end if
    call print_ok
end subroutine test_factory_netcdf

end program test_field_netcdf
