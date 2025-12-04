program test_babin_coordinates
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_coordinates, only: coordinate_system_t, &
        make_babin_coordinate_system, babin_coordinate_system_t
    use math_constants, only: TWOPI
    use nctools_module, only: nc_open, nc_close, nc_get
    implicit none

    integer :: nerrors
    logical :: all_passed
    class(coordinate_system_t), allocatable :: cs
    type(babin_coordinate_system_t), pointer :: bcs
    real(dp) :: u(3), x(3), u_back(3)
    integer :: ierr
    character(len=*), parameter :: volume_file = "babin_circular.nc"

    nerrors = 0
    all_passed = .true.

    print *, "Testing Babin coordinate system (map2disc volume)..."

    call make_babin_coordinate_system(cs, volume_file)

    if (.not. allocated(cs)) then
        print *, "  FAIL: make_babin_coordinate_system did not allocate cs"
        nerrors = nerrors + 1
    end if

    select type (bcs => cs)
    type is (babin_coordinate_system_t)
        call run_roundtrip_check(bcs, nerrors)
        call run_boundary_check(bcs, nerrors)
        call run_metric_check(bcs, nerrors)
    class default
        print *, "  FAIL: coordinate system is not Babin type"
        nerrors = nerrors + 1
    end select

    if (nerrors > 0) then
        all_passed = .false.
        print *, "FAILED: ", nerrors, " error(s) detected in Babin coordinate tests"
    else
        print *, "All Babin coordinate tests passed!"
    end if

    if (.not. all_passed) error stop 1

contains

    subroutine run_roundtrip_check(bcs, nerrors)
        type(babin_coordinate_system_t), intent(in) :: bcs
        integer, intent(inout) :: nerrors
        real(dp), allocatable :: R_ref(:, :, :), Z_ref(:, :, :)
        real(dp) :: rho_val, theta_val
        real(dp) :: u(3), x(3), u_back(3)
        integer :: ierr, ncid
        integer :: i_rho, i_theta
        integer, parameter :: nrho = 63, ntheta = 64, nzeta = 65
        real(dp), parameter :: tol_u = 1.0e-8_dp
        real(dp), parameter :: tol_x = 1.0e-10_dp

        call nc_open(volume_file, ncid)
        allocate(R_ref(nrho, ntheta, nzeta))
        allocate(Z_ref(nrho, ntheta, nzeta))
        call nc_get(ncid, "R", R_ref)
        call nc_get(ncid, "Z", Z_ref)
        call nc_close(ncid)

        i_rho = 12
        i_theta = 17
        rho_val = real(i_rho - 1, dp) / real(nrho - 1, dp)
        theta_val = TWOPI * real(i_theta - 1, dp) / real(ntheta, dp)

        u = [rho_val, theta_val, 0.0_dp]
        call bcs%evaluate_point(u, x)

        if (abs(x(1) - R_ref(i_rho, i_theta, 1)) > tol_x .or. &
            abs(x(3) - Z_ref(i_rho, i_theta, 1)) > tol_x) then
            print *, "  FAIL: forward map does not reproduce reference R,Z"
            nerrors = nerrors + 1
        end if

        call bcs%from_cyl([x(1), u(3), x(3)], u_back, ierr)

        if (ierr /= 0) then
            print *, "  FAIL: inverse mapping reported error code ", ierr
            nerrors = nerrors + 1
        else if (abs(u_back(1) - rho_val) > tol_u .or. &
                 abs(modulo(u_back(2) - theta_val, TWOPI)) > tol_u) then
            print *, "  FAIL: roundtrip mismatch in Babin coordinates"
            nerrors = nerrors + 1
        else
            print *, "  PASS: forward/inverse roundtrip against map2disc volume"
        end if

    end subroutine run_roundtrip_check

    subroutine run_boundary_check(bcs, nerrors)
        type(babin_coordinate_system_t), intent(in) :: bcs
        integer, intent(inout) :: nerrors
        real(dp), allocatable :: R_ref(:, :, :), Z_ref(:, :, :)
        real(dp) :: u(3), x(3)
        integer :: ncid
        integer :: i_theta
        integer, parameter :: nrho = 63, ntheta = 64, nzeta = 65
        real(dp), parameter :: tol = 1.0e-8_dp

        call nc_open(volume_file, ncid)
        allocate(R_ref(nrho, ntheta, nzeta))
        allocate(Z_ref(nrho, ntheta, nzeta))
        call nc_get(ncid, "R", R_ref)
        call nc_get(ncid, "Z", Z_ref)
        call nc_close(ncid)

        do i_theta = 1, 4
            u(1) = 1.0_dp
            u(2) = TWOPI * real(i_theta - 1, dp) / real(ntheta, dp)
            u(3) = 0.0_dp

            call bcs%evaluate_point(u, x)

            if (abs(x(1) - R_ref(nrho, i_theta, 1)) > tol .or. &
                abs(x(3) - Z_ref(nrho, i_theta, 1)) > tol) then
                print *, "  FAIL: boundary point mismatch at theta index=", &
                    i_theta
                nerrors = nerrors + 1
            end if
        end do

        if (nerrors == 0) then
            print *, "  PASS: boundary recovery matches map2disc volume"
        end if
    end subroutine run_boundary_check

    subroutine run_metric_check(bcs, nerrors)
        type(babin_coordinate_system_t), intent(in) :: bcs
        integer, intent(inout) :: nerrors
        real(dp) :: u(3), g(3,3), ginv(3,3), sqrtg
        real(dp) :: identity(3,3), prod(3,3)
        real(dp), parameter :: tol = 1.0e-10_dp
        integer :: i, j, k
        logical :: sym_ok, id_ok

        u = [0.2_dp, 0.7_dp, 0.9_dp]
        call bcs%metric_tensor(u, g, ginv, sqrtg)

        if (sqrtg <= 0.0_dp) then
            print *, "  FAIL: Jacobian determinant not positive"
            nerrors = nerrors + 1
        else
            print *, "  PASS: Jacobian determinant positive"
        end if

        sym_ok = .true.
        do i = 1, 3
            do j = i + 1, 3
                if (abs(g(i, j) - g(j, i)) > tol) sym_ok = .false.
                if (abs(ginv(i, j) - ginv(j, i)) > tol) sym_ok = .false.
            end do
        end do
        if (.not. sym_ok) then
            print *, "  FAIL: metric tensor not symmetric"
            nerrors = nerrors + 1
        else
            print *, "  PASS: metric tensor symmetry"
        end if

        identity = 0.0_dp
        do i = 1, 3
            identity(i, i) = 1.0_dp
        end do

        prod = 0.0_dp
        do i = 1, 3
            do j = 1, 3
                do k = 1, 3
                    prod(i, j) = prod(i, j) + g(i, k) * ginv(k, j)
                end do
            end do
        end do

        id_ok = .true.
        do i = 1, 3
            do j = 1, 3
                if (abs(prod(i, j) - identity(i, j)) > tol) id_ok = .false.
            end do
        end do
        if (.not. id_ok) then
            print *, "  FAIL: g * ginv /= identity"
            nerrors = nerrors + 1
        else
            print *, "  PASS: g * ginv = identity"
        end if
    end subroutine run_metric_check

end program test_babin_coordinates
