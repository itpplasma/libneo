program test_chartmap_coordinates
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_coordinates, only: coordinate_system_t, &
                                  make_chartmap_coordinate_system, &
                                  chartmap_coordinate_system_t
    use math_constants, only: TWOPI
    use nctools_module, only: nc_close, nc_get, nc_inq_dim, nc_open
    use netcdf
    implicit none

    integer :: nerrors
    logical :: all_passed
    class(coordinate_system_t), allocatable :: cs
    character(len=*), parameter :: volume_file = "chartmap.nc"

    nerrors = 0
    all_passed = .true.

    print *, "Testing chartmap coordinate system..."

    call make_chartmap_coordinate_system(cs, volume_file)

    if (.not. allocated(cs)) then
        print *, "  FAIL: make_chartmap_coordinate_system did not allocate cs"
        nerrors = nerrors + 1
    end if

    select type (ccs => cs)
    type is (chartmap_coordinate_system_t)
        call run_roundtrip_u_check(ccs, nerrors)
        call run_roundtrip_check(ccs, nerrors)
        call run_boundary_check(ccs, nerrors)
        call run_metric_check(ccs, nerrors)
        call run_covariant_basis_fd_check(ccs, nerrors)
    class default
        print *, "  FAIL: coordinate system is not chartmap type"
        nerrors = nerrors + 1
    end select

    if (nerrors > 0) then
        all_passed = .false.
        print *, "FAILED: ", nerrors, " error(s) detected in chartmap tests"
    else
        print *, "All chartmap coordinate tests passed!"
    end if

    if (.not. all_passed) error stop 1

contains

    subroutine nc_check(status)
        integer, intent(in) :: status
        if (status /= NF90_NOERR) then
            print *, trim(nf90_strerror(status))
            error stop 2
        end if
    end subroutine nc_check

    subroutine run_roundtrip_check(ccs, nerrors)
        type(chartmap_coordinate_system_t), intent(in) :: ccs
        integer, intent(inout) :: nerrors
        real(dp), allocatable :: x_ref(:, :, :), y_ref(:, :, :), z_ref(:, :, :)
        real(dp) :: rho_val, theta_val
        real(dp) :: u(3), x(3), u_back(3), xcyl(3), x_round(3)
        real(dp) :: diff_x
        integer :: ierr, ncid
        integer :: i_rho, i_theta
        integer, parameter :: nrho = 63, ntheta = 64, nzeta = 65
        real(dp), parameter :: tol_x = 1.0e-10_dp

        call nc_open(volume_file, ncid)
        allocate (x_ref(nrho, ntheta, nzeta))
        allocate (y_ref(nrho, ntheta, nzeta))
        allocate (z_ref(nrho, ntheta, nzeta))
        call nc_get(ncid, "x", x_ref)
        call nc_get(ncid, "y", y_ref)
        call nc_get(ncid, "z", z_ref)
        call nc_close(ncid)

        i_rho = 12
        i_theta = 17
        rho_val = real(i_rho - 1, dp)/real(nrho - 1, dp)
        theta_val = TWOPI*real(i_theta - 1, dp)/real(ntheta, dp)

        u = [rho_val, theta_val, 0.0_dp]
        call ccs%evaluate_cart(u, x)

        if (abs(x(1) - x_ref(i_rho, i_theta, 1)) > tol_x .or. &
            abs(x(2) - y_ref(i_rho, i_theta, 1)) > tol_x .or. &
            abs(x(3) - z_ref(i_rho, i_theta, 1)) > tol_x) then
            print *, "  FAIL: forward map does not reproduce reference X,Y,Z"
            nerrors = nerrors + 1
        end if

        xcyl(1) = sqrt(x(1)**2 + x(2)**2)
        xcyl(2) = atan2(x(2), x(1))
        xcyl(3) = x(3)

        call ccs%from_cyl(xcyl, u_back, ierr)

        if (ierr /= 0) then
            print *, "  FAIL: inverse mapping reported error code ", ierr
            nerrors = nerrors + 1
        else
            call ccs%evaluate_cart(u_back, x_round)
            diff_x = maxval(abs(x_round - x))

            if (diff_x > tol_x) then
                print *, "  FAIL: roundtrip mismatch in Cartesian coordinates, "// &
                    "max |x_round - x| = ", diff_x
                nerrors = nerrors + 1
            else
                print *, "  PASS: forward/inverse roundtrip against chartmap volume"
            end if
        end if

    end subroutine run_roundtrip_check

    subroutine run_roundtrip_u_check(ccs, nerrors)
        type(chartmap_coordinate_system_t), intent(in) :: ccs
        integer, intent(inout) :: nerrors

        real(dp) :: u(3), x(3), u_back(3), xcyl(3), x_round(3)
        real(dp) :: dtheta, dzeta
        integer :: ierr
        integer :: ir, it, iz
        real(dp), parameter :: rho_vals(3) = [0.05_dp, 0.55_dp, 0.95_dp]
        real(dp), parameter :: theta_vals(4) = &
                               [0.0_dp, TWOPI*0.5_dp, TWOPI*0.25_dp, TWOPI*0.75_dp]
        real(dp), parameter :: zeta_vals(4) = &
                               [TWOPI*0.5_dp, 0.0_dp, TWOPI*0.25_dp, TWOPI*0.75_dp]
        real(dp), parameter :: tol_u = 1.0e-8_dp
        real(dp), parameter :: tol_x = 2.0e-8_dp
        integer :: nerrors_start

        nerrors_start = nerrors

        do iz = 1, 4
            do it = 1, 4
                do ir = 0, 2
                    u(1) = rho_vals(ir + 1)
                    u(2) = theta_vals(it)
                    u(3) = zeta_vals(iz)

                    call ccs%evaluate_cart(u, x)

                    xcyl(1) = sqrt(x(1)**2 + x(2)**2)
                    xcyl(2) = atan2(x(2), x(1))
                    xcyl(3) = x(3)

                    call ccs%from_cyl(xcyl, u_back, ierr)
                    if (ierr /= 0) then
                        print *, "  FAIL: from_cyl ierr=", ierr, " for u=", u
                        nerrors = nerrors + 1
                        cycle
                    end if

                    dtheta = abs(modulo(u_back(2) - u(2) + 0.5_dp*TWOPI, TWOPI) - &
                                 0.5_dp*TWOPI)
                    dzeta = abs(modulo(u_back(3) - u(3) + 0.5_dp*TWOPI, TWOPI) - &
                                0.5_dp*TWOPI)

                    if (abs(u_back(1) - u(1)) > tol_u .or. dtheta > tol_u .or. &
                        dzeta > tol_u) then
                        print *, "  FAIL: u roundtrip mismatch u_back-u=", u_back - u
                        nerrors = nerrors + 1
                        cycle
                    end if

                    call ccs%evaluate_cart(u_back, x_round)
                    if (maxval(abs(x_round - x)) > tol_x) then
                        print *, "  FAIL: x(u_back) mismatch max|dx|=", &
                            maxval(abs(x_round - x))
                        nerrors = nerrors + 1
                        cycle
                    end if
                end do
            end do
        end do

        if (nerrors == nerrors_start) then
            print *, &
                "  PASS: u->x->u_back roundtrip across zeta slices and near-axis points"
        end if
    end subroutine run_roundtrip_u_check

    subroutine run_boundary_check(ccs, nerrors)
        type(chartmap_coordinate_system_t), intent(in) :: ccs
        integer, intent(inout) :: nerrors
        real(dp), allocatable :: x_ref(:, :, :), y_ref(:, :, :), z_ref(:, :, :)
        real(dp) :: u(3), x(3)
        integer :: ncid
        integer :: i_theta
        integer, parameter :: nrho = 63, ntheta = 64, nzeta = 65
        real(dp), parameter :: tol = 1.0e-8_dp
        integer :: nerrors_local

        nerrors_local = 0

        call nc_open(volume_file, ncid)
        allocate (x_ref(nrho, ntheta, nzeta))
        allocate (y_ref(nrho, ntheta, nzeta))
        allocate (z_ref(nrho, ntheta, nzeta))
        call nc_get(ncid, "x", x_ref)
        call nc_get(ncid, "y", y_ref)
        call nc_get(ncid, "z", z_ref)
        call nc_close(ncid)

        do i_theta = 1, 4
            u(1) = 1.0_dp
            u(2) = TWOPI*real(i_theta - 1, dp)/real(ntheta, dp)
            u(3) = 0.0_dp

            call ccs%evaluate_cart(u, x)

            if (abs(x(1) - x_ref(nrho, i_theta, 1)) > tol .or. &
                abs(x(2) - y_ref(nrho, i_theta, 1)) > tol .or. &
                abs(x(3) - z_ref(nrho, i_theta, 1)) > tol) then
                print *, "  FAIL: boundary point mismatch at theta index=", &
                    i_theta
                nerrors_local = nerrors_local + 1
            end if
        end do

        nerrors = nerrors + nerrors_local
        if (nerrors_local == 0) then
            print *, "  PASS: boundary recovery matches chartmap volume"
        end if
    end subroutine run_boundary_check

    subroutine run_metric_check(ccs, nerrors)
        type(chartmap_coordinate_system_t), intent(in) :: ccs
        integer, intent(inout) :: nerrors
        real(dp) :: u(3), g(3, 3), ginv(3, 3), sqrtg
        real(dp) :: identity(3, 3), prod(3, 3)
        real(dp), parameter :: tol = 1.0e-10_dp
        integer :: i, j, k
        logical :: sym_ok, id_ok

        u = [0.2_dp, 0.7_dp, 0.9_dp]
        call ccs%metric_tensor(u, g, ginv, sqrtg)

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
                    prod(i, j) = prod(i, j) + g(i, k)*ginv(k, j)
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

    subroutine run_covariant_basis_fd_check(ccs, nerrors)
        type(chartmap_coordinate_system_t), intent(in) :: ccs
        integer, intent(inout) :: nerrors

        real(dp) :: u(3), e_cov(3, 3), e_cov_fd(3, 3)
        real(dp) :: u_plus(3), u_minus(3), x_plus(3), x_minus(3)
        real(dp) :: rel_err, max_err
        real(dp), parameter :: h = 1.0e-6_dp
        real(dp), parameter :: tol = 1.0e-3_dp
        integer :: k
        integer :: n_tests, n_passed
        real(dp) :: rho_vals(3), theta_vals(3), zeta_vals(2)
        integer :: ir, it, iz

        rho_vals = [0.3_dp, 0.5_dp, 0.8_dp]
        theta_vals = [0.5_dp, 2.0_dp, 4.5_dp]
        zeta_vals = [0.3_dp, 1.5_dp]

        n_tests = 0
        n_passed = 0
        max_err = 0.0_dp

        do iz = 1, 2
            do it = 1, 3
                do ir = 1, 3
                    u = [rho_vals(ir), theta_vals(it), zeta_vals(iz)]
                    n_tests = n_tests + 1

                    call ccs%covariant_basis(u, e_cov)

                    do k = 1, 3
                        u_plus = u
                        u_minus = u
                        u_plus(k) = u(k) + h
                        u_minus(k) = u(k) - h
                        call ccs%evaluate_cart(u_plus, x_plus)
                        call ccs%evaluate_cart(u_minus, x_minus)
                        e_cov_fd(:, k) = (x_plus - x_minus) / (2.0_dp * h)
                    end do

                    rel_err = maxval(abs(e_cov - e_cov_fd)) / &
                              max(maxval(abs(e_cov_fd)), 1.0e-10_dp)

                    if (rel_err > max_err) max_err = rel_err

                    if (rel_err < tol) then
                        n_passed = n_passed + 1
                    else
                        print *, "  covariant_basis FD mismatch at u=", u
                        print *, "    rel_err=", rel_err
                        print *, "    e_cov(1,:)=", e_cov(1, :)
                        print *, "    e_cov_fd(1,:)=", e_cov_fd(1, :)
                    end if
                end do
            end do
        end do

        if (n_passed == n_tests) then
            print *, "  PASS: covariant_basis matches FD (", n_tests, &
                " points, max_err=", max_err, ")"
        else
            print *, "  FAIL: covariant_basis FD check failed ", &
                n_tests - n_passed, " of ", n_tests, " points"
            print *, "    max relative error: ", max_err
            nerrors = nerrors + 1
        end if
    end subroutine run_covariant_basis_fd_check

end program test_chartmap_coordinates
