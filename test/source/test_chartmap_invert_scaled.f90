! A minimal uniformly-scaled chartmap subclass that mirrors SIMPLE's
! scaled_chartmap_coordinate_system_t (RZ scale): it overrides evaluate_cart,
! evaluate_cyl, covariant_basis, and metric_tensor to multiply the parent map by a
! constant, and INHERITS invert_cart from chartmap_coordinate_system_t. The test below
! checks that the inherited inverse inverts the SCALED map -- i.e. that the chartmap
! override routes its forward map and basis through the polymorphic bindings, not the
! internal unscaled spline. Without that fix the inverse would solve the unscaled map
! against the scaled target and miss the root.
module scaled_chartmap_test_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_coordinates, only: chartmap_coordinate_system_t
    implicit none
    private

    public :: scaled_chartmap_test_t

    type, extends(chartmap_coordinate_system_t) :: scaled_chartmap_test_t
        real(dp) :: cart_scale = 1.0_dp
    contains
        procedure :: evaluate_cart => scaled_test_evaluate_cart
        procedure :: evaluate_cyl => scaled_test_evaluate_cyl
        procedure :: covariant_basis => scaled_test_covariant_basis
        procedure :: metric_tensor => scaled_test_metric_tensor
    end type scaled_chartmap_test_t

contains

    subroutine scaled_test_evaluate_cart(self, u, x)
        class(scaled_chartmap_test_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: x(3)

        call self%chartmap_coordinate_system_t%evaluate_cart(u, x)
        x = self%cart_scale*x
    end subroutine scaled_test_evaluate_cart

    subroutine scaled_test_evaluate_cyl(self, u, x)
        class(scaled_chartmap_test_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: x(3)

        call self%chartmap_coordinate_system_t%evaluate_cyl(u, x)
        x(1) = self%cart_scale*x(1)
        x(3) = self%cart_scale*x(3)
    end subroutine scaled_test_evaluate_cyl

    subroutine scaled_test_covariant_basis(self, u, e_cov)
        class(scaled_chartmap_test_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: e_cov(3, 3)

        call self%chartmap_coordinate_system_t%covariant_basis(u, e_cov)
        e_cov = self%cart_scale*e_cov
    end subroutine scaled_test_covariant_basis

    subroutine scaled_test_metric_tensor(self, u, g, ginv, sqrtg)
        class(scaled_chartmap_test_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: g(3, 3), ginv(3, 3), sqrtg

        real(dp) :: scale2

        call self%chartmap_coordinate_system_t%metric_tensor(u, g, ginv, sqrtg)
        scale2 = self%cart_scale*self%cart_scale
        g = scale2*g
        ginv = ginv/scale2
        sqrtg = self%cart_scale*scale2*sqrtg
    end subroutine scaled_test_metric_tensor

end module scaled_chartmap_test_mod

program test_chartmap_invert_scaled
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_coordinates, only: coordinate_system_t, chartmap_coordinate_system_t, &
        make_chartmap_coordinate_system, CHARTMAP_LOCATED
    use scaled_chartmap_test_mod, only: scaled_chartmap_test_t
    use math_constants, only: TWOPI
    use netcdf
    implicit none

    character(len=*), parameter :: filename = "chartmap_invert_scaled.nc"
    real(dp), parameter :: cart_scale = 3.0_dp
    integer :: nerrors

    nerrors = 0

    call write_boozer_file(filename)
    call run_scaled_checks(filename, nerrors)

    if (nerrors > 0) then
        print *, "FAILED: ", nerrors, " error(s) in scaled invert_cart tests"
        error stop 1
    end if

    print *, "All scaled chartmap invert_cart tests passed!"

contains

    subroutine nc_check(status)
        integer, intent(in) :: status
        if (status /= NF90_NOERR) then
            print *, trim(nf90_strerror(status))
            error stop 2
        end if
    end subroutine nc_check

    pure real(dp) function angular_diff(a, b, period)
        real(dp), intent(in) :: a, b, period
        angular_diff = abs(modulo(a - b + 0.5_dp*period, period) - 0.5_dp*period)
    end function angular_diff

    ! Build the scaled subclass from a freshly loaded plain chartmap, then for a grid of
    ! interior points assert invert_cart(evaluate_cart(u0)) == u0 to tolerance with
    ! LOCATED. The unscaled chartmap with the SAME x_target would fail to locate, which
    ! is what makes this a genuine scaling test.
    subroutine run_scaled_checks(path, nerrors)
        character(len=*), intent(in) :: path
        integer, intent(inout) :: nerrors

        class(coordinate_system_t), allocatable :: cs
        type(scaled_chartmap_test_t) :: scaled
        real(dp) :: zeta_period, u0(3), x(3), u(3), u_guess(3), dth, dze
        integer :: ir, it, iz, status, nfail
        real(dp), parameter :: rho_vals(3) = [0.2_dp, 0.55_dp, 0.85_dp]
        real(dp), parameter :: theta_vals(3) = [0.6_dp, 2.7_dp, 4.9_dp]
        real(dp), parameter :: zeta_fracs(2) = [0.2_dp, 0.7_dp]
        real(dp), parameter :: tol = 1.0e-6_dp

        call make_chartmap_coordinate_system(cs, path)
        select type (ccs => cs)
            type is (chartmap_coordinate_system_t)
            scaled%chartmap_coordinate_system_t = ccs
        class default
            print *, "  FAIL: chartmap file did not load as chartmap type"
            nerrors = nerrors + 1
            return
        end select
        scaled%cart_scale = cart_scale

        zeta_period = TWOPI/real(scaled%num_field_periods, dp)
        nfail = 0
        do iz = 1, size(zeta_fracs)
            do it = 1, size(theta_vals)
                do ir = 1, size(rho_vals)
                    u0 = [rho_vals(ir), theta_vals(it), zeta_fracs(iz)*zeta_period]
                    call scaled%evaluate_cart(u0, x)
                    u_guess = [u0(1) + 0.06_dp, u0(2) + 0.12_dp, u0(3) - 0.04_dp]
                    call scaled%invert_cart(x, u_guess, u, status)

                    dth = angular_diff(u(2), u0(2), TWOPI)
                    dze = angular_diff(u(3), u0(3), zeta_period)
                    if (status /= CHARTMAP_LOCATED .or. abs(u(1) - u0(1)) > tol .or. &
                        dth > tol .or. dze > tol) then
                        print *, "  FAIL: scaled roundtrip u0=", u0
                        print *, "        u=", u, " status=", status
                        nfail = nfail + 1
                    end if
                end do
            end do
        end do
        nerrors = nerrors + nfail
        if (nfail == 0) print *, "  PASS: scaled invert_cart inverts the scaled map"

        call check_scale_actually_matters(scaled, zeta_period, nerrors)
    end subroutine run_scaled_checks

    ! Guard against a false pass: confirm the SCALED target is genuinely off the
    ! unscaled map. Inverting the scaled target on the PARENT (unscaled) inverse must
    ! NOT recover u0; if it did, the scale would be a no-op and the test above would be
    ! vacuous.
    subroutine check_scale_actually_matters(scaled, zeta_period, nerrors)
        type(scaled_chartmap_test_t), intent(in), target :: scaled
        real(dp), intent(in) :: zeta_period
        integer, intent(inout) :: nerrors

        real(dp) :: u0(3), x(3), u(3), u_guess(3)
        integer :: status

        u0 = [0.55_dp, 2.7_dp, 0.2_dp*zeta_period]
        call scaled%evaluate_cart(u0, x)
        u_guess = [u0(1) + 0.06_dp, u0(2) + 0.12_dp, u0(3) - 0.04_dp]
        call scaled%chartmap_coordinate_system_t%invert_cart(x, u_guess, u, status)
        if (abs(u(1) - u0(1)) < 1.0e-3_dp .and. status == CHARTMAP_LOCATED) then
            print *, "  FAIL: unscaled inverse recovered u0 from a scaled target;", &
                " scale is a no-op, scaled test is vacuous"
            nerrors = nerrors + 1
        else
            print *, "  PASS: scaled target is genuinely off the unscaled map"
        end if
    end subroutine check_scale_actually_matters

    ! Reactor-scale synthetic Boozer chartmap, nfp = 3, identical construction to
    ! test_chartmap_invert_cart so the fixture needs no external generator.
    subroutine write_boozer_file(path)
        character(len=*), intent(in) :: path

        integer :: ncid
        integer :: dim_rho, dim_theta, dim_zeta
        integer :: var_rho, var_theta, var_zeta, var_x, var_y, var_z, var_nfp
        integer, parameter :: nrho = 9, ntheta = 11, nzeta = 10, nfp = 3
        real(dp), parameter :: r0 = 180.0_dp
        real(dp), parameter :: a = 45.0_dp
        real(dp) :: rho(nrho), theta(ntheta), zeta(nzeta)
        real(dp) :: x(nrho, ntheta, nzeta), y(nrho, ntheta, nzeta), &
            z(nrho, ntheta, nzeta)
        real(dp) :: rho_val, theta_val, zeta_val
        real(dp) :: phi_geom, rmaj, zpos, period
        integer :: ir, it, iz

        period = TWOPI/real(nfp, dp)

        do ir = 1, nrho
            rho(ir) = real(ir - 1, dp)/real(nrho - 1, dp)
        end do
        do it = 1, ntheta
            theta(it) = TWOPI*real(it - 1, dp)/real(ntheta, dp)
        end do
        do iz = 1, nzeta
            zeta(iz) = period*real(iz - 1, dp)/real(nzeta, dp)
        end do

        do iz = 1, nzeta
            zeta_val = zeta(iz)
            do it = 1, ntheta
                theta_val = theta(it)
                do ir = 1, nrho
                    rho_val = rho(ir)
                    phi_geom = zeta_val + 0.08_dp*rho_val*sin(theta_val) + &
                        0.02_dp*sin(2.0_dp*zeta_val)
                    rmaj = r0 + a*rho_val*cos(theta_val) + &
                        0.04_dp*a*rho_val*cos(theta_val - 2.0_dp*zeta_val)
                    zpos = a*rho_val*sin(theta_val) + &
                        0.02_dp*a*rho_val*sin(zeta_val)

                    x(ir, it, iz) = rmaj*cos(phi_geom)
                    y(ir, it, iz) = rmaj*sin(phi_geom)
                    z(ir, it, iz) = zpos
                end do
            end do
        end do

        call nc_check(nf90_create(trim(path), NF90_NETCDF4, ncid))
        call nc_check(nf90_put_att(ncid, NF90_GLOBAL, "zeta_convention", "boozer"))
        call nc_check(nf90_def_dim(ncid, "rho", nrho, dim_rho))
        call nc_check(nf90_def_dim(ncid, "theta", ntheta, dim_theta))
        call nc_check(nf90_def_dim(ncid, "zeta", nzeta, dim_zeta))
        call nc_check(nf90_def_var(ncid, "rho", NF90_DOUBLE, [dim_rho], var_rho))
        call nc_check(nf90_def_var(ncid, "theta", NF90_DOUBLE, [dim_theta], var_theta))
        call nc_check(nf90_def_var(ncid, "zeta", NF90_DOUBLE, [dim_zeta], var_zeta))
        call nc_check(nf90_def_var(ncid, "x", NF90_DOUBLE, &
            [dim_rho, dim_theta, dim_zeta], var_x))
        call nc_check(nf90_def_var(ncid, "y", NF90_DOUBLE, &
            [dim_rho, dim_theta, dim_zeta], var_y))
        call nc_check(nf90_def_var(ncid, "z", NF90_DOUBLE, &
            [dim_rho, dim_theta, dim_zeta], var_z))
        call nc_check(nf90_def_var(ncid, "num_field_periods", NF90_INT, var_nfp))
        call nc_check(nf90_put_att(ncid, var_x, "units", "cm"))
        call nc_check(nf90_put_att(ncid, var_y, "units", "cm"))
        call nc_check(nf90_put_att(ncid, var_z, "units", "cm"))
        call nc_check(nf90_enddef(ncid))

        call nc_check(nf90_put_var(ncid, var_rho, rho))
        call nc_check(nf90_put_var(ncid, var_theta, theta))
        call nc_check(nf90_put_var(ncid, var_zeta, zeta))
        call nc_check(nf90_put_var(ncid, var_x, x))
        call nc_check(nf90_put_var(ncid, var_y, y))
        call nc_check(nf90_put_var(ncid, var_z, z))
        call nc_check(nf90_put_var(ncid, var_nfp, nfp))
        call nc_check(nf90_close(ncid))
    end subroutine write_boozer_file

end program test_chartmap_invert_scaled
