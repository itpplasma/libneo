program test_chartmap_matches_vmec
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_chartmap_vmec_generator, only: write_chartmap_from_vmec
    use libneo_coordinates, only: coordinate_system_t, chartmap_coordinate_system_t, &
                                  make_chartmap_coordinate_system, &
                                  vmec_coordinate_system_t
    use math_constants, only: TWOPI
    use nctools_module, only: nc_open, nc_close, nc_get
    use new_vmec_stuff_mod, only: netcdffile
    use spline_vmec_sub, only: spline_vmec_data
    implicit none

    character(len=*), parameter :: wout_file = "wout.nc"
    character(len=*), parameter :: chartmap_file = "wout_vmec.chartmap.nc"

    class(coordinate_system_t), allocatable :: cs
    type(vmec_coordinate_system_t) :: vmec
    integer :: ierr
    integer :: i
    integer, parameter :: ntest = 25
    real(dp) :: u(3), u_shift(3), u_vmec(3)
    real(dp) :: x_chart(3), x_vmec_cart(3), x_vmec_cyl(3)
    real(dp) :: dx(3)
    real(dp) :: max_dx
    real(dp) :: max_dx_grid
    real(dp) :: zeta_period
    character(len=2048) :: message
    integer :: ncid
    real(dp), allocatable :: rho(:), theta(:), zeta(:)

    max_dx_grid = 0.0_dp
    max_dx = 0.0_dp

    call write_chartmap_from_vmec(wout_file, chartmap_file, 63, 64, 65, ierr, message)
    if (ierr /= 0) then
        print *, "  FAIL: write_chartmap_from_vmec: ", trim(message)
        error stop 1
    end if

    netcdffile = wout_file
    call spline_vmec_data

    call make_chartmap_coordinate_system(cs, chartmap_file)
    select type (ccs => cs)
    type is (chartmap_coordinate_system_t)
        zeta_period = TWOPI/real(ccs%num_field_periods, dp)
        call nc_open(chartmap_file, ncid)
        allocate (rho(63), theta(64), zeta(65))
        call nc_get(ncid, "rho", rho)
        call nc_get(ncid, "theta", theta)
        call nc_get(ncid, "zeta", zeta)
        call nc_close(ncid)

        do i = 1, min(ntest, 20)
            u(1) = rho(1 + 3*i)
            u(2) = theta(1 + 2*i)
            u(3) = zeta(1 + i)
            call ccs%evaluate_point(u, x_chart)
            u_vmec = [u(1)**2, u(2), u(3)]
            call vmec%evaluate_point(u_vmec, x_vmec_cyl)
            x_vmec_cart(1) = x_vmec_cyl(1)*cos(x_vmec_cyl(2))
            x_vmec_cart(2) = x_vmec_cyl(1)*sin(x_vmec_cyl(2))
            x_vmec_cart(3) = x_vmec_cyl(3)
            max_dx_grid = max(max_dx_grid, maxval(abs(x_chart - x_vmec_cart)))
        end do

        do i = 1, ntest
            u(1) = real(i, dp)/real(ntest + 1, dp)
            u(2) = modulo(real(37*i, dp)*0.017_dp*TWOPI, TWOPI)
            u(3) = modulo(real(29*i, dp)*0.013_dp*zeta_period, zeta_period)

            call ccs%evaluate_point(u, x_chart)

            u_vmec = [u(1)**2, u(2), u(3)]
            call vmec%evaluate_point(u_vmec, x_vmec_cyl)
            x_vmec_cart(1) = x_vmec_cyl(1)*cos(x_vmec_cyl(2))
            x_vmec_cart(2) = x_vmec_cyl(1)*sin(x_vmec_cyl(2))
            x_vmec_cart(3) = x_vmec_cyl(3)

            dx = x_chart - x_vmec_cart
            max_dx = max(max_dx, maxval(abs(dx)))

            u_shift = [u(1), u(2) + TWOPI, u(3)]
            call ccs%evaluate_point(u_shift, x_vmec_cart)
            dx = x_chart - x_vmec_cart
            max_dx = max(max_dx, maxval(abs(dx)))

            u_shift = [u(1), u(2), u(3) + zeta_period]
            call ccs%evaluate_point(u_shift, x_vmec_cart)
            dx = x_chart - x_vmec_cart
            max_dx = max(max_dx, maxval(abs(dx)))
        end do
    class default
        print *, "  FAIL: make_chartmap_coordinate_system did not return chartmap type"
        error stop 1
    end select

    if (max_dx > 2.0e-2_dp) then
        print *, "  FAIL: chartmap evaluate_point differs from VMEC"
        print *, "    max|dx| at grid nodes      =", max_dx_grid
        print *, "    max|dx| off grid points    =", max_dx
        error stop 1
    end if

    print *, "  PASS: chartmap matches VMEC max|dx|=", max_dx
end program test_chartmap_matches_vmec
