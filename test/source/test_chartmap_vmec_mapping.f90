program test_chartmap_vmec_mapping
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_chartmap_vmec_generator, only: write_chartmap_from_vmec
    use libneo_coordinates, only: coordinate_system_t, chartmap_coordinate_system_t, &
                                  make_chartmap_coordinate_system, &
                                  vmec_coordinate_system_t
    use libneo_coordinates_mapping, only: map_vmec_u_to_chartmap_u
    use math_constants, only: TWOPI
    use new_vmec_stuff_mod, only: netcdffile
    use spline_vmec_sub, only: spline_vmec_data
    implicit none

    character(len=*), parameter :: wout_file = "wout.nc"
    character(len=*), parameter :: chartmap_file = "wout_vmec.mapping.chartmap.nc"

    type(vmec_coordinate_system_t) :: vmec
    class(coordinate_system_t), allocatable :: cs
    integer :: ierr
    integer :: i
    integer, parameter :: ntest = 25
    real(dp) :: u_vmec(3)
    real(dp) :: u_chart(3)
    real(dp) :: x_chart(3)
    real(dp) :: x_vmec_cyl(3), x_vmec_cart(3)
    real(dp) :: zeta_period
    real(dp) :: max_dx
    character(len=2048) :: message

    max_dx = 0.0_dp

    call write_chartmap_from_vmec(wout_file, chartmap_file, 32, 33, 34, ierr, message)
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
        do i = 1, ntest
            call build_test_u_vmec(i, ntest, zeta_period, u_vmec)

            call vmec%evaluate_point(u_vmec, x_vmec_cyl)
            call cyl_to_cart(x_vmec_cyl, x_vmec_cart)

            call map_vmec_u_to_chartmap_u(vmec, ccs, u_vmec, u_chart, ierr, message)
            if (ierr /= 0) then
                print *, "  FAIL: map_vmec_u_to_chartmap_u: ", trim(message)
                error stop 1
            end if

            call ccs%evaluate_point(u_chart, x_chart)
            max_dx = max(max_dx, maxval(abs(x_chart - x_vmec_cart)))
        end do
    class default
        print *, "  FAIL: chartmap file did not load as chartmap type"
        error stop 1
    end select

    if (max_dx > 2.0e-2_dp) then
        print *, "  FAIL: vmec->chartmap mapping mismatch"
        print *, "    max|dx| =", max_dx
        error stop 1
    end if

    print *, "  PASS: vmec->chartmap mapping max|dx|=", max_dx

contains

    subroutine build_test_u_vmec(i, n, zeta_period, u)
        integer, intent(in) :: i
        integer, intent(in) :: n
        real(dp), intent(in) :: zeta_period
        real(dp), intent(out) :: u(3)

        u(1) = real(i, dp)/real(n + 1, dp)
        u(2) = modulo(real(37*i, dp)*0.017_dp*TWOPI, TWOPI)
        u(3) = modulo(real(29*i, dp)*0.013_dp*zeta_period, zeta_period)
    end subroutine build_test_u_vmec

    pure subroutine cyl_to_cart(xcyl, x)
        real(dp), intent(in) :: xcyl(3)
        real(dp), intent(out) :: x(3)

        x(1) = xcyl(1)*cos(xcyl(2))
        x(2) = xcyl(1)*sin(xcyl(2))
        x(3) = xcyl(3)
    end subroutine cyl_to_cart

end program test_chartmap_vmec_mapping

