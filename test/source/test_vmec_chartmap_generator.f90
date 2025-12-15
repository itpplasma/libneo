program test_vmec_chartmap_generator
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_chartmap_vmec_generator, only: write_chartmap_from_vmec
    use libneo_coordinates, only: validate_chartmap_file, &
                                  make_chartmap_coordinate_system, &
                                  coordinate_system_t, chartmap_coordinate_system_t
    use math_constants, only: TWOPI
    implicit none

    integer :: ierr
    integer :: nerrors
    character(len=2048) :: message
    character(len=*), parameter :: wout_file = "wout.nc"
    character(len=*), parameter :: chartmap_file = "wout_vmec.chartmap.nc"

    nerrors = 0

    call write_chartmap_from_vmec(wout_file, chartmap_file, 32, 33, 34, ierr, message)
    if (ierr /= 0) then
        print *, "  FAIL: write_chartmap_from_vmec: ", trim(message)
        nerrors = nerrors + 1
    end if

    call validate_chartmap_file(chartmap_file, ierr, message)
    if (ierr /= 0) then
        print *, "  FAIL: validate_chartmap_file: ", trim(message)
        nerrors = nerrors + 1
    else
        print *, "  PASS: generated chartmap validated"
    end if

    call run_roundtrip(chartmap_file, nerrors)

    if (nerrors > 0) then
        print *, "FAILED: ", nerrors, " error(s) in VMEC->chartmap generator test"
        error stop 1
    end if

contains

    subroutine run_roundtrip(filename, nerrors)
        character(len=*), intent(in) :: filename
        integer, intent(inout) :: nerrors

        class(coordinate_system_t), allocatable :: cs
        real(dp) :: u(3), x(3), u_back(3), xcyl(3), x_round(3)
        integer :: ierr_local
        real(dp) :: zeta_period
        real(dp), parameter :: tol_x = 1.0e-2_dp

        call make_chartmap_coordinate_system(cs, filename)
        select type (ccs => cs)
        type is (chartmap_coordinate_system_t)
            zeta_period = TWOPI/real(ccs%num_field_periods, dp)
            u = [0.55_dp, 0.25_dp*TWOPI, 0.4_dp*zeta_period]
            call ccs%evaluate_cart(u, x)
            xcyl(1) = sqrt(x(1)**2 + x(2)**2)
            xcyl(2) = atan2(x(2), x(1))
            xcyl(3) = x(3)

            call ccs%from_cyl(xcyl, u_back, ierr_local)
            if (ierr_local /= 0) then
                print *, "  FAIL: from_cyl ierr=", ierr_local
                nerrors = nerrors + 1
                return
            end if

            call ccs%evaluate_cart(u_back, x_round)
            if (maxval(abs(x_round - x)) > tol_x) then
                print *, "  FAIL: x roundtrip mismatch max|dx|=", &
                    maxval(abs(x_round - x))
                nerrors = nerrors + 1
                return
            end if

            print *, "  PASS: generated chartmap roundtrip u->x->u_back"
        class default
            print *, &
                "  FAIL: make_chartmap_coordinate_system did not return chartmap type"
            nerrors = nerrors + 1
        end select
    end subroutine run_roundtrip

end program test_vmec_chartmap_generator
