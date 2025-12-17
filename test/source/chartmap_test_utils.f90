module chartmap_test_utils
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_coordinates, only: make_chartmap_coordinate_system, &
                                  coordinate_system_t, chartmap_coordinate_system_t
    use math_constants, only: TWOPI
    implicit none

    private
    public :: chartmap_roundtrip_check

contains

    subroutine chartmap_roundtrip_check(filename, tol_x, nerrors)
        character(len=*), intent(in) :: filename
        real(dp), intent(in) :: tol_x
        integer, intent(inout) :: nerrors

        class(coordinate_system_t), allocatable :: cs
        real(dp) :: u(3), x(3), u_back(3), xcyl(3), x_round(3)
        integer :: ierr_local
        real(dp) :: zeta_period

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
    end subroutine chartmap_roundtrip_check

end module chartmap_test_utils

