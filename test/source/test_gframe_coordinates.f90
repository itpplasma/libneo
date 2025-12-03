program test_gframe_coordinates
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_coordinates, only: coordinate_system_t, &
        make_gframe_coordinate_system, gframe_coordinate_system_t
    use gframe_boundary, only: gframe_boundary_t, make_circular_boundary
    implicit none

    integer :: nerrors
    logical :: all_passed
    class(coordinate_system_t), allocatable :: cs
    type(gframe_boundary_t) :: boundary
    type(gframe_coordinate_system_t), pointer :: gcs
    real(dp) :: u(3), x(3), u_back(3)
    integer :: ierr

    nerrors = 0
    all_passed = .true.

    print *, "Testing Babin/G-frame coordinate system..."
    call make_circular_boundary(boundary, ntheta = 64, nzeta = 8, r0 = 1.7_dp, &
        a = 0.35_dp)

    call make_gframe_coordinate_system(cs, boundary)

    if (.not. allocated(cs)) then
        print *, "  FAIL: make_gframe_coordinate_system did not allocate cs"
        nerrors = nerrors + 1
    end if

    select type (gcs => cs)
    type is (gframe_coordinate_system_t)
        call run_roundtrip_check(gcs, nerrors)
        call run_boundary_check(gcs, nerrors)
        call run_metric_check(gcs, nerrors)
    class default
        print *, "  FAIL: coordinate system is not gframe type"
        nerrors = nerrors + 1
    end select

    if (nerrors > 0) then
        all_passed = .false.
        print *, "FAILED: ", nerrors, " error(s) detected in G-frame tests"
    else
        print *, "All G-frame coordinate tests passed!"
    end if

    if (.not. all_passed) error stop 1

contains

    subroutine run_roundtrip_check(gcs, nerrors)
        type(gframe_coordinate_system_t), intent(in) :: gcs
        integer, intent(inout) :: nerrors
        real(dp) :: u(3), x(3), u_back(3)
        integer :: ierr

        u = [0.35_dp, 1.1_dp, 0.6_dp]
        call gcs%evaluate_point(u, x)
        call gcs%from_cyl(x, u_back, ierr)

        if (ierr /= 0) then
            print *, "  FAIL: inverse mapping reported error code ", ierr
            nerrors = nerrors + 1
        else if (maxval(abs(u_back - u)) > 1.0e-10_dp) then
            print *, "  FAIL: roundtrip mismatch: |u_back - u| = ", &
                maxval(abs(u_back - u))
            nerrors = nerrors + 1
        else
            print *, "  PASS: forward/inverse roundtrip within tolerance"
        end if

        call verify_against_analytic(u, x, nerrors)
    end subroutine run_roundtrip_check

    subroutine run_boundary_check(gcs, nerrors)
        type(gframe_coordinate_system_t), intent(in) :: gcs
        integer, intent(inout) :: nerrors
        real(dp) :: u(3), x(3), expected(3)
        real(dp), parameter :: tol = 5.0e-13_dp
        integer :: k

        do k = 0, 3
            u = [1.0_dp, 0.5_dp*pi_dp()*k, 1.1_dp]
            call gcs%evaluate_point(u, x)
            call analytic_map(u, expected)
            if (maxval(abs(x - expected)) > tol) then
                print *, "  FAIL: boundary point mismatch at k=", k
                nerrors = nerrors + 1
            end if
        end do
        if (k == 4) print *, "  PASS: boundary recovery matches analytic curve"
    end subroutine run_boundary_check

    subroutine run_metric_check(gcs, nerrors)
        type(gframe_coordinate_system_t), intent(in) :: gcs
        integer, intent(inout) :: nerrors
        real(dp) :: u(3), g(3,3), ginv(3,3), sqrtg

        u = [0.2_dp, 0.7_dp, 0.9_dp]
        call gcs%metric_tensor(u, g, ginv, sqrtg)
        if (sqrtg <= 0.0_dp) then
            print *, "  FAIL: Jacobian determinant not positive"
            nerrors = nerrors + 1
        else
            print *, "  PASS: Jacobian determinant positive"
        end if
    end subroutine run_metric_check

    subroutine verify_against_analytic(u, x, nerrors)
        real(dp), intent(in) :: u(3)
        real(dp), intent(in) :: x(3)
        integer, intent(inout) :: nerrors
        real(dp) :: expected(3)
        real(dp), parameter :: tol = 5.0e-13_dp

        call analytic_map(u, expected)
        if (maxval(abs(expected - x)) > tol) then
            print *, "  FAIL: forward map differs from analytic by ", &
                maxval(abs(expected - x))
            nerrors = nerrors + 1
        else
            print *, "  PASS: forward map matches analytic reference"
        end if
    end subroutine verify_against_analytic

    subroutine analytic_map(u, x)
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: x(3)
        real(dp), parameter :: r0 = 1.7_dp, a = 0.35_dp

        x(1) = r0 + a*u(1)*cos(u(2))
        x(2) = u(3)
        x(3) = a*u(1)*sin(u(2))
    end subroutine analytic_map

    pure real(dp) function pi_dp()
        pi_dp = acos(-1.0_dp)
    end function pi_dp

end program test_gframe_coordinates
