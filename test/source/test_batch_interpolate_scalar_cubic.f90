program test_batch_interpolate_scalar_cubic
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use batch_interpolate_types, only: BatchSplineData3D
    use batch_interpolate_3d, only: construct_batch_splines_3d, &
                                    destroy_batch_splines_3d, &
                                    evaluate_batch_splines_3d_der2, &
                                    evaluate_batch_spline_3d_scalar_cubic_der2

    implicit none

    integer, parameter :: N1 = 17, N2 = 19, N3 = 21
    integer, parameter :: NEVAL = 128
    real(dp), parameter :: TOL = 3.0e-12_dp
    real(dp), parameter :: XMIN(3) = [-0.4_dp, 0.0_dp, -1.2_dp]
    real(dp), parameter :: XMAX(3) = [1.1_dp, 2.0_dp, 0.8_dp]
    integer, parameter :: ORDER(3) = [3, 3, 3]
    logical, parameter :: PERIODIC(3) = [.false., .true., .true.]

    type(BatchSplineData3D) :: spl
    real(dp) :: samples(N1, N2, N3, 1)
    real(dp) :: x(3), y_ref(1), dy_ref(3, 1), d2y_ref(6, 1)
    real(dp) :: y_got, dy_got(3), d2y_got(6)
    real(dp) :: max_error
    integer :: i1, i2, i3, ie

    do i3 = 1, N3
        do i2 = 1, N2
            do i1 = 1, N1
                samples(i1, i2, i3, 1) = &
                    exp(0.17_dp*real(i1 - 1, dp))* &
                    (1.4_dp + 0.2_dp*sin(0.31_dp*real(i2 - 1, dp)))* &
                    (1.2_dp + 0.3_dp*cos(0.23_dp*real(i3 - 1, dp)))
            end do
        end do
    end do

    call construct_batch_splines_3d(XMIN, XMAX, samples, ORDER, PERIODIC, spl)

    max_error = 0.0_dp
    do ie = 1, NEVAL
        x(1) = XMIN(1) + (XMAX(1) - XMIN(1))* &
               real(mod(37*ie, 127), dp)/127.0_dp
        ! Exercise periodic wrapping on both sides of the tabulated interval.
        x(2) = XMIN(2) + (XMAX(2) - XMIN(2))* &
               (3.0_dp*real(mod(53*ie, 127), dp)/127.0_dp - 1.0_dp)
        x(3) = XMIN(3) + (XMAX(3) - XMIN(3))* &
               (3.0_dp*real(mod(71*ie, 127), dp)/127.0_dp - 1.0_dp)

        call evaluate_batch_splines_3d_der2(spl, x, y_ref, dy_ref, d2y_ref)
        call evaluate_batch_spline_3d_scalar_cubic_der2(spl, x, y_got, dy_got, d2y_got)

        max_error = max(max_error, relative_error(y_got, y_ref(1)))
        max_error = max(max_error, maxval(relative_error_vec(dy_got, dy_ref(:, 1))))
        max_error = max(max_error, maxval(relative_error_vec(d2y_got, d2y_ref(:, 1))))
    end do

    call destroy_batch_splines_3d(spl)

    if (max_error > TOL) then
        print *, "scalar cubic max relative error:", max_error
        error stop "scalar cubic evaluator disagrees with independent general evaluator"
    end if

    print *, "PASSED: scalar cubic evaluator max relative error", max_error

contains

    pure real(dp) function relative_error(got, expected)
        real(dp), intent(in) :: got, expected
        relative_error = abs(got - expected)/max(1.0_dp, abs(expected))
    end function relative_error

    pure function relative_error_vec(got, expected) result(error)
        real(dp), intent(in) :: got(:), expected(:)
        real(dp) :: error(size(got))
        error = abs(got - expected)/max(1.0_dp, abs(expected))
    end function relative_error_vec

end program test_batch_interpolate_scalar_cubic
