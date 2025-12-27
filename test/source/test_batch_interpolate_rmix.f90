program test_batch_interpolate_rmix
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use interpolate, only: BatchSplineData3D, construct_batch_splines_3d, &
                           destroy_batch_splines_3d, evaluate_batch_splines_3d_der2
    use batch_interpolate_3d, only: evaluate_batch_splines_3d_der2_rmix

    implicit none

    real(dp), parameter :: TOL = 1.0d-12

    call test_rmix_matches_full_der2_nq1_o555
    call test_rmix_matches_full_der2_nq1_general_periodic
    call test_rmix_matches_full_der2_nq2_general
    call test_rmix_matches_full_der2_nq2_o555

    print *, "PASSED: evaluate_batch_splines_3d_der2_rmix matches full der2"

contains

    subroutine test_rmix_matches_full_der2_nq1_o555
        integer, parameter :: N_POINTS(3) = [24, 26, 28]
        integer, parameter :: ORDER(3) = [5, 5, 5]
        logical, parameter :: PERIODIC(3) = [.false., .false., .false.]

        type(BatchSplineData3D) :: spl
        real(dp) :: y_data(N_POINTS(1), N_POINTS(2), N_POINTS(3), 1)
        real(dp) :: x_min(3), x_max(3)
        real(dp) :: x_eval(3)
        real(dp) :: y1(1), y2(1)
        real(dp) :: dy1(3, 1), dy2(3, 1)
        real(dp) :: d2_full(6, 1), d2_rmix(3, 1)
        integer :: i1, i2, i3, n_test

        x_min = [0.1d0, -0.2d0, 0.3d0]
        x_max = [1.3d0, 0.9d0, 2.1d0]

        do i3 = 1, N_POINTS(3)
            do i2 = 1, N_POINTS(2)
                do i1 = 1, N_POINTS(1)
                    y_data(i1, i2, i3, 1) = sin(real(i1 - 1, dp)*0.07d0) + &
                                            cos(real(i2 - 1, dp)*0.05d0) + &
                                            exp(-real(i3 - 1, dp)*0.03d0)
                end do
            end do
        end do

        call construct_batch_splines_3d(x_min, x_max, y_data, ORDER, PERIODIC, spl)

        do n_test = 1, 7
            x_eval(1) = x_min(1) + (x_max(1) - x_min(1)) * real(1 + n_test, dp) / 9.0d0
            x_eval(2) = x_min(2) + (x_max(2) - x_min(2)) * real(2 + n_test, dp) / 9.0d0
            x_eval(3) = x_min(3) + (x_max(3) - x_min(3)) * real(3 + n_test, dp) / 9.0d0

            call evaluate_batch_splines_3d_der2(spl, x_eval, y1, dy1, d2_full)
            call evaluate_batch_splines_3d_der2_rmix(spl, x_eval, y2, dy2, d2_rmix)

            call assert_close_vec("value", y2, y1, TOL)
            call assert_close_mat("dy", dy2, dy1, TOL)
            call assert_close_vec("d2_rmix(1)", d2_rmix(1, :), d2_full(1, :), TOL)
            call assert_close_vec("d2_rmix(2)", d2_rmix(2, :), d2_full(2, :), TOL)
            call assert_close_vec("d2_rmix(3)", d2_rmix(3, :), d2_full(3, :), TOL)
        end do

        call destroy_batch_splines_3d(spl)
    end subroutine test_rmix_matches_full_der2_nq1_o555

    subroutine test_rmix_matches_full_der2_nq2_general
        integer, parameter :: N_POINTS(3) = [21, 23, 25]
        integer, parameter :: ORDER(3) = [3, 5, 3]
        logical, parameter :: PERIODIC(3) = [.false., .true., .false.]
        integer, parameter :: NQ = 2

        type(BatchSplineData3D) :: spl
        real(dp) :: y_data(N_POINTS(1), N_POINTS(2), N_POINTS(3), NQ)
        real(dp) :: x_min(3), x_max(3)
        real(dp) :: x_eval(3)
        real(dp) :: y1(NQ), y2(NQ)
        real(dp) :: dy1(3, NQ), dy2(3, NQ)
        real(dp) :: d2_full(6, NQ), d2_rmix(3, NQ)
        integer :: i1, i2, i3, iq, n_test

        x_min = [-0.5d0, 0.0d0, -1.0d0]
        x_max = [0.7d0, 2.0d0, 1.2d0]

        do iq = 1, NQ
            do i3 = 1, N_POINTS(3)
                do i2 = 1, N_POINTS(2)
                    do i1 = 1, N_POINTS(1)
                        y_data(i1, i2, i3, iq) = real(iq, dp) * &
                                                 (0.1d0*real(i1, dp)**2 + &
                                                  0.2d0*sin(real(i2 - 1, dp)*0.07d0) + &
                                                  0.3d0*cos(real(i3 - 1, dp)*0.05d0))
                    end do
                end do
            end do
        end do

        call construct_batch_splines_3d(x_min, x_max, y_data, ORDER, PERIODIC, spl)

        do n_test = 1, 7
            x_eval(1) = x_min(1) + (x_max(1) - x_min(1)) * real(1 + n_test, dp) / 9.0d0
            x_eval(2) = x_min(2) + (x_max(2) - x_min(2)) * real(2 + n_test, dp) / 9.0d0
            x_eval(3) = x_min(3) + (x_max(3) - x_min(3)) * real(3 + n_test, dp) / 9.0d0

            call evaluate_batch_splines_3d_der2(spl, x_eval, y1, dy1, d2_full)
            call evaluate_batch_splines_3d_der2_rmix(spl, x_eval, y2, dy2, d2_rmix)

            call assert_close_vec("value", y2, y1, TOL)
            call assert_close_mat("dy", dy2, dy1, TOL)
            call assert_close_vec("d2_rmix(1)", d2_rmix(1, :), d2_full(1, :), TOL)
            call assert_close_vec("d2_rmix(2)", d2_rmix(2, :), d2_full(2, :), TOL)
            call assert_close_vec("d2_rmix(3)", d2_rmix(3, :), d2_full(3, :), TOL)
        end do

        call destroy_batch_splines_3d(spl)
    end subroutine test_rmix_matches_full_der2_nq2_general

    subroutine test_rmix_matches_full_der2_nq1_general_periodic
        integer, parameter :: N_POINTS(3) = [23, 25, 27]
        integer, parameter :: ORDER(3) = [3, 5, 3]
        logical, parameter :: PERIODIC(3) = [.false., .true., .false.]

        type(BatchSplineData3D) :: spl
        real(dp) :: y_data(N_POINTS(1), N_POINTS(2), N_POINTS(3), 1)
        real(dp) :: x_min(3), x_max(3)
        real(dp) :: x_eval(3)
        real(dp) :: y1(1), y2(1)
        real(dp) :: dy1(3, 1), dy2(3, 1)
        real(dp) :: d2_full(6, 1), d2_rmix(3, 1)
        integer :: i1, i2, i3, n_test

        x_min = [-0.4d0, 0.0d0, -0.9d0]
        x_max = [0.8d0, 2.0d0, 1.1d0]

        do i3 = 1, N_POINTS(3)
            do i2 = 1, N_POINTS(2)
                do i1 = 1, N_POINTS(1)
                    y_data(i1, i2, i3, 1) = 0.1d0*real(i1, dp)**2 + &
                                            0.2d0*sin(real(i2 - 1, dp)*0.07d0) + &
                                            0.3d0*cos(real(i3 - 1, dp)*0.05d0)
                end do
            end do
        end do

        call construct_batch_splines_3d(x_min, x_max, y_data, ORDER, PERIODIC, spl)

        do n_test = 1, 7
            x_eval(1) = x_min(1) + (x_max(1) - x_min(1)) * real(1 + n_test, dp) / 9.0d0
            x_eval(2) = x_min(2) + (x_max(2) - x_min(2)) * real(2 + n_test, dp) / 9.0d0
            x_eval(3) = x_min(3) + (x_max(3) - x_min(3)) * real(3 + n_test, dp) / 9.0d0

            call evaluate_batch_splines_3d_der2(spl, x_eval, y1, dy1, d2_full)
            call evaluate_batch_splines_3d_der2_rmix(spl, x_eval, y2, dy2, d2_rmix)

            call assert_close_vec("value", y2, y1, TOL)
            call assert_close_mat("dy", dy2, dy1, TOL)
            call assert_close_vec("d2_rmix(1)", d2_rmix(1, :), d2_full(1, :), TOL)
            call assert_close_vec("d2_rmix(2)", d2_rmix(2, :), d2_full(2, :), TOL)
            call assert_close_vec("d2_rmix(3)", d2_rmix(3, :), d2_full(3, :), TOL)
        end do

        call destroy_batch_splines_3d(spl)
    end subroutine test_rmix_matches_full_der2_nq1_general_periodic

    subroutine test_rmix_matches_full_der2_nq2_o555
        integer, parameter :: N_POINTS(3) = [24, 26, 28]
        integer, parameter :: ORDER(3) = [5, 5, 5]
        logical, parameter :: PERIODIC(3) = [.false., .true., .false.]
        integer, parameter :: NQ = 2

        type(BatchSplineData3D) :: spl
        real(dp) :: y_data(N_POINTS(1), N_POINTS(2), N_POINTS(3), NQ)
        real(dp) :: x_min(3), x_max(3)
        real(dp) :: x_eval(3)
        real(dp) :: y1(NQ), y2(NQ)
        real(dp) :: dy1(3, NQ), dy2(3, NQ)
        real(dp) :: d2_full(6, NQ), d2_rmix(3, NQ)
        integer :: i1, i2, i3, iq, n_test

        x_min = [0.2d0, 0.0d0, -0.6d0]
        x_max = [1.4d0, 2.0d0, 0.8d0]

        do iq = 1, NQ
            do i3 = 1, N_POINTS(3)
                do i2 = 1, N_POINTS(2)
                    do i1 = 1, N_POINTS(1)
                        y_data(i1, i2, i3, iq) = real(iq, dp) * &
                                                 (sin(real(i1 - 1, dp)*0.05d0) + &
                                                  cos(real(i2 - 1, dp)*0.04d0) + &
                                                  exp(-real(i3 - 1, dp)*0.03d0))
                    end do
                end do
            end do
        end do

        call construct_batch_splines_3d(x_min, x_max, y_data, ORDER, PERIODIC, spl)

        do n_test = 1, 7
            x_eval(1) = x_min(1) + (x_max(1) - x_min(1)) * real(1 + n_test, dp) / 9.0d0
            x_eval(2) = x_min(2) + (x_max(2) - x_min(2)) * real(2 + n_test, dp) / 9.0d0
            x_eval(3) = x_min(3) + (x_max(3) - x_min(3)) * real(3 + n_test, dp) / 9.0d0

            call evaluate_batch_splines_3d_der2(spl, x_eval, y1, dy1, d2_full)
            call evaluate_batch_splines_3d_der2_rmix(spl, x_eval, y2, dy2, d2_rmix)

            call assert_close_vec("value", y2, y1, TOL)
            call assert_close_mat("dy", dy2, dy1, TOL)
            call assert_close_vec("d2_rmix(1)", d2_rmix(1, :), d2_full(1, :), TOL)
            call assert_close_vec("d2_rmix(2)", d2_rmix(2, :), d2_full(2, :), TOL)
            call assert_close_vec("d2_rmix(3)", d2_rmix(3, :), d2_full(3, :), TOL)
        end do

        call destroy_batch_splines_3d(spl)
    end subroutine test_rmix_matches_full_der2_nq2_o555

    subroutine assert_close_vec(label, got, ref, tol)
        character(*), intent(in) :: label
        real(dp), intent(in) :: got(:), ref(:)
        real(dp), intent(in) :: tol
        integer :: i

        if (size(got) /= size(ref)) error stop "assert_close_vec: size mismatch"

        do i = 1, size(got)
            if (abs(got(i) - ref(i)) > tol) then
                print *, "Mismatch: ", trim(label), " i=", i
                print *, "  got=", got(i)
                print *, "  ref=", ref(i)
                print *, "  absdiff=", abs(got(i) - ref(i))
                error stop "assert_close_vec failed"
            end if
        end do
    end subroutine assert_close_vec

    subroutine assert_close_mat(label, got, ref, tol)
        character(*), intent(in) :: label
        real(dp), intent(in) :: got(:, :), ref(:, :)
        real(dp), intent(in) :: tol
        integer :: i, j

        if (any(shape(got) /= shape(ref))) error stop "assert_close_mat: shape mismatch"

        do j = 1, size(got, 2)
            do i = 1, size(got, 1)
                if (abs(got(i, j) - ref(i, j)) > tol) then
                    print *, "Mismatch: ", trim(label), " (", i, ",", j, ")"
                    print *, "  got=", got(i, j)
                    print *, "  ref=", ref(i, j)
                    print *, "  absdiff=", abs(got(i, j) - ref(i, j))
                    error stop "assert_close_mat failed"
                end if
            end do
        end do
    end subroutine assert_close_mat

end program test_batch_interpolate_rmix
