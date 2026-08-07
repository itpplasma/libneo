program test_batch_interpolate_scalar_cubic
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use batch_interpolate_types, only: BatchSplineData3D
    use batch_interpolate_3d, only: construct_batch_splines_3d, &
                                    destroy_batch_splines_3d, &
                                    evaluate_batch_splines_3d_der2, &
                                    evaluate_batch_spline_3d_scalar_cubic_der, &
                                    evaluate_batch_spline_3d_scalar_quintic_der, &
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

    call test_value_gradient_oracle
    call test_quintic_value_gradient_oracle

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

    subroutine test_value_gradient_oracle
        type(BatchSplineData3D) :: oracle_spl
        real(dp) :: basis(0:3, 3), dbasis(0:3, 3)
        real(dp) :: coeff, expected_y, expected_dy(3), x_local(3)
        real(dp) :: oracle_x(3), oracle_y, oracle_dy(3)
        integer :: j, k1, k2, k3

        oracle_spl%order = 3
        oracle_spl%num_points = 2
        oracle_spl%periodic = .false.
        oracle_spl%x_min = [-0.5_dp, 0.25_dp, -1.0_dp]
        oracle_spl%h_step = [2.0_dp, 1.5_dp, 0.75_dp]
        oracle_spl%inv_h_step = 1.0_dp/oracle_spl%h_step
        oracle_spl%period = oracle_spl%h_step
        oracle_spl%num_quantities = 1
        allocate (oracle_spl%coeff(1, 0:3, 0:3, 0:3, 1, 1, 1))

        do k3 = 0, 3
            do k2 = 0, 3
                do k1 = 0, 3
                    oracle_spl%coeff(1, k1, k2, k3, 1, 1, 1) = &
                        0.125_dp*real(1 + 2*k1 - 3*k2 + 5*k3, dp)
                end do
            end do
        end do

        x_local = [0.37_dp, 0.61_dp, 0.29_dp]*oracle_spl%h_step
        oracle_x = oracle_spl%x_min + x_local
        do j = 1, 3
            basis(0, j) = 1.0_dp
            dbasis(0, j) = 0.0_dp
            do k1 = 1, 3
                basis(k1, j) = x_local(j)**k1
                dbasis(k1, j) = real(k1, dp)*x_local(j)**(k1 - 1)
            end do
        end do
        expected_y = 0.0_dp
        expected_dy = 0.0_dp
        do k3 = 0, 3
            do k2 = 0, 3
                do k1 = 0, 3
                    coeff = oracle_spl%coeff(1, k1, k2, k3, 1, 1, 1)
                    expected_y = expected_y + coeff*basis(k1, 1)* &
                                 basis(k2, 2)*basis(k3, 3)
                    expected_dy(1) = expected_dy(1) + coeff*dbasis(k1, 1)* &
                                     basis(k2, 2)*basis(k3, 3)
                    expected_dy(2) = expected_dy(2) + coeff*basis(k1, 1)* &
                                     dbasis(k2, 2)*basis(k3, 3)
                    expected_dy(3) = expected_dy(3) + coeff*basis(k1, 1)* &
                                     basis(k2, 2)*dbasis(k3, 3)
                end do
            end do
        end do

        call evaluate_batch_spline_3d_scalar_cubic_der(oracle_spl, oracle_x, &
                                                       oracle_y, oracle_dy)
        if (relative_error(oracle_y, expected_y) > TOL .or. &
            maxval(relative_error_vec(oracle_dy, expected_dy)) > TOL) then
            error stop &
                "scalar cubic value/gradient failed analytic polynomial oracle"
        end if
        deallocate (oracle_spl%coeff)
    end subroutine test_value_gradient_oracle

    subroutine test_quintic_value_gradient_oracle
        type(BatchSplineData3D) :: oracle_spl
        real(dp) :: basis(0:5, 3), dbasis(0:5, 3)
        real(dp) :: coeff, expected_y, expected_dy(3), x_local(3)
        real(dp) :: oracle_x(3), oracle_y, oracle_dy(3)
        integer :: j, k1, k2, k3

        oracle_spl%order = 5
        oracle_spl%num_points = 2
        oracle_spl%periodic = .false.
        oracle_spl%x_min = [-0.5_dp, 0.25_dp, -1.0_dp]
        oracle_spl%h_step = [2.0_dp, 1.5_dp, 0.75_dp]
        oracle_spl%inv_h_step = 1.0_dp/oracle_spl%h_step
        oracle_spl%period = oracle_spl%h_step
        oracle_spl%num_quantities = 1
        allocate (oracle_spl%coeff(1, 0:5, 0:5, 0:5, 1, 1, 1))

        do k3 = 0, 5
            do k2 = 0, 5
                do k1 = 0, 5
                    oracle_spl%coeff(1, k1, k2, k3, 1, 1, 1) = &
                        0.03125_dp*real(1 + 2*k1 - 3*k2 + 5*k3, dp)
                end do
            end do
        end do

        x_local = [0.37_dp, 0.61_dp, 0.29_dp]*oracle_spl%h_step
        oracle_x = oracle_spl%x_min + x_local
        do j = 1, 3
            basis(0, j) = 1.0_dp
            dbasis(0, j) = 0.0_dp
            do k1 = 1, 5
                basis(k1, j) = x_local(j)**k1
                dbasis(k1, j) = real(k1, dp)*x_local(j)**(k1 - 1)
            end do
        end do
        expected_y = 0.0_dp
        expected_dy = 0.0_dp
        do k3 = 0, 5
            do k2 = 0, 5
                do k1 = 0, 5
                    coeff = oracle_spl%coeff(1, k1, k2, k3, 1, 1, 1)
                    expected_y = expected_y + coeff*basis(k1, 1)* &
                        basis(k2, 2)*basis(k3, 3)
                    expected_dy(1) = expected_dy(1) + coeff*dbasis(k1, 1)* &
                        basis(k2, 2)*basis(k3, 3)
                    expected_dy(2) = expected_dy(2) + coeff*basis(k1, 1)* &
                        dbasis(k2, 2)*basis(k3, 3)
                    expected_dy(3) = expected_dy(3) + coeff*basis(k1, 1)* &
                        basis(k2, 2)*dbasis(k3, 3)
                end do
            end do
        end do

        call evaluate_batch_spline_3d_scalar_quintic_der(oracle_spl, oracle_x, &
            oracle_y, oracle_dy)
        if (relative_error(oracle_y, expected_y) > TOL .or. &
                maxval(relative_error_vec(oracle_dy, expected_dy)) > TOL) then
            error stop &
                "scalar quintic value/gradient failed analytic polynomial oracle"
        end if
        deallocate (oracle_spl%coeff)
    end subroutine test_quintic_value_gradient_oracle

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
