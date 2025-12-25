module spline_eval_reference
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use interpolate, only: SplineData1D, SplineData2D, SplineData3D

    implicit none
    private

    public :: evaluate_splines_1d_ref, evaluate_splines_1d_der_ref, &
              evaluate_splines_1d_der2_ref
    public :: evaluate_splines_2d_ref, evaluate_splines_2d_der_ref
    public :: evaluate_splines_3d_ref, evaluate_splines_3d_der_ref, &
              evaluate_splines_3d_der2_ref

contains

    subroutine evaluate_splines_1d_ref(spl, x, y)
        type(SplineData1D), intent(in) :: spl
        real(dp), intent(in) :: x
        real(dp), intent(out) :: y

        real(dp) :: x_norm, x_local, xj
        integer :: interval_index, k_power

        if (spl%periodic) then
            xj = modulo(x - spl%x_min, spl%h_step*(spl%num_points-1)) + spl%x_min
        else
            xj = x
        end if
        x_norm = (xj - spl%x_min) / spl%h_step
        interval_index = max(0, min(spl%num_points-1, int(x_norm)))
        x_local = (x_norm - dble(interval_index))*spl%h_step

        y = spl%coeff(spl%order, interval_index+1)
        do k_power = spl%order-1, 0, -1
            y = spl%coeff(k_power, interval_index+1) + x_local*y
        enddo
    end subroutine evaluate_splines_1d_ref


    subroutine evaluate_splines_1d_der_ref(spl, x, y, dy)
        type(SplineData1D), intent(in) :: spl
        real(dp), intent(in) :: x
        real(dp), intent(out) :: y, dy

        real(dp) :: x_norm, x_local, xj
        integer :: interval_index, k_power

        if (spl%periodic) then
            xj = modulo(x - spl%x_min, spl%h_step*(spl%num_points-1)) + spl%x_min
        else
            xj = x
        end if
        x_norm = (xj - spl%x_min) / spl%h_step
        interval_index = max(0, min(spl%num_points-1, int(x_norm)))
        x_local = (x_norm - dble(interval_index))*spl%h_step

        y = spl%coeff(spl%order, interval_index+1)
        do k_power = spl%order-1, 0, -1
            y = spl%coeff(k_power, interval_index+1) + x_local*y
        enddo
        dy = spl%coeff(spl%order, interval_index+1)*spl%order
        do k_power = spl%order-1, 1, -1
            dy = spl%coeff(k_power, interval_index+1)*k_power + x_local*dy
        enddo
    end subroutine evaluate_splines_1d_der_ref


    subroutine evaluate_splines_1d_der2_ref(spl, x, y, dy, d2y)
        type(SplineData1D), intent(in) :: spl
        real(dp), intent(in) :: x
        real(dp), intent(out) :: y, dy, d2y

        real(dp) :: x_norm, x_local, xj
        integer :: interval_index, k_power

        if (spl%periodic) then
            xj = modulo(x - spl%x_min, spl%h_step*(spl%num_points-1)) + spl%x_min
        else
            xj = x
        end if
        x_norm = (xj - spl%x_min) / spl%h_step
        interval_index = max(0, min(spl%num_points-1, int(x_norm)))
        x_local = (x_norm - dble(interval_index))*spl%h_step

        y = spl%coeff(spl%order, interval_index+1)
        do k_power = spl%order-1, 0, -1
            y = spl%coeff(k_power, interval_index+1) + x_local*y
        enddo
        dy = spl%coeff(spl%order, interval_index+1)*spl%order
        do k_power = spl%order-1, 1, -1
            dy = spl%coeff(k_power, interval_index+1)*k_power + x_local*dy
        enddo
        d2y = spl%coeff(spl%order, interval_index+1)*spl%order*(spl%order-1)
        do k_power = spl%order-1, 2, -1
            d2y = spl%coeff(k_power, interval_index+1)*k_power*(k_power-1) + x_local*d2y
        enddo
    end subroutine evaluate_splines_1d_der2_ref


    subroutine evaluate_splines_2d_ref(spl, x, y)
        type(SplineData2D), intent(in) :: spl
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: y

        real(dp) :: x_norm(2), x_local(2), xj
        real(dp) :: coeff_2(0:spl%order(2))
        integer :: interval_index(2), k1, k2, j

        do j=1,2
            if (spl%periodic(j)) then
                xj = modulo(x(j) - spl%x_min(j), &
                    spl%h_step(j)*(spl%num_points(j)-1)) + spl%x_min(j)
            else
                xj = x(j)
            end if
            x_norm(j) = (xj - spl%x_min(j))/spl%h_step(j)
            interval_index(j) = max(0, min(spl%num_points(j)-1, int(x_norm(j))))
            x_local(j) = (x_norm(j) - dble(interval_index(j)))*spl%h_step(j)
        end do

        coeff_2(:) = spl%coeff(spl%order(1), 0:spl%order(2), &
                               interval_index(1) + 1, interval_index(2) + 1)
        do k1 = spl%order(1)-1, 0, -1
            coeff_2(:) = spl%coeff(k1, :, interval_index(1) + 1, interval_index(2) + 1) &
                         + x_local(1)*coeff_2
        enddo

        y = coeff_2(spl%order(2))
        do k2 = spl%order(2)-1, 0, -1
            y = coeff_2(k2) + x_local(2)*y
        enddo
    end subroutine evaluate_splines_2d_ref


    subroutine evaluate_splines_2d_der_ref(spl, x, y, dy)
        type(SplineData2D), intent(in) :: spl
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: y
        real(dp), intent(out) :: dy(2)

        real(dp) :: x_norm(2), x_local(2), xj
        real(dp) :: coeff_2(0:spl%order(2))
        real(dp) :: coeff_2_dx1(0:spl%order(2))
        integer :: interval_index(2), k1, k2, j

        do j=1,2
            if (spl%periodic(j)) then
                xj = modulo(x(j) - spl%x_min(j), &
                    spl%h_step(j)*(spl%num_points(j)-1)) + spl%x_min(j)
            else
                xj = x(j)
            end if
            x_norm(j) = (xj - spl%x_min(j))/spl%h_step(j)
            interval_index(j) = max(0, min(spl%num_points(j)-1, int(x_norm(j))))
            x_local(j) = (x_norm(j) - dble(interval_index(j)))*spl%h_step(j)
        end do

        coeff_2(:) = spl%coeff(spl%order(1), 0:spl%order(2), &
                               interval_index(1) + 1, interval_index(2) + 1)
        do k1 = spl%order(1)-1, 0, -1
            coeff_2(:) = spl%coeff(k1, :, interval_index(1) + 1, interval_index(2) + 1) &
                         + x_local(1)*coeff_2
        enddo

        coeff_2_dx1(:) = spl%coeff(spl%order(1), 0:spl%order(2), &
                                   interval_index(1) + 1, interval_index(2) + 1)*spl%order(1)
        do k1 = spl%order(1)-1, 1, -1
            coeff_2_dx1(:) = spl%coeff(k1, :, interval_index(1) + 1, interval_index(2) + 1)*k1 &
                             + x_local(1)*coeff_2_dx1
        enddo

        y = coeff_2(spl%order(2))
        do k2 = spl%order(2)-1, 0, -1
            y = coeff_2(k2) + x_local(2)*y
        enddo

        dy(1) = coeff_2_dx1(spl%order(2))
        do k2 = spl%order(2)-1, 0, -1
            dy(1) = coeff_2_dx1(k2) + x_local(2)*dy(1)
        enddo

        dy(2) = coeff_2(spl%order(2))*spl%order(2)
        do k2 = spl%order(2)-1, 1, -1
            dy(2) = coeff_2(k2)*k2 + x_local(2)*dy(2)
        enddo
    end subroutine evaluate_splines_2d_der_ref


    subroutine evaluate_splines_3d_ref(spl, x, y)
        type(SplineData3D), intent(in) :: spl
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: y

        real(dp) :: x_norm(3), x_local(3), xj
        real(dp) :: coeff_3(0:spl%order(3))
        real(dp) :: coeff_23(0:spl%order(2),0:spl%order(3))
        integer :: interval_index(3), k1, k2, k3, j

        do j=1,3
            if (spl%periodic(j)) then
                xj = modulo(x(j) - spl%x_min(j), &
                    spl%h_step(j)*(spl%num_points(j)-1)) + spl%x_min(j)
            else
                xj = x(j)
            end if
            x_norm(j) = (xj - spl%x_min(j))/spl%h_step(j)
            interval_index(j) = max(0, min(spl%num_points(j)-1, int(x_norm(j))))
            x_local(j) = (x_norm(j) - dble(interval_index(j)))*spl%h_step(j)
        end do

        coeff_23(:, :) = spl%coeff(0, :, :, interval_index(1) + 1, &
                                   interval_index(2) + 1, interval_index(3) + 1)
        do k1 = 1, spl%order(1)
            coeff_23(:, :) = spl%coeff(k1, :, :, interval_index(1) + 1, &
                                       interval_index(2) + 1, interval_index(3) + 1) &
                             + x_local(1)*coeff_23(:, :)
        enddo

        coeff_3(:) = coeff_23(spl%order(2), :)
        do k2 = spl%order(2)-1, 0, -1
            coeff_3(:) = coeff_23(k2, :) + x_local(2)*coeff_3
        enddo

        y = coeff_3(spl%order(3))
        do k3 = spl%order(3)-1, 0, -1
            y = coeff_3(k3) + x_local(3)*y
        enddo
    end subroutine evaluate_splines_3d_ref


    subroutine evaluate_splines_3d_der_ref(spl, x, y, dy)
        type(SplineData3D), intent(in) :: spl
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: y, dy(3)

        real(dp) :: x_norm(3), x_local(3), xj
        real(dp) :: coeff_23(0:spl%order(2),0:spl%order(3))
        real(dp) :: coeff_23_dx1(0:spl%order(2),0:spl%order(3))
        real(dp) :: coeff_3(0:spl%order(3))
        real(dp) :: coeff_3_dx1(0:spl%order(3))
        real(dp) :: coeff_3_dx2(0:spl%order(3))
        integer :: interval_index(3), k1, k2, k3, j

        dy = 0d0

        do j=1,3
            if (spl%periodic(j)) then
                xj = modulo(x(j) - spl%x_min(j), &
                    spl%h_step(j)*(spl%num_points(j)-1)) + spl%x_min(j)
            else
                xj = x(j)
            end if
            x_norm(j) = (xj - spl%x_min(j))/spl%h_step(j)
            interval_index(j) = max(0, min(spl%num_points(j)-1, int(x_norm(j))))
            x_local(j) = (x_norm(j) - dble(interval_index(j)))*spl%h_step(j)
        end do

        associate(N1 => spl%order(1), N2 => spl%order(2), N3 => spl%order(3))

            do k3 = 0, N3
                do k2 = 0, N2
                    coeff_23(k2, k3) = spl%coeff(0, k2, k3, interval_index(1) + 1, &
                        interval_index(2) + 1, interval_index(3) + 1)
                    do k1 = 1, N1
                        coeff_23(k2, k3) = spl%coeff(k1, k2, k3, &
                                                     interval_index(1) + 1, &
                                                     interval_index(2) + 1, &
                                                     interval_index(3) + 1) &
                                           + x_local(1)*coeff_23(k2, k3)
                    enddo
                enddo
            enddo

            do k3 = 0, N3
                do k2 = 0, N2
                    coeff_23_dx1(k2, k3) = spl%coeff(0, k2, k3, interval_index(1) + 1, &
                        interval_index(2) + 1, interval_index(3) + 1)*N1
                    do k1 = 1, N1-1
                        coeff_23_dx1(k2, k3) = spl%coeff(k1, k2, k3, &
                                                         interval_index(1) + 1, &
                                                         interval_index(2) + 1, &
                                                         interval_index(3) + 1) * &
                                               (N1-k1) + x_local(1)*coeff_23_dx1(k2, k3)
                    enddo
                enddo
            enddo

            coeff_3(:) = coeff_23(N2, :)
            coeff_3_dx1(:) = coeff_23_dx1(N2, :)
            do k2 = N2-1, 0, -1
                coeff_3(:) = coeff_23(k2, :) + x_local(2)*coeff_3
                coeff_3_dx1(:) = coeff_23_dx1(k2, :) + x_local(2)*coeff_3_dx1
            enddo
            coeff_3_dx2(:) = coeff_23(N2, :)*N2
            do k2 = N2-1, 1, -1
                coeff_3_dx2(:) = coeff_23(k2, :)*k2 + x_local(2)*coeff_3_dx2(:)
            enddo

            y = coeff_3(N3)
            dy(1) = coeff_3_dx1(N3)
            dy(2) = coeff_3_dx2(N3)
            do k3 = N3-1, 0, -1
                y = coeff_3(k3) + x_local(3)*y
                dy(1) = coeff_3_dx1(k3) + x_local(3)*dy(1)
                dy(2) = coeff_3_dx2(k3) + x_local(3)*dy(2)
            enddo
            dy(3) = coeff_3(N3)*N3
            do k3 = N3-1, 1, -1
                dy(3) = coeff_3(k3)*k3 + x_local(3)*dy(3)
            enddo

        end associate
    end subroutine evaluate_splines_3d_der_ref


    subroutine evaluate_splines_3d_der2_ref(spl, x, y, dy, d2y)
        type(SplineData3D), intent(in) :: spl
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: y, dy(3), d2y(6)

        real(dp) :: x_norm(3), x_local(3), xj
        real(dp) :: coeff_23(0:spl%order(2),0:spl%order(3))
        real(dp) :: coeff_23_dx1(0:spl%order(2),0:spl%order(3))
        real(dp) :: coeff_23_dx1x1(0:spl%order(2),0:spl%order(3))
        real(dp) :: coeff_3(0:spl%order(3))
        real(dp) :: coeff_3_dx1(0:spl%order(3))
        real(dp) :: coeff_3_dx2(0:spl%order(3))
        real(dp) :: coeff_3_dx1x1(0:spl%order(3))
        real(dp) :: coeff_3_dx1x2(0:spl%order(3))
        real(dp) :: coeff_3_dx2x2(0:spl%order(3))
        integer :: interval_index(3), k1, k2, k3, j

        dy = 0d0
        d2y = 0d0

        do j=1,3
            if (spl%periodic(j)) then
                xj = modulo(x(j) - spl%x_min(j), &
                    spl%h_step(j)*(spl%num_points(j)-1)) + spl%x_min(j)
            else
                xj = x(j)
            end if
            x_norm(j) = (xj - spl%x_min(j))/spl%h_step(j)
            interval_index(j) = max(0, min(spl%num_points(j)-1, int(x_norm(j))))
            x_local(j) = (x_norm(j) - dble(interval_index(j)))*spl%h_step(j)
        end do

        associate(N1 => spl%order(1), N2 => spl%order(2), N3 => spl%order(3))

            do k3 = 0, N3
                do k2 = 0, N2
                    coeff_23(k2, k3) = spl%coeff(0, k2, k3, interval_index(1) + 1, &
                        interval_index(2) + 1, interval_index(3) + 1)
                    do k1 = 1, N1
                        coeff_23(k2, k3) = spl%coeff(k1, k2, k3, &
                                                     interval_index(1) + 1, &
                                                     interval_index(2) + 1, &
                                                     interval_index(3) + 1) &
                                           + x_local(1)*coeff_23(k2, k3)
                    enddo
                enddo
            enddo

            do k3 = 0, N3
                do k2 = 0, N2
                    coeff_23_dx1(k2, k3) = spl%coeff(0, k2, k3, interval_index(1) + 1, &
                        interval_index(2) + 1, interval_index(3) + 1)*N1
                    do k1 = 1, N1-1
                        coeff_23_dx1(k2, k3) = spl%coeff(k1, k2, k3, &
                                                         interval_index(1) + 1, &
                                                         interval_index(2) + 1, &
                                                         interval_index(3) + 1) * &
                                               (N1-k1) + x_local(1)*coeff_23_dx1(k2, k3)
                    enddo
                enddo
            enddo

            do k3 = 0, N3
                do k2 = 0, N2
                    coeff_23_dx1x1(k2, k3) = spl%coeff(0, k2, k3, interval_index(1) + 1, &
                        interval_index(2) + 1, interval_index(3) + 1)*N1*(N1-1)
                    do k1 = 1, N1-2
                        coeff_23_dx1x1(k2, k3) = spl%coeff(k1, k2, k3, &
                                                           interval_index(1) + 1, &
                                                           interval_index(2) + 1, &
                                                           interval_index(3) + 1) * &
                                                 (N1-k1)*(N1-k1-1) &
                                                 + x_local(1)*coeff_23_dx1x1(k2, k3)
                    enddo
                enddo
            enddo

            coeff_3(:) = coeff_23(N2, :)
            coeff_3_dx1(:) = coeff_23_dx1(N2, :)
            coeff_3_dx1x1(:) = coeff_23_dx1x1(N2, :)
            do k2 = N2-1, 0, -1
                coeff_3(:) = coeff_23(k2, :) + x_local(2)*coeff_3
                coeff_3_dx1(:) = coeff_23_dx1(k2, :) + x_local(2)*coeff_3_dx1
                coeff_3_dx1x1(:) = coeff_23_dx1x1(k2, :) + x_local(2)*coeff_3_dx1x1
            enddo
            coeff_3_dx2(:) = coeff_23(N2, :)*N2
            coeff_3_dx1x2(:) = coeff_23_dx1(N2, :)*N2
            do k2 = N2-1, 1, -1
                coeff_3_dx2(:) = coeff_23(k2, :)*k2 + x_local(2)*coeff_3_dx2(:)
                coeff_3_dx1x2(:) = coeff_23_dx1(k2, :)*k2 + x_local(2)*coeff_3_dx1x2(:)
            enddo
            coeff_3_dx2x2(:) = coeff_23(N2, :)*N2*(N2-1)
            do k2 = N2-1, 2, -1
                coeff_3_dx2x2(:) = coeff_23(k2, :)*k2*(k2-1) + x_local(2)*coeff_3_dx2x2(:)
            enddo

            y = coeff_3(N3)
            dy(1) = coeff_3_dx1(N3)
            dy(2) = coeff_3_dx2(N3)
            d2y(1) = coeff_3_dx1x1(N3)
            d2y(2) = coeff_3_dx1x2(N3)
            d2y(4) = coeff_3_dx2x2(N3)
            do k3 = N3-1, 0, -1
                y = coeff_3(k3) + x_local(3)*y
                dy(1) = coeff_3_dx1(k3) + x_local(3)*dy(1)
                dy(2) = coeff_3_dx2(k3) + x_local(3)*dy(2)
                d2y(1) = coeff_3_dx1x1(k3) + x_local(3)*d2y(1)
                d2y(2) = coeff_3_dx1x2(k3) + x_local(3)*d2y(2)
                d2y(4) = coeff_3_dx2x2(k3) + x_local(3)*d2y(4)
            enddo
            dy(3) = coeff_3(N3)*N3
            d2y(3) = coeff_3_dx1(N3)*N3
            d2y(5) = coeff_3_dx2(N3)*N3
            do k3 = N3-1, 1, -1
                dy(3) = coeff_3(k3)*k3 + x_local(3)*dy(3)
                d2y(3) = coeff_3_dx1(k3)*k3 + x_local(3)*d2y(3)
                d2y(5) = coeff_3_dx2(k3)*k3 + x_local(3)*d2y(5)
            enddo
            d2y(6) = coeff_3(N3)*N3*(N3-1)
            do k3 = N3-1, 2, -1
                d2y(6) = coeff_3(k3)*k3*(k3-1) + x_local(3)*d2y(6)
            enddo

        end associate
    end subroutine evaluate_splines_3d_der2_ref

end module spline_eval_reference
