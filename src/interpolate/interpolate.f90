module interpolate
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use spl_three_to_five_sub, only: spl_per, spl_reg
    use batch_interpolate, only: BatchSplineData1D, BatchSplineData2D, BatchSplineData3D
    use batch_interpolate, only: construct_batch_splines_1d, &
                                 construct_batch_splines_2d, &
                                 construct_batch_splines_3d, &
                                 construct_batch_splines_1d_resident, &
                                 construct_batch_splines_1d_resident_device, &
                                 construct_batch_splines_2d_resident, &
                                 construct_batch_splines_2d_resident_device, &
                                 construct_batch_splines_3d_resident, &
                                 construct_batch_splines_3d_resident_device
    use batch_interpolate, only: destroy_batch_splines_1d, destroy_batch_splines_2d, &
                                 destroy_batch_splines_3d
    use batch_interpolate, only: evaluate_batch_splines_1d, &
                                 evaluate_batch_splines_1d_single, &
                                 evaluate_batch_splines_1d_many, &
                                 evaluate_batch_splines_1d_many_resident, &
                                 evaluate_batch_splines_1d_many_der, &
                                 evaluate_batch_splines_1d_many_der2, &
                                 evaluate_batch_splines_1d_many_der3, &
                                 evaluate_batch_splines_1d_many_der2_mask, &
                                 evaluate_batch_splines_1d_many_der3_mask, &
                                 evaluate_batch_splines_1d_der, &
                                 evaluate_batch_splines_1d_der2, &
                                 evaluate_batch_splines_1d_der3
    use batch_interpolate, only: evaluate_batch_splines_2d, &
                                 evaluate_batch_splines_2d_der, &
                                 evaluate_batch_splines_2d_many, &
                                 evaluate_batch_splines_2d_many_mask, &
                                 evaluate_batch_splines_2d_many_resident, &
                                 evaluate_batch_splines_2d_many_der, &
                                 evaluate_batch_splines_2d_many_der_mask
    use batch_interpolate, only: evaluate_batch_splines_3d, &
                                 evaluate_batch_splines_3d_der, &
                                 evaluate_batch_splines_3d_der2, &
                                 evaluate_batch_splines_3d_many, &
                                 evaluate_batch_splines_3d_many_resident, &
                                 evaluate_batch_splines_3d_many_der, &
                                 evaluate_batch_splines_3d_many_der2, &
                                 evaluate_batch_splines_3d_many_der2_mask

    implicit none

    ! Re-export batch interpolate types and procedures
    public :: BatchSplineData1D, BatchSplineData2D, BatchSplineData3D
    public :: construct_batch_splines_1d, construct_batch_splines_2d, &
              construct_batch_splines_3d
    public :: construct_batch_splines_1d_resident, construct_batch_splines_2d_resident
    public :: construct_batch_splines_1d_resident_device
    public :: construct_batch_splines_2d_resident_device
    public :: construct_batch_splines_3d_resident
    public :: construct_batch_splines_3d_resident_device
    public :: destroy_batch_splines_1d, destroy_batch_splines_2d, &
              destroy_batch_splines_3d
    public :: evaluate_batch_splines_1d, evaluate_batch_splines_1d_single
    public :: evaluate_batch_splines_1d_many, evaluate_batch_splines_1d_many_resident
    public :: evaluate_batch_splines_1d_many_der
    public :: evaluate_batch_splines_1d_many_der2, evaluate_batch_splines_1d_many_der3
    public :: evaluate_batch_splines_1d_many_der2_mask
    public :: evaluate_batch_splines_1d_many_der3_mask
    public :: evaluate_batch_splines_1d_der, evaluate_batch_splines_1d_der2
    public :: evaluate_batch_splines_1d_der3
    public :: evaluate_batch_splines_2d, evaluate_batch_splines_2d_der
    public :: evaluate_batch_splines_2d_many, evaluate_batch_splines_2d_many_mask
    public :: evaluate_batch_splines_2d_many_resident
    public :: evaluate_batch_splines_2d_many_der, &
              evaluate_batch_splines_2d_many_der_mask
    public :: evaluate_batch_splines_3d, evaluate_batch_splines_3d_der, &
              evaluate_batch_splines_3d_der2
    public :: evaluate_batch_splines_3d_many, evaluate_batch_splines_3d_many_resident
    public :: evaluate_batch_splines_3d_many_der, evaluate_batch_splines_3d_many_der2
    public :: evaluate_batch_splines_3d_many_der2_mask

    ! Single-quantity non-batch spline types and routines
    public :: SplineData1D, SplineData2D, SplineData3D
    public :: construct_splines_1d, construct_splines_2d, construct_splines_3d
    public :: destroy_splines_1d, destroy_splines_2d, destroy_splines_3d
    public :: evaluate_splines_1d, evaluate_splines_1d_der, evaluate_splines_1d_der2
    public :: evaluate_splines_1d_many, evaluate_splines_1d_many_der, &
              evaluate_splines_1d_many_der2
    public :: evaluate_splines_2d, evaluate_splines_2d_der
    public :: evaluate_splines_2d_many, evaluate_splines_2d_many_der
    public :: evaluate_splines_3d, evaluate_splines_3d_der, evaluate_splines_3d_der2
    public :: evaluate_splines_3d_many, evaluate_splines_3d_many_der, &
              evaluate_splines_3d_many_der2

    type :: SplineData1D
        integer :: order
        integer :: num_points
        logical :: periodic
        real(dp) :: x_min
        real(dp) :: h_step
        real(dp), dimension(:, :), allocatable :: coeff
    end type SplineData1D

    type :: SplineData2D
        integer :: order(2)
        integer :: num_points(2)
        logical :: periodic(2)
        real(dp) :: h_step(2)
        real(dp) :: x_min(2)
        real(dp), dimension(:, :, :, :), allocatable :: coeff
    end type SplineData2D

    type :: SplineData3D
        integer :: order(3)
        integer :: num_points(3)
        logical :: periodic(3)
        real(dp) :: h_step(3)
        real(dp) :: x_min(3)
        real(dp), dimension(:, :, :, :, :, :), allocatable :: coeff
    end type SplineData3D

    interface disp
        module procedure disp_3d
    end interface disp

contains

    subroutine construct_splines_1d(x_min, x_max, y, order, periodic, spl)
        real(dp), intent(in) :: x_min, x_max, y(:)
        integer, intent(in) :: order
        logical, intent(in) :: periodic

        type(SplineData1D), intent(out) :: spl

        spl%x_min = x_min
        spl%order = order
        spl%periodic = periodic
        spl%num_points = size(y)
        spl%h_step = (x_max - x_min)/(spl%num_points - 1)

        if (allocated(spl%coeff)) deallocate (spl%coeff)
        allocate (spl%coeff(0:order, spl%num_points))
        spl%coeff(0, :) = y

        if (periodic) then
            call spl_per(spl%order, spl%num_points, spl%h_step, spl%coeff)
        else
            call spl_reg(spl%order, spl%num_points, spl%h_step, spl%coeff)
        end if
    end subroutine construct_splines_1d

    subroutine destroy_splines_1d(spl)
        type(SplineData1D), intent(inout) :: spl

        if (allocated(spl%coeff)) deallocate (spl%coeff)
    end subroutine destroy_splines_1d

    subroutine evaluate_splines_1d(spl, x, y)
        type(SplineData1D), intent(in) :: spl
        real(dp), intent(in) :: x
        real(dp), intent(out) :: y

        real(dp) :: x_arr(1), y_arr(1)

        x_arr(1) = x
        call evaluate_splines_1d_many(spl, x_arr, y_arr)
        y = y_arr(1)
    end subroutine evaluate_splines_1d

    subroutine evaluate_splines_1d_der(spl, x, y, dy)
        type(SplineData1D), intent(in) :: spl
        real(dp), intent(in) :: x
        real(dp), intent(out) :: y, dy

        real(dp) :: x_arr(1), y_arr(1), dy_arr(1)

        x_arr(1) = x
        call evaluate_splines_1d_many_der(spl, x_arr, y_arr, dy_arr)
        y = y_arr(1)
        dy = dy_arr(1)
    end subroutine evaluate_splines_1d_der

    subroutine evaluate_splines_1d_der2(spl, x, y, dy, d2y)
        type(SplineData1D), intent(in) :: spl
        real(dp), intent(in) :: x
        real(dp), intent(out) :: y, dy, d2y

        real(dp) :: x_arr(1), y_arr(1), dy_arr(1), d2y_arr(1)

        x_arr(1) = x
        call evaluate_splines_1d_many_der2(spl, x_arr, y_arr, dy_arr, d2y_arr)
        y = y_arr(1)
        dy = dy_arr(1)
        d2y = d2y_arr(1)
    end subroutine evaluate_splines_1d_der2

    subroutine evaluate_splines_1d_many(spl, x, y)
        type(SplineData1D), intent(in) :: spl
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)

        real(dp) :: x_norm, x_local, xj
        integer :: interval_index, k_power, ie, n_eval

        n_eval = size(x)
        do ie = 1, n_eval
            if (spl%periodic) then
                xj = modulo(x(ie) - spl%x_min, spl%h_step*(spl%num_points - 1)) &
                     + spl%x_min
            else
                xj = x(ie)
            end if
            x_norm = (xj - spl%x_min)/spl%h_step
            interval_index = max(0, min(spl%num_points - 2, int(x_norm)))
            x_local = (x_norm - dble(interval_index))*spl%h_step

            y(ie) = spl%coeff(spl%order, interval_index + 1)
            do k_power = spl%order - 1, 0, -1
                y(ie) = spl%coeff(k_power, interval_index + 1) + x_local*y(ie)
            end do
        end do
    end subroutine evaluate_splines_1d_many

    subroutine evaluate_splines_1d_many_der(spl, x, y, dy)
        type(SplineData1D), intent(in) :: spl
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:), dy(:)

        real(dp) :: x_norm, x_local, xj
        integer :: interval_index, k_power, ie, n_eval

        n_eval = size(x)
        do ie = 1, n_eval
            if (spl%periodic) then
                xj = modulo(x(ie) - spl%x_min, spl%h_step*(spl%num_points - 1)) &
                     + spl%x_min
            else
                xj = x(ie)
            end if
            x_norm = (xj - spl%x_min)/spl%h_step
            interval_index = max(0, min(spl%num_points - 2, int(x_norm)))
            x_local = (x_norm - dble(interval_index))*spl%h_step

            y(ie) = spl%coeff(spl%order, interval_index + 1)
            do k_power = spl%order - 1, 0, -1
                y(ie) = spl%coeff(k_power, interval_index + 1) + x_local*y(ie)
            end do
            dy(ie) = spl%coeff(spl%order, interval_index + 1)*spl%order
            do k_power = spl%order - 1, 1, -1
                dy(ie) = spl%coeff(k_power, interval_index + 1)*k_power + x_local*dy(ie)
            end do
        end do
    end subroutine evaluate_splines_1d_many_der

    subroutine evaluate_splines_1d_many_der2(spl, x, y, dy, d2y)
        type(SplineData1D), intent(in) :: spl
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:), dy(:), d2y(:)

        real(dp) :: x_norm, x_local, xj
        integer :: interval_index, k_power, ie, n_eval

        n_eval = size(x)
        do ie = 1, n_eval
            if (spl%periodic) then
                xj = modulo(x(ie) - spl%x_min, spl%h_step*(spl%num_points - 1)) &
                     + spl%x_min
            else
                xj = x(ie)
            end if
            x_norm = (xj - spl%x_min)/spl%h_step
            interval_index = max(0, min(spl%num_points - 2, int(x_norm)))
            x_local = (x_norm - dble(interval_index))*spl%h_step

            y(ie) = spl%coeff(spl%order, interval_index + 1)
            do k_power = spl%order - 1, 0, -1
                y(ie) = spl%coeff(k_power, interval_index + 1) + x_local*y(ie)
            end do
            dy(ie) = spl%coeff(spl%order, interval_index + 1)*spl%order
            do k_power = spl%order - 1, 1, -1
                dy(ie) = spl%coeff(k_power, interval_index + 1)*k_power + x_local*dy(ie)
            end do
            d2y(ie) = spl%coeff(spl%order, interval_index + 1)*spl%order*(spl%order - 1)
            do k_power = spl%order - 1, 2, -1
                d2y(ie) = spl%coeff(k_power, interval_index + 1)*k_power*(k_power - 1) &
                          + x_local*d2y(ie)
            end do
        end do
    end subroutine evaluate_splines_1d_many_der2

    subroutine construct_splines_2d(x_min, x_max, y, order, periodic, spl)
        real(dp), intent(in) :: x_min(:), x_max(:), y(:, :)
        integer, intent(in) :: order(:)
        logical, intent(in) :: periodic(:)

        type(SplineData2D), intent(out) :: spl

        real(dp), dimension(:, :), allocatable  :: splcoe

        integer :: i1, i2  ! Loop indices for points (1 ... num_points)
        integer :: k2      ! Loop indices for polynomial order (0 ... order)

        spl%x_min = x_min
        spl%order = order
        spl%periodic = periodic
        spl%num_points = shape(y)
        spl%h_step = (x_max - x_min)/(spl%num_points - 1)

        if (allocated(spl%coeff)) deallocate (spl%coeff)
        allocate (spl%coeff(0:order(1), 0:order(2), &
                            spl%num_points(1), spl%num_points(2)))

        ! Spline over x2
        allocate (splcoe(0:spl%order(2), spl%num_points(2)))
        do i1 = 1, spl%num_points(1)
            splcoe(0, :) = y(i1, :)
            if (spl%periodic(2)) then
                call spl_per(spl%order(2), spl%num_points(2), spl%h_step(2), splcoe)
            else
                call spl_reg(spl%order(2), spl%num_points(2), spl%h_step(2), splcoe)
            end if
            spl%coeff(0, :, i1, :) = splcoe
        end do
        deallocate (splcoe)

        ! Spline over x1
        allocate (splcoe(0:spl%order(1), spl%num_points(1)))
        do i2 = 1, spl%num_points(2)
            do k2 = 0, spl%order(2)
                splcoe(0, :) = spl%coeff(0, k2, :, i2)
                if (spl%periodic(1)) then
                    call spl_per(spl%order(1), spl%num_points(1), spl%h_step(1), splcoe)
                else
                    call spl_reg(spl%order(1), spl%num_points(1), spl%h_step(1), splcoe)
                end if
                spl%coeff(:, k2, :, i2) = splcoe
            end do
        end do
        deallocate (splcoe)

    end subroutine construct_splines_2d

    subroutine evaluate_splines_2d(spl, x, y)
        type(SplineData2D), intent(in) :: spl
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: y

        real(dp) :: x_arr(2, 1), y_arr(1)

        x_arr(:, 1) = x
        call evaluate_splines_2d_many(spl, x_arr, y_arr)
        y = y_arr(1)
    end subroutine evaluate_splines_2d

    subroutine evaluate_splines_2d_der(spl, x, y, dy)
        type(SplineData2D), intent(in) :: spl
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: y
        real(dp), intent(out) :: dy(2)

        real(dp) :: x_arr(2, 1), y_arr(1), dy_arr(2, 1)

        x_arr(:, 1) = x
        call evaluate_splines_2d_many_der(spl, x_arr, y_arr, dy_arr)
        y = y_arr(1)
        dy = dy_arr(:, 1)
    end subroutine evaluate_splines_2d_der

    subroutine evaluate_splines_2d_many(spl, x, y)
        type(SplineData2D), intent(in) :: spl
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: y(:)

        real(dp) :: x_norm(2), x_local(2), xj
        real(dp) :: coeff_2(0:spl%order(2))
        integer :: interval_index(2), k1, k2, j, ie, n_eval

        n_eval = size(x, 2)
        do ie = 1, n_eval
            do j = 1, 2
                if (spl%periodic(j)) then
                    xj = modulo(x(j, ie) - spl%x_min(j), &
                                spl%h_step(j)*(spl%num_points(j) - 1)) + spl%x_min(j)
                else
                    xj = x(j, ie)
                end if
                x_norm(j) = (xj - spl%x_min(j))/spl%h_step(j)
                interval_index(j) = max(0, min(spl%num_points(j) - 2, int(x_norm(j))))
                x_local(j) = (x_norm(j) - dble(interval_index(j)))*spl%h_step(j)
            end do

            coeff_2(:) = spl%coeff(spl%order(1), 0:spl%order(2), &
                                   interval_index(1) + 1, interval_index(2) + 1)
            do k1 = spl%order(1) - 1, 0, -1
                coeff_2(:) = spl%coeff(k1, :, interval_index(1) + 1, &
                                       interval_index(2) + 1) + x_local(1)*coeff_2
            end do

            y(ie) = coeff_2(spl%order(2))
            do k2 = spl%order(2) - 1, 0, -1
                y(ie) = coeff_2(k2) + x_local(2)*y(ie)
            end do
        end do
    end subroutine evaluate_splines_2d_many

    subroutine evaluate_splines_2d_many_der(spl, x, y, dy)
        type(SplineData2D), intent(in) :: spl
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: y(:)
        real(dp), intent(out) :: dy(:, :)

        real(dp) :: x_norm(2), x_local(2), xj
        real(dp) :: coeff_2(0:spl%order(2))
        real(dp) :: coeff_2_dx1(0:spl%order(2))
        integer :: interval_index(2), k1, k2, j, ie, n_eval

        n_eval = size(x, 2)
        do ie = 1, n_eval
            do j = 1, 2
                if (spl%periodic(j)) then
                    xj = modulo(x(j, ie) - spl%x_min(j), &
                                spl%h_step(j)*(spl%num_points(j) - 1)) + spl%x_min(j)
                else
                    xj = x(j, ie)
                end if
                x_norm(j) = (xj - spl%x_min(j))/spl%h_step(j)
                interval_index(j) = max(0, min(spl%num_points(j) - 2, int(x_norm(j))))
                x_local(j) = (x_norm(j) - dble(interval_index(j)))*spl%h_step(j)
            end do

            coeff_2(:) = spl%coeff(spl%order(1), 0:spl%order(2), &
                                   interval_index(1) + 1, interval_index(2) + 1)
            do k1 = spl%order(1) - 1, 0, -1
                coeff_2(:) = spl%coeff(k1, :, interval_index(1) + 1, &
                                       interval_index(2) + 1) + x_local(1)*coeff_2
            end do

            coeff_2_dx1(:) = spl%coeff(spl%order(1), 0:spl%order(2), &
                                       interval_index(1) + 1, interval_index(2) + 1) &
                             *spl%order(1)
            do k1 = spl%order(1) - 1, 1, -1
                coeff_2_dx1(:) = spl%coeff(k1, :, interval_index(1) + 1, &
                                           interval_index(2) + 1)*k1 &
                                 + x_local(1)*coeff_2_dx1
            end do

            y(ie) = coeff_2(spl%order(2))
            do k2 = spl%order(2) - 1, 0, -1
                y(ie) = coeff_2(k2) + x_local(2)*y(ie)
            end do

            dy(1, ie) = coeff_2_dx1(spl%order(2))
            do k2 = spl%order(2) - 1, 0, -1
                dy(1, ie) = coeff_2_dx1(k2) + x_local(2)*dy(1, ie)
            end do

            dy(2, ie) = coeff_2(spl%order(2))*spl%order(2)
            do k2 = spl%order(2) - 1, 1, -1
                dy(2, ie) = coeff_2(k2)*k2 + x_local(2)*dy(2, ie)
            end do
        end do
    end subroutine evaluate_splines_2d_many_der

    subroutine destroy_splines_2d(spl)
        type(SplineData2D), intent(inout) :: spl

        if (allocated(spl%coeff)) deallocate (spl%coeff)
    end subroutine destroy_splines_2d

    subroutine construct_splines_3d(x_min, x_max, y, order, periodic, spl)
        real(dp), intent(in) :: x_min(:), x_max(:), y(:, :, :)
        integer, intent(in) :: order(3)
        logical, intent(in) :: periodic(3)

        type(SplineData3D), intent(out) :: spl

        real(dp), dimension(:, :), allocatable  :: splcoe

        integer :: i1, i2, i3  ! Loop indices for points (1 ... num_points)
        integer :: k2, k3      ! Loop indices for polynomial order (0 ... order)

        spl%x_min = x_min
        spl%order = order
        spl%periodic = periodic
        spl%num_points = shape(y)
        spl%h_step = (x_max - x_min)/(spl%num_points - 1)

        if (allocated(spl%coeff)) deallocate (spl%coeff)
        allocate (spl%coeff(0:order(1), 0:order(2), 0:order(3), &
                            spl%num_points(1), spl%num_points(2), spl%num_points(3)))

        ! Spline over x3
        allocate (splcoe(0:spl%order(3), spl%num_points(3)))
        do i2 = 1, spl%num_points(2)
        do i1 = 1, spl%num_points(1)
            splcoe(0, :) = y(i1, i2, :)
            if (spl%periodic(3)) then
                call spl_per(spl%order(3), spl%num_points(3), spl%h_step(3), splcoe)
            else
                call spl_reg(spl%order(3), spl%num_points(3), spl%h_step(3), splcoe)
            end if
            spl%coeff(order(1), 0, :, i1, i2, :) = splcoe
        end do
        end do
        deallocate (splcoe)

        ! Spline over x2
        allocate (splcoe(0:spl%order(2), spl%num_points(2)))
        do i3 = 1, spl%num_points(3)
        do i1 = 1, spl%num_points(1)
            do k3 = 0, spl%order(3)
                splcoe(0, :) = spl%coeff(order(1), 0, k3, i1, :, i3)
                if (spl%periodic(2)) then
                    call spl_per( &
                        spl%order(2), spl%num_points(2), spl%h_step(2), splcoe)
                else
                    call spl_reg( &
                        spl%order(2), spl%num_points(2), spl%h_step(2), splcoe)
                end if
                spl%coeff(order(1), :, k3, i1, :, i3) = splcoe
            end do
        end do
        end do
        deallocate (splcoe)

        ! Spline over x1
        allocate (splcoe(0:spl%order(1), spl%num_points(1)))
        do i3 = 1, spl%num_points(3)
        do i2 = 1, spl%num_points(2)
            do k3 = 0, spl%order(3)
            do k2 = 0, spl%order(2)
                splcoe(0, :) = spl%coeff(order(1), k2, k3, :, i2, i3)
                if (spl%periodic(1)) then
                    call spl_per( &
                        spl%order(1), spl%num_points(1), spl%h_step(1), splcoe)
                else
                    call spl_reg( &
                        spl%order(1), spl%num_points(1), spl%h_step(1), splcoe)
                end if
                spl%coeff(order(1):0:-1, k2, k3, :, i2, i3) = splcoe
            end do
            end do
        end do
        end do
        deallocate (splcoe)
    end subroutine construct_splines_3d

    subroutine evaluate_splines_3d(spl, x, y)
        type(SplineData3D), intent(in) :: spl
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: y

        real(dp) :: x_arr(3, 1), y_arr(1)

        x_arr(:, 1) = x
        call evaluate_splines_3d_many(spl, x_arr, y_arr)
        y = y_arr(1)
    end subroutine evaluate_splines_3d

    subroutine evaluate_splines_3d_der(spl, x, y, dy)
        type(SplineData3D), intent(in) :: spl
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: y, dy(3)

        real(dp) :: x_arr(3, 1), y_arr(1), dy_arr(3, 1)

        x_arr(:, 1) = x
        call evaluate_splines_3d_many_der(spl, x_arr, y_arr, dy_arr)
        y = y_arr(1)
        dy = dy_arr(:, 1)
    end subroutine evaluate_splines_3d_der

    subroutine evaluate_splines_3d_der2(spl, x, y, dy, d2y)
        type(SplineData3D), intent(in) :: spl
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: y, dy(3), d2y(6)

        real(dp) :: x_arr(3, 1), y_arr(1), dy_arr(3, 1), d2y_arr(6, 1)

        x_arr(:, 1) = x
        call evaluate_splines_3d_many_der2(spl, x_arr, y_arr, dy_arr, d2y_arr)
        y = y_arr(1)
        dy = dy_arr(:, 1)
        d2y = d2y_arr(:, 1)
    end subroutine evaluate_splines_3d_der2

    subroutine evaluate_splines_3d_many(spl, x, y)
        type(SplineData3D), intent(in) :: spl
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: y(:)

        real(dp) :: x_norm(3), x_local(3), xj
        real(dp) :: coeff_3(0:spl%order(3))
        real(dp) :: coeff_23(0:spl%order(2), 0:spl%order(3))
        integer :: interval_index(3), k1, k2, k3, j, ie, n_eval

        n_eval = size(x, 2)
        do ie = 1, n_eval
            do j = 1, 3
                if (spl%periodic(j)) then
                    xj = modulo(x(j, ie) - spl%x_min(j), &
                                spl%h_step(j)*(spl%num_points(j) - 1)) + spl%x_min(j)
                else
                    xj = x(j, ie)
                end if
                x_norm(j) = (xj - spl%x_min(j))/spl%h_step(j)
                interval_index(j) = max(0, min(spl%num_points(j) - 2, int(x_norm(j))))
                x_local(j) = (x_norm(j) - dble(interval_index(j)))*spl%h_step(j)
            end do

            coeff_23(:, :) = spl%coeff(0, :, :, interval_index(1) + 1, &
                                       interval_index(2) + 1, interval_index(3) + 1)
            do k1 = 1, spl%order(1)
                coeff_23(:, :) = spl%coeff(k1, :, :, interval_index(1) + 1, &
                                           interval_index(2) + 1, &
                                           interval_index(3) + 1) &
                                 + x_local(1)*coeff_23(:, :)
            end do

            coeff_3(:) = coeff_23(spl%order(2), :)
            do k2 = spl%order(2) - 1, 0, -1
                coeff_3(:) = coeff_23(k2, :) + x_local(2)*coeff_3
            end do

            y(ie) = coeff_3(spl%order(3))
            do k3 = spl%order(3) - 1, 0, -1
                y(ie) = coeff_3(k3) + x_local(3)*y(ie)
            end do
        end do
    end subroutine evaluate_splines_3d_many

    subroutine evaluate_splines_3d_many_der(spl, x, y, dy)
        type(SplineData3D), intent(in) :: spl
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: y(:), dy(:, :)

        real(dp) :: x_norm(3), x_local(3), xj
        real(dp) :: coeff_23(0:spl%order(2), 0:spl%order(3))
        real(dp) :: coeff_23_dx1(0:spl%order(2), 0:spl%order(3))
        real(dp) :: coeff_3(0:spl%order(3))
        real(dp) :: coeff_3_dx1(0:spl%order(3))
        real(dp) :: coeff_3_dx2(0:spl%order(3))
        integer :: interval_index(3), k1, k2, k3, j, ie, n_eval

        n_eval = size(x, 2)
        do ie = 1, n_eval
            do j = 1, 3
                if (spl%periodic(j)) then
                    xj = modulo(x(j, ie) - spl%x_min(j), &
                                spl%h_step(j)*(spl%num_points(j) - 1)) + spl%x_min(j)
                else
                    xj = x(j, ie)
                end if
                x_norm(j) = (xj - spl%x_min(j))/spl%h_step(j)
                interval_index(j) = max(0, min(spl%num_points(j) - 2, int(x_norm(j))))
                x_local(j) = (x_norm(j) - dble(interval_index(j)))*spl%h_step(j)
            end do

            associate (N1 => spl%order(1), N2 => spl%order(2), N3 => spl%order(3))
                do k3 = 0, N3
                    do k2 = 0, N2
                        coeff_23(k2, k3) = spl%coeff(0, k2, k3, interval_index(1) + 1, &
                                                     interval_index(2) + 1, &
                                                         interval_index(3) + 1)
                        do k1 = 1, N1
                            coeff_23(k2, k3) = spl%coeff(k1, k2, k3, &
                                                         interval_index(1) + 1, &
                                                             interval_index(2) + 1, &
                                                         interval_index(3) + 1) + &
                                                             x_local(1)*coeff_23(k2, &
                                                             k3)
                        end do
                    end do
                end do

                do k3 = 0, N3
                    do k2 = 0, N2
                        coeff_23_dx1(k2, k3) = spl%coeff(0, k2, k3, &
                                                         interval_index(1) + 1, &
                                                         interval_index(2) + 1, &
                                                             interval_index(3) + 1)*N1
                        do k1 = 1, N1 - 1
                            coeff_23_dx1(k2, k3) = spl%coeff(k1, k2, k3, &
                                         interval_index(1) + 1, interval_index(2) + 1, &
                                                             interval_index(3) + &
                                                                 1)*(N1 - k1) + &
                                                   x_local(1)*coeff_23_dx1(k2, k3)
                        end do
                    end do
                end do

                coeff_3(:) = coeff_23(N2, :)
                coeff_3_dx1(:) = coeff_23_dx1(N2, :)
                do k2 = N2 - 1, 0, -1
                    coeff_3(:) = coeff_23(k2, :) + x_local(2)*coeff_3
                    coeff_3_dx1(:) = coeff_23_dx1(k2, :) + x_local(2)*coeff_3_dx1
                end do
                coeff_3_dx2(:) = coeff_23(N2, :)*N2
                do k2 = N2 - 1, 1, -1
                    coeff_3_dx2(:) = coeff_23(k2, :)*k2 + x_local(2)*coeff_3_dx2(:)
                end do

                y(ie) = coeff_3(N3)
                dy(1, ie) = coeff_3_dx1(N3)
                dy(2, ie) = coeff_3_dx2(N3)
                do k3 = N3 - 1, 0, -1
                    y(ie) = coeff_3(k3) + x_local(3)*y(ie)
                    dy(1, ie) = coeff_3_dx1(k3) + x_local(3)*dy(1, ie)
                    dy(2, ie) = coeff_3_dx2(k3) + x_local(3)*dy(2, ie)
                end do
                dy(3, ie) = coeff_3(N3)*N3
                do k3 = N3 - 1, 1, -1
                    dy(3, ie) = coeff_3(k3)*k3 + x_local(3)*dy(3, ie)
                end do
            end associate
        end do
    end subroutine evaluate_splines_3d_many_der

    subroutine evaluate_splines_3d_many_der2(spl, x, y, dy, d2y)
        type(SplineData3D), intent(in) :: spl
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: y(:), dy(:, :), d2y(:, :)

        real(dp) :: x_norm(3), x_local(3), xj
        real(dp) :: coeff_23(0:spl%order(2), 0:spl%order(3))
        real(dp) :: coeff_23_dx1(0:spl%order(2), 0:spl%order(3))
        real(dp) :: coeff_23_dx1x1(0:spl%order(2), 0:spl%order(3))
        real(dp) :: coeff_3(0:spl%order(3))
        real(dp) :: coeff_3_dx1(0:spl%order(3))
        real(dp) :: coeff_3_dx2(0:spl%order(3))
        real(dp) :: coeff_3_dx1x1(0:spl%order(3))
        real(dp) :: coeff_3_dx1x2(0:spl%order(3))
        real(dp) :: coeff_3_dx2x2(0:spl%order(3))
        integer :: interval_index(3), k1, k2, k3, j, ie, n_eval

        n_eval = size(x, 2)
        do ie = 1, n_eval
            do j = 1, 3
                if (spl%periodic(j)) then
                    xj = modulo(x(j, ie) - spl%x_min(j), &
                                spl%h_step(j)*(spl%num_points(j) - 1)) + spl%x_min(j)
                else
                    xj = x(j, ie)
                end if
                x_norm(j) = (xj - spl%x_min(j))/spl%h_step(j)
                interval_index(j) = max(0, min(spl%num_points(j) - 2, int(x_norm(j))))
                x_local(j) = (x_norm(j) - dble(interval_index(j)))*spl%h_step(j)
            end do

            associate (N1 => spl%order(1), N2 => spl%order(2), N3 => spl%order(3))
                do k3 = 0, N3
                    do k2 = 0, N2
                        coeff_23(k2, k3) = spl%coeff(0, k2, k3, interval_index(1) + 1, &
                                                     interval_index(2) + 1, &
                                                         interval_index(3) + 1)
                        do k1 = 1, N1
                            coeff_23(k2, k3) = spl%coeff(k1, k2, k3, &
                                                         interval_index(1) + 1, &
                                                             interval_index(2) + 1, &
                                                         interval_index(3) + 1) + &
                                                             x_local(1)*coeff_23(k2, &
                                                             k3)
                        end do
                    end do
                end do

                do k3 = 0, N3
                    do k2 = 0, N2
                        coeff_23_dx1(k2, k3) = spl%coeff(0, k2, k3, &
                                                         interval_index(1) + 1, &
                                                         interval_index(2) + 1, &
                                                             interval_index(3) + 1)*N1
                        do k1 = 1, N1 - 1
                            coeff_23_dx1(k2, k3) = spl%coeff(k1, k2, k3, &
                                         interval_index(1) + 1, interval_index(2) + 1, &
                                                             interval_index(3) + &
                                                                 1)*(N1 - k1) + &
                                                   x_local(1)*coeff_23_dx1(k2, k3)
                        end do
                    end do
                end do

                do k3 = 0, N3
                    do k2 = 0, N2
                        coeff_23_dx1x1(k2, k3) = spl%coeff(0, k2, k3, &
                                                           interval_index(1) + 1, &
                                                           interval_index(2) + 1, &
                                                               interval_index(3) + &
                                                               1)*N1*(N1 - 1)
                        do k1 = 1, N1 - 2
                            coeff_23_dx1x1(k2, k3) = spl%coeff(k1, k2, k3, &
                                         interval_index(1) + 1, interval_index(2) + 1, &
                                                               interval_index(3) + &
                                                                   1)*(N1 - k1)*(N1 &
                                                                   - k1 - 1) &
                                                     + x_local(1)*coeff_23_dx1x1(k2, k3)
                        end do
                    end do
                end do

                coeff_3(:) = coeff_23(N2, :)
                coeff_3_dx1(:) = coeff_23_dx1(N2, :)
                coeff_3_dx1x1(:) = coeff_23_dx1x1(N2, :)
                do k2 = N2 - 1, 0, -1
                    coeff_3(:) = coeff_23(k2, :) + x_local(2)*coeff_3
                    coeff_3_dx1(:) = coeff_23_dx1(k2, :) + x_local(2)*coeff_3_dx1
                    coeff_3_dx1x1(:) = coeff_23_dx1x1(k2, :) + x_local(2)*coeff_3_dx1x1
                end do
                coeff_3_dx2(:) = coeff_23(N2, :)*N2
                coeff_3_dx1x2(:) = coeff_23_dx1(N2, :)*N2
                do k2 = N2 - 1, 1, -1
                    coeff_3_dx2(:) = coeff_23(k2, :)*k2 + x_local(2)*coeff_3_dx2(:)
                    coeff_3_dx1x2(:) = coeff_23_dx1(k2, :)*k2 + x_local(2)*coeff_3_dx1x2
                end do
                coeff_3_dx2x2(:) = coeff_23(N2, :)*N2*(N2 - 1)
                do k2 = N2 - 1, 2, -1
                    coeff_3_dx2x2(:) = coeff_23(k2, :)*k2*(k2 - 1) + &
                                       x_local(2)*coeff_3_dx2x2
                end do

                y(ie) = coeff_3(N3)
                dy(1, ie) = coeff_3_dx1(N3)
                dy(2, ie) = coeff_3_dx2(N3)
                d2y(1, ie) = coeff_3_dx1x1(N3)
                d2y(2, ie) = coeff_3_dx1x2(N3)
                d2y(4, ie) = coeff_3_dx2x2(N3)
                do k3 = N3 - 1, 0, -1
                    y(ie) = coeff_3(k3) + x_local(3)*y(ie)
                    dy(1, ie) = coeff_3_dx1(k3) + x_local(3)*dy(1, ie)
                    dy(2, ie) = coeff_3_dx2(k3) + x_local(3)*dy(2, ie)
                    d2y(1, ie) = coeff_3_dx1x1(k3) + x_local(3)*d2y(1, ie)
                    d2y(2, ie) = coeff_3_dx1x2(k3) + x_local(3)*d2y(2, ie)
                    d2y(4, ie) = coeff_3_dx2x2(k3) + x_local(3)*d2y(4, ie)
                end do
                dy(3, ie) = coeff_3(N3)*N3
                d2y(3, ie) = coeff_3_dx1(N3)*N3
                d2y(5, ie) = coeff_3_dx2(N3)*N3
                do k3 = N3 - 1, 1, -1
                    dy(3, ie) = coeff_3(k3)*k3 + x_local(3)*dy(3, ie)
                    d2y(3, ie) = coeff_3_dx1(k3)*k3 + x_local(3)*d2y(3, ie)
                    d2y(5, ie) = coeff_3_dx2(k3)*k3 + x_local(3)*d2y(5, ie)
                end do
                d2y(6, ie) = coeff_3(N3)*N3*(N3 - 1)
                do k3 = N3 - 1, 2, -1
                    d2y(6, ie) = coeff_3(k3)*k3*(k3 - 1) + x_local(3)*d2y(6, ie)
                end do
            end associate
        end do
    end subroutine evaluate_splines_3d_many_der2

    subroutine destroy_splines_3d(spl)
        type(SplineData3D), intent(inout) :: spl

        if (allocated(spl%coeff)) deallocate (spl%coeff)
    end subroutine destroy_splines_3d

    subroutine disp_3d(spl)
        type(SplineData3D), intent(in) :: spl
        print *, "SplineData3D"
        print *, "  order = ", spl%order
        print *, "  num_points = ", spl%num_points
        print *, "  periodic = ", spl%periodic
        print *, "  x_min = ", spl%x_min
        print *, "  x_max = ", spl%x_min + (spl%num_points - 1)*spl%h_step
        print *, "  h_step = ", spl%h_step
    end subroutine disp_3d

end module interpolate
