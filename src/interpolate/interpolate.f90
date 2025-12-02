module interpolate
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use spl_three_to_five_sub, only: spl_reg, spl_per
    use batch_interpolate, only: BatchSplineData1D, BatchSplineData2D, &
        BatchSplineData3D, construct_batch_splines_1d, &
        construct_batch_splines_2d, construct_batch_splines_3d, &
        destroy_batch_splines_1d, destroy_batch_splines_2d, &
        destroy_batch_splines_3d, evaluate_batch_splines_1d, &
        evaluate_batch_splines_1d_single, evaluate_batch_splines_1d_der, &
        evaluate_batch_splines_1d_der2, evaluate_batch_splines_2d, &
        evaluate_batch_splines_2d_der, evaluate_batch_splines_3d, &
        evaluate_batch_splines_3d_der, evaluate_batch_splines_3d_der2
    
    implicit none
    
    ! Re-export batch interpolate types and procedures
    public :: BatchSplineData1D, BatchSplineData2D, BatchSplineData3D
    public :: construct_batch_splines_1d, construct_batch_splines_2d, construct_batch_splines_3d
    public :: destroy_batch_splines_1d, destroy_batch_splines_2d, destroy_batch_splines_3d
    public :: evaluate_batch_splines_1d, evaluate_batch_splines_1d_single
    public :: evaluate_batch_splines_1d_der, evaluate_batch_splines_1d_der2
    public :: evaluate_batch_splines_2d, evaluate_batch_splines_2d_der
    public :: evaluate_batch_splines_3d, evaluate_batch_splines_3d_der, evaluate_batch_splines_3d_der2

    type :: SplineData1D
        integer :: order
        integer :: num_points
        logical :: periodic
        logical :: has_mask = .false.
        real(dp) :: x_min
        real(dp) :: h_step
        real(dp), dimension(:,:), allocatable :: coeff
        logical, dimension(:), allocatable :: valid_interval
    end type SplineData1D

    type :: SplineData2D
        integer :: order(2)
        integer :: num_points(2)
        logical :: periodic(2)
        real(dp) :: h_step(2)
        real(dp) :: x_min(2)
        real(dp), dimension(:,:,:,:), allocatable :: coeff
        logical :: has_mask = .false.
        logical, dimension(:,:), allocatable :: valid_cell
    end type SplineData2D

    type :: SplineData3D
        integer :: order(3)
        integer :: num_points(3)
        logical :: periodic(3)
        real(dp) :: h_step(3)
        real(dp) :: x_min(3)
        real(dp), dimension(:,:,:,:,:,:), allocatable :: coeff
        logical :: has_mask = .false.
        logical, dimension(:,:,:), allocatable :: valid_cell
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
        spl%h_step = (x_max - x_min) / (spl%num_points - 1)

        if(allocated(spl%coeff)) deallocate(spl%coeff)
        allocate(spl%coeff(0:order, spl%num_points))
        spl%coeff(0,:) = y

        if (periodic) then
            call spl_per(spl%order, spl%num_points, spl%h_step, spl%coeff)
        else
            call spl_reg(spl%order, spl%num_points, spl%h_step, spl%coeff)
        endif
    end subroutine construct_splines_1d


    subroutine construct_splines_1d_masked(x_min, x_max, y, order, periodic, &
                                           mask_nodes, spl)
        real(dp), intent(in) :: x_min, x_max, y(:)
        integer, intent(in) :: order
        logical, intent(in) :: periodic
        logical, intent(in) :: mask_nodes(:)
        type(SplineData1D), intent(out) :: spl

        integer :: i

        call construct_splines_1d(x_min, x_max, y, order, periodic, spl)

        if (size(mask_nodes) /= spl%num_points) then
            error stop "construct_splines_1d_masked: mask size mismatch"
        end if

        if (allocated(spl%valid_interval)) deallocate(spl%valid_interval)
        allocate(spl%valid_interval(spl%num_points - 1))

        do i = 1, spl%num_points - 1
            spl%valid_interval(i) = mask_nodes(i) .and. mask_nodes(i + 1)
        end do
        spl%has_mask = .true.
    end subroutine construct_splines_1d_masked


    subroutine destroy_splines_1d(spl)
        type(SplineData1D), intent(inout) :: spl

        if(allocated(spl%coeff)) deallocate(spl%coeff)
        if(allocated(spl%valid_interval)) deallocate(spl%valid_interval)
        spl%has_mask = .false.
    end subroutine destroy_splines_1d


    subroutine evaluate_splines_1d(spl, x, y)
        type(SplineData1D), intent(in) :: spl
        real(dp), intent(in) :: x
        real(dp), intent(out) :: y

        real(dp) :: x_norm, x_local, coeff_local(0:spl%order), xj
        integer :: interval_index, k_power

        if (spl%periodic) then
            xj = modulo(x, spl%h_step*(spl%num_points-1))
        else
            xj = x
        end if
        x_norm = (xj - spl%x_min) / spl%h_step
        interval_index = max(0, min(spl%num_points-1, int(x_norm)))
        x_local = (x_norm - dble(interval_index))*spl%h_step  ! Distance to grid point

        coeff_local(:) = spl%coeff(:, interval_index+1)

        ! Start with largest power and then multiply recursively
        y = coeff_local(spl%order)
        do k_power = spl%order-1, 0, -1
            y = coeff_local(k_power) + x_local*y
        enddo
    end subroutine evaluate_splines_1d


    subroutine evaluate_splines_1d_masked(spl, x, y, is_valid)
        type(SplineData1D), intent(in) :: spl
        real(dp), intent(in) :: x
        real(dp), intent(out) :: y
        logical, intent(out), optional :: is_valid

        real(dp) :: x_norm, xj
        integer :: interval_index, cell_index

        if (.not. spl%has_mask) then
            call evaluate_splines_1d(spl, x, y)
            if (present(is_valid)) is_valid = .true.
            return
        end if

        if (spl%periodic) then
            xj = modulo(x, spl%h_step*(spl%num_points - 1))
        else
            xj = x
        end if
        x_norm = (xj - spl%x_min)/spl%h_step
        interval_index = max(0, min(spl%num_points - 1, int(x_norm)))

        cell_index = max(1, min(spl%num_points - 1, interval_index + 1))

        if (.not. spl%valid_interval(cell_index)) then
            y = 0.0_dp
            if (present(is_valid)) is_valid = .false.
            return
        end if

        call evaluate_splines_1d(spl, x, y)
        if (present(is_valid)) is_valid = .true.
    end subroutine evaluate_splines_1d_masked


    subroutine evaluate_splines_1d_der(spl, x, y, dy)
        type(SplineData1D), intent(in) :: spl
        real(dp), intent(in) :: x
        real(dp), intent(out) :: y, dy

        real(dp) :: x_norm, x_local, coeff_local(0:spl%order)
        integer :: interval_index, k_power

        x_norm = (x - spl%x_min) / spl%h_step
        interval_index = max(0, min(spl%num_points-1, int(x_norm)))
        x_local = (x_norm - dble(interval_index))*spl%h_step

        coeff_local(:) = spl%coeff(:, interval_index+1)

        y = coeff_local(spl%order)
        do k_power = spl%order-1, 0, -1
            y = coeff_local(k_power) + x_local*y
        enddo
        dy = coeff_local(spl%order)*spl%order
        do k_power = spl%order-1, 1, -1
            dy = coeff_local(k_power)*k_power + x_local*dy
        enddo
    end subroutine evaluate_splines_1d_der


    subroutine evaluate_splines_1d_der2(spl, x, y, dy, d2y)
        type(SplineData1D), intent(in) :: spl
        real(dp), intent(in) :: x
        real(dp), intent(out) :: y, dy, d2y

        real(dp) :: x_norm, x_local, coeff_local(0:spl%order)
        integer :: interval_index, k_power

        x_norm = (x - spl%x_min) / spl%h_step
        interval_index = max(0, min(spl%num_points-1, int(x_norm)))
        x_local = (x_norm - dble(interval_index))*spl%h_step

        coeff_local(:) = spl%coeff(:, interval_index+1)

        y = coeff_local(spl%order)
        do k_power = spl%order-1, 0, -1
            y = coeff_local(k_power) + x_local*y
        enddo
        dy = coeff_local(spl%order)*spl%order
        do k_power = spl%order-1, 1, -1
            dy = coeff_local(k_power)*k_power + x_local*dy
        enddo
        d2y = coeff_local(spl%order)*spl%order*(spl%order-1)
        do k_power = spl%order-1, 2, -1
            d2y = coeff_local(k_power)*k_power*(k_power-1) + x_local*d2y
        enddo
    end subroutine evaluate_splines_1d_der2


    subroutine construct_splines_2d(x_min, x_max, y, order, periodic, spl)
        real(dp), intent(in) :: x_min(:), x_max(:), y(:,:)
        integer, intent(in) :: order(:)
        logical, intent(in) :: periodic(:)

        type(SplineData2D), intent(out) :: spl

        real(dp), dimension(:,:), allocatable  :: splcoe

        integer :: i1, i2  ! Loop indices for points (1 ... num_points)
        integer :: k2      ! Loop indices for polynomial order (0 ... order)

        spl%x_min = x_min
        spl%order = order
        spl%periodic = periodic
        spl%num_points = shape(y)
        spl%h_step = (x_max - x_min) / (spl%num_points - 1)

        if(allocated(spl%coeff)) deallocate(spl%coeff)
        allocate(spl%coeff(0:order(1), 0:order(2), &
                 spl%num_points(1), spl%num_points(2)))

        ! Spline over x2
        allocate(splcoe(0:spl%order(2), spl%num_points(2)))
        do i1=1,spl%num_points(1)
            splcoe(0,:) = y(i1, :)
            if (spl%periodic(2)) then
                call spl_per(spl%order(2), spl%num_points(2), spl%h_step(2), splcoe)
            else
                call spl_reg(spl%order(2), spl%num_points(2), spl%h_step(2), splcoe)
            endif
            spl%coeff(0, :, i1, :) = splcoe
        enddo
        deallocate(splcoe)

        ! Spline over x1
        allocate(splcoe(0:spl%order(1), spl%num_points(1)))
        do i2=1,spl%num_points(2)
            do k2=0,spl%order(2)
                splcoe(0,:) = spl%coeff(0, k2, :, i2)
                if(spl%periodic(1)) then
                    call spl_per(spl%order(1), spl%num_points(1), spl%h_step(1), splcoe)
                else
                    call spl_reg(spl%order(1), spl%num_points(1), spl%h_step(1), splcoe)
                endif
                spl%coeff(:, k2, :, i2) = splcoe
            enddo
        enddo
        deallocate(splcoe)

    end subroutine construct_splines_2d


    subroutine construct_splines_2d_masked(x_min, x_max, y, order, periodic, &
                                           mask_nodes, spl)
        real(dp), intent(in) :: x_min(:), x_max(:), y(:,:)
        integer, intent(in) :: order(:)
        logical, intent(in) :: periodic(:)
        logical, intent(in) :: mask_nodes(:,:)
        type(SplineData2D), intent(out) :: spl

        integer :: n1, n2, i1, i2

        call construct_splines_2d(x_min, x_max, y, order, periodic, spl)

        if (any(shape(mask_nodes) /= spl%num_points)) then
            error stop "construct_splines_2d_masked: mask size mismatch"
        end if

        n1 = spl%num_points(1)
        n2 = spl%num_points(2)

        if (allocated(spl%valid_cell)) deallocate(spl%valid_cell)
        if (n1 > 1 .and. n2 > 1) then
            allocate(spl%valid_cell(n1 - 1, n2 - 1))
            do i2 = 1, n2 - 1
                do i1 = 1, n1 - 1
                    spl%valid_cell(i1, i2) = mask_nodes(i1, i2) .and. &
                        mask_nodes(i1 + 1, i2) .and. &
                        mask_nodes(i1, i2 + 1) .and. &
                        mask_nodes(i1 + 1, i2 + 1)
                end do
            end do
            spl%has_mask = .true.
        else
            spl%has_mask = .false.
        end if
    end subroutine construct_splines_2d_masked


    subroutine evaluate_splines_2d(spl, x, y)
        type(SplineData2D), intent(in) :: spl
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: y

        real(dp) :: x_norm(2), x_local(2), xj
        real(dp) :: coeff_2(0:spl%order(2)), &
                    coeff_local(0:spl%order(1),0:spl%order(2))
        integer :: interval_index(2), k1, k2, j

        do j=1,2
            if (spl%periodic(j)) then
                xj = modulo(x(j), spl%h_step(j)*(spl%num_points(j)-1))
            else
                xj = x(j)
            end if
            x_norm(j) = (xj - spl%x_min(j))/spl%h_step(j)
            interval_index(j) = max(0, min(spl%num_points(j)-1, int(x_norm(j))))
            x_local(j) = (x_norm(j) - dble(interval_index(j)))*spl%h_step(j)
        end do

        coeff_local(:,:) = &
            spl%coeff(:, :, interval_index(1) + 1, interval_index(2) + 1)

        coeff_2(:) = coeff_local(spl%order(1), 0:spl%order(2))
        do k1 = spl%order(1)-1, 0, -1
            coeff_2(:) = coeff_local(k1, :) + x_local(1)*coeff_2
        enddo

        y = coeff_2(spl%order(2))
        do k2 = spl%order(2)-1, 0, -1
            y = coeff_2(k2) + x_local(2)*y
        enddo
    end subroutine evaluate_splines_2d


    subroutine evaluate_splines_2d_masked(spl, x, y, is_valid)
        type(SplineData2D), intent(in) :: spl
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: y
        logical, intent(out), optional :: is_valid

        real(dp) :: x_norm(2), xj
        integer :: interval_index(2), cell_index(2), j

        if (.not. spl%has_mask) then
            call evaluate_splines_2d(spl, x, y)
            if (present(is_valid)) is_valid = .true.
            return
        end if

        do j = 1, 2
            if (spl%periodic(j)) then
                xj = modulo(x(j), spl%h_step(j)*(spl%num_points(j) - 1))
            else
                xj = x(j)
            end if
            x_norm(j) = (xj - spl%x_min(j))/spl%h_step(j)
            interval_index(j) = max(0, min(spl%num_points(j) - 1, &
                int(x_norm(j))))
            cell_index(j) = max(1, min(spl%num_points(j) - 1, &
                interval_index(j) + 1))
        end do

        if (.not. spl%valid_cell(cell_index(1), cell_index(2))) then
            y = 0.0_dp
            if (present(is_valid)) is_valid = .false.
            return
        end if

        call evaluate_splines_2d(spl, x, y)
        if (present(is_valid)) is_valid = .true.
    end subroutine evaluate_splines_2d_masked


    subroutine evaluate_splines_2d_der(spl, x, y, dy)
        type(SplineData2D), intent(in) :: spl
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: y
        real(dp), intent(out) :: dy(2)
        
        real(dp) :: x_norm(2), x_local(2), xj
        real(dp) :: coeff_2(0:spl%order(2))
        real(dp) :: coeff_2_dx1(0:spl%order(2))
        real(dp) :: coeff_local(0:spl%order(1),0:spl%order(2))
        integer :: interval_index(2), k1, k2, j
        
        do j=1,2
            if (spl%periodic(j)) then
                xj = modulo(x(j), spl%h_step(j)*(spl%num_points(j)-1))
            else
                xj = x(j)
            end if
            x_norm(j) = (xj - spl%x_min(j))/spl%h_step(j)
            interval_index(j) = max(0, min(spl%num_points(j)-1, int(x_norm(j))))
            x_local(j) = (x_norm(j) - dble(interval_index(j)))*spl%h_step(j)
        end do
        
        coeff_local(:,:) = &
            spl%coeff(:, :, interval_index(1) + 1, interval_index(2) + 1)
        
        ! Interpolation over x1
        coeff_2(:) = coeff_local(spl%order(1), 0:spl%order(2))
        do k1 = spl%order(1)-1, 0, -1
            coeff_2(:) = coeff_local(k1, :) + x_local(1)*coeff_2
        enddo
        
        ! Derivative over x1
        coeff_2_dx1(:) = coeff_local(spl%order(1), 0:spl%order(2))*spl%order(1)
        do k1 = spl%order(1)-1, 1, -1
            coeff_2_dx1(:) = coeff_local(k1, :)*k1 + x_local(1)*coeff_2_dx1
        enddo
        
        ! Interpolation over x2
        y = coeff_2(spl%order(2))
        do k2 = spl%order(2)-1, 0, -1
            y = coeff_2(k2) + x_local(2)*y
        enddo
        
        ! Derivative w.r.t. x1
        dy(1) = coeff_2_dx1(spl%order(2))
        do k2 = spl%order(2)-1, 0, -1
            dy(1) = coeff_2_dx1(k2) + x_local(2)*dy(1)
        enddo
        
        ! Derivative w.r.t. x2
        dy(2) = coeff_2(spl%order(2))*spl%order(2)
        do k2 = spl%order(2)-1, 1, -1
            dy(2) = coeff_2(k2)*k2 + x_local(2)*dy(2)
        enddo
        
    end subroutine evaluate_splines_2d_der
    
    
    subroutine destroy_splines_2d(spl)
        type(SplineData2D), intent(inout) :: spl

        if(allocated(spl%coeff)) deallocate(spl%coeff)
        if(allocated(spl%valid_cell)) deallocate(spl%valid_cell)
        spl%has_mask = .false.
    end subroutine destroy_splines_2d


    subroutine construct_splines_3d(x_min, x_max, y, order, periodic, spl)
        real(dp), intent(in) :: x_min(:), x_max(:), y(:,:,:)
        integer, intent(in) :: order(3)
        logical, intent(in) :: periodic(3)

        type(SplineData3D), intent(out) :: spl

        real(dp), dimension(:,:), allocatable  :: splcoe

        integer :: i1, i2, i3  ! Loop indices for points (1 ... num_points)
        integer :: k2, k3      ! Loop indices for polynomial order (0 ... order)

        spl%x_min = x_min
        spl%order = order
        spl%periodic = periodic
        spl%num_points = shape(y)
        spl%h_step = (x_max - x_min) / (spl%num_points - 1)

        if(allocated(spl%coeff)) deallocate(spl%coeff)
        allocate(spl%coeff(0:order(1), 0:order(2), 0:order(3), &
                 spl%num_points(1), spl%num_points(2), spl%num_points(3)))

        ! Spline over x3
        allocate(splcoe(0:spl%order(3), spl%num_points(3)))
        do i2=1,spl%num_points(2)
        do i1=1,spl%num_points(1)
            splcoe(0,:) = y(i1, i2, :)
            if (spl%periodic(3)) then
                call spl_per(spl%order(3), spl%num_points(3), spl%h_step(3), splcoe)
            else
                call spl_reg(spl%order(3), spl%num_points(3), spl%h_step(3), splcoe)
            endif
            spl%coeff(order(1), 0, :, i1, i2, :) = splcoe
        enddo
        enddo
        deallocate(splcoe)

        ! Spline over x2
        allocate(splcoe(0:spl%order(2), spl%num_points(2)))
        do i3=1,spl%num_points(3)
        do i1=1,spl%num_points(1)
            do k3=0,spl%order(3)
                splcoe(0,:) = spl%coeff(order(1), 0, k3, i1, :, i3)
                if(spl%periodic(2)) then
                    call spl_per( &
                        spl%order(2), spl%num_points(2), spl%h_step(2), splcoe)
                else
                    call spl_reg( &
                        spl%order(2), spl%num_points(2), spl%h_step(2), splcoe)
                endif
                spl%coeff(order(1), :, k3, i1, :, i3) = splcoe
            enddo
        enddo
        enddo
        deallocate(splcoe)

        ! Spline over x1
        allocate(splcoe(0:spl%order(1), spl%num_points(1)))
        do i3=1,spl%num_points(3)
        do i2=1,spl%num_points(2)
            do k3=0,spl%order(3)
            do k2=0,spl%order(2)
                splcoe(0,:) = spl%coeff(order(1), k2, k3, :, i2, i3)
                if(spl%periodic(1)) then
                    call spl_per( &
                        spl%order(1), spl%num_points(1), spl%h_step(1), splcoe)
                else
                    call spl_reg( &
                        spl%order(1), spl%num_points(1), spl%h_step(1), splcoe)
                endif
                spl%coeff(order(1):0:-1, k2, k3, :, i2, i3) = splcoe
            enddo
            enddo
        enddo
        enddo
        deallocate(splcoe)
    end subroutine construct_splines_3d


    subroutine construct_splines_3d_masked(x_min, x_max, y, order, periodic, &
                                           mask_nodes, spl)
        real(dp), intent(in) :: x_min(:), x_max(:), y(:,:,:)
        integer, intent(in) :: order(3)
        logical, intent(in) :: periodic(3)
        logical, intent(in) :: mask_nodes(:,:,:)
        type(SplineData3D), intent(out) :: spl

        integer :: n1, n2, n3, i1, i2, i3

        call construct_splines_3d(x_min, x_max, y, order, periodic, spl)

        if (any(shape(mask_nodes) /= spl%num_points)) then
            error stop "construct_splines_3d_masked: mask size mismatch"
        end if

        n1 = spl%num_points(1)
        n2 = spl%num_points(2)
        n3 = spl%num_points(3)

        if (allocated(spl%valid_cell)) deallocate(spl%valid_cell)
        if (n1 > 1 .and. n2 > 1 .and. n3 > 1) then
            allocate(spl%valid_cell(n1 - 1, n2 - 1, n3 - 1))
            do i3 = 1, n3 - 1
                do i2 = 1, n2 - 1
                    do i1 = 1, n1 - 1
                        spl%valid_cell(i1, i2, i3) = &
                            mask_nodes(i1, i2, i3) .and. &
                            mask_nodes(i1 + 1, i2, i3) .and. &
                            mask_nodes(i1, i2 + 1, i3) .and. &
                            mask_nodes(i1 + 1, i2 + 1, i3) .and. &
                            mask_nodes(i1, i2, i3 + 1) .and. &
                            mask_nodes(i1 + 1, i2, i3 + 1) .and. &
                            mask_nodes(i1, i2 + 1, i3 + 1) .and. &
                            mask_nodes(i1 + 1, i2 + 1, i3 + 1)
                    end do
                end do
            end do
            spl%has_mask = .true.
        else
            spl%has_mask = .false.
        end if
    end subroutine construct_splines_3d_masked


    subroutine evaluate_splines_3d(spl, x, y)
        type(SplineData3D), intent(in) :: spl
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: y

        real(dp) :: x_norm(3), x_local(3), xj
        real(dp) :: coeff_3(0:spl%order(3)), &
                    coeff_23(0:spl%order(2),0:spl%order(3)), &
                    coeff_local(0:spl%order(1),0:spl%order(2),0:spl%order(3))
        integer :: interval_index(3), k1, k2, k3, j

        do j=1,3
            if (spl%periodic(j)) then
                xj = modulo(x(j), spl%h_step(j)*(spl%num_points(j)-1))
            else
                xj = x(j)
            end if
            x_norm(j) = (xj - spl%x_min(j))/spl%h_step(j)
            interval_index(j) = max(0, min(spl%num_points(j)-1, int(x_norm(j))))
            x_local(j) = (x_norm(j) - dble(interval_index(j)))*spl%h_step(j)
        end do

        coeff_local(:, :, :) = spl%coeff(:, :, :, &
            interval_index(1) + 1, interval_index(2) + 1, interval_index(3) + 1)

        ! Interpolation over x1
        coeff_23(:, :) = coeff_local(0, :, :)
        do k1 = 1, spl%order(1)
            coeff_23(:, :) = coeff_local(k1, :, :) + x_local(1)*coeff_23(:, :)
        enddo

        ! Interpolation over x2
        coeff_3(:) = coeff_23(spl%order(2), :)
        do k2 = spl%order(2)-1, 0, -1
            coeff_3(:) = coeff_23(k2, :) + x_local(2)*coeff_3
        enddo

        ! Interpolation over x3
        y = coeff_3(spl%order(3))
        do k3 = spl%order(3)-1, 0, -1
            y = coeff_3(k3) + x_local(3)*y
        enddo

    end subroutine evaluate_splines_3d


    subroutine evaluate_splines_3d_masked(spl, x, y, is_valid)
        type(SplineData3D), intent(in) :: spl
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: y
        logical, intent(out), optional :: is_valid

        real(dp) :: x_norm(3), xj
        integer :: interval_index(3), cell_index(3), j

        if (.not. spl%has_mask) then
            call evaluate_splines_3d(spl, x, y)
            if (present(is_valid)) is_valid = .true.
            return
        end if

        do j = 1, 3
            if (spl%periodic(j)) then
                xj = modulo(x(j), spl%h_step(j)*(spl%num_points(j) - 1))
            else
                xj = x(j)
            end if
            x_norm(j) = (xj - spl%x_min(j))/spl%h_step(j)
            interval_index(j) = max(0, min(spl%num_points(j) - 1, &
                int(x_norm(j))))
            cell_index(j) = max(1, min(spl%num_points(j) - 1, &
                interval_index(j) + 1))
        end do

        if (.not. spl%valid_cell(cell_index(1), cell_index(2), &
            cell_index(3))) then
            y = 0.0_dp
            if (present(is_valid)) is_valid = .false.
            return
        end if

        call evaluate_splines_3d(spl, x, y)
        if (present(is_valid)) is_valid = .true.
    end subroutine evaluate_splines_3d_masked


    subroutine evaluate_splines_3d_der(spl, x, y, dy)

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
                xj = modulo(x(j), spl%h_step(j)*(spl%num_points(j)-1))
            else
                xj = x(j)
            end if
            x_norm(j) = (xj - spl%x_min(j))/spl%h_step(j)
            interval_index(j) = max(0, min(spl%num_points(j)-1, int(x_norm(j))))
            x_local(j) = (x_norm(j) - dble(interval_index(j)))*spl%h_step(j)
        end do

        associate(N1 => spl%order(1), N2 => spl%order(2), N3 => spl%order(3))

            ! Interpolation over x1
            do k3 = 0, N3
                do k2 = 0, N2
                    coeff_23(k2, k3) = spl%coeff(0, k2, k3, interval_index(1) + 1, &
                        interval_index(2) + 1, interval_index(3) + 1)
                    do k1 = 1, N1
                        coeff_23(k2, k3) = spl%coeff(k1, k2, k3, interval_index(1) + 1, &
                            interval_index(2) + 1, interval_index(3) + 1) &
                            + x_local(1)*coeff_23(k2, k3)
                    enddo
                enddo
            enddo

            ! First derivitative over x1
            do k3 = 0, N3
                do k2 = 0, N2
                    coeff_23_dx1(k2, k3) = spl%coeff(0, k2, k3, interval_index(1) + 1, &
                        interval_index(2) + 1, interval_index(3) + 1)*N1
                    do k1 = 1, N1-1
                        coeff_23_dx1(k2, k3) = spl%coeff(k1, k2, k3, interval_index(1) + 1, &
                            interval_index(2) + 1, interval_index(3) + 1)*(N1-k1) &
                            + x_local(1)*coeff_23_dx1(k2, k3)
                    enddo
                enddo
            enddo

            ! Interpolation over x2 and pure derivatives over x1
            coeff_3(:) = coeff_23(N2, :)
            coeff_3_dx1(:) = coeff_23_dx1(N2, :)
            do k2 = N2-1, 0, -1
                coeff_3(:) = coeff_23(k2, :) + x_local(2)*coeff_3
                coeff_3_dx1(:) = coeff_23_dx1(k2, :) &
                    + x_local(2)*coeff_3_dx1
            enddo
            ! First derivitatives over x2
            coeff_3_dx2(:) = coeff_23(N2, :)*N2
            do k2 = N2-1, 1, -1
                coeff_3_dx2(:) = coeff_23(k2, :)*k2 &
                    + x_local(2)*coeff_3_dx2(:)
            enddo

            ! Interpolation over x3
            y = coeff_3(N3)
            dy(1) = coeff_3_dx1(N3)
            dy(2) = coeff_3_dx2(N3)
            do k3 = N3-1, 0, -1
                y = coeff_3(k3) + x_local(3)*y
                dy(1) = coeff_3_dx1(k3) + x_local(3)*dy(1)
                dy(2) = coeff_3_dx2(k3) + x_local(3)*dy(2)
            enddo
            ! First derivitatives over x3
            dy(3) = coeff_3(N3)*N3
            do k3 = N3-1, 1, -1
                dy(3) = coeff_3(k3)*k3 + x_local(3)*dy(3)
            enddo

        end associate

    end subroutine evaluate_splines_3d_der


    subroutine evaluate_splines_3d_der2(spl, x, y, dy, d2y)

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
                xj = modulo(x(j), spl%h_step(j)*(spl%num_points(j)-1))
            else
                xj = x(j)
            end if
            x_norm(j) = (xj - spl%x_min(j))/spl%h_step(j)
            interval_index(j) = max(0, min(spl%num_points(j)-1, int(x_norm(j))))
            x_local(j) = (x_norm(j) - dble(interval_index(j)))*spl%h_step(j)
        end do

        associate(N1 => spl%order(1), N2 => spl%order(2), N3 => spl%order(3))

            ! Interpolation over x1
            do k3 = 0, N3
                do k2 = 0, N2
                    coeff_23(k2, k3) = spl%coeff(0, k2, k3, interval_index(1) + 1, &
                        interval_index(2) + 1, interval_index(3) + 1)
                    do k1 = 1, N1
                        coeff_23(k2, k3) = spl%coeff(k1, k2, k3, interval_index(1) + 1, &
                            interval_index(2) + 1, interval_index(3) + 1) &
                            + x_local(1)*coeff_23(k2, k3)
                    enddo
                enddo
            enddo

            ! First derivitative over x1
            do k3 = 0, N3
                do k2 = 0, N2
                    coeff_23_dx1(k2, k3) = spl%coeff(0, k2, k3, interval_index(1) + 1, &
                        interval_index(2) + 1, interval_index(3) + 1)*N1
                    do k1 = 1, N1-1
                        coeff_23_dx1(k2, k3) = spl%coeff(k1, k2, k3, interval_index(1) + 1, &
                            interval_index(2) + 1, interval_index(3) + 1)*(N1-k1) &
                            + x_local(1)*coeff_23_dx1(k2, k3)
                    enddo
                enddo
            enddo

            ! Second derivitative over x1
            do k3 = 0, N3
                do k2 = 0, N2
                    coeff_23_dx1x1(k2,k3) = spl%coeff(0, k2, k3, interval_index(1) + 1, &
                        interval_index(2) + 1, interval_index(3) + 1)*N1*(N1-1)
                    do k1 = 1, N1-2
                        coeff_23_dx1x1(k2, k3)=spl%coeff(k1,k2,k3, interval_index(1) + 1, &
                            interval_index(2) + 1, interval_index(3) + 1)*(N1-k1)*(N1-k1-1) &
                            + x_local(1)*coeff_23_dx1x1(k2, k3)
                    enddo
                enddo
            enddo

            ! Interpolation over x2 and pure derivatives over x1
            coeff_3(:) = coeff_23(N2, :)
            coeff_3_dx1(:) = coeff_23_dx1(N2, :)
            coeff_3_dx1x1(:) = coeff_23_dx1x1(N2, :)
            do k2 = N2-1, 0, -1
                coeff_3(:) = coeff_23(k2, :) + x_local(2)*coeff_3
                coeff_3_dx1(:) = coeff_23_dx1(k2, :) &
                    + x_local(2)*coeff_3_dx1
                coeff_3_dx1x1(:) = coeff_23_dx1x1(k2, :) &
                    + x_local(2)*coeff_3_dx1x1
            enddo
            ! First derivitatives over x2
            coeff_3_dx2(:) = coeff_23(N2, :)*N2
            coeff_3_dx1x2(:) = coeff_23_dx1(N2, :)*N2
            do k2 = N2-1, 1, -1
                coeff_3_dx2(:) = coeff_23(k2, :)*k2 &
                    + x_local(2)*coeff_3_dx2(:)
                coeff_3_dx1x2(:) = coeff_23_dx1(k2, :)*k2 &
                    + x_local(2)*coeff_3_dx1x2
            enddo
            ! Second derivitative over x2
            coeff_3_dx2x2(:) = coeff_23(N2, :)*N2*(N2-1)
            do k2 = N2-1, 2, -1
                coeff_3_dx2x2(:) = coeff_23(k2, :)*k2*(k2-1) &
                    + x_local(2)*coeff_3_dx2x2
            enddo

            ! Interpolation over x3
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
            ! First derivitatives over x3
            dy(3) = coeff_3(N3)*N3
            d2y(3) = coeff_3_dx1(N3)*N3
            d2y(5) = coeff_3_dx2(N3)*N3
            do k3 = N3-1, 1, -1
                dy(3) = coeff_3(k3)*k3 + x_local(3)*dy(3)
                d2y(3) = coeff_3_dx1(k3)*k3 + x_local(3)*d2y(3)
                d2y(5) = coeff_3_dx2(k3)*k3 + x_local(3)*d2y(5)
            enddo
            ! Second derivitative over x3
            d2y(6) = coeff_3(N3)*N3*(N3-1)
            do k3 = N3-1, 2, -1
                d2y(6) = coeff_3(k3)*k3*(k3-1) + x_local(3)*d2y(6)
            enddo

        end associate

    end subroutine evaluate_splines_3d_der2


    subroutine destroy_splines_3d(spl)
        type(SplineData3D), intent(inout) :: spl

        if(allocated(spl%coeff)) deallocate(spl%coeff)
        if(allocated(spl%valid_cell)) deallocate(spl%valid_cell)
        spl%has_mask = .false.
    end subroutine destroy_splines_3d


    subroutine disp_3d(spl)
        type(SplineData3D), intent(in) :: spl
        print *, "SplineData3D"
        print *, "  order = ", spl%order
        print *, "  num_points = ", spl%num_points
        print *, "  periodic = ", spl%periodic
        print *, "  x_min = ", spl%x_min
        print *, "  x_max = ", spl%x_min + (spl%num_points-1)*spl%h_step
        print *, "  h_step = ", spl%h_step
    end subroutine disp_3d

end module interpolate
