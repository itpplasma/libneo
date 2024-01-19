module interpolate
    implicit none

    integer, parameter :: MAX_ORDER = 5

    type :: SplineData1D
        integer :: order
        integer :: num_points
        logical :: periodic
        real(8) :: h_step
        real(8), dimension(:,:), allocatable :: coeff
    end type SplineData1D

    type :: SplineData2D
        integer :: order(2)
        integer :: num_points(2)
        logical :: periodic(2)
        real(8) :: h_step(2)
        real(8), dimension(:,:,:,:), allocatable :: coeff
    end type SplineData2D

    type :: SplineData3D
        integer :: order(3)
        integer :: num_points(3)
        logical :: periodic(3)
        real(8) :: h_step(3)
        real(8), dimension(:,:,:,:,:,:), allocatable :: coeff
    end type SplineData3D

contains

    subroutine construct_splines_1d(x, y, order, periodic, spl)
        real(8), intent(in) :: x(:), y(:)
        integer, intent(in) :: order
        logical, intent(in) :: periodic

        type(SplineData1D), intent(out) :: spl

        spl%order = order
        spl%periodic = periodic
        spl%num_points = size(x)
        spl%h_step = x(2) - x(1)

        allocate(spl%coeff(0:order, spl%num_points))
        spl%coeff(0,:) = y

        if (periodic) then
            call spl_per(spl%order, spl%num_points, spl%h_step, spl%coeff)
        else
            call spl_reg(spl%order, spl%num_points, spl%h_step, spl%coeff)
        endif
    end subroutine construct_splines_1d


    subroutine destroy_splines_1d(spl)
        type(SplineData1D), intent(inout) :: spl

        deallocate(spl%coeff)
    end subroutine destroy_splines_1d


    subroutine evaluate_splines_1d(x, spl, y)
        real(8), intent(in) :: x
        type(SplineData1D), intent(in) :: spl
        real(8), intent(out) :: y

        real(8) :: x_norm, x_local, local_coeff(0:MAX_ORDER)
        integer :: interval_index, k_power

        x_norm = x/spl%h_step
        interval_index = max(0, min(spl%num_points-1, int(x_norm)))
        x_local = (x_norm - dble(interval_index))*spl%h_step  ! Distance to grid point

        local_coeff(0:spl%order) = spl%coeff(:, interval_index+1)

        ! Start with largest power and then multiply recursively
        y = local_coeff(spl%order)
        do k_power = spl%order, 0, -1
            y = local_coeff(k_power) + x_local*y
        enddo
    end subroutine evaluate_splines_1d


    subroutine construct_splines_2d(x1, x2, y, order, periodic, spl)
        real(8), intent(in) :: x1(:), x2(:), y(:,:)
        integer, intent(in) :: order(2)
        logical, intent(in) :: periodic(2)

        type(SplineData2D), intent(out) :: spl

        real(8), dimension(:,:), allocatable  :: splcoe

        integer :: i1, i2  ! Loop indices for points (1 ... num_points)
        integer :: k1, k2  ! Loop indices for polynomial order (0 ... order)

        spl%order = order
        spl%periodic = periodic
        spl%num_points(1) = size(x1)
        spl%num_points(2) = size(x2)
        spl%h_step(1) = x1(2) - x1(1)
        spl%h_step(2) = x2(2) - x2(1)

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
            do k2 = 0, spl%order(2)
                spl%coeff(1, k2, i1, :) = splcoe(k2, :)
            enddo
        enddo
        deallocate(splcoe)

        ! Spline over x1
        allocate(splcoe(0:spl%order(1), spl%num_points(1)))
        do i2=1,spl%num_points(2)
            do k2=0,spl%order(2)
                splcoe(0,:) = spl%coeff(1, k2, :, i2)
                if(spl%periodic(1)) then
                    call spl_per(spl%order(1), spl%num_points(1), spl%h_step(1), splcoe)
                else
                    call spl_reg(spl%order(1), spl%num_points(1), spl%h_step(1), splcoe)
                endif
                do k1 = 0, spl%order(1)
                    spl%coeff(k1, k2, :, i2) = splcoe(k1, :)
                enddo
            enddo
        enddo
        deallocate(splcoe)

    end subroutine construct_splines_2d


    subroutine evaluate_splines_2d(x, spl, y)
        real(8), intent(in) :: x(2)
        type(SplineData2D), intent(in) :: spl
        real(8), intent(out) :: y

        real(8) :: x_norm(2), x_local(2)
        real(8) :: coeff_2(0:MAX_ORDER), coeff_12(0:MAX_ORDER,0:MAX_ORDER)
        integer :: interval_index(2), k1, k2, j

        do j=1,2
            x_norm(j) = x(j)/spl%h_step(j)
            interval_index(j) = max(0, min(spl%num_points(j)-1, int(x_norm(j))))
            x_local(j) = (x_norm(j) - dble(interval_index(j)))*spl%h_step(j)
        end do

        coeff_12(0:spl%order(1),0:spl%order(2)) = &
            spl%coeff(:, :, interval_index(1) + 1, interval_index(2) + 1)

        coeff_2(0:spl%order(2)) = coeff_12(spl%order(1), 0:spl%order(2))
        do k1 = spl%order(1)-1, 0, -1
            coeff_2(0:spl%order(1)) = &
                coeff_12(k1, 0:spl%order(2)) + x_local(1)*coeff_2(0:spl%order(1))
        enddo

        y = coeff_2(spl%order(2))
        do k2 = spl%order(2)-1, 0, -1
            y = coeff_2(k2) + x_local(2)*y
        enddo
    end subroutine evaluate_splines_2d


    subroutine destroy_splines_2d(spl)
        type(SplineData2D), intent(inout) :: spl

        deallocate(spl%coeff)
    end subroutine destroy_splines_2d


    subroutine construct_splines_3d(x1, x2, x3, y, order, periodic, spl)
        real(8), intent(in) :: x1(:), x2(:), x3(:), y(:,:,:)
        integer, intent(in) :: order(3)
        logical, intent(in) :: periodic(3)

        type(SplineData3D), intent(out) :: spl

        real(8), dimension(:,:), allocatable  :: splcoe

        integer :: i1, i2, i3  ! Loop indices for points (1 ... num_points)
        integer :: k1, k2, k3  ! Loop indices for polynomial order (0 ... order)

        spl%order = order
        spl%periodic = periodic
        spl%num_points(1) = size(x1)
        spl%num_points(2) = size(x2)
        spl%num_points(3) = size(x3)
        spl%h_step(1) = x1(2) - x1(1)
        spl%h_step(2) = x2(2) - x2(1)
        spl%h_step(3) = x3(2) - x3(1)

        allocate(spl%coeff(0:order(1), 0:order(2), 0:order(3), &
                 spl%num_points(1), spl%num_points(2), spl%num_points(3)))

        ! Spline over x3
        allocate(splcoe(0:spl%order(3), spl%num_points(3)))
        do i1=1,spl%num_points(1)
        do i2=1,spl%num_points(2)
            splcoe(0,:) = y(i1, i2, :)
            if (spl%periodic(3)) then
                call spl_per(spl%order(3), spl%num_points(3), spl%h_step(3), splcoe)
            else
                call spl_reg(spl%order(3), spl%num_points(3), spl%h_step(3), splcoe)
            endif
            do k3 = 0, spl%order(3)
                spl%coeff(1, 1, k3, i1, i2, :) = splcoe(k3, :)
            enddo
        enddo
        enddo
        deallocate(splcoe)

        ! Spline over x2
        allocate(splcoe(0:spl%order(2), spl%num_points(2)))
        do i1=1,spl%num_points(1)
        do i3=1,spl%num_points(3)
            do k3=0,spl%order(3)
                splcoe(0,:) = spl%coeff(1, 1, k3, i1, :, i3)
                if(spl%periodic(1)) then
                    call spl_per(spl%order(2), spl%num_points(2), spl%h_step(2), splcoe)
                else
                    call spl_reg(spl%order(2), spl%num_points(2), spl%h_step(2), splcoe)
                endif
                do k2 = 0, spl%order(2)
                    spl%coeff(1, k2, k3, i1, :, i3) = splcoe(k2, :)
                enddo
            enddo
        enddo
        enddo
        deallocate(splcoe)

        ! Spline over x1
        allocate(splcoe(0:spl%order(1), spl%num_points(1)))
        do i2=1,spl%num_points(2)
        do i3=1,spl%num_points(3)
            do k2=0,spl%order(2)
            do k3=0,spl%order(3)
                splcoe(0,:) = spl%coeff(1, k2, k3, i1, :, i3)
                if(spl%periodic(1)) then
                    call spl_per(spl%order(1), spl%num_points(1), spl%h_step(1), splcoe)
                else
                    call spl_reg(spl%order(1), spl%num_points(1), spl%h_step(1), splcoe)
                endif
                do k1 = 0, spl%order(1)
                    spl%coeff(k1, k2, k3, i1, :, i3) = splcoe(k1, :)
                enddo
            enddo
            enddo
        enddo
        enddo
        deallocate(splcoe)
    end subroutine construct_splines_3d


    subroutine evaluate_splines_3d(x, spl, y)
        real(8), intent(in) :: x(3)
        type(SplineData3D), intent(in) :: spl
        real(8), intent(out) :: y

        real(8) :: x_norm(3), x_local(3)
        real(8) :: coeff_2(0:MAX_ORDER), coeff_12(0:MAX_ORDER,0:MAX_ORDER), &
            coeff_123(0:MAX_ORDER,0:MAX_ORDER,0:MAX_ORDER)
        integer :: interval_index(3), k1, k2, k3, j

        do j=1,3
            x_norm(j) = x(j)/spl%h_step(j)
            interval_index(j) = max(0, min(spl%num_points(j)-1, int(x_norm(j))))
            x_local(j) = (x_norm(j) - dble(interval_index(j)))*spl%h_step(j)
        end do

        coeff_123(0:spl%order(1), 0:spl%order(2), 0:spl%order(3)) = &
            spl%coeff(:, :, :, &
            interval_index(1) + 1, interval_index(2) + 1, interval_index(3) + 1)

        coeff_12(0:spl%order(1), 0:spl%order(2)) = &
            coeff_123(1, 0:spl%order(2), 0:spl%order(3))
        do k1 = spl%order(1)-1, 0, -1
            coeff_12(0:spl%order(1), 0:spl%order(2)) = &
            coeff_123(k1, 0:spl%order(2), 0:spl%order(3)) &
                + x_local(1)*coeff_12(0:spl%order(1), 0:spl%order(2))
        enddo

        coeff_2(0:spl%order(2)) = coeff_12(spl%order(1), 0:spl%order(2))
        do k2 = spl%order(1)-1, 0, -1
            coeff_2(0:spl%order(1)) = &
                coeff_12(k1, 0:spl%order(2)) + x_local(2)*coeff_2(0:spl%order(2))
        enddo

        y = coeff_2(spl%order(3))
        do k3 = spl%order(3)-1, 0, -1
            y = coeff_2(k3) + x_local(3)*y
        enddo
    end subroutine evaluate_splines_3d


    subroutine destroy_splines_3d(spl)
        type(SplineData3D), intent(inout) :: spl

        deallocate(spl%coeff)
    end subroutine destroy_splines_3d
end module interpolate
