module interpolate
    use spl_three_to_five_sub

    implicit none
    integer, parameter :: dp = kind(1.0d0)

    type :: SplineData1D
        integer :: order
        integer :: num_points
        logical :: periodic
        real(dp) :: x_min
        real(dp) :: h_step
        real(dp), dimension(:,:), allocatable :: coeff
    end type SplineData1D

    type :: SplineData2D
        integer :: order(2)
        integer :: num_points(2)
        logical :: periodic(2)
        real(dp) :: h_step(2)
        real(dp) :: x_min(2)
        real(dp), dimension(:,:,:,:), allocatable :: coeff
    end type SplineData2D

    type :: SplineData3D
        integer :: order(3)
        integer :: num_points(3)
        logical :: periodic(3)
        real(dp) :: h_step(3)
        real(dp) :: x_min(3)
        real(dp), dimension(:,:,:,:,:,:), allocatable :: coeff
    end type SplineData3D
    
    ! Batch spline types for multiple quantities on shared grid
    type :: BatchSplineData1D
        ! Shared grid data
        integer :: order
        integer :: num_points  
        logical :: periodic
        real(dp) :: x_min
        real(dp) :: h_step
        
        ! Batched coefficients: (num_quantities, order+1, num_points)
        ! Optimized layout for cache-efficient evaluation at single points
        integer :: num_quantities
        real(dp), dimension(:,:,:), allocatable :: coeff
    end type BatchSplineData1D
    
    type :: BatchSplineData2D
        ! Shared grid data
        integer :: order(2)
        integer :: num_points(2)
        logical :: periodic(2)
        real(dp) :: h_step(2)
        real(dp) :: x_min(2)
        
        ! Batched coefficients: (num_quantities, order1+1, order2+1, n1, n2)
        ! Optimized layout for cache-efficient evaluation at single points
        integer :: num_quantities
        real(dp), dimension(:,:,:,:,:), allocatable :: coeff
    end type BatchSplineData2D
    
    type :: BatchSplineData3D
        ! Shared grid data
        integer :: order(3)
        integer :: num_points(3)
        logical :: periodic(3)
        real(dp) :: h_step(3)
        real(dp) :: x_min(3)
        
        ! Batched coefficients: (num_quantities, order1+1, order2+1, order3+1, n1, n2, n3)
        ! Optimized layout for cache-efficient evaluation at single points
        integer :: num_quantities
        real(dp), dimension(:,:,:,:,:,:,:), allocatable :: coeff
    end type BatchSplineData3D


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


    subroutine destroy_splines_1d(spl)
        type(SplineData1D), intent(inout) :: spl

        if(allocated(spl%coeff)) deallocate(spl%coeff)
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
    
    
    ! ============================================================================
    ! Batch Spline Routines - 1D
    ! ============================================================================
    
    subroutine construct_batch_splines_1d(x_min, x_max, y_batch, order, periodic, spl)
        real(dp), intent(in) :: x_min, x_max
        real(dp), intent(in) :: y_batch(:,:)  ! (n_points, n_quantities)
        integer, intent(in) :: order
        logical, intent(in) :: periodic
        type(BatchSplineData1D), intent(out) :: spl
        
        integer :: iq, n_points, n_quantities
        real(dp), dimension(:,:), allocatable :: splcoe_temp
        
        n_points = size(y_batch, 1)
        n_quantities = size(y_batch, 2)
        
        ! Set metadata
        spl%x_min = x_min
        spl%order = order
        spl%periodic = periodic
        spl%num_points = n_points
        spl%num_quantities = n_quantities
        spl%h_step = (x_max - x_min) / (n_points - 1)
        
        ! Allocate coefficient array with optimized layout
        if(allocated(spl%coeff)) deallocate(spl%coeff)
        allocate(spl%coeff(n_quantities, 0:order, n_points))
        
        ! Temporary array for spline coefficient computation
        allocate(splcoe_temp(0:order, n_points))
        
        ! Compute spline coefficients for each quantity
        do iq = 1, n_quantities
            ! Copy y values for this quantity
            splcoe_temp(0,:) = y_batch(:, iq)
            
            ! Compute spline coefficients
            if (periodic) then
                call spl_per(order, n_points, spl%h_step, splcoe_temp)
            else
                call spl_reg(order, n_points, spl%h_step, splcoe_temp)
            endif
            
            ! Store in batch array with new layout
            spl%coeff(iq,:,:) = splcoe_temp
        end do
        
        deallocate(splcoe_temp)
        
    end subroutine construct_batch_splines_1d
    
    
    subroutine destroy_batch_splines_1d(spl)
        type(BatchSplineData1D), intent(inout) :: spl
        
        if(allocated(spl%coeff)) deallocate(spl%coeff)
    end subroutine destroy_batch_splines_1d
    
    
    subroutine evaluate_batch_splines_1d(spl, x, y_batch)
        type(BatchSplineData1D), intent(in) :: spl
        real(dp), intent(in) :: x
        real(dp), intent(out) :: y_batch(:)  ! (n_quantities)
        
        real(dp) :: x_norm, x_local, xj
        integer :: interval_index, k_power, iq
        
        ! Handle periodic boundary
        if (spl%periodic) then
            xj = modulo(x, spl%h_step*(spl%num_points-1))
        else
            xj = x
        end if
        
        ! Find interval and local coordinate
        x_norm = (xj - spl%x_min) / spl%h_step
        interval_index = max(0, min(spl%num_points-1, int(x_norm)))
        x_local = (x_norm - dble(interval_index))*spl%h_step
        
        ! Initialize with highest order coefficients
        ! Memory access: spl%coeff(:, spl%order, interval_index+1) is contiguous!
        !$omp simd
        do iq = 1, spl%num_quantities
            y_batch(iq) = spl%coeff(iq, spl%order, interval_index+1)
        end do
        
        ! Horner method with SIMD for all quantities at once
        do k_power = spl%order-1, 0, -1
            !$omp simd
            do iq = 1, spl%num_quantities
                y_batch(iq) = spl%coeff(iq, k_power, interval_index+1) + &
                              x_local*y_batch(iq)
            end do
        end do
        
    end subroutine evaluate_batch_splines_1d
    
    
    subroutine evaluate_batch_splines_1d_single(spl, x, iq, y)
        type(BatchSplineData1D), intent(in) :: spl
        real(dp), intent(in) :: x
        integer, intent(in) :: iq  ! quantity index
        real(dp), intent(out) :: y
        
        real(dp) :: x_norm, x_local, xj
        integer :: interval_index, k_power
        
        ! Handle periodic boundary
        if (spl%periodic) then
            xj = modulo(x, spl%h_step*(spl%num_points-1))
        else
            xj = x
        end if
        
        ! Find interval and local coordinate
        x_norm = (xj - spl%x_min) / spl%h_step
        interval_index = max(0, min(spl%num_points-1, int(x_norm)))
        x_local = (x_norm - dble(interval_index))*spl%h_step
        
        ! Horner method for single quantity evaluation
        y = spl%coeff(iq, spl%order, interval_index+1)
        do k_power = spl%order-1, 0, -1
            y = spl%coeff(iq, k_power, interval_index+1) + x_local*y
        enddo
        
    end subroutine evaluate_batch_splines_1d_single
    
    
    subroutine evaluate_batch_splines_1d_der(spl, x, y_batch, dy_batch)
        type(BatchSplineData1D), intent(in) :: spl
        real(dp), intent(in) :: x
        real(dp), intent(out) :: y_batch(:)   ! (n_quantities)
        real(dp), intent(out) :: dy_batch(:)  ! (n_quantities)
        
        real(dp) :: x_norm, x_local, xj
        integer :: interval_index, k_power, iq
        
        ! Handle periodic boundary
        if (spl%periodic) then
            xj = modulo(x, spl%h_step*(spl%num_points-1))
        else
            xj = x
        end if
        
        ! Find interval and local coordinate
        x_norm = (xj - spl%x_min) / spl%h_step
        interval_index = max(0, min(spl%num_points-1, int(x_norm)))
        x_local = (x_norm - dble(interval_index))*spl%h_step
        
        ! Initialize with highest order coefficients
        !$omp simd
        do iq = 1, spl%num_quantities
            y_batch(iq) = spl%coeff(iq, spl%order, interval_index+1)
            dy_batch(iq) = spl%coeff(iq, spl%order, interval_index+1)*spl%order
        end do
        
        ! Compute value and derivative with SIMD
        do k_power = spl%order-1, 1, -1
            !$omp simd
            do iq = 1, spl%num_quantities
                y_batch(iq) = spl%coeff(iq, k_power, interval_index+1) + &
                              x_local*y_batch(iq)
                dy_batch(iq) = spl%coeff(iq, k_power, interval_index+1)*k_power + &
                               x_local*dy_batch(iq)
            end do
        end do
        
        ! Final step for y_batch only
        !$omp simd
        do iq = 1, spl%num_quantities
            y_batch(iq) = spl%coeff(iq, 0, interval_index+1) + x_local*y_batch(iq)
        end do
        
    end subroutine evaluate_batch_splines_1d_der
    
    
    subroutine evaluate_batch_splines_1d_der2(spl, x, y_batch, dy_batch, d2y_batch)
        type(BatchSplineData1D), intent(in) :: spl
        real(dp), intent(in) :: x
        real(dp), intent(out) :: y_batch(:)    ! (n_quantities)
        real(dp), intent(out) :: dy_batch(:)   ! (n_quantities)
        real(dp), intent(out) :: d2y_batch(:)  ! (n_quantities)
        
        real(dp) :: x_norm, x_local, xj
        integer :: interval_index, k_power, iq
        
        ! Handle periodic boundary
        if (spl%periodic) then
            xj = modulo(x, spl%h_step*(spl%num_points-1))
        else
            xj = x
        end if
        
        ! Find interval and local coordinate
        x_norm = (xj - spl%x_min) / spl%h_step
        interval_index = max(0, min(spl%num_points-1, int(x_norm)))
        x_local = (x_norm - dble(interval_index))*spl%h_step
        
        ! Initialize with highest order coefficients
        !$omp simd
        do iq = 1, spl%num_quantities
            y_batch(iq) = spl%coeff(iq, spl%order, interval_index+1)
            dy_batch(iq) = spl%coeff(iq, spl%order, interval_index+1)*spl%order
            d2y_batch(iq) = spl%coeff(iq, spl%order, interval_index+1)*spl%order* &
                            (spl%order-1)
        end do
        
        ! Compute value and derivatives with SIMD
        do k_power = spl%order-1, 2, -1
            !$omp simd
            do iq = 1, spl%num_quantities
                y_batch(iq) = spl%coeff(iq, k_power, interval_index+1) + &
                              x_local*y_batch(iq)
                dy_batch(iq) = spl%coeff(iq, k_power, interval_index+1)*k_power + &
                               x_local*dy_batch(iq)
                d2y_batch(iq) = spl%coeff(iq, k_power, interval_index+1)*k_power* &
                                (k_power-1) + x_local*d2y_batch(iq)
            end do
        end do
        
        ! Handle k_power = 1
        if (spl%order >= 1) then
            !$omp simd
            do iq = 1, spl%num_quantities
                y_batch(iq) = spl%coeff(iq, 1, interval_index+1) + x_local*y_batch(iq)
                dy_batch(iq) = spl%coeff(iq, 1, interval_index+1) + x_local*dy_batch(iq)
            end do
        endif
        
        ! Final step for y_batch only
        !$omp simd
        do iq = 1, spl%num_quantities
            y_batch(iq) = spl%coeff(iq, 0, interval_index+1) + x_local*y_batch(iq)
        end do
        
    end subroutine evaluate_batch_splines_1d_der2
    
    
    ! ============================================================================
    ! Batch Spline Routines - 2D
    ! ============================================================================
    
    subroutine construct_batch_splines_2d(x_min, x_max, y_batch, order, periodic, spl)
        real(dp), intent(in) :: x_min(:), x_max(:)
        real(dp), intent(in) :: y_batch(:,:,:)  ! (n1, n2, n_quantities)
        integer, intent(in) :: order(:)
        logical, intent(in) :: periodic(:)
        type(BatchSplineData2D), intent(out) :: spl
        
        real(dp), dimension(:,:), allocatable  :: splcoe
        real(dp), dimension(:,:,:,:), allocatable :: temp_coeff
        integer :: i1, i2, iq  ! Loop indices
        integer :: k2          ! Loop index for polynomial order
        
        spl%x_min = x_min
        spl%order = order
        spl%periodic = periodic
        spl%num_points(1) = size(y_batch, 1)
        spl%num_points(2) = size(y_batch, 2)
        spl%num_quantities = size(y_batch, 3)
        spl%h_step = (x_max - x_min) / (spl%num_points - 1)
        
        if(allocated(spl%coeff)) deallocate(spl%coeff)
        allocate(spl%coeff(spl%num_quantities, 0:order(1), 0:order(2), &
                 spl%num_points(1), spl%num_points(2)))
        
        ! Allocate temporary array for intermediate results
        allocate(temp_coeff(0:order(1), 0:order(2), &
                 spl%num_points(1), spl%num_points(2)))
        
        ! Process each quantity
        do iq = 1, spl%num_quantities
            ! Spline over x2
            allocate(splcoe(0:spl%order(2), spl%num_points(2)))
            do i1=1,spl%num_points(1)
                splcoe(0,:) = y_batch(i1, :, iq)
                if (spl%periodic(2)) then
                    call spl_per(spl%order(2), spl%num_points(2), spl%h_step(2), splcoe)
                else
                    call spl_reg(spl%order(2), spl%num_points(2), spl%h_step(2), splcoe)
                endif
                temp_coeff(0, :, i1, :) = splcoe
            enddo
            deallocate(splcoe)
            
            ! Spline over x1
            allocate(splcoe(0:spl%order(1), spl%num_points(1)))
            do i2=1,spl%num_points(2)
                do k2=0,spl%order(2)
                    splcoe(0,:) = temp_coeff(0, k2, :, i2)
                    if(spl%periodic(1)) then
                        call spl_per(spl%order(1), spl%num_points(1), spl%h_step(1), splcoe)
                    else
                        call spl_reg(spl%order(1), spl%num_points(1), spl%h_step(1), splcoe)
                    endif
                    temp_coeff(:, k2, :, i2) = splcoe
                enddo
            enddo
            deallocate(splcoe)
            
            ! Store with new memory layout
            spl%coeff(iq, :, :, :, :) = temp_coeff
        end do
        
        deallocate(temp_coeff)
        
    end subroutine construct_batch_splines_2d
    
    
    subroutine destroy_batch_splines_2d(spl)
        type(BatchSplineData2D), intent(inout) :: spl
        
        if(allocated(spl%coeff)) deallocate(spl%coeff)
    end subroutine destroy_batch_splines_2d
    
    
    subroutine evaluate_batch_splines_2d(spl, x, y_batch)
        type(BatchSplineData2D), intent(in) :: spl
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: y_batch(:)  ! (n_quantities)
        
        real(dp) :: x_norm(2), x_local(2), xj
        real(dp) :: coeff_2(spl%num_quantities, 0:spl%order(2))
        integer :: interval_index(2), k1, k2, j, iq
        
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
        
        ! First reduction: evaluate along x1 dimension
        ! Initialize with highest order in x1
        !$omp simd
        do iq = 1, spl%num_quantities
            do k2 = 0, spl%order(2)
                coeff_2(iq, k2) = spl%coeff(iq, spl%order(1), k2, &
                                   interval_index(1)+1, interval_index(2)+1)
            end do
        end do
        
        ! Apply Horner method along x1
        do k1 = spl%order(1)-1, 0, -1
            !$omp simd
            do iq = 1, spl%num_quantities
                do k2 = 0, spl%order(2)
                    coeff_2(iq, k2) = spl%coeff(iq, k1, k2, interval_index(1)+1, &
                                       interval_index(2)+1) + x_local(1)*coeff_2(iq, k2)
                end do
            end do
        end do
        
        ! Second reduction: evaluate along x2 dimension
        ! Initialize with highest order in x2
        !$omp simd
        do iq = 1, spl%num_quantities
            y_batch(iq) = coeff_2(iq, spl%order(2))
        end do
        
        ! Apply Horner method along x2
        do k2 = spl%order(2)-1, 0, -1
            !$omp simd
            do iq = 1, spl%num_quantities
                y_batch(iq) = coeff_2(iq, k2) + x_local(2)*y_batch(iq)
            end do
        end do
        
    end subroutine evaluate_batch_splines_2d
    
    
    subroutine evaluate_batch_splines_2d_der(spl, x, y_batch, dy_batch)
        type(BatchSplineData2D), intent(in) :: spl
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: y_batch(:)     ! (n_quantities)
        real(dp), intent(out) :: dy_batch(:,:)  ! (2, n_quantities)
        
        real(dp) :: x_norm(2), x_local(2), xj
        real(dp) :: coeff_2(0:spl%order(2))
        real(dp) :: coeff_2_dx1(0:spl%order(2))
        real(dp) :: coeff_local(0:spl%order(1),0:spl%order(2))
        integer :: interval_index(2), k1, k2, j, iq
        
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
        
        ! Evaluate for each quantity
        do iq = 1, spl%num_quantities
            coeff_local(:,:) = &
                spl%coeff(iq, :, :, interval_index(1) + 1, interval_index(2) + 1)
            
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
            y_batch(iq) = coeff_2(spl%order(2))
            do k2 = spl%order(2)-1, 0, -1
                y_batch(iq) = coeff_2(k2) + x_local(2)*y_batch(iq)
            enddo
            
            ! Derivative w.r.t. x1
            dy_batch(1, iq) = coeff_2_dx1(spl%order(2))
            do k2 = spl%order(2)-1, 0, -1
                dy_batch(1, iq) = coeff_2_dx1(k2) + x_local(2)*dy_batch(1, iq)
            enddo
            
            ! Derivative w.r.t. x2
            dy_batch(2, iq) = coeff_2(spl%order(2))*spl%order(2)
            do k2 = spl%order(2)-1, 1, -1
                dy_batch(2, iq) = coeff_2(k2)*k2 + x_local(2)*dy_batch(2, iq)
            enddo
        end do
        
    end subroutine evaluate_batch_splines_2d_der
    
    
    ! ============================================================================
    ! Batch Spline Routines - 3D
    ! ============================================================================
    
    subroutine construct_batch_splines_3d(x_min, x_max, y_batch, order, periodic, spl)
        real(dp), intent(in) :: x_min(:), x_max(:)
        real(dp), intent(in) :: y_batch(:,:,:,:)  ! (n1, n2, n3, n_quantities)
        integer, intent(in) :: order(3)
        logical, intent(in) :: periodic(3)
        type(BatchSplineData3D), intent(out) :: spl
        
        real(dp), dimension(:,:), allocatable  :: splcoe
        real(dp), dimension(:,:,:,:,:,:), allocatable :: temp_coeff
        integer :: i1, i2, i3, iq  ! Loop indices
        integer :: k2, k3          ! Loop indices for polynomial order
        
        spl%x_min = x_min
        spl%order = order
        spl%periodic = periodic
        spl%num_points(1) = size(y_batch, 1)
        spl%num_points(2) = size(y_batch, 2)
        spl%num_points(3) = size(y_batch, 3)
        spl%num_quantities = size(y_batch, 4)
        spl%h_step = (x_max - x_min) / (spl%num_points - 1)
        
        if(allocated(spl%coeff)) deallocate(spl%coeff)
        allocate(spl%coeff(spl%num_quantities, 0:order(1), 0:order(2), 0:order(3), &
                 spl%num_points(1), spl%num_points(2), spl%num_points(3)))
        
        ! Allocate temporary array for intermediate results
        allocate(temp_coeff(0:order(1), 0:order(2), 0:order(3), &
                 spl%num_points(1), spl%num_points(2), spl%num_points(3)))
        
        ! Process each quantity
        do iq = 1, spl%num_quantities
            ! Spline over x3
            allocate(splcoe(0:spl%order(3), spl%num_points(3)))
            do i2=1,spl%num_points(2)
            do i1=1,spl%num_points(1)
                splcoe(0,:) = y_batch(i1, i2, :, iq)
                if (spl%periodic(3)) then
                    call spl_per(spl%order(3), spl%num_points(3), spl%h_step(3), splcoe)
                else
                    call spl_reg(spl%order(3), spl%num_points(3), spl%h_step(3), splcoe)
                endif
                temp_coeff(0, 0, :, i1, i2, :) = splcoe
            enddo
            enddo
            deallocate(splcoe)
            
            ! Spline over x2
            allocate(splcoe(0:spl%order(2), spl%num_points(2)))
            do i3=1,spl%num_points(3)
            do i1=1,spl%num_points(1)
                do k3=0,spl%order(3)
                    splcoe(0,:) = temp_coeff(0, 0, k3, i1, :, i3)
                    if(spl%periodic(2)) then
                        call spl_per( &
                            spl%order(2), spl%num_points(2), spl%h_step(2), splcoe)
                    else
                        call spl_reg( &
                            spl%order(2), spl%num_points(2), spl%h_step(2), splcoe)
                    endif
                    temp_coeff(0, :, k3, i1, :, i3) = splcoe
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
                    splcoe(0,:) = temp_coeff(0, k2, k3, :, i2, i3)
                    if(spl%periodic(1)) then
                        call spl_per( &
                            spl%order(1), spl%num_points(1), spl%h_step(1), splcoe)
                    else
                        call spl_reg( &
                            spl%order(1), spl%num_points(1), spl%h_step(1), splcoe)
                    endif
                    temp_coeff(:, k2, k3, :, i2, i3) = splcoe
                enddo
                enddo
            enddo
            enddo
            deallocate(splcoe)
            
            ! Store with new memory layout
            spl%coeff(iq, :, :, :, :, :, :) = temp_coeff
        end do
        
        deallocate(temp_coeff)
        
    end subroutine construct_batch_splines_3d
    
    
    subroutine destroy_batch_splines_3d(spl)
        type(BatchSplineData3D), intent(inout) :: spl
        
        if(allocated(spl%coeff)) deallocate(spl%coeff)
    end subroutine destroy_batch_splines_3d
    
    
    subroutine evaluate_batch_splines_3d(spl, x, y_batch)
        type(BatchSplineData3D), intent(in) :: spl
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: y_batch(:)  ! (n_quantities)
        
        real(dp) :: x_norm(3), x_local(3), xj
        real(dp) :: coeff_3(spl%num_quantities, 0:spl%order(3))
        real(dp) :: coeff_23(spl%num_quantities, 0:spl%order(2), 0:spl%order(3))
        integer :: interval_index(3), k1, k2, k3, j, iq
        
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
        
        ! First reduction: evaluate along x1 dimension
        ! Initialize with highest order in x1
        !$omp simd
        do iq = 1, spl%num_quantities
            do k3 = 0, spl%order(3)
                do k2 = 0, spl%order(2)
                    coeff_23(iq, k2, k3) = spl%coeff(iq, spl%order(1), k2, k3, &
                        interval_index(1)+1, interval_index(2)+1, interval_index(3)+1)
                end do
            end do
        end do
        
        ! Apply Horner method along x1
        do k1 = spl%order(1)-1, 0, -1
            !$omp simd
            do iq = 1, spl%num_quantities
                do k3 = 0, spl%order(3)
                    do k2 = 0, spl%order(2)
                        coeff_23(iq, k2, k3) = spl%coeff(iq, k1, k2, k3, &
                            interval_index(1)+1, interval_index(2)+1, interval_index(3)+1) &
                            + x_local(1)*coeff_23(iq, k2, k3)
                    end do
                end do
            end do
        end do
        
        ! Second reduction: evaluate along x2 dimension
        ! Initialize with highest order in x2
        !$omp simd
        do iq = 1, spl%num_quantities
            do k3 = 0, spl%order(3)
                coeff_3(iq, k3) = coeff_23(iq, spl%order(2), k3)
            end do
        end do
        
        ! Apply Horner method along x2
        do k2 = spl%order(2)-1, 0, -1
            !$omp simd
            do iq = 1, spl%num_quantities
                do k3 = 0, spl%order(3)
                    coeff_3(iq, k3) = coeff_23(iq, k2, k3) + x_local(2)*coeff_3(iq, k3)
                end do
            end do
        end do
        
        ! Third reduction: evaluate along x3 dimension
        ! Initialize with highest order in x3
        !$omp simd
        do iq = 1, spl%num_quantities
            y_batch(iq) = coeff_3(iq, spl%order(3))
        end do
        
        ! Apply Horner method along x3
        do k3 = spl%order(3)-1, 0, -1
            !$omp simd
            do iq = 1, spl%num_quantities
                y_batch(iq) = coeff_3(iq, k3) + x_local(3)*y_batch(iq)
            end do
        end do
        
    end subroutine evaluate_batch_splines_3d
    
    
    subroutine evaluate_batch_splines_3d_der(spl, x, y_batch, dy_batch)
        type(BatchSplineData3D), intent(in) :: spl
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: y_batch(:)     ! (n_quantities)
        real(dp), intent(out) :: dy_batch(:,:)  ! (3, n_quantities)
        
        real(dp) :: x_norm(3), x_local(3), xj
        ! Arrays with quantities first for SIMD
        real(dp) :: coeff_23(spl%num_quantities, 0:spl%order(2), 0:spl%order(3))
        real(dp) :: coeff_23_dx1(spl%num_quantities, 0:spl%order(2), 0:spl%order(3))
        real(dp) :: coeff_3(spl%num_quantities, 0:spl%order(3))
        real(dp) :: coeff_3_dx1(spl%num_quantities, 0:spl%order(3))
        real(dp) :: coeff_3_dx2(spl%num_quantities, 0:spl%order(3))
        
        integer :: interval_index(3), k1, k2, k3, j, iq
        integer :: N1, N2, N3
        
        N1 = spl%order(1)
        N2 = spl%order(2)
        N3 = spl%order(3)
        
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
        
        ! First reduction: interpolation over x1 for value
        ! Initialize with highest order coefficient
        !$omp simd
        do iq = 1, spl%num_quantities
            do k3 = 0, N3
                do k2 = 0, N2
                    coeff_23(iq, k2, k3) = spl%coeff(iq, N1, k2, k3, &
                        interval_index(1)+1, interval_index(2)+1, interval_index(3)+1)
                end do
            end do
        end do
        
        ! Apply Horner's method for value
        do k1 = N1-1, 0, -1
            !$omp simd
            do iq = 1, spl%num_quantities
                do k3 = 0, N3
                    do k2 = 0, N2
                        coeff_23(iq, k2, k3) = spl%coeff(iq, k1, k2, k3, &
                            interval_index(1)+1, interval_index(2)+1, interval_index(3)+1) &
                            + x_local(1)*coeff_23(iq, k2, k3)
                    end do
                end do
            end do
        end do
        
        ! First derivative over x1
        !$omp simd
        do iq = 1, spl%num_quantities
            do k3 = 0, N3
                do k2 = 0, N2
                    coeff_23_dx1(iq, k2, k3) = N1 * spl%coeff(iq, N1, k2, k3, &
                        interval_index(1)+1, interval_index(2)+1, interval_index(3)+1)
                end do
            end do
        end do
        
        do k1 = N1-1, 1, -1
            !$omp simd
            do iq = 1, spl%num_quantities
                do k3 = 0, N3
                    do k2 = 0, N2
                        coeff_23_dx1(iq, k2, k3) = k1 * spl%coeff(iq, k1, k2, k3, &
                            interval_index(1)+1, interval_index(2)+1, interval_index(3)+1) &
                            + x_local(1)*coeff_23_dx1(iq, k2, k3)
                    end do
                end do
            end do
        end do
        
        ! Second reduction: interpolation over x2
        ! Initialize with highest order
        !$omp simd
        do iq = 1, spl%num_quantities
            do k3 = 0, N3
                coeff_3(iq, k3) = coeff_23(iq, N2, k3)
                coeff_3_dx1(iq, k3) = coeff_23_dx1(iq, N2, k3)
                coeff_3_dx2(iq, k3) = N2 * coeff_23(iq, N2, k3)
            end do
        end do
        
        ! Apply Horner for value and dx1
        do k2 = N2-1, 0, -1
            !$omp simd
            do iq = 1, spl%num_quantities
                do k3 = 0, N3
                    coeff_3(iq, k3) = coeff_23(iq, k2, k3) + x_local(2)*coeff_3(iq, k3)
                    coeff_3_dx1(iq, k3) = coeff_23_dx1(iq, k2, k3) + &
                                           x_local(2)*coeff_3_dx1(iq, k3)
                end do
            end do
        end do
        
        ! Derivative over x2
        do k2 = N2-1, 1, -1
            !$omp simd
            do iq = 1, spl%num_quantities
                do k3 = 0, N3
                    coeff_3_dx2(iq, k3) = k2 * coeff_23(iq, k2, k3) + &
                                           x_local(2)*coeff_3_dx2(iq, k3)
                end do
            end do
        end do
        
        ! Third reduction: interpolation over x3
        ! Initialize with highest order
        !$omp simd
        do iq = 1, spl%num_quantities
            y_batch(iq) = coeff_3(iq, N3)
            dy_batch(1, iq) = coeff_3_dx1(iq, N3)
            dy_batch(2, iq) = coeff_3_dx2(iq, N3)
            dy_batch(3, iq) = N3 * coeff_3(iq, N3)
        end do
        
        ! Apply Horner for all
        do k3 = N3-1, 0, -1
            !$omp simd
            do iq = 1, spl%num_quantities
                y_batch(iq) = coeff_3(iq, k3) + x_local(3)*y_batch(iq)
                dy_batch(1, iq) = coeff_3_dx1(iq, k3) + x_local(3)*dy_batch(1, iq)
                dy_batch(2, iq) = coeff_3_dx2(iq, k3) + x_local(3)*dy_batch(2, iq)
            end do
        end do
        
        ! Derivative over x3
        do k3 = N3-1, 1, -1
            !$omp simd
            do iq = 1, spl%num_quantities
                dy_batch(3, iq) = k3 * coeff_3(iq, k3) + x_local(3)*dy_batch(3, iq)
            end do
        end do
        
    end subroutine evaluate_batch_splines_3d_der
    
    
    subroutine evaluate_batch_splines_3d_der2(spl, x, y_batch, dy_batch, d2y_batch)
        type(BatchSplineData3D), intent(in) :: spl
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: y_batch(:)      ! (n_quantities)
        real(dp), intent(out) :: dy_batch(:,:)   ! (3, n_quantities)
        real(dp), intent(out) :: d2y_batch(:,:)  ! (6, n_quantities)
        
        real(dp) :: x_norm(3), x_local(3), xj
        ! Arrays with quantities first for SIMD
        real(dp) :: coeff_23(spl%num_quantities, 0:spl%order(2), 0:spl%order(3))
        real(dp) :: coeff_23_dx1(spl%num_quantities, 0:spl%order(2), 0:spl%order(3))
        real(dp) :: coeff_23_dx1x1(spl%num_quantities, 0:spl%order(2), 0:spl%order(3))
        
        real(dp) :: coeff_3(spl%num_quantities, 0:spl%order(3))
        real(dp) :: coeff_3_dx1(spl%num_quantities, 0:spl%order(3))
        real(dp) :: coeff_3_dx2(spl%num_quantities, 0:spl%order(3))
        real(dp) :: coeff_3_dx1x1(spl%num_quantities, 0:spl%order(3))
        real(dp) :: coeff_3_dx1x2(spl%num_quantities, 0:spl%order(3))
        real(dp) :: coeff_3_dx2x2(spl%num_quantities, 0:spl%order(3))
        
        integer :: interval_index(3), k1, k2, k3, j, iq
        integer :: N1, N2, N3
        
        N1 = spl%order(1)
        N2 = spl%order(2)
        N3 = spl%order(3)
        
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
        
        ! First reduction: interpolation over x1
        ! Initialize with k1=0 (highest degree coefficient in the original convention)
        do iq = 1, spl%num_quantities
            do k3 = 0, N3
                do k2 = 0, N2
                    coeff_23(iq, k2, k3) = spl%coeff(iq, 0, k2, k3, &
                        interval_index(1)+1, interval_index(2)+1, interval_index(3)+1)
                end do
            end do
        end do
        
        ! Apply forward Horner evaluation (matches original convention)
        do k1 = 1, N1
            !$omp simd
            do iq = 1, spl%num_quantities
                do k3 = 0, N3
                    do k2 = 0, N2
                        coeff_23(iq, k2, k3) = spl%coeff(iq, k1, k2, k3, &
                            interval_index(1)+1, interval_index(2)+1, interval_index(3)+1) &
                            + x_local(1)*coeff_23(iq, k2, k3)
                    end do
                end do
            end do
        end do
        
        ! First derivative over x1
        !$omp simd
        do iq = 1, spl%num_quantities
            do k3 = 0, N3
                do k2 = 0, N2
                    coeff_23_dx1(iq, k2, k3) = N1 * spl%coeff(iq, 0, k2, k3, &
                        interval_index(1)+1, interval_index(2)+1, interval_index(3)+1)
                end do
            end do
        end do
        
        do k1 = 1, N1-1
            !$omp simd
            do iq = 1, spl%num_quantities
                do k3 = 0, N3
                    do k2 = 0, N2
                        coeff_23_dx1(iq, k2, k3) = (N1-k1) * spl%coeff(iq, k1, k2, k3, &
                            interval_index(1)+1, interval_index(2)+1, interval_index(3)+1) &
                            + x_local(1)*coeff_23_dx1(iq, k2, k3)
                    end do
                end do
            end do
        end do
        
        ! Second derivative over x1
        !$omp simd
        do iq = 1, spl%num_quantities
            do k3 = 0, N3
                do k2 = 0, N2
                    coeff_23_dx1x1(iq, k2, k3) = N1*(N1-1) * spl%coeff(iq, 0, k2, k3, &
                        interval_index(1)+1, interval_index(2)+1, interval_index(3)+1)
                end do
            end do
        end do
        
        do k1 = 1, N1-2
            !$omp simd
            do iq = 1, spl%num_quantities
                do k3 = 0, N3
                    do k2 = 0, N2
                        coeff_23_dx1x1(iq, k2, k3) = (N1-k1)*(N1-k1-1) * spl%coeff(iq, k1, k2, k3, &
                            interval_index(1)+1, interval_index(2)+1, interval_index(3)+1) &
                            + x_local(1)*coeff_23_dx1x1(iq, k2, k3)
                    end do
                end do
            end do
        end do
        
        ! Second reduction: interpolation over x2
        ! Initialize with k2=0
        !$omp simd
        do iq = 1, spl%num_quantities
            do k3 = 0, N3
                coeff_3(iq, k3) = coeff_23(iq, 0, k3)
                coeff_3_dx1(iq, k3) = coeff_23_dx1(iq, 0, k3)
                coeff_3_dx1x1(iq, k3) = coeff_23_dx1x1(iq, 0, k3)
            end do
        end do
        
        ! Apply forward Horner for value and dx1 derivatives
        do k2 = 1, N2
            !$omp simd
            do iq = 1, spl%num_quantities
                do k3 = 0, N3
                    coeff_3(iq, k3) = coeff_23(iq, k2, k3) + x_local(2)*coeff_3(iq, k3)
                    coeff_3_dx1(iq, k3) = coeff_23_dx1(iq, k2, k3) + x_local(2)*coeff_3_dx1(iq, k3)
                    coeff_3_dx1x1(iq, k3) = coeff_23_dx1x1(iq, k2, k3) + x_local(2)*coeff_3_dx1x1(iq, k3)
                end do
            end do
        end do
        
        ! First derivative over x2 and mixed derivative
        !$omp simd
        do iq = 1, spl%num_quantities
            do k3 = 0, N3
                coeff_3_dx2(iq, k3) = N2 * coeff_23(iq, 0, k3)
                coeff_3_dx1x2(iq, k3) = N2 * coeff_23_dx1(iq, 0, k3)
                coeff_3_dx2x2(iq, k3) = N2*(N2-1) * coeff_23(iq, 0, k3)
            end do
        end do
        
        do k2 = 1, N2-1
            !$omp simd
            do iq = 1, spl%num_quantities
                do k3 = 0, N3
                    coeff_3_dx2(iq, k3) = (N2-k2) * coeff_23(iq, k2, k3) + x_local(2)*coeff_3_dx2(iq, k3)
                    coeff_3_dx1x2(iq, k3) = (N2-k2) * coeff_23_dx1(iq, k2, k3) + x_local(2)*coeff_3_dx1x2(iq, k3)
                end do
            end do
        end do
        
        ! Second derivative over x2
        do k2 = 1, N2-2
            !$omp simd
            do iq = 1, spl%num_quantities
                do k3 = 0, N3
                    coeff_3_dx2x2(iq, k3) = (N2-k2)*(N2-k2-1) * coeff_23(iq, k2, k3) &
                        + x_local(2)*coeff_3_dx2x2(iq, k3)
                end do
            end do
        end do
        
        ! Third reduction: interpolation over x3
        ! Initialize with k3=0
        !$omp simd
        do iq = 1, spl%num_quantities
            y_batch(iq) = coeff_3(iq, 0)
            dy_batch(1, iq) = coeff_3_dx1(iq, 0)
            dy_batch(2, iq) = coeff_3_dx2(iq, 0)
            d2y_batch(1, iq) = coeff_3_dx1x1(iq, 0)
            d2y_batch(2, iq) = coeff_3_dx1x2(iq, 0)
            d2y_batch(4, iq) = coeff_3_dx2x2(iq, 0)
        end do
        
        ! Apply forward Horner for value and all derivatives
        do k3 = 1, N3
            !$omp simd
            do iq = 1, spl%num_quantities
                y_batch(iq) = coeff_3(iq, k3) + x_local(3)*y_batch(iq)
                dy_batch(1, iq) = coeff_3_dx1(iq, k3) + x_local(3)*dy_batch(1, iq)
                dy_batch(2, iq) = coeff_3_dx2(iq, k3) + x_local(3)*dy_batch(2, iq)
                d2y_batch(1, iq) = coeff_3_dx1x1(iq, k3) + x_local(3)*d2y_batch(1, iq)
                d2y_batch(2, iq) = coeff_3_dx1x2(iq, k3) + x_local(3)*d2y_batch(2, iq)
                d2y_batch(4, iq) = coeff_3_dx2x2(iq, k3) + x_local(3)*d2y_batch(4, iq)
            end do
        end do
        
        ! First derivatives over x3 and mixed derivatives
        !$omp simd
        do iq = 1, spl%num_quantities
            dy_batch(3, iq) = N3 * coeff_3(iq, 0)
            d2y_batch(3, iq) = N3 * coeff_3_dx1(iq, 0)
            d2y_batch(5, iq) = N3 * coeff_3_dx2(iq, 0)
            d2y_batch(6, iq) = N3*(N3-1) * coeff_3(iq, 0)
        end do
        
        do k3 = 1, N3-1
            !$omp simd
            do iq = 1, spl%num_quantities
                dy_batch(3, iq) = (N3-k3) * coeff_3(iq, k3) + x_local(3)*dy_batch(3, iq)
                d2y_batch(3, iq) = (N3-k3) * coeff_3_dx1(iq, k3) + x_local(3)*d2y_batch(3, iq)
                d2y_batch(5, iq) = (N3-k3) * coeff_3_dx2(iq, k3) + x_local(3)*d2y_batch(5, iq)
            end do
        end do
        
        ! Second derivative over x3
        do k3 = 1, N3-2
            !$omp simd
            do iq = 1, spl%num_quantities
                d2y_batch(6, iq) = (N3-k3)*(N3-k3-1) * coeff_3(iq, k3) + x_local(3)*d2y_batch(6, iq)
            end do
        end do
        
    end subroutine evaluate_batch_splines_3d_der2
    
end module interpolate
