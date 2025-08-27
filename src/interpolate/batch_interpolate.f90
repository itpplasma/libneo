module batch_interpolate
    use spl_three_to_five_sub
    
    implicit none
    private
    
    integer, parameter :: dp = kind(1.0d0)
    
    ! Export batch spline types
    public :: BatchSplineData1D, BatchSplineData2D, BatchSplineData3D
    
    ! Export batch spline construction routines
    public :: construct_batch_splines_1d, construct_batch_splines_2d, construct_batch_splines_3d
    public :: destroy_batch_splines_1d, destroy_batch_splines_2d, destroy_batch_splines_3d
    
    ! Export batch spline evaluation routines
    public :: evaluate_batch_splines_1d, evaluate_batch_splines_1d_single
    public :: evaluate_batch_splines_1d_der, evaluate_batch_splines_1d_der2
    public :: evaluate_batch_splines_2d, evaluate_batch_splines_2d_der
    public :: evaluate_batch_splines_3d, evaluate_batch_splines_3d_der, evaluate_batch_splines_3d_der2
    
    ! Batch spline types for multiple quantities on shared grid
    type :: BatchSplineData1D
        ! Shared grid data
        integer :: order
        integer :: num_points  
        logical :: periodic
        real(dp) :: x_min
        real(dp) :: h_step
        
        ! Batch data
        integer :: num_quantities
        ! Memory layout: (n_quantities, 0:order, n_points) for cache efficiency
        real(dp), dimension(:,:,:), allocatable :: coeff  
    end type BatchSplineData1D
    
    type :: BatchSplineData2D
        ! Shared grid data
        integer :: order(2)
        integer :: num_points(2)
        logical :: periodic(2)
        real(dp) :: h_step(2)
        real(dp) :: x_min(2)
        
        ! Batch data
        integer :: num_quantities
        ! Memory layout: (n_quantities, 0:order1, 0:order2, n1, n2) for cache efficiency
        real(dp), dimension(:,:,:,:,:), allocatable :: coeff
    end type BatchSplineData2D
    
    type :: BatchSplineData3D
        ! Shared grid data
        integer :: order(3)
        integer :: num_points(3)
        logical :: periodic(3)
        real(dp) :: h_step(3)
        real(dp) :: x_min(3)
        
        ! Batch data
        integer :: num_quantities
        ! Memory layout: (n_quantities, 0:order1, 0:order2, 0:order3, n1, n2, n3)
        real(dp), dimension(:,:,:,:,:,:,:), allocatable :: coeff
    end type BatchSplineData3D
    
contains
    
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
        
        ! Store metadata
        spl%order = order
        spl%num_points = n_points
        spl%periodic = periodic
        spl%x_min = x_min
        spl%h_step = (x_max - x_min) / dble(n_points - 1)
        spl%num_quantities = n_quantities
        
        ! Allocate batch coefficients with cache-friendly layout
        allocate(spl%coeff(n_quantities, 0:order, n_points))
        
        ! Temporary array for single spline construction
        allocate(splcoe_temp(0:order, n_points))
        
        ! Process each quantity
        do iq = 1, n_quantities
            ! Set the data values
            splcoe_temp(0, :) = y_batch(:, iq)
            
            ! Construct spline coefficients
            if (periodic) then
                call spl_per(order, n_points, spl%h_step, splcoe_temp)
            else
                call spl_reg(order, n_points, spl%h_step, splcoe_temp)
            end if
            
            ! Store coefficients in batch array with quantities first
            spl%coeff(iq, :, :) = splcoe_temp(:, :)
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
        
        ! Handle periodic boundary conditions
        if (spl%periodic) then
            xj = modulo(x, spl%h_step*(spl%num_points-1))
        else
            xj = x
        end if
        
        ! Find interval and local coordinate
        x_norm = (xj - spl%x_min)/spl%h_step
        interval_index = max(0, min(spl%num_points-1, int(x_norm)))
        x_local = (x_norm - dble(interval_index))*spl%h_step
        
        ! Initialize with highest order coefficient
        !$omp simd
        do iq = 1, spl%num_quantities
            y_batch(iq) = spl%coeff(iq, spl%order, interval_index+1)
        end do
        
        ! Apply Horner method (backward)
        do k_power = spl%order-1, 0, -1
            !$omp simd
            do iq = 1, spl%num_quantities
                y_batch(iq) = spl%coeff(iq, k_power, interval_index+1) + &
                              x_local * y_batch(iq)
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
        
        ! Handle periodic boundary conditions
        if (spl%periodic) then
            xj = modulo(x, spl%h_step*(spl%num_points-1))
        else
            xj = x
        end if
        
        ! Find interval and local coordinate
        x_norm = (xj - spl%x_min)/spl%h_step
        interval_index = max(0, min(spl%num_points-1, int(x_norm)))
        x_local = (x_norm - dble(interval_index))*spl%h_step
        
        ! Horner method for single quantity
        y = spl%coeff(iq, spl%order, interval_index+1)
        do k_power = spl%order-1, 0, -1
            y = spl%coeff(iq, k_power, interval_index+1) + x_local * y
        end do
    end subroutine evaluate_batch_splines_1d_single
    
    
    subroutine evaluate_batch_splines_1d_der(spl, x, y_batch, dy_batch)
        type(BatchSplineData1D), intent(in) :: spl
        real(dp), intent(in) :: x
        real(dp), intent(out) :: y_batch(:)   ! (n_quantities)
        real(dp), intent(out) :: dy_batch(:)  ! (n_quantities)
        
        real(dp) :: x_norm, x_local, xj
        integer :: interval_index, k_power, iq
        integer :: N
        
        N = spl%order
        
        ! Handle periodic boundary conditions
        if (spl%periodic) then
            xj = modulo(x, spl%h_step*(spl%num_points-1))
        else
            xj = x
        end if
        
        ! Find interval and local coordinate
        x_norm = (xj - spl%x_min)/spl%h_step
        interval_index = max(0, min(spl%num_points-1, int(x_norm)))
        x_local = (x_norm - dble(interval_index))*spl%h_step
        
        ! Initialize value with highest order coefficient
        !$omp simd
        do iq = 1, spl%num_quantities
            y_batch(iq) = spl%coeff(iq, N, interval_index+1)
        end do
        
        ! Apply Horner method for value
        do k_power = N-1, 0, -1
            !$omp simd
            do iq = 1, spl%num_quantities
                y_batch(iq) = spl%coeff(iq, k_power, interval_index+1) + &
                              x_local * y_batch(iq)
            end do
        end do
        
        ! Initialize derivative with highest order coefficient
        !$omp simd
        do iq = 1, spl%num_quantities
            dy_batch(iq) = N * spl%coeff(iq, N, interval_index+1)
        end do
        
        ! Apply Horner method for derivative
        do k_power = N-1, 1, -1
            !$omp simd
            do iq = 1, spl%num_quantities
                dy_batch(iq) = k_power * spl%coeff(iq, k_power, interval_index+1) + &
                               x_local * dy_batch(iq)
            end do
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
        integer :: N
        
        N = spl%order
        
        ! Handle periodic boundary conditions
        if (spl%periodic) then
            xj = modulo(x, spl%h_step*(spl%num_points-1))
        else
            xj = x
        end if
        
        ! Find interval and local coordinate
        x_norm = (xj - spl%x_min)/spl%h_step
        interval_index = max(0, min(spl%num_points-1, int(x_norm)))
        x_local = (x_norm - dble(interval_index))*spl%h_step
        
        ! Initialize value with highest order coefficient
        !$omp simd
        do iq = 1, spl%num_quantities
            y_batch(iq) = spl%coeff(iq, N, interval_index+1)
        end do
        
        ! Apply Horner method for value
        do k_power = N-1, 0, -1
            !$omp simd
            do iq = 1, spl%num_quantities
                y_batch(iq) = spl%coeff(iq, k_power, interval_index+1) + &
                              x_local * y_batch(iq)
            end do
        end do
        
        ! Initialize first derivative
        !$omp simd
        do iq = 1, spl%num_quantities
            dy_batch(iq) = N * spl%coeff(iq, N, interval_index+1)
        end do
        
        ! Apply Horner method for first derivative
        do k_power = N-1, 1, -1
            !$omp simd
            do iq = 1, spl%num_quantities
                dy_batch(iq) = k_power * spl%coeff(iq, k_power, interval_index+1) + &
                               x_local * dy_batch(iq)
            end do
        end do
        
        ! Initialize second derivative
        !$omp simd
        do iq = 1, spl%num_quantities
            d2y_batch(iq) = N*(N-1) * spl%coeff(iq, N, interval_index+1)
        end do
        
        ! Apply Horner method for second derivative
        do k_power = N-1, 2, -1
            !$omp simd
            do iq = 1, spl%num_quantities
                d2y_batch(iq) = k_power*(k_power-1) * spl%coeff(iq, k_power, interval_index+1) + &
                                x_local * d2y_batch(iq)
            end do
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
        integer :: n1, n2, n_quantities
        integer :: N1_order, N2_order
        
        n1 = size(y_batch, 1)
        n2 = size(y_batch, 2)
        n_quantities = size(y_batch, 3)
        N1_order = order(1)
        N2_order = order(2)
        
        ! Store metadata
        spl%order = order
        spl%num_points = [n1, n2]
        spl%periodic = periodic
        spl%h_step = [(x_max(1) - x_min(1))/dble(n1-1), &
                      (x_max(2) - x_min(2))/dble(n2-1)]
        spl%x_min = x_min
        spl%num_quantities = n_quantities
        
        ! Allocate batch coefficients with cache-friendly layout
        allocate(spl%coeff(n_quantities, 0:N1_order, 0:N2_order, n1, n2))
        allocate(temp_coeff(0:N1_order, 0:N2_order, n1, n2))
        
        ! Process each quantity
        do iq = 1, n_quantities
            ! Spline over x2 (second dimension)
            allocate(splcoe(0:N2_order, n2))
            do i1 = 1, n1
                splcoe(0,:) = y_batch(i1, :, iq)
                if (periodic(2)) then
                    call spl_per(N2_order, n2, spl%h_step(2), splcoe)
                else
                    call spl_reg(N2_order, n2, spl%h_step(2), splcoe)
                end if
                temp_coeff(0, :, i1, :) = splcoe
            end do
            
            deallocate(splcoe)
            allocate(splcoe(0:N1_order, n1))
            
            ! Spline over x1 (first dimension)
            do i2 = 1, n2
                do k2 = 0, N2_order
                    splcoe(0,:) = temp_coeff(0, k2, :, i2)
                    if (periodic(1)) then
                        call spl_per(N1_order, n1, spl%h_step(1), splcoe)
                    else
                        call spl_reg(N1_order, n1, spl%h_step(1), splcoe)
                    end if
                    temp_coeff(:, k2, :, i2) = splcoe
                end do
            end do
            
            ! Store in batch array with quantities first
            spl%coeff(iq, :, :, :, :) = temp_coeff
            
            deallocate(splcoe)
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
        integer :: N1, N2
        
        N1 = spl%order(1)
        N2 = spl%order(2)
        
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
        
        ! First reduction: interpolation over x1
        ! Initialize with highest order
        !$omp simd
        do iq = 1, spl%num_quantities
            do k2 = 0, N2
                coeff_2(iq, k2) = spl%coeff(iq, N1, k2, &
                    interval_index(1)+1, interval_index(2)+1)
            end do
        end do
        
        ! Apply Horner method
        do k1 = N1-1, 0, -1
            !$omp simd
            do iq = 1, spl%num_quantities
                do k2 = 0, N2
                    coeff_2(iq, k2) = spl%coeff(iq, k1, k2, &
                        interval_index(1)+1, interval_index(2)+1) + &
                        x_local(1) * coeff_2(iq, k2)
                end do
            end do
        end do
        
        ! Second reduction: interpolation over x2
        ! Initialize with highest order
        !$omp simd
        do iq = 1, spl%num_quantities
            y_batch(iq) = coeff_2(iq, N2)
        end do
        
        ! Apply Horner method
        do k2 = N2-1, 0, -1
            !$omp simd
            do iq = 1, spl%num_quantities
                y_batch(iq) = coeff_2(iq, k2) + x_local(2) * y_batch(iq)
            end do
        end do
    end subroutine evaluate_batch_splines_2d
    
    
    subroutine evaluate_batch_splines_2d_der(spl, x, y_batch, dy_batch)
        type(BatchSplineData2D), intent(in) :: spl
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: y_batch(:)     ! (n_quantities)
        real(dp), intent(out) :: dy_batch(:,:)  ! (2, n_quantities)
        
        real(dp) :: x_norm(2), x_local(2), xj
        real(dp) :: coeff_2(spl%num_quantities, 0:spl%order(2))
        real(dp) :: coeff_2_dx1(spl%num_quantities, 0:spl%order(2))
        integer :: interval_index(2), k1, k2, j, iq
        integer :: N1, N2
        
        N1 = spl%order(1)
        N2 = spl%order(2)
        
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
        
        ! First reduction: interpolation over x1 for value
        !$omp simd
        do iq = 1, spl%num_quantities
            do k2 = 0, N2
                coeff_2(iq, k2) = spl%coeff(iq, N1, k2, &
                    interval_index(1)+1, interval_index(2)+1)
            end do
        end do
        
        do k1 = N1-1, 0, -1
            !$omp simd
            do iq = 1, spl%num_quantities
                do k2 = 0, N2
                    coeff_2(iq, k2) = spl%coeff(iq, k1, k2, &
                        interval_index(1)+1, interval_index(2)+1) + &
                        x_local(1) * coeff_2(iq, k2)
                end do
            end do
        end do
        
        ! First derivative over x1
        !$omp simd
        do iq = 1, spl%num_quantities
            do k2 = 0, N2
                coeff_2_dx1(iq, k2) = N1 * spl%coeff(iq, N1, k2, &
                    interval_index(1)+1, interval_index(2)+1)
            end do
        end do
        
        do k1 = N1-1, 1, -1
            !$omp simd
            do iq = 1, spl%num_quantities
                do k2 = 0, N2
                    coeff_2_dx1(iq, k2) = k1 * spl%coeff(iq, k1, k2, &
                        interval_index(1)+1, interval_index(2)+1) + &
                        x_local(1) * coeff_2_dx1(iq, k2)
                end do
            end do
        end do
        
        ! Second reduction: interpolation over x2
        !$omp simd
        do iq = 1, spl%num_quantities
            y_batch(iq) = coeff_2(iq, N2)
            dy_batch(1, iq) = coeff_2_dx1(iq, N2)
            dy_batch(2, iq) = N2 * coeff_2(iq, N2)
        end do
        
        do k2 = N2-1, 0, -1
            !$omp simd
            do iq = 1, spl%num_quantities
                y_batch(iq) = coeff_2(iq, k2) + x_local(2) * y_batch(iq)
                dy_batch(1, iq) = coeff_2_dx1(iq, k2) + x_local(2) * dy_batch(1, iq)
            end do
        end do
        
        do k2 = N2-1, 1, -1
            !$omp simd
            do iq = 1, spl%num_quantities
                dy_batch(2, iq) = k2 * coeff_2(iq, k2) + x_local(2) * dy_batch(2, iq)
            end do
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
        integer :: n1, n2, n3, n_quantities
        integer :: N1_order, N2_order, N3_order
        
        n1 = size(y_batch, 1)
        n2 = size(y_batch, 2)
        n3 = size(y_batch, 3)
        n_quantities = size(y_batch, 4)
        N1_order = order(1)
        N2_order = order(2)
        N3_order = order(3)
        
        ! Store metadata
        spl%order = order
        spl%num_points = [n1, n2, n3]
        spl%periodic = periodic
        spl%h_step = [(x_max(1) - x_min(1))/dble(n1-1), &
                      (x_max(2) - x_min(2))/dble(n2-1), &
                      (x_max(3) - x_min(3))/dble(n3-1)]
        spl%x_min = x_min
        spl%num_quantities = n_quantities
        
        ! Allocate batch coefficients with cache-friendly layout
        allocate(spl%coeff(n_quantities, 0:N1_order, 0:N2_order, 0:N3_order, n1, n2, n3))
        allocate(temp_coeff(0:N1_order, 0:N2_order, 0:N3_order, n1, n2, n3))
        
        ! Process each quantity
        do iq = 1, n_quantities
            ! Spline over x3 first (like original)
            allocate(splcoe(0:N3_order, n3))
            do i2 = 1, n2
                do i1 = 1, n1
                    splcoe(0,:) = y_batch(i1, i2, :, iq)
                    if (periodic(3)) then
                        call spl_per(N3_order, n3, spl%h_step(3), splcoe)
                    else
                        call spl_reg(N3_order, n3, spl%h_step(3), splcoe)
                    end if
                    temp_coeff(N1_order, 0, :, i1, i2, :) = splcoe
                end do
            end do
            deallocate(splcoe)
            
            ! Spline over x2 second
            allocate(splcoe(0:N2_order, n2))
            do i3 = 1, n3
                do i1 = 1, n1
                    do k3 = 0, N3_order
                        splcoe(0,:) = temp_coeff(N1_order, 0, k3, i1, :, i3)
                        if (periodic(2)) then
                            call spl_per(N2_order, n2, spl%h_step(2), splcoe)
                        else
                            call spl_reg(N2_order, n2, spl%h_step(2), splcoe)
                        end if
                        temp_coeff(N1_order, :, k3, i1, :, i3) = splcoe
                    end do
                end do
            end do
            deallocate(splcoe)
            
            ! Spline over x1 last
            allocate(splcoe(0:N1_order, n1))
            do i3 = 1, n3
                do i2 = 1, n2
                    do k3 = 0, N3_order
                        do k2 = 0, N2_order
                            splcoe(0,:) = temp_coeff(N1_order, k2, k3, :, i2, i3)
                            if (periodic(1)) then
                                call spl_per(N1_order, n1, spl%h_step(1), splcoe)
                            else
                                call spl_reg(N1_order, n1, spl%h_step(1), splcoe)
                            end if
                            ! Store reversed for consistency with original 3D (reversed indexing)
                            temp_coeff(N1_order:0:-1, k2, k3, :, i2, i3) = splcoe
                        end do
                    end do
                end do
            end do
            deallocate(splcoe)
            
            ! Copy to final coefficient array
            do i3 = 1, n3
                do i2 = 1, n2
                    do i1 = 1, n1
                        spl%coeff(iq, :, :, :, i1, i2, i3) = &
                            temp_coeff(:, :, :, i1, i2, i3)
                    end do
                end do
            end do
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
        real(dp) :: coeff_23(spl%num_quantities, 0:spl%order(2), 0:spl%order(3))
        real(dp) :: coeff_3(spl%num_quantities, 0:spl%order(3))
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
        
        ! First reduction: interpolation over x1 (matching original exactly)
        ! Initialize with index 0 (which is highest order due to reversed storage)
        !$omp simd
        do iq = 1, spl%num_quantities
            do k3 = 0, N3
                do k2 = 0, N2
                    coeff_23(iq, k2, k3) = spl%coeff(iq, 0, k2, k3, &
                        interval_index(1)+1, interval_index(2)+1, interval_index(3)+1)
                end do
            end do
        end do
        
        ! Apply forward iteration (descending powers due to reversed storage)
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
        
        ! Second reduction: interpolation over x2
        ! Initialize with highest order
        !$omp simd
        do iq = 1, spl%num_quantities
            do k3 = 0, N3
                coeff_3(iq, k3) = coeff_23(iq, N2, k3)
            end do
        end do
        
        ! Apply Horner method
        do k2 = N2-1, 0, -1
            !$omp simd
            do iq = 1, spl%num_quantities
                do k3 = 0, N3
                    coeff_3(iq, k3) = coeff_23(iq, k2, k3) + x_local(2)*coeff_3(iq, k3)
                end do
            end do
        end do
        
        ! Third reduction: interpolation over x3
        ! Initialize with highest order
        !$omp simd
        do iq = 1, spl%num_quantities
            y_batch(iq) = coeff_3(iq, N3)
        end do
        
        ! Apply Horner method
        do k3 = N3-1, 0, -1
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
        
        ! First reduction: interpolation over x1 for value (matching original exactly)
        ! Initialize with index 0 (highest order due to reversed storage)
        !$omp simd
        do iq = 1, spl%num_quantities
            do k3 = 0, N3
                do k2 = 0, N2
                    coeff_23(iq, k2, k3) = spl%coeff(iq, 0, k2, k3, &
                        interval_index(1)+1, interval_index(2)+1, interval_index(3)+1)
                end do
            end do
        end do
        
        ! Apply forward iteration for value (descending powers due to reversed storage)
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
        
        ! First derivative over x1 (matching original exactly)
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
        
        ! First reduction: interpolation over x1 (matching original exactly)
        ! Initialize with index 0 (highest order due to reversed storage)
        !$omp simd
        do iq = 1, spl%num_quantities
            do k3 = 0, N3
                do k2 = 0, N2
                    coeff_23(iq, k2, k3) = spl%coeff(iq, 0, k2, k3, &
                        interval_index(1)+1, interval_index(2)+1, interval_index(3)+1)
                end do
            end do
        end do
        
        ! Apply forward iteration for value (descending powers due to reversed storage)
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
        
        ! First derivative over x1 (matching original exactly)
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
        
        ! Second derivative over x1 (matching original pattern)
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
        
        ! Second reduction: interpolation over x2 using backward Horner
        ! Initialize with highest degree
        !$omp simd
        do iq = 1, spl%num_quantities
            do k3 = 0, N3
                coeff_3(iq, k3) = coeff_23(iq, N2, k3)
                coeff_3_dx1(iq, k3) = coeff_23_dx1(iq, N2, k3)
                coeff_3_dx1x1(iq, k3) = coeff_23_dx1x1(iq, N2, k3)
            end do
        end do
        
        ! Apply backward Horner for value and dx1 derivatives
        do k2 = N2-1, 0, -1
            !$omp simd
            do iq = 1, spl%num_quantities
                do k3 = 0, N3
                    coeff_3(iq, k3) = coeff_23(iq, k2, k3) + x_local(2)*coeff_3(iq, k3)
                    coeff_3_dx1(iq, k3) = coeff_23_dx1(iq, k2, k3) + x_local(2)*coeff_3_dx1(iq, k3)
                    coeff_3_dx1x1(iq, k3) = coeff_23_dx1x1(iq, k2, k3) + x_local(2)*coeff_3_dx1x1(iq, k3)
                end do
            end do
        end do
        
        ! First derivative over x2 and mixed derivative using backward Horner
        !$omp simd
        do iq = 1, spl%num_quantities
            do k3 = 0, N3
                coeff_3_dx2(iq, k3) = N2 * coeff_23(iq, N2, k3)
                coeff_3_dx1x2(iq, k3) = N2 * coeff_23_dx1(iq, N2, k3)
                coeff_3_dx2x2(iq, k3) = N2*(N2-1) * coeff_23(iq, N2, k3)
            end do
        end do
        
        do k2 = N2-1, 1, -1
            !$omp simd
            do iq = 1, spl%num_quantities
                do k3 = 0, N3
                    coeff_3_dx2(iq, k3) = k2 * coeff_23(iq, k2, k3) + x_local(2)*coeff_3_dx2(iq, k3)
                    coeff_3_dx1x2(iq, k3) = k2 * coeff_23_dx1(iq, k2, k3) + x_local(2)*coeff_3_dx1x2(iq, k3)
                end do
            end do
        end do
        
        ! Second derivative over x2 using backward Horner
        do k2 = N2-1, 2, -1
            !$omp simd
            do iq = 1, spl%num_quantities
                do k3 = 0, N3
                    coeff_3_dx2x2(iq, k3) = k2*(k2-1) * coeff_23(iq, k2, k3) &
                        + x_local(2)*coeff_3_dx2x2(iq, k3)
                end do
            end do
        end do
        
        ! Third reduction: interpolation over x3 using backward Horner
        ! Initialize with highest degree
        !$omp simd
        do iq = 1, spl%num_quantities
            y_batch(iq) = coeff_3(iq, N3)
            dy_batch(1, iq) = coeff_3_dx1(iq, N3)
            dy_batch(2, iq) = coeff_3_dx2(iq, N3)
            d2y_batch(1, iq) = coeff_3_dx1x1(iq, N3)
            d2y_batch(2, iq) = coeff_3_dx1x2(iq, N3)
            d2y_batch(4, iq) = coeff_3_dx2x2(iq, N3)
        end do
        
        ! Apply backward Horner for value and all derivatives
        do k3 = N3-1, 0, -1
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
        
        ! First derivatives over x3 and mixed derivatives using backward Horner
        !$omp simd
        do iq = 1, spl%num_quantities
            dy_batch(3, iq) = N3 * coeff_3(iq, N3)
            d2y_batch(3, iq) = N3 * coeff_3_dx1(iq, N3)
            d2y_batch(5, iq) = N3 * coeff_3_dx2(iq, N3)
            d2y_batch(6, iq) = N3*(N3-1) * coeff_3(iq, N3)
        end do
        
        do k3 = N3-1, 1, -1
            !$omp simd
            do iq = 1, spl%num_quantities
                dy_batch(3, iq) = k3 * coeff_3(iq, k3) + x_local(3)*dy_batch(3, iq)
                d2y_batch(3, iq) = k3 * coeff_3_dx1(iq, k3) + x_local(3)*d2y_batch(3, iq)
                d2y_batch(5, iq) = k3 * coeff_3_dx2(iq, k3) + x_local(3)*d2y_batch(5, iq)
            end do
        end do
        
        ! Second derivative over x3 using backward Horner
        do k3 = N3-1, 2, -1
            !$omp simd
            do iq = 1, spl%num_quantities
                d2y_batch(6, iq) = k3*(k3-1) * coeff_3(iq, k3) + x_local(3)*d2y_batch(6, iq)
            end do
        end do
        
    end subroutine evaluate_batch_splines_3d_der2
    
end module batch_interpolate