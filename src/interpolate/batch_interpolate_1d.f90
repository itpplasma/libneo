module batch_interpolate_1d
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use batch_interpolate_types, only: BatchSplineData1D
    use spl_three_to_five_sub, only: spl_per, spl_reg
#ifdef _OPENACC
    use openacc, only: acc_is_present
#endif

    implicit none
    private
    
    ! Export batch spline construction/destruction routines
    public :: construct_batch_splines_1d
    public :: construct_batch_splines_1d_resident
    public :: destroy_batch_splines_1d
    
    ! Export batch spline evaluation routines
    public :: evaluate_batch_splines_1d
    public :: evaluate_batch_splines_1d_single
    public :: evaluate_batch_splines_1d_many
    public :: evaluate_batch_splines_1d_many_resident
    public :: evaluate_batch_splines_1d_der
    public :: evaluate_batch_splines_1d_der2
    
contains
    
    subroutine construct_batch_splines_1d(x_min, x_max, y_batch, order, periodic, spl)
        real(dp), intent(in) :: x_min, x_max
        real(dp), intent(in) :: y_batch(:,:)  ! (n_points, n_quantities)
        integer, intent(in) :: order
        logical, intent(in) :: periodic
        type(BatchSplineData1D), intent(out) :: spl
        
        integer :: iq, n_points, n_quantities, istat
        real(dp), dimension(:,:), allocatable :: splcoe_temp
        
        n_points = size(y_batch, 1)
        n_quantities = size(y_batch, 2)
        
        ! Validate input
        if (n_points < 2) then
            error stop "construct_batch_splines_1d: Need at least 2 points"
        end if
        if (n_quantities < 1) then
            error stop "construct_batch_splines_1d: Need at least 1 quantity"
        end if
        if (order < 3 .or. order > 5) then
            error stop "construct_batch_splines_1d: Order must be between 3 and 5"
        end if
        
        ! Store metadata
        spl%order = order
        spl%num_points = n_points
        spl%periodic = periodic
        spl%x_min = x_min
        spl%h_step = (x_max - x_min) / dble(n_points - 1)
        spl%num_quantities = n_quantities
        
        ! Allocate batch coefficients with cache-friendly layout
        allocate(spl%coeff(n_quantities, 0:order, n_points), stat=istat)
        if (istat /= 0) then
            error stop "construct_batch_splines_1d: Allocation failed for coeff"
        end if
        
        ! Temporary array for single spline construction
        allocate(splcoe_temp(0:order, n_points), stat=istat)
        if (istat /= 0) then
            error stop "construct_batch_splines_1d: Allocation failed for splcoe_temp"
        end if
        
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

    subroutine construct_batch_splines_1d_resident(x_min, x_max, y_batch, order, &
                                                   periodic, spl)
        real(dp), intent(in) :: x_min, x_max
        real(dp), intent(in) :: y_batch(:, :)  ! (n_points, n_quantities)
        integer, intent(in) :: order
        logical, intent(in) :: periodic
        type(BatchSplineData1D), intent(out) :: spl

        call construct_batch_splines_1d(x_min, x_max, y_batch, order, periodic, spl)
#ifdef _OPENACC
        !$acc enter data copyin(spl%coeff)
#endif
    end subroutine construct_batch_splines_1d_resident


    subroutine destroy_batch_splines_1d(spl)
        type(BatchSplineData1D), intent(inout) :: spl

#ifdef _OPENACC
        if (allocated(spl%coeff)) then
            if (acc_is_present(spl%coeff)) then
                !$acc exit data delete(spl%coeff)
            end if
        end if
#endif
        if (allocated(spl%coeff)) deallocate(spl%coeff)
    end subroutine destroy_batch_splines_1d
    
    
    subroutine evaluate_batch_splines_1d(spl, x, y_batch)
        type(BatchSplineData1D), intent(in) :: spl
        real(dp), intent(in) :: x
        real(dp), intent(out) :: y_batch(:)  ! (n_quantities)
        
        real(dp) :: x_norm, x_local, xj
        integer :: interval_index, k_power, iq
        
        ! Validate output size
        if (size(y_batch) < spl%num_quantities) then
            error stop "evaluate_batch_splines_1d: Output array too small"
        end if
        
        ! Handle periodic boundary conditions
        if (spl%periodic) then
            xj = modulo(x - spl%x_min, spl%h_step*(spl%num_points-1)) + spl%x_min
        else
            xj = x
        end if
        
        ! Find interval and local coordinate
        x_norm = (xj - spl%x_min)/spl%h_step
        interval_index = max(0, min(spl%num_points-2, int(x_norm)))
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
        
        ! Validate input
        if (iq < 1 .or. iq > spl%num_quantities) then
            error stop "evaluate_batch_splines_1d_single: Invalid quantity index"
        end if
        
        ! Handle periodic boundary conditions
        if (spl%periodic) then
            xj = modulo(x - spl%x_min, spl%h_step*(spl%num_points-1)) + spl%x_min
        else
            xj = x
        end if
        
        ! Find interval and local coordinate
        x_norm = (xj - spl%x_min)/spl%h_step
        interval_index = max(0, min(spl%num_points-2, int(x_norm)))
        x_local = (x_norm - dble(interval_index))*spl%h_step
        
        ! Horner method for single quantity
        y = spl%coeff(iq, spl%order, interval_index+1)
        do k_power = spl%order-1, 0, -1
            y = spl%coeff(iq, k_power, interval_index+1) + x_local * y
        end do
    end subroutine evaluate_batch_splines_1d_single

    subroutine evaluate_batch_splines_1d_many(spl, x, y_batch)
        type(BatchSplineData1D), intent(in) :: spl
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y_batch(:, :)  ! (n_quantities, n_points)

        integer :: ipt, iq, k_power, idx, k_wrap
        integer :: num_points, nq, order
        real(dp) :: xj, x_norm, x_local, x_min, h_step, period, t, w

        if (size(y_batch, 1) < spl%num_quantities) then
            error stop "evaluate_batch_splines_1d_many: First dimension too small"
        end if
        if (size(y_batch, 2) /= size(x)) then
            error stop "evaluate_batch_splines_1d_many: y_batch second dim mismatch"
        end if

        order = spl%order
        num_points = spl%num_points
        nq = spl%num_quantities
        x_min = spl%x_min
        h_step = spl%h_step
        period = h_step*real(num_points - 1, dp)

        do ipt = 1, size(x)
            include "spline1d_many_point_body.inc"
        end do
    end subroutine evaluate_batch_splines_1d_many

    subroutine evaluate_batch_splines_1d_many_resident(spl, x, y_batch)
        type(BatchSplineData1D), intent(in) :: spl
        real(dp), intent(in) :: x(:)
        real(dp), intent(inout) :: y_batch(:, :)  ! (n_quantities, n_points)

        integer :: ipt, iq, k_power, idx, k_wrap
        integer :: num_points, nq, order
        real(dp) :: xj, x_norm, x_local, x_min, h_step, period, t, w

        if (size(y_batch, 1) < spl%num_quantities) then
            error stop "evaluate_batch_splines_1d_many_resident: y_batch dim1 too small"
        end if
        if (size(y_batch, 2) /= size(x)) then
            error stop "evaluate_batch_splines_1d_many_resident: y_batch dim2 mismatch"
        end if

        order = spl%order
        num_points = spl%num_points
        nq = spl%num_quantities
        x_min = spl%x_min
        h_step = spl%h_step
        period = h_step*real(num_points - 1, dp)

        !$acc parallel loop present(spl%coeff, x, y_batch) &
        !$acc& private(ipt, iq, k_power, idx, k_wrap, xj, x_norm, x_local, t, w) &
        !$acc& gang vector vector_length(256)
        do ipt = 1, size(x)
            include "spline1d_many_point_body.inc"
        end do
        !$acc end parallel loop
    end subroutine evaluate_batch_splines_1d_many_resident
    
    
    subroutine evaluate_batch_splines_1d_der(spl, x, y_batch, dy_batch)
        type(BatchSplineData1D), intent(in) :: spl
        real(dp), intent(in) :: x
        real(dp), intent(out) :: y_batch(:)   ! (n_quantities)
        real(dp), intent(out) :: dy_batch(:)  ! (n_quantities)
        
        real(dp) :: x_norm, x_local, xj
        integer :: interval_index, k_power, iq
        integer :: N
        
        N = spl%order
        
        ! Validate output sizes
        if (size(y_batch) < spl%num_quantities .or. &
            size(dy_batch) < spl%num_quantities) then
            error stop "evaluate_batch_splines_1d_der: Output arrays too small"
        end if
        
        ! Handle periodic boundary conditions
        if (spl%periodic) then
            xj = modulo(x - spl%x_min, spl%h_step*(spl%num_points-1)) + spl%x_min
        else
            xj = x
        end if
        
        ! Find interval and local coordinate
        x_norm = (xj - spl%x_min)/spl%h_step
        interval_index = max(0, min(spl%num_points-2, int(x_norm)))
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
        
        ! Validate output sizes
        if (size(y_batch) < spl%num_quantities .or. &
            size(dy_batch) < spl%num_quantities .or. &
            size(d2y_batch) < spl%num_quantities) then
            error stop "evaluate_batch_splines_1d_der2: Output arrays too small"
        end if
        
        ! Handle periodic boundary conditions
        if (spl%periodic) then
            xj = modulo(x - spl%x_min, spl%h_step*(spl%num_points-1)) + spl%x_min
        else
            xj = x
        end if
        
        ! Find interval and local coordinate
        x_norm = (xj - spl%x_min)/spl%h_step
        interval_index = max(0, min(spl%num_points-2, int(x_norm)))
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
                d2y_batch(iq) = k_power*(k_power-1) * &
                    spl%coeff(iq, k_power, interval_index+1) + &
                    x_local * d2y_batch(iq)
            end do
        end do
    end subroutine evaluate_batch_splines_1d_der2
    
end module batch_interpolate_1d
