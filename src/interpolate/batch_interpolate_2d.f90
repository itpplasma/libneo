module batch_interpolate_2d
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use batch_interpolate_types, only: BatchSplineData2D
    use spl_three_to_five_sub, only: spl_per, spl_reg

    implicit none
    private
    
    ! Export batch spline construction/destruction routines
    public :: construct_batch_splines_2d
    public :: destroy_batch_splines_2d
    
    ! Export batch spline evaluation routines
    public :: evaluate_batch_splines_2d
    public :: evaluate_batch_splines_2d_der
    public :: evaluate_batch_splines_2d_many
    public :: evaluate_batch_splines_2d_many_resident
    
contains
    
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
        integer :: n1, n2, n_quantities, istat
        integer :: N1_order, N2_order
        
        n1 = size(y_batch, 1)
        n2 = size(y_batch, 2)
        n_quantities = size(y_batch, 3)
        N1_order = order(1)
        N2_order = order(2)
        
        ! Validate input
        if (n1 < 2 .or. n2 < 2) then
            error stop "construct_batch_splines_2d: Need at least 2 points in each dimension"
        end if
        if (n_quantities < 1) then
            error stop "construct_batch_splines_2d: Need at least 1 quantity"
        end if
        if (N1_order < 3 .or. N1_order > 5 .or. N2_order < 3 .or. N2_order > 5) then
            error stop "construct_batch_splines_2d: Orders must be between 3 and 5"
        end if
        
        ! Store metadata
        spl%order = order
        spl%num_points = [n1, n2]
        spl%periodic = periodic
        spl%h_step = [(x_max(1) - x_min(1))/dble(n1-1), &
                      (x_max(2) - x_min(2))/dble(n2-1)]
        spl%x_min = x_min
        spl%num_quantities = n_quantities
        
        ! Allocate batch coefficients with cache-friendly layout
        allocate(spl%coeff(n_quantities, 0:N1_order, 0:N2_order, n1, n2), stat=istat)
        if (istat /= 0) then
            error stop "construct_batch_splines_2d: Allocation failed for coeff"
        end if
        
        allocate(temp_coeff(0:N1_order, 0:N2_order, n1, n2), stat=istat)
        if (istat /= 0) then
            error stop "construct_batch_splines_2d: Allocation failed for temp_coeff"
        end if
        
        ! Process each quantity
        do iq = 1, n_quantities
            ! Spline over x2 (second dimension)
            allocate(splcoe(0:N2_order, n2), stat=istat)
            if (istat /= 0) then
                error stop "construct_batch_splines_2d: Allocation failed for splcoe"
            end if
            
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
            allocate(splcoe(0:N1_order, n1), stat=istat)
            if (istat /= 0) then
                error stop "construct_batch_splines_2d: Allocation failed for splcoe"
            end if
            
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
        
        ! Validate output size
        if (size(y_batch) < spl%num_quantities) then
            error stop "evaluate_batch_splines_2d: Output array too small"
        end if
        
        do j=1,2
            if (spl%periodic(j)) then
                xj = modulo(x(j) - spl%x_min(j), &
                    spl%h_step(j)*(spl%num_points(j)-1)) + spl%x_min(j)
            else
                xj = x(j)
            end if
            x_norm(j) = (xj - spl%x_min(j))/spl%h_step(j)
            interval_index(j) = max(0, min(spl%num_points(j)-2, int(x_norm(j))))
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

    subroutine evaluate_batch_splines_2d_many(spl, x, y_batch)
        type(BatchSplineData2D), intent(in) :: spl
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: y_batch(:, :)  ! (n_quantities, n_points)

        integer :: ipt, iq, k1, k2, i1, i2, k_wrap
        integer :: nq, order1, order2
        integer :: num_points(2)
        real(dp) :: xj1, xj2, x_norm1, x_norm2, x_local1, x_local2
        real(dp) :: x_min(2), h_step(2), period(2)
        real(dp) :: t, w, v, yq

        if (size(x, 1) /= 2) then
            error stop "evaluate_batch_splines_2d_many: First dimension of x must be 2"
        end if
        if (size(y_batch, 1) < spl%num_quantities) then
            error stop "evaluate_batch_splines_2d_many: First dimension too small"
        end if
        if (size(y_batch, 2) /= size(x, 2)) then
            error stop "evaluate_batch_splines_2d_many: Second dimension must equal size(x,2)"
        end if

        nq = spl%num_quantities
        num_points = spl%num_points
        order1 = spl%order(1)
        order2 = spl%order(2)
        x_min = spl%x_min
        h_step = spl%h_step
        period(1) = h_step(1)*real(num_points(1) - 1, dp)
        period(2) = h_step(2)*real(num_points(2) - 1, dp)

        do ipt = 1, size(x, 2)
            include "spline2d_many_point_body.inc"
        end do
    end subroutine evaluate_batch_splines_2d_many

    subroutine evaluate_batch_splines_2d_many_resident(spl, x, y_batch)
        type(BatchSplineData2D), intent(in) :: spl
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(inout) :: y_batch(:, :)  ! (n_quantities, n_points)

        integer :: ipt, iq, k1, k2, i1, i2, k_wrap
        integer :: nq, order1, order2
        integer :: num_points(2)
        real(dp) :: xj1, xj2, x_norm1, x_norm2, x_local1, x_local2
        real(dp) :: x_min(2), h_step(2), period(2)
        real(dp) :: t, w, v, yq

        if (size(x, 1) /= 2) then
            error stop "evaluate_batch_splines_2d_many_resident: First dimension of x must be 2"
        end if
        if (size(y_batch, 1) < spl%num_quantities) then
            error stop "evaluate_batch_splines_2d_many_resident: First dimension too small"
        end if
        if (size(y_batch, 2) /= size(x, 2)) then
            error stop "evaluate_batch_splines_2d_many_resident: Second dimension must equal size(x,2)"
        end if

        nq = spl%num_quantities
        num_points = spl%num_points
        order1 = spl%order(1)
        order2 = spl%order(2)
        x_min = spl%x_min
        h_step = spl%h_step
        period(1) = h_step(1)*real(num_points(1) - 1, dp)
        period(2) = h_step(2)*real(num_points(2) - 1, dp)

        !$acc parallel loop present(spl%coeff, x, y_batch) &
        !$acc& private(ipt, iq, k1, k2, i1, i2, k_wrap) &
        !$acc& private(xj1, xj2, x_norm1, x_norm2, x_local1, x_local2) &
        !$acc& private(t, w, v, yq) gang vector vector_length(256)
        do ipt = 1, size(x, 2)
            include "spline2d_many_point_body.inc"
        end do
        !$acc end parallel loop
    end subroutine evaluate_batch_splines_2d_many_resident
    
    
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
        
        ! Validate output sizes
        if (size(y_batch) < spl%num_quantities) then
            error stop "evaluate_batch_splines_2d_der: y_batch array too small"
        end if
        if (size(dy_batch, 1) < 2 .or. size(dy_batch, 2) < spl%num_quantities) then
            error stop "evaluate_batch_splines_2d_der: dy_batch array too small"
        end if
        
        do j=1,2
            if (spl%periodic(j)) then
                xj = modulo(x(j) - spl%x_min(j), &
                    spl%h_step(j)*(spl%num_points(j)-1)) + spl%x_min(j)
            else
                xj = x(j)
            end if
            x_norm(j) = (xj - spl%x_min(j))/spl%h_step(j)
            interval_index(j) = max(0, min(spl%num_points(j)-2, int(x_norm(j))))
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
    
end module batch_interpolate_2d
