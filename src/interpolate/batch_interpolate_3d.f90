module batch_interpolate_3d
    use batch_interpolate_types
    use spl_three_to_five_sub
    
    implicit none
    private
    
    integer, parameter :: dp = kind(1.0d0)
    
    ! Export batch spline construction/destruction routines
    public :: construct_batch_splines_3d
    public :: destroy_batch_splines_3d
    
    ! Export batch spline evaluation routines
    public :: evaluate_batch_splines_3d
    public :: evaluate_batch_splines_3d_der
    public :: evaluate_batch_splines_3d_der2
    
contains
    
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
        integer :: n1, n2, n3, n_quantities, istat
        integer :: N1_order, N2_order, N3_order
        
        n1 = size(y_batch, 1)
        n2 = size(y_batch, 2)
        n3 = size(y_batch, 3)
        n_quantities = size(y_batch, 4)
        N1_order = order(1)
        N2_order = order(2)
        N3_order = order(3)
        
        ! Validate input
        if (n1 < 2 .or. n2 < 2 .or. n3 < 2) then
            error stop &
                "construct_batch_splines_3d: Need at least 2 points in each dimension"
        end if
        if (n_quantities < 1) then
            error stop "construct_batch_splines_3d: Need at least 1 quantity"
        end if
        if (N1_order < 3 .or. N1_order > 5 .or. &
            N2_order < 3 .or. N2_order > 5 .or. &
            N3_order < 3 .or. N3_order > 5) then
            error stop "construct_batch_splines_3d: Orders must be between 3 and 5"
        end if
        
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
        allocate(spl%coeff(n_quantities, 0:N1_order, 0:N2_order, 0:N3_order, &
            n1, n2, n3), stat=istat)
        if (istat /= 0) then
            error stop "construct_batch_splines_3d: Allocation failed for coeff"
        end if
        
        allocate(temp_coeff(0:N1_order, 0:N2_order, 0:N3_order, n1, n2, n3), &
            stat=istat)
        if (istat /= 0) then
            error stop "construct_batch_splines_3d: Allocation failed for temp_coeff"
        end if
        
        ! Process each quantity
        do iq = 1, n_quantities
            ! Spline over x3 first
            allocate(splcoe(0:N3_order, n3), stat=istat)
            if (istat /= 0) then
                error stop "construct_batch_splines_3d: Allocation failed for splcoe"
            end if
            
            do i2 = 1, n2
                do i1 = 1, n1
                    splcoe(0,:) = y_batch(i1, i2, :, iq)
                    if (periodic(3)) then
                        call spl_per(N3_order, n3, spl%h_step(3), splcoe)
                    else
                        call spl_reg(N3_order, n3, spl%h_step(3), splcoe)
                    end if
                    temp_coeff(0, 0, :, i1, i2, :) = splcoe
                end do
            end do
            deallocate(splcoe)
            
            ! Spline over x2 second
            allocate(splcoe(0:N2_order, n2), stat=istat)
            if (istat /= 0) then
                error stop "construct_batch_splines_3d: Allocation failed for splcoe"
            end if
            
            do i3 = 1, n3
                do i1 = 1, n1
                    do k3 = 0, N3_order
                        splcoe(0,:) = temp_coeff(0, 0, k3, i1, :, i3)
                        if (periodic(2)) then
                            call spl_per(N2_order, n2, spl%h_step(2), splcoe)
                        else
                            call spl_reg(N2_order, n2, spl%h_step(2), splcoe)
                        end if
                        temp_coeff(0, :, k3, i1, :, i3) = splcoe
                    end do
                end do
            end do
            deallocate(splcoe)
            
            ! Spline over x1 last - store in STANDARD order for consistency
            allocate(splcoe(0:N1_order, n1), stat=istat)
            if (istat /= 0) then
                error stop "construct_batch_splines_3d: Allocation failed for splcoe"
            end if
            
            do i3 = 1, n3
                do i2 = 1, n2
                    do k3 = 0, N3_order
                        do k2 = 0, N2_order
                            splcoe(0,:) = temp_coeff(0, k2, k3, :, i2, i3)
                            if (periodic(1)) then
                                call spl_per(N1_order, n1, spl%h_step(1), splcoe)
                            else
                                call spl_reg(N1_order, n1, spl%h_step(1), splcoe)
                            end if
                            ! Store in STANDARD order for consistency and clarity
                            temp_coeff(:, k2, k3, :, i2, i3) = splcoe
                        end do
                    end do
                end do
            end do
            deallocate(splcoe)
            
            ! Copy to final coefficient array
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
        real(dp) :: coeff_23(spl%num_quantities, 0:spl%order(2), 0:spl%order(3))
        real(dp) :: coeff_3(spl%num_quantities, 0:spl%order(3))
        integer :: interval_index(3), k1, k2, k3, j, iq
        integer :: N1, N2, N3
        
        N1 = spl%order(1)
        N2 = spl%order(2)
        N3 = spl%order(3)
        
        ! Validate output size
        if (size(y_batch) < spl%num_quantities) then
            error stop "evaluate_batch_splines_3d: Output array too small"
        end if
        
        do j=1,3
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
        
        ! First reduction: interpolation over x1 using backward Horner
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
        
        ! Apply backward Horner method (standard polynomial evaluation)
        do k1 = N1-1, 0, -1
            !$omp simd
            do iq = 1, spl%num_quantities
                do k3 = 0, N3
                    do k2 = 0, N2
                        coeff_23(iq, k2, k3) = spl%coeff(iq, k1, k2, k3, &
                            interval_index(1)+1, interval_index(2)+1, &
                            interval_index(3)+1) + x_local(1)*coeff_23(iq, k2, k3)
                    end do
                end do
            end do
        end do
        
        ! Second reduction: interpolation over x2 using backward Horner
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
        
        ! Third reduction: interpolation over x3 using backward Horner
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
        real(dp) :: coeff_23_dx1(spl%num_quantities, 0:spl%order(2), &
            0:spl%order(3))
        real(dp) :: coeff_3(spl%num_quantities, 0:spl%order(3))
        real(dp) :: coeff_3_dx1(spl%num_quantities, 0:spl%order(3))
        real(dp) :: coeff_3_dx2(spl%num_quantities, 0:spl%order(3))
        
        integer :: interval_index(3), k1, k2, k3, j, iq
        integer :: N1, N2, N3
        
        N1 = spl%order(1)
        N2 = spl%order(2)
        N3 = spl%order(3)
        
        ! Validate output sizes
        if (size(y_batch) < spl%num_quantities) then
            error stop "evaluate_batch_splines_3d_der: y_batch array too small"
        end if
        if (size(dy_batch, 1) < 3 .or. size(dy_batch, 2) < spl%num_quantities) then
            error stop "evaluate_batch_splines_3d_der: dy_batch array too small"
        end if
        
        do j=1,3
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
        
        ! First reduction: interpolation over x1 for value using backward Horner
        ! Initialize with highest order
        !$omp simd
        do iq = 1, spl%num_quantities
            do k3 = 0, N3
                do k2 = 0, N2
                    coeff_23(iq, k2, k3) = spl%coeff(iq, N1, k2, k3, &
                        interval_index(1)+1, interval_index(2)+1, &
                        interval_index(3)+1)
                end do
            end do
        end do
        
        ! Apply backward Horner for value
        do k1 = N1-1, 0, -1
            !$omp simd
            do iq = 1, spl%num_quantities
                do k3 = 0, N3
                    do k2 = 0, N2
                        coeff_23(iq, k2, k3) = spl%coeff(iq, k1, k2, k3, &
                            interval_index(1)+1, interval_index(2)+1, &
                            interval_index(3)+1) + x_local(1)*coeff_23(iq, k2, k3)
                    end do
                end do
            end do
        end do
        
        ! First derivative over x1 using backward Horner for derivative
        !$omp simd
        do iq = 1, spl%num_quantities
            do k3 = 0, N3
                do k2 = 0, N2
                    coeff_23_dx1(iq, k2, k3) = N1 * spl%coeff(iq, N1, k2, k3, &
                        interval_index(1)+1, interval_index(2)+1, &
                        interval_index(3)+1)
                end do
            end do
        end do
        
        do k1 = N1-1, 1, -1
            !$omp simd
            do iq = 1, spl%num_quantities
                do k3 = 0, N3
                    do k2 = 0, N2
                        coeff_23_dx1(iq, k2, k3) = k1 * spl%coeff(iq, k1, k2, k3, &
                            interval_index(1)+1, interval_index(2)+1, &
                            interval_index(3)+1) + x_local(1)*coeff_23_dx1(iq, k2, k3)
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
                    coeff_3(iq, k3) = coeff_23(iq, k2, k3) + &
                        x_local(2)*coeff_3(iq, k3)
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
        real(dp) :: coeff_23_dx1(spl%num_quantities, 0:spl%order(2), &
            0:spl%order(3))
        real(dp) :: coeff_23_dx1x1(spl%num_quantities, 0:spl%order(2), &
            0:spl%order(3))
        
        real(dp) :: coeff_3(spl%num_quantities, 0:spl%order(3))
        real(dp) :: coeff_3_dx1(spl%num_quantities, 0:spl%order(3))
        real(dp) :: coeff_3_dx2(spl%num_quantities, 0:spl%order(3))
        real(dp) :: coeff_3_dx1x1(spl%num_quantities, 0:spl%order(3))
        real(dp) :: coeff_3_dx1x2(spl%num_quantities, 0:spl%order(3))
        real(dp) :: coeff_3_dx2x2(spl%num_quantities, 0:spl%order(3))
        
        integer :: interval_index(3), k1, k2, k3, j, iq
        integer :: N1, N2, N3
        real(dp) :: c  ! Temporary variable for coefficient reuse
        
        N1 = spl%order(1)
        N2 = spl%order(2)
        N3 = spl%order(3)
        
        ! Validate output sizes
        if (size(y_batch) < spl%num_quantities) then
            error stop "evaluate_batch_splines_3d_der2: y_batch array too small"
        end if
        if (size(dy_batch, 1) < 3 .or. size(dy_batch, 2) < spl%num_quantities) then
            error stop "evaluate_batch_splines_3d_der2: dy_batch array too small"
        end if
        if (size(d2y_batch, 1) < 6 .or. size(d2y_batch, 2) < spl%num_quantities) then
            error stop "evaluate_batch_splines_3d_der2: d2y_batch array too small"
        end if
        
        do j=1,3
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
        
        ! First reduction: x1 interpolation with fused derivative computation
        ! Initialize all polynomials with highest order coefficient (k1 = N1)
        !$omp simd
        do iq = 1, spl%num_quantities
            do k3 = 0, N3
                do k2 = 0, N2
                    c = spl%coeff(iq, N1, k2, k3, interval_index(1)+1, &
                        interval_index(2)+1, interval_index(3)+1)
                    coeff_23(iq, k2, k3) = c
                    coeff_23_dx1(iq, k2, k3) = N1 * c
                    coeff_23_dx1x1(iq, k2, k3) = N1*(N1-1) * c
                end do
            end do
        end do
        
        ! Fused backward Horner: compute all derivatives simultaneously
        ! Range k1 = N1-1 down to 2: update all three polynomials
        do k1 = N1-1, 2, -1
            !$omp simd
            do iq = 1, spl%num_quantities
                do k3 = 0, N3
                    do k2 = 0, N2
                        c = spl%coeff(iq, k1, k2, k3, interval_index(1)+1, &
                            interval_index(2)+1, interval_index(3)+1)
                        coeff_23(iq, k2, k3) = c + x_local(1)*coeff_23(iq, k2, k3)
                        coeff_23_dx1(iq, k2, k3) = k1*c + &
                            x_local(1)*coeff_23_dx1(iq, k2, k3)
                        coeff_23_dx1x1(iq, k2, k3) = k1*(k1-1)*c + &
                            x_local(1)*coeff_23_dx1x1(iq, k2, k3)
                    end do
                end do
            end do
        end do
        
        ! k1 = 1: update value and first derivative only
        if (N1 > 1) then
            k1 = 1
            !$omp simd
            do iq = 1, spl%num_quantities
                do k3 = 0, N3
                    do k2 = 0, N2
                        c = spl%coeff(iq, k1, k2, k3, interval_index(1)+1, &
                            interval_index(2)+1, interval_index(3)+1)
                        coeff_23(iq, k2, k3) = c + x_local(1)*coeff_23(iq, k2, k3)
                        coeff_23_dx1(iq, k2, k3) = k1*c + &
                            x_local(1)*coeff_23_dx1(iq, k2, k3)
                    end do
                end do
            end do
        end if
        
        ! k1 = 0: update value only
        k1 = 0
        !$omp simd
        do iq = 1, spl%num_quantities
            do k3 = 0, N3
                do k2 = 0, N2
                    c = spl%coeff(iq, k1, k2, k3, interval_index(1)+1, &
                        interval_index(2)+1, interval_index(3)+1)
                    coeff_23(iq, k2, k3) = c + x_local(1)*coeff_23(iq, k2, k3)
                end do
            end do
        end do
        
        ! Second reduction: x2 interpolation with fused derivative computation
        ! Initialize all polynomials with highest order coefficient (k2 = N2)
        !$omp simd
        do iq = 1, spl%num_quantities
            do k3 = 0, N3
                coeff_3(iq, k3) = coeff_23(iq, N2, k3)
                coeff_3_dx1(iq, k3) = coeff_23_dx1(iq, N2, k3)
                coeff_3_dx1x1(iq, k3) = coeff_23_dx1x1(iq, N2, k3)
                coeff_3_dx2(iq, k3) = N2 * coeff_23(iq, N2, k3)
                coeff_3_dx1x2(iq, k3) = N2 * coeff_23_dx1(iq, N2, k3)
                coeff_3_dx2x2(iq, k3) = N2*(N2-1) * coeff_23(iq, N2, k3)
            end do
        end do
        
        ! Fused backward Horner: update all polynomials simultaneously
        ! Range k2 = N2-1 down to 2: update all six polynomials
        do k2 = N2-1, 2, -1
            !$omp simd
            do iq = 1, spl%num_quantities
                do k3 = 0, N3
                    coeff_3(iq, k3) = coeff_23(iq, k2, k3) + &
                        x_local(2)*coeff_3(iq, k3)
                    coeff_3_dx1(iq, k3) = coeff_23_dx1(iq, k2, k3) + &
                        x_local(2)*coeff_3_dx1(iq, k3)
                    coeff_3_dx1x1(iq, k3) = coeff_23_dx1x1(iq, k2, k3) + &
                        x_local(2)*coeff_3_dx1x1(iq, k3)
                    coeff_3_dx2(iq, k3) = k2 * coeff_23(iq, k2, k3) + &
                        x_local(2)*coeff_3_dx2(iq, k3)
                    coeff_3_dx1x2(iq, k3) = k2 * coeff_23_dx1(iq, k2, k3) + &
                        x_local(2)*coeff_3_dx1x2(iq, k3)
                    coeff_3_dx2x2(iq, k3) = k2*(k2-1) * coeff_23(iq, k2, k3) + &
                        x_local(2)*coeff_3_dx2x2(iq, k3)
                end do
            end do
        end do
        
        ! k2 = 1: update value, dx1, dx1x1, dx2, dx1x2 only
        if (N2 > 1) then
            k2 = 1
            !$omp simd
            do iq = 1, spl%num_quantities
                do k3 = 0, N3
                    coeff_3(iq, k3) = coeff_23(iq, k2, k3) + &
                        x_local(2)*coeff_3(iq, k3)
                    coeff_3_dx1(iq, k3) = coeff_23_dx1(iq, k2, k3) + &
                        x_local(2)*coeff_3_dx1(iq, k3)
                    coeff_3_dx1x1(iq, k3) = coeff_23_dx1x1(iq, k2, k3) + &
                        x_local(2)*coeff_3_dx1x1(iq, k3)
                    coeff_3_dx2(iq, k3) = k2 * coeff_23(iq, k2, k3) + &
                        x_local(2)*coeff_3_dx2(iq, k3)
                    coeff_3_dx1x2(iq, k3) = k2 * coeff_23_dx1(iq, k2, k3) + &
                        x_local(2)*coeff_3_dx1x2(iq, k3)
                end do
            end do
        end if
        
        ! k2 = 0: update value, dx1, dx1x1 only
        k2 = 0
        !$omp simd
        do iq = 1, spl%num_quantities
            do k3 = 0, N3
                coeff_3(iq, k3) = coeff_23(iq, k2, k3) + &
                    x_local(2)*coeff_3(iq, k3)
                coeff_3_dx1(iq, k3) = coeff_23_dx1(iq, k2, k3) + &
                    x_local(2)*coeff_3_dx1(iq, k3)
                coeff_3_dx1x1(iq, k3) = coeff_23_dx1x1(iq, k2, k3) + &
                    x_local(2)*coeff_3_dx1x1(iq, k3)
            end do
        end do
        
        ! Third reduction: x3 interpolation with fused derivative computation
        ! Initialize all outputs with highest order coefficient (k3 = N3)
        !$omp simd
        do iq = 1, spl%num_quantities
            y_batch(iq) = coeff_3(iq, N3)
            dy_batch(1, iq) = coeff_3_dx1(iq, N3)
            dy_batch(2, iq) = coeff_3_dx2(iq, N3)
            dy_batch(3, iq) = N3 * coeff_3(iq, N3)
            d2y_batch(1, iq) = coeff_3_dx1x1(iq, N3)
            d2y_batch(2, iq) = coeff_3_dx1x2(iq, N3)
            d2y_batch(3, iq) = N3 * coeff_3_dx1(iq, N3)
            d2y_batch(4, iq) = coeff_3_dx2x2(iq, N3)
            d2y_batch(5, iq) = N3 * coeff_3_dx2(iq, N3)
            d2y_batch(6, iq) = N3*(N3-1) * coeff_3(iq, N3)
        end do
        
        ! Fused backward Horner: update all outputs simultaneously
        ! Range k3 = N3-1 down to 2: update all outputs
        do k3 = N3-1, 2, -1
            !$omp simd
            do iq = 1, spl%num_quantities
                y_batch(iq) = coeff_3(iq, k3) + x_local(3)*y_batch(iq)
                dy_batch(1, iq) = coeff_3_dx1(iq, k3) + x_local(3)*dy_batch(1, iq)
                dy_batch(2, iq) = coeff_3_dx2(iq, k3) + x_local(3)*dy_batch(2, iq)
                dy_batch(3, iq) = k3 * coeff_3(iq, k3) + x_local(3)*dy_batch(3, iq)
                d2y_batch(1, iq) = coeff_3_dx1x1(iq, k3) + &
                    x_local(3)*d2y_batch(1, iq)
                d2y_batch(2, iq) = coeff_3_dx1x2(iq, k3) + &
                    x_local(3)*d2y_batch(2, iq)
                d2y_batch(3, iq) = k3 * coeff_3_dx1(iq, k3) + &
                    x_local(3)*d2y_batch(3, iq)
                d2y_batch(4, iq) = coeff_3_dx2x2(iq, k3) + &
                    x_local(3)*d2y_batch(4, iq)
                d2y_batch(5, iq) = k3 * coeff_3_dx2(iq, k3) + &
                    x_local(3)*d2y_batch(5, iq)
                d2y_batch(6, iq) = k3*(k3-1) * coeff_3(iq, k3) + &
                    x_local(3)*d2y_batch(6, iq)
            end do
        end do
        
        ! k3 = 1: update all except d2y_batch(6)
        if (N3 > 1) then
            k3 = 1
            !$omp simd
            do iq = 1, spl%num_quantities
                y_batch(iq) = coeff_3(iq, k3) + x_local(3)*y_batch(iq)
                dy_batch(1, iq) = coeff_3_dx1(iq, k3) + x_local(3)*dy_batch(1, iq)
                dy_batch(2, iq) = coeff_3_dx2(iq, k3) + x_local(3)*dy_batch(2, iq)
                dy_batch(3, iq) = k3 * coeff_3(iq, k3) + x_local(3)*dy_batch(3, iq)
                d2y_batch(1, iq) = coeff_3_dx1x1(iq, k3) + &
                    x_local(3)*d2y_batch(1, iq)
                d2y_batch(2, iq) = coeff_3_dx1x2(iq, k3) + &
                    x_local(3)*d2y_batch(2, iq)
                d2y_batch(3, iq) = k3 * coeff_3_dx1(iq, k3) + &
                    x_local(3)*d2y_batch(3, iq)
                d2y_batch(4, iq) = coeff_3_dx2x2(iq, k3) + &
                    x_local(3)*d2y_batch(4, iq)
                d2y_batch(5, iq) = k3 * coeff_3_dx2(iq, k3) + &
                    x_local(3)*d2y_batch(5, iq)
            end do
        end if
        
        ! k3 = 0: update value and first-order derivatives only
        k3 = 0
        !$omp simd
        do iq = 1, spl%num_quantities
            y_batch(iq) = coeff_3(iq, k3) + x_local(3)*y_batch(iq)
            dy_batch(1, iq) = coeff_3_dx1(iq, k3) + x_local(3)*dy_batch(1, iq)
            dy_batch(2, iq) = coeff_3_dx2(iq, k3) + x_local(3)*dy_batch(2, iq)
            d2y_batch(1, iq) = coeff_3_dx1x1(iq, k3) + &
                x_local(3)*d2y_batch(1, iq)
            d2y_batch(2, iq) = coeff_3_dx1x2(iq, k3) + &
                x_local(3)*d2y_batch(2, iq)
            d2y_batch(4, iq) = coeff_3_dx2x2(iq, k3) + &
                x_local(3)*d2y_batch(4, iq)
        end do
        
    end subroutine evaluate_batch_splines_3d_der2
    
end module batch_interpolate_3d