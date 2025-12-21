module batch_interpolate_1d
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use batch_interpolate_types, only: BatchSplineData1D
    use acc_spline_build_order5, only: spl_five_per_line, spl_five_reg_line
    use acc_tridiag_pcr, only: tridiag_pcr_solve_inplace, tridiag_pcr_solve_cyclic
    use spl_three_to_five_sub, only: spl_per, spl_reg
#ifdef _OPENACC
    use openacc, only: acc_is_present
#endif

    implicit none
    private
    
    ! Export batch spline construction/destruction routines
    public :: construct_batch_splines_1d
    public :: construct_batch_splines_1d_resident
    public :: construct_batch_splines_1d_resident_device
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

    subroutine construct_batch_splines_1d_resident_device(x_min, x_max, y_batch, &
                                                          order, periodic, spl, &
                                                          update_host, &
                                                          assume_y_present)
        real(dp), intent(in) :: x_min, x_max
        real(dp), intent(in) :: y_batch(:, :)  ! (n_points, n_quantities)
        integer, intent(in) :: order
        logical, intent(in) :: periodic
        type(BatchSplineData1D), intent(out) :: spl
        logical, intent(in), optional :: update_host
        logical, intent(in), optional :: assume_y_present

        integer :: n_points, n_quantities
        logical :: do_update
        logical :: do_assume_present

        do_update = .true.
        if (present(update_host)) do_update = update_host
        do_assume_present = .false.
        if (present(assume_y_present)) do_assume_present = assume_y_present

        n_points = size(y_batch, 1)
        n_quantities = size(y_batch, 2)

#ifdef _OPENACC
        call construct_batch_splines_1d_resident_device_impl(x_min, x_max, n_points, &
                                                            n_quantities, y_batch, &
                                                            order, periodic, spl, &
                                                            do_update, &
                                                            do_assume_present)
#else
        call construct_batch_splines_1d_resident(x_min, x_max, y_batch, order, &
                                                 periodic, spl)
#endif
    end subroutine construct_batch_splines_1d_resident_device

    subroutine construct_batch_splines_1d_resident_device_impl(x_min, x_max, n_points, &
                                                               n_quantities, y_batch, &
                                                               order, periodic, spl, &
                                                               do_update, &
                                                               do_assume_present)
        real(dp), intent(in) :: x_min, x_max
        integer, intent(in) :: n_points, n_quantities
        real(dp), intent(in) :: y_batch(n_points, n_quantities)
        integer, intent(in) :: order
        logical, intent(in) :: periodic
        type(BatchSplineData1D), intent(out) :: spl
        logical, intent(in) :: do_update
        logical, intent(in) :: do_assume_present

        integer :: istat
        integer :: ip, iq, k
        real(dp), allocatable :: work(:, :, :)
        real(dp), allocatable :: a_tri(:, :), b_tri(:, :), c_tri(:, :)
        real(dp), allocatable :: rhs(:, :), a_tmp(:, :), b_tmp(:, :)
        real(dp), allocatable :: c_tmp(:, :), rhs_tmp(:, :), x(:, :)
        real(dp), allocatable :: bdiag(:, :), y_sol(:, :), z_sol(:, :)
        real(dp), allocatable :: u_vec(:, :), factor(:)
        integer :: m
        real(dp) :: psi, inv_h
        real(dp) :: h_step

        if (n_points < 2) then
            error stop "construct_batch_splines_1d_resident_device:" // &
                " Need at least 2 points"
        end if
        if (n_quantities < 1) then
            error stop "construct_batch_splines_1d_resident_device:" // &
                " Need at least 1 quantity"
        end if
        if (order < 3 .or. order > 5) then
            error stop "construct_batch_splines_1d_resident_device:" // &
                " Order must be between 3 and 5"
        end if

        spl%order = order
        spl%num_points = n_points
        spl%periodic = periodic
        spl%x_min = x_min
        spl%h_step = (x_max - x_min) / dble(n_points - 1)
        spl%num_quantities = n_quantities

        allocate(spl%coeff(n_quantities, 0:order, n_points), stat=istat)
        if (istat /= 0) then
            error stop "construct_batch_splines_1d_resident_device:" // &
                " Allocation failed for coeff"
        end if

        if (.not. do_assume_present) then
            !$acc enter data copyin(y_batch)
        end if
        !$acc enter data create(spl%coeff)

        allocate(work(n_points, n_quantities, 0:order), stat=istat)
        if (istat /= 0) then
            error stop "construct_batch_splines_1d_resident_device:" // &
                " Allocation failed for work"
        end if

        inv_h = 1.0d0/spl%h_step
        h_step = spl%h_step
        psi = 3.0d0*inv_h*inv_h

        !$acc data present(y_batch, spl%coeff) create(work)
        !$acc parallel loop collapse(2) gang
        do iq = 1, n_quantities
            do ip = 1, n_points
                work(ip, iq, 0) = y_batch(ip, iq)
            end do
        end do

        if (order == 5) then
            !$acc parallel loop gang
            do iq = 1, n_quantities
                if (periodic) then
                    call spl_five_per_line(n_points, h_step, work(:, iq, 0), &
                                           work(:, iq, 1), work(:, iq, 2), &
                                           work(:, iq, 3), work(:, iq, 4), &
                                           work(:, iq, 5))
                else
                    call spl_five_reg_line(n_points, h_step, work(:, iq, 0), &
                                           work(:, iq, 1), work(:, iq, 2), &
                                           work(:, iq, 3), work(:, iq, 4), &
                                           work(:, iq, 5))
                end if
            end do
        else if (order == 3) then
            if (.not. periodic) then
                allocate(a_tri(n_points, n_quantities), b_tri(n_points, n_quantities), &
                         c_tri(n_points, n_quantities), rhs(n_points, n_quantities), &
                         a_tmp(n_points, n_quantities), b_tmp(n_points, n_quantities), &
                         c_tmp(n_points, n_quantities), rhs_tmp(n_points, n_quantities), &
                         x(n_points, n_quantities), stat=istat)
                if (istat /= 0) then
                    error stop "construct_batch_splines_1d_resident_device:" // &
                        " Allocation failed for cubic work arrays"
                end if

                !$acc data present(work) create(a_tri, b_tri, c_tri, rhs, a_tmp, b_tmp, &
                !$acc& c_tmp, rhs_tmp, x)
                !$acc parallel loop collapse(2) gang vector
                do iq = 1, n_quantities
                    do ip = 1, n_points
                        a_tri(ip, iq) = 1.0d0
                        b_tri(ip, iq) = 4.0d0
                        c_tri(ip, iq) = 1.0d0
                        rhs(ip, iq) = 0.0d0
                    end do
                end do

                !$acc parallel loop gang vector
                do iq = 1, n_quantities
                    a_tri(1, iq) = 0.0d0
                    c_tri(1, iq) = 0.0d0
                    b_tri(1, iq) = 1.0d0
                    a_tri(n_points, iq) = 0.0d0
                    c_tri(n_points, iq) = 0.0d0
                    b_tri(n_points, iq) = 1.0d0
                end do

                !$acc parallel loop collapse(2) gang vector
                do iq = 1, n_quantities
                    do ip = 2, n_points - 1
                        rhs(ip, iq) = (work(ip + 1, iq, 0) - 2.0d0*work(ip, iq, 0) + &
                                       work(ip - 1, iq, 0))*psi
                    end do
                end do

                call tridiag_pcr_solve_inplace(a_tri, b_tri, c_tri, rhs, a_tmp, b_tmp, &
                                               c_tmp, rhs_tmp, x)

                !$acc parallel loop collapse(2) gang vector
                do iq = 1, n_quantities
                    do ip = 1, n_points
                        work(ip, iq, 2) = x(ip, iq)
                    end do
                end do

                !$acc parallel loop collapse(2) gang vector
                do iq = 1, n_quantities
                    do ip = 1, n_points - 1
                        work(ip, iq, 1) = (work(ip + 1, iq, 0) - work(ip, iq, 0))*inv_h - &
                                          h_step*(work(ip + 1, iq, 2) + &
                                          2.0d0*work(ip, iq, 2))/3.0d0
                        work(ip, iq, 3) = (work(ip + 1, iq, 2) - work(ip, iq, 2))*inv_h / &
                                          3.0d0
                    end do
                end do
                !$acc parallel loop gang vector
                do iq = 1, n_quantities
                    work(n_points, iq, 1) = 0.0d0
                    work(n_points, iq, 3) = 0.0d0
                end do
                !$acc end data

                deallocate(a_tri, b_tri, c_tri, rhs, a_tmp, b_tmp, c_tmp, rhs_tmp, x)
            else
                if (n_points < 3) then
                    error stop "construct_batch_splines_1d_resident_device:" // &
                        " periodic cubic requires at least 3 points"
                end if

                m = n_points - 1
                allocate(a_tri(m, n_quantities), b_tri(m, n_quantities), c_tri(m, n_quantities), &
                         rhs(m, n_quantities), a_tmp(m, n_quantities), b_tmp(m, n_quantities), &
                         c_tmp(m, n_quantities), rhs_tmp(m, n_quantities), x(m, n_quantities), &
                         bdiag(m, n_quantities), y_sol(m, n_quantities), z_sol(m, n_quantities), &
                         u_vec(m, n_quantities), factor(n_quantities), stat=istat)
                if (istat /= 0) then
                    error stop "construct_batch_splines_1d_resident_device:" // &
                        " Allocation failed for periodic cubic work arrays"
                end if

                !$acc data present(work) create(a_tri, b_tri, c_tri, rhs, a_tmp, b_tmp, &
                !$acc& c_tmp, rhs_tmp, x, bdiag, y_sol, z_sol, u_vec, factor)
                !$acc parallel loop collapse(2) gang vector
                do iq = 1, n_quantities
                    do ip = 1, m
                        bdiag(ip, iq) = 4.0d0
                        rhs(ip, iq) = 0.0d0
                    end do
                end do

                !$acc parallel loop gang vector
                do iq = 1, n_quantities
                    rhs(1, iq) = (work(2, iq, 0) - work(1, iq, 0) - work(n_points, iq, 0) + &
                                  work(n_points - 1, iq, 0))*psi
                    rhs(m, iq) = (work(n_points, iq, 0) - 2.0d0*work(n_points - 1, iq, 0) + &
                                  work(n_points - 2, iq, 0))*psi
                end do
                !$acc parallel loop collapse(2) gang vector
                do iq = 1, n_quantities
                    do ip = 2, m - 1
                        rhs(ip, iq) = (work(ip + 1, iq, 0) - 2.0d0*work(ip, iq, 0) + &
                                       work(ip - 1, iq, 0))*psi
                    end do
                end do

                call tridiag_pcr_solve_cyclic(bdiag, rhs, a_tri, b_tri, c_tri, rhs, a_tmp, &
                                              b_tmp, c_tmp, rhs_tmp, y_sol, z_sol, u_vec, &
                                              factor, x)

                !$acc parallel loop collapse(2) gang vector
                do iq = 1, n_quantities
                    do ip = 1, m
                        work(ip, iq, 2) = x(ip, iq)
                    end do
                end do
                !$acc parallel loop gang vector
                do iq = 1, n_quantities
                    work(n_points, iq, 2) = work(1, iq, 2)
                end do

                !$acc parallel loop collapse(2) gang vector
                do iq = 1, n_quantities
                    do ip = 1, m - 1
                        work(ip, iq, 1) = (work(ip + 1, iq, 0) - work(ip, iq, 0))*inv_h - &
                                          h_step*(work(ip + 1, iq, 2) + &
                                          2.0d0*work(ip, iq, 2))/3.0d0
                        work(ip, iq, 3) = (work(ip + 1, iq, 2) - work(ip, iq, 2))*inv_h / &
                                          3.0d0
                    end do
                end do
                !$acc parallel loop gang vector
                do iq = 1, n_quantities
                    work(m, iq, 1) = (work(n_points, iq, 0) - work(n_points - 1, iq, 0))*inv_h - &
                                     h_step*(work(1, iq, 2) + 2.0d0*work(m, iq, 2))/3.0d0
                    work(m, iq, 3) = (work(1, iq, 2) - work(m, iq, 2))*inv_h / 3.0d0
                    work(n_points, iq, 1) = work(1, iq, 1)
                    work(n_points, iq, 3) = work(1, iq, 3)
                end do
                !$acc end data

                deallocate(a_tri, b_tri, c_tri, rhs, a_tmp, b_tmp, c_tmp, rhs_tmp, x)
                deallocate(bdiag, y_sol, z_sol, u_vec, factor)
            end if
        else
            error stop "construct_batch_splines_1d_resident_device:" // &
                " order=4 not yet supported on device"
        end if

        !$acc parallel loop collapse(2) gang
        do ip = 1, n_points
            do k = 0, order
                !$acc loop vector
                do iq = 1, n_quantities
                    spl%coeff(iq, k, ip) = work(ip, iq, k)
                end do
            end do
        end do
        !$acc end data

        deallocate(work)
        if (.not. do_assume_present) then
            !$acc exit data delete(y_batch)
        end if
        if (do_update) then
            !$acc update self(spl%coeff)
        end if
    end subroutine construct_batch_splines_1d_resident_device_impl


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
