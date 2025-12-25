module batch_interpolate_2d
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use batch_interpolate_types, only: BatchSplineData2D
    use spline_build_lines, only: spl_build_line_inplace, &
                                  spl_three_per_line, spl_three_reg_line, &
                                  spl_four_per_line, spl_four_reg_line, &
                                  spl_five_per_line, spl_five_reg_line
    use spl_three_to_five_sub, only: spl_per, spl_reg
#ifdef _OPENACC
    use openacc, only: acc_is_present
#endif

    implicit none
    private
    
    ! Export batch spline construction/destruction routines
    public :: construct_batch_splines_2d
    public :: construct_batch_splines_2d_lines
    public :: construct_batch_splines_2d_resident
    public :: construct_batch_splines_2d_resident_device
    public :: destroy_batch_splines_2d

#ifdef LIBNEO_ENABLE_SPLINE_ORACLE
    public :: construct_batch_splines_2d_legacy
#endif
    
    ! Export batch spline evaluation routines
    public :: evaluate_batch_splines_2d
    public :: evaluate_batch_splines_2d_der
    public :: evaluate_batch_splines_2d_many
    public :: evaluate_batch_splines_2d_many_resident
    public :: evaluate_batch_splines_2d_many_der
    
contains
    
    subroutine construct_batch_splines_2d(x_min, x_max, y_batch, order, periodic, spl)
        real(dp), intent(in) :: x_min(:), x_max(:)
        real(dp), intent(in) :: y_batch(:,:,:)  ! (n1, n2, n_quantities)
        integer, intent(in) :: order(:)
        logical, intent(in) :: periodic(:)
        type(BatchSplineData2D), intent(inout) :: spl

#ifdef _OPENACC
        ! Use GPU construction by default (1.5x faster than CPU)
        call construct_batch_splines_2d_resident_device(x_min, x_max, y_batch, &
                                                        order, periodic, spl, &
                                                        update_host=.true., &
                                                        assume_y_present=.false.)
#else
        call construct_batch_splines_2d_legacy(x_min, x_max, y_batch, order, periodic, spl)
#endif
    end subroutine construct_batch_splines_2d

    subroutine construct_batch_splines_2d_legacy(x_min, x_max, y_batch, order, periodic, spl)
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
            error stop "construct_batch_splines_2d: Need at least 2 points in " // &
                       "each dimension"
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
    end subroutine construct_batch_splines_2d_legacy

    subroutine construct_batch_splines_2d_lines(x_min, x_max, y_batch, order, periodic, spl)
        real(dp), intent(in) :: x_min(2), x_max(2)
        real(dp), intent(in) :: y_batch(:, :, :)  ! (n1, n2, n_quantities)
        integer, intent(in) :: order(2)
        logical, intent(in) :: periodic(2)
        type(BatchSplineData2D), intent(out) :: spl

        integer :: n1, n2, n_quantities
        integer :: N1_order, N2_order
        integer :: istat
        integer :: i1, i2, iq, k1, k2
        integer :: line, line2
        real(dp), allocatable :: work2(:, :, :)
        real(dp), allocatable :: work1(:, :, :)

        n1 = size(y_batch, 1)
        n2 = size(y_batch, 2)
        n_quantities = size(y_batch, 3)
        N1_order = order(1)
        N2_order = order(2)

        if (n1 < 2 .or. n2 < 2) then
            error stop "construct_batch_splines_2d_lines: Need at least 2 points"
        end if
        if (n_quantities < 1) then
            error stop "construct_batch_splines_2d_lines: Need at least 1 quantity"
        end if
        if (any(order < 3) .or. any(order > 5)) then
            error stop "construct_batch_splines_2d_lines: Order must be between 3 and 5"
        end if
        if (N1_order == 3) then
            if (periodic(1) .and. n1 < 3) then
                error stop "construct_batch_splines_2d_lines: Need at least 3 points " // &
                           "for periodic order(1)=3"
            end if
        else if (N1_order == 4) then
            if (n1 < 5) then
                error stop "construct_batch_splines_2d_lines: Need at least 5 points for order(1)=4"
            end if
        else
            if (n1 < 6) then
                error stop "construct_batch_splines_2d_lines: Need at least 6 points for order(1)=5"
            end if
        end if
        if (N2_order == 3) then
            if (periodic(2) .and. n2 < 3) then
                error stop "construct_batch_splines_2d_lines: Need at least 3 points " // &
                           "for periodic order(2)=3"
            end if
        else if (N2_order == 4) then
            if (n2 < 5) then
                error stop "construct_batch_splines_2d_lines: Need at least 5 points for order(2)=4"
            end if
        else
            if (n2 < 6) then
                error stop "construct_batch_splines_2d_lines: Need at least 6 points for order(2)=5"
            end if
        end if

        spl%order = order
        spl%num_points = [n1, n2]
        spl%periodic = periodic
        spl%h_step = [(x_max(1) - x_min(1))/dble(n1 - 1), &
                      (x_max(2) - x_min(2))/dble(n2 - 1)]
        spl%x_min = x_min
        spl%num_quantities = n_quantities

        allocate(spl%coeff(n_quantities, 0:N1_order, 0:N2_order, n1, n2), stat=istat)
        if (istat /= 0) then
            error stop "construct_batch_splines_2d_lines: Allocation failed for coeff"
        end if

        allocate(work2(n2, n1*n_quantities, 0:N2_order), stat=istat)
        if (istat /= 0) then
            error stop "construct_batch_splines_2d_lines: Allocation failed for work2"
        end if
        allocate(work1(n1, n2*n_quantities, 0:N1_order), stat=istat)
        if (istat /= 0) then
            error stop "construct_batch_splines_2d_lines: Allocation failed for work1"
        end if

        do iq = 1, n_quantities
            do i1 = 1, n1
                line = (iq - 1)*n1 + i1
                work2(:, line, 0) = y_batch(i1, :, iq)
            end do
        end do

#ifdef _OPENMP
        !$omp parallel do default(none) shared(work2, n1, n2, n_quantities, N2_order, periodic, &
        !$omp& spl) private(line)
#endif
        do line = 1, n1*n_quantities
            if (N2_order == 3) then
                if (periodic(2)) then
                    call spl_three_per_line(n2, spl%h_step(2), work2(:, line, 0), &
                                            work2(:, line, 1), work2(:, line, 2), &
                                            work2(:, line, 3))
                else
                    call spl_three_reg_line(n2, spl%h_step(2), work2(:, line, 0), &
                                            work2(:, line, 1), work2(:, line, 2), &
                                            work2(:, line, 3))
                end if
            else if (N2_order == 4) then
                if (periodic(2)) then
                    call spl_four_per_line(n2, spl%h_step(2), work2(:, line, 0), &
                                           work2(:, line, 1), work2(:, line, 2), &
                                           work2(:, line, 3), work2(:, line, 4))
                else
                    call spl_four_reg_line(n2, spl%h_step(2), work2(:, line, 0), &
                                           work2(:, line, 1), work2(:, line, 2), &
                                           work2(:, line, 3), work2(:, line, 4))
                end if
            else
                if (periodic(2)) then
                    call spl_five_per_line(n2, spl%h_step(2), work2(:, line, 0), &
                                           work2(:, line, 1), work2(:, line, 2), &
                                           work2(:, line, 3), work2(:, line, 4), &
                                           work2(:, line, 5))
                else
                    call spl_five_reg_line(n2, spl%h_step(2), work2(:, line, 0), &
                                           work2(:, line, 1), work2(:, line, 2), &
                                           work2(:, line, 3), work2(:, line, 4), &
                                           work2(:, line, 5))
                end if
            end if
        end do
#ifdef _OPENMP
        !$omp end parallel do
#endif

        do k2 = 0, N2_order
            do iq = 1, n_quantities
                do i2 = 1, n2
                    line = (iq - 1)*n2 + i2
                    do i1 = 1, n1
                        line2 = (iq - 1)*n1 + i1
                        work1(i1, line, 0) = work2(i2, line2, k2)
                    end do
                end do
            end do

#ifdef _OPENMP
            !$omp parallel do default(none) shared(work1, n1, n2, n_quantities, N1_order, periodic, &
            !$omp& spl) private(line)
#endif
            do line = 1, n2*n_quantities
                if (N1_order == 3) then
                    if (periodic(1)) then
                        call spl_three_per_line(n1, spl%h_step(1), work1(:, line, 0), &
                                                work1(:, line, 1), work1(:, line, 2), &
                                                work1(:, line, 3))
                    else
                        call spl_three_reg_line(n1, spl%h_step(1), work1(:, line, 0), &
                                                work1(:, line, 1), work1(:, line, 2), &
                                                work1(:, line, 3))
                    end if
                else if (N1_order == 4) then
                    if (periodic(1)) then
                        call spl_four_per_line(n1, spl%h_step(1), work1(:, line, 0), &
                                               work1(:, line, 1), work1(:, line, 2), &
                                               work1(:, line, 3), work1(:, line, 4))
                    else
                        call spl_four_reg_line(n1, spl%h_step(1), work1(:, line, 0), &
                                               work1(:, line, 1), work1(:, line, 2), &
                                               work1(:, line, 3), work1(:, line, 4))
                    end if
                else
                    if (periodic(1)) then
                        call spl_five_per_line(n1, spl%h_step(1), work1(:, line, 0), &
                                               work1(:, line, 1), work1(:, line, 2), &
                                               work1(:, line, 3), work1(:, line, 4), &
                                               work1(:, line, 5))
                    else
                        call spl_five_reg_line(n1, spl%h_step(1), work1(:, line, 0), &
                                               work1(:, line, 1), work1(:, line, 2), &
                                               work1(:, line, 3), work1(:, line, 4), &
                                               work1(:, line, 5))
                    end if
                end if
            end do
#ifdef _OPENMP
            !$omp end parallel do
#endif

            do iq = 1, n_quantities
                do i2 = 1, n2
                    line = (iq - 1)*n2 + i2
                    do i1 = 1, n1
                        do k1 = 0, N1_order
                            spl%coeff(iq, k1, k2, i1, i2) = work1(i1, line, k1)
                        end do
                    end do
                end do
            end do
        end do

        deallocate(work1)
        deallocate(work2)
    end subroutine construct_batch_splines_2d_lines

    subroutine construct_batch_splines_2d_resident(x_min, x_max, y_batch, order, &
                                                   periodic, spl)
        real(dp), intent(in) :: x_min(:), x_max(:)
        real(dp), intent(in) :: y_batch(:, :, :)  ! (n1, n2, n_quantities)
        integer, intent(in) :: order(:)
        logical, intent(in) :: periodic(:)
        type(BatchSplineData2D), intent(out) :: spl

        call construct_batch_splines_2d(x_min, x_max, y_batch, order, periodic, spl)
    end subroutine construct_batch_splines_2d_resident

    subroutine construct_batch_splines_2d_resident_device(x_min, x_max, y_batch, &
                                                          order, periodic, spl, &
                                                          update_host, &
                                                          assume_y_present)
        real(dp), intent(in) :: x_min(:), x_max(:)
        real(dp), intent(in) :: y_batch(:, :, :)  ! (n1, n2, n_quantities)
        integer, intent(in) :: order(:)
        logical, intent(in) :: periodic(:)
        type(BatchSplineData2D), intent(out) :: spl
        logical, intent(in), optional :: update_host
        logical, intent(in), optional :: assume_y_present

        integer :: n1, n2, n_quantities
        logical :: do_update, do_assume_present

        do_update = .true.
        if (present(update_host)) do_update = update_host
        do_assume_present = .false.
        if (present(assume_y_present)) do_assume_present = assume_y_present

        n1 = size(y_batch, 1)
        n2 = size(y_batch, 2)
        n_quantities = size(y_batch, 3)

#ifdef _OPENACC
        call construct_batch_splines_2d_resident_device_impl(x_min, x_max, n1, n2, &
                                                            n_quantities, y_batch, &
                                                            order, periodic, spl, &
                                                            do_update, do_assume_present)
#else
        call construct_batch_splines_2d_resident(x_min, x_max, y_batch, order, &
                                                 periodic, spl)
#endif
    end subroutine construct_batch_splines_2d_resident_device

    subroutine construct_batch_splines_2d_resident_device_impl(x_min, x_max, n1, n2, &
                                                               n_quantities, y_batch, &
                                                               order, periodic, spl, &
                                                               do_update, &
                                                               do_assume_present)
        real(dp), intent(in) :: x_min(2), x_max(2)
        integer, intent(in) :: n1, n2, n_quantities
        real(dp), intent(in) :: y_batch(n1, n2, n_quantities)
        integer, intent(in) :: order(2)
        logical, intent(in) :: periodic(2)
        type(BatchSplineData2D), intent(inout) :: spl
        logical, intent(in) :: do_update
        logical, intent(in) :: do_assume_present

        integer :: istat
        integer :: N1_order, N2_order
        integer :: order1, order2
        integer :: i1, i2, iq, k1, k2
        integer :: line, line2
        ! GCC OpenACC bug workaround: keep work arrays persistent across calls
        real(dp), allocatable, save :: work2(:, :, :)
        real(dp), allocatable, save :: work1(:, :, :)
        real(dp) :: h1, h2
        logical :: periodic1, periodic2
        logical :: work_needs_alloc

	        N1_order = order(1)
	        N2_order = order(2)
	        order1 = N1_order
	        order2 = N2_order
	        periodic1 = periodic(1)
	        periodic2 = periodic(2)

        if (n1 < 2 .or. n2 < 2) then
            error stop "construct_batch_splines_2d_resident_device:" // &
                " Need at least 2 points"
        end if
        if (n_quantities < 1) then
            error stop "construct_batch_splines_2d_resident_device:" // &
                " Need at least 1 quantity"
        end if
        if (any(order < 3) .or. any(order > 5)) then
            error stop "construct_batch_splines_2d_resident_device:" // &
                " Order must be between 3 and 5"
        end if
        if (N1_order == 3) then
            if (periodic(1) .and. n1 < 3) then
                error stop "construct_batch_splines_2d_resident_device:" // &
                    " Need at least 3 points for periodic order(1)=3"
            end if
        else if (N1_order == 4) then
            if (n1 < 5) then
                error stop "construct_batch_splines_2d_resident_device:" // &
                    " Need at least 5 points for order(1)=4"
            end if
        else
            if (n1 < 6) then
                error stop "construct_batch_splines_2d_resident_device:" // &
                    " Need at least 6 points for order(1)=5"
            end if
        end if
        if (N2_order == 3) then
            if (periodic(2) .and. n2 < 3) then
                error stop "construct_batch_splines_2d_resident_device:" // &
                    " Need at least 3 points for periodic order(2)=3"
            end if
        else if (N2_order == 4) then
            if (n2 < 5) then
                error stop "construct_batch_splines_2d_resident_device:" // &
                    " Need at least 5 points for order(2)=4"
            end if
        else
            if (n2 < 6) then
                error stop "construct_batch_splines_2d_resident_device:" // &
                    " Need at least 6 points for order(2)=5"
            end if
        end if

        spl%order = order
        spl%num_points = [n1, n2]
        spl%periodic = periodic
        spl%h_step = [(x_max(1) - x_min(1))/dble(n1 - 1), &
                      (x_max(2) - x_min(2))/dble(n2 - 1)]
        spl%x_min = x_min
        spl%num_quantities = n_quantities

        ! GCC OpenACC bug workaround: reuse existing allocation if possible
        if (.not. allocated(spl%coeff)) then
            allocate(spl%coeff(n_quantities, 0:N1_order, 0:N2_order, n1, n2), stat=istat)
            if (istat /= 0) then
                error stop "construct_batch_splines_2d_resident_device:" // &
                    " Allocation failed for coeff"
            end if
            !$acc enter data create(spl%coeff)
        else if (size(spl%coeff, 1) /= n_quantities .or. &
                 size(spl%coeff, 2) /= N1_order + 1 .or. &
                 size(spl%coeff, 3) /= N2_order + 1 .or. &
                 size(spl%coeff, 4) /= n1 .or. &
                 size(spl%coeff, 5) /= n2) then
#ifdef _OPENACC
            if (acc_is_present(spl%coeff)) then
                !$acc exit data delete(spl%coeff)
            end if
#endif
            deallocate(spl%coeff)
            allocate(spl%coeff(n_quantities, 0:N1_order, 0:N2_order, n1, n2), stat=istat)
            if (istat /= 0) then
                error stop "construct_batch_splines_2d_resident_device:" // &
                    " Allocation failed for coeff"
            end if
            !$acc enter data create(spl%coeff)
        end if

        ! GCC OpenACC bug workaround: check if work arrays need (re)allocation
        ! Must check BOTH work arrays since their sizes depend on different orders
        work_needs_alloc = .not. allocated(work2)
        if (allocated(work2)) then
            if (size(work2, 1) /= n2 .or. &
                size(work2, 2) /= n1*n_quantities .or. &
                size(work2, 3) /= N2_order + 1 .or. &
                size(work1, 1) /= n1 .or. &
                size(work1, 2) /= n2*n_quantities .or. &
                size(work1, 3) /= N1_order + 1) then
#ifdef _OPENACC
                if (acc_is_present(work2)) then
                    !$acc exit data delete(work2, work1)
                end if
#endif
                deallocate(work2, work1)
                work_needs_alloc = .true.
            end if
        end if

        if (work_needs_alloc) then
            allocate(work2(n2, n1*n_quantities, 0:N2_order), stat=istat)
            if (istat /= 0) then
                error stop "construct_batch_splines_2d_resident_device:" // &
                    " Allocation failed for work2"
            end if
            allocate(work1(n1, n2*n_quantities, 0:N1_order), stat=istat)
            if (istat /= 0) then
                error stop "construct_batch_splines_2d_resident_device:" // &
                    " Allocation failed for work1"
            end if
        end if

        h1 = spl%h_step(1)
        h2 = spl%h_step(2)

#ifdef _OPENACC
        block
            if (work_needs_alloc) then
                !$acc enter data create(work2, work1)
            end if

            ! GCC OpenACC bug workaround: first parallel loop MUST write to spl%coeff
            !$acc parallel loop collapse(3) gang present(work2, spl%coeff) &
            !$acc&    present_or_copyin(y_batch(1:n1, 1:n2, 1:n_quantities))
            do iq = 1, n_quantities
                do i1 = 1, n1
                    do i2 = 1, n2
                        line = (iq - 1)*n1 + i1
                        work2(i2, line, 0) = y_batch(i1, i2, iq)
                        ! GCC bug workaround: write to coeff in first loop
                        spl%coeff(iq, 0, 0, i1, i2) = 0.0d0
                    end do
                end do
            end do

            !$acc parallel loop gang present(work2)
            do line = 1, n1*n_quantities
                call spl_build_line_inplace(order2, periodic2, n2, h2, work2, line)
            end do

            do k2 = 0, N2_order
                !$acc parallel loop collapse(2) gang present(work1, work2)
                do iq = 1, n_quantities
                    do i2 = 1, n2
                        line = (iq - 1)*n2 + i2
                        !$acc loop vector
                        do i1 = 1, n1
                            line2 = (iq - 1)*n1 + i1
                            work1(i1, line, 0) = work2(i2, line2, k2)
                        end do
                    end do
                end do

                !$acc parallel loop gang present(work1)
                do line = 1, n2*n_quantities
                    call spl_build_line_inplace(order1, periodic1, n1, h1, work1, &
                                                line)
                end do

                !$acc parallel loop collapse(3) gang present(work1, spl%coeff)
                do iq = 1, n_quantities
                    do i2 = 1, n2
                        do i1 = 1, n1
                            line = (iq - 1)*n2 + i2
                            !$acc loop vector
                            do k1 = 0, N1_order
                                spl%coeff(iq, k1, k2, i1, i2) = work1(i1, line, k1)
                            end do
                        end do
                    end do
                end do
            end do

            ! GCC bug workaround: keep work arrays mapped for next call
        end block
#endif

        ! GCC bug workaround: keep work arrays allocated for next call
        if (do_update) then
            !$acc update self(spl%coeff(1:n_quantities, 0:N1_order, 0:N2_order, 1:n1, 1:n2))
        end if
    end subroutine construct_batch_splines_2d_resident_device_impl
    
    
    subroutine destroy_batch_splines_2d(spl)
        type(BatchSplineData2D), intent(inout) :: spl

#ifdef _OPENACC
        if (allocated(spl%coeff)) then
            if (acc_is_present(spl%coeff)) then
                !$acc exit data delete(spl%coeff)
            end if
        end if
#endif
        if (allocated(spl%coeff)) deallocate(spl%coeff)
    end subroutine destroy_batch_splines_2d
    
    
    subroutine evaluate_batch_splines_2d(spl, x, y_batch)
        type(BatchSplineData2D), intent(in) :: spl
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: y_batch(:)  ! (n_quantities)

        real(dp) :: x_arr(2, 1)
        real(dp) :: y_arr(spl%num_quantities, 1)

        x_arr(:, 1) = x
        call evaluate_batch_splines_2d_many(spl, x_arr, y_arr)
        y_batch(1:spl%num_quantities) = y_arr(:, 1)
    end subroutine evaluate_batch_splines_2d

    subroutine evaluate_batch_splines_2d_many(spl, x, y_batch)
        type(BatchSplineData2D), intent(in) :: spl
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: y_batch(:, :)  ! (n_quantities, n_points)

        integer :: ipt, iq, k1, k2, i1, i2, k_wrap
        integer :: nq, order1, order2
        integer :: num_points(2)
        logical :: periodic(2)
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
            error stop "evaluate_batch_splines_2d_many: y_batch second dim mismatch"
        end if

#ifdef _OPENACC
        if (acc_is_present(spl%coeff)) then
            !$acc data copyin(x) copy(y_batch)
            call evaluate_batch_splines_2d_many_resident(spl, x, y_batch)
            !$acc end data
            return
        end if
#endif

        nq = spl%num_quantities
        num_points = spl%num_points
        order1 = spl%order(1)
        order2 = spl%order(2)
        periodic = spl%periodic
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
        real(dp), intent(inout) :: y_batch(:, :)

        integer :: ipt, iq, k1, k2, i1, i2, k_wrap
        integer :: nq, order1, order2
        integer :: num_points(2)
        logical :: periodic(2)
        real(dp) :: xj1, xj2, x_norm1, x_norm2, x_local1, x_local2
        real(dp) :: x_min(2), h_step(2), period(2)
        real(dp) :: t, w, v, yq

        if (size(x, 1) /= 2) then
            error stop "evaluate_batch_splines_2d_many_resident: x first dim must be 2"
        end if
        if (size(y_batch, 1) < spl%num_quantities) then
            error stop "evaluate_batch_splines_2d_many_resident: y_batch dim1 too small"
        end if
        if (size(y_batch, 2) /= size(x, 2)) then
            error stop "evaluate_batch_splines_2d_many_resident: y_batch dim2 mismatch"
        end if

        nq = spl%num_quantities
        num_points = spl%num_points
        order1 = spl%order(1)
        order2 = spl%order(2)
        periodic = spl%periodic
        x_min = spl%x_min
        h_step = spl%h_step
        period(1) = h_step(1)*real(num_points(1) - 1, dp)
        period(2) = h_step(2)*real(num_points(2) - 1, dp)

        !$acc parallel loop present(spl%coeff, x, y_batch) &
        !$acc& private(ipt, iq, k1, k2, i1, i2, k_wrap) &
        !$acc& private(xj1, xj2, x_norm1, x_norm2, x_local1, x_local2, t, w, v, yq)
        do ipt = 1, size(x, 2)
            include "spline2d_many_point_body.inc"
        end do
        !$acc end parallel loop
    end subroutine evaluate_batch_splines_2d_many_resident

    subroutine evaluate_batch_splines_2d_many_der(spl, x, y_batch, dy_batch)
        type(BatchSplineData2D), intent(in) :: spl
        real(dp), intent(in) :: x(:,:)           ! (2, npts)
        real(dp), intent(out) :: y_batch(:,:)    ! (nq, npts)
        real(dp), intent(out) :: dy_batch(:,:,:) ! (2, nq, npts)

        integer :: ipt, iq, k1, k2, j, idx(2), npts, nq, N1, N2
        real(dp) :: xj, x_norm(2), x_local(2), period(2)
        real(dp), allocatable :: coeff_2(:,:), coeff_2_dx1(:,:)

        npts = size(x, 2)
        nq = spl%num_quantities
        N1 = spl%order(1)
        N2 = spl%order(2)
        period = spl%h_step * real(spl%num_points - 1, dp)

        allocate(coeff_2(nq, 0:N2), coeff_2_dx1(nq, 0:N2))

        do ipt = 1, npts
            do j = 1, 2
                if (spl%periodic(j)) then
                    xj = modulo(x(j, ipt) - spl%x_min(j), period(j)) + spl%x_min(j)
                else
                    xj = x(j, ipt)
                end if
                x_norm(j) = (xj - spl%x_min(j)) / spl%h_step(j)
                idx(j) = max(0, min(spl%num_points(j) - 2, int(x_norm(j)))) + 1
                x_local(j) = (x_norm(j) - real(idx(j) - 1, dp)) * spl%h_step(j)
            end do

            do iq = 1, nq
                do k2 = 0, N2
                    coeff_2(iq, k2) = spl%coeff(iq, N1, k2, idx(1), idx(2))
                    coeff_2_dx1(iq, k2) = N1 * spl%coeff(iq, N1, k2, idx(1), idx(2))
                end do
            end do

            do k1 = N1 - 1, 0, -1
                do iq = 1, nq
                    do k2 = 0, N2
                        coeff_2(iq, k2) = spl%coeff(iq, k1, k2, idx(1), idx(2)) + &
                                          x_local(1) * coeff_2(iq, k2)
                    end do
                end do
            end do

            do k1 = N1 - 1, 1, -1
                do iq = 1, nq
                    do k2 = 0, N2
                        coeff_2_dx1(iq, k2) = k1 * spl%coeff(iq, k1, k2, idx(1), idx(2)) + &
                                              x_local(1) * coeff_2_dx1(iq, k2)
                    end do
                end do
            end do

            do iq = 1, nq
                y_batch(iq, ipt) = coeff_2(iq, N2)
                dy_batch(1, iq, ipt) = coeff_2_dx1(iq, N2)
                dy_batch(2, iq, ipt) = N2 * coeff_2(iq, N2)
            end do

            do k2 = N2 - 1, 0, -1
                do iq = 1, nq
                    y_batch(iq, ipt) = coeff_2(iq, k2) + x_local(2) * y_batch(iq, ipt)
                    dy_batch(1, iq, ipt) = coeff_2_dx1(iq, k2) + &
                                           x_local(2) * dy_batch(1, iq, ipt)
                end do
            end do

            do k2 = N2 - 1, 1, -1
                do iq = 1, nq
                    dy_batch(2, iq, ipt) = k2 * coeff_2(iq, k2) + &
                                           x_local(2) * dy_batch(2, iq, ipt)
                end do
            end do
        end do

        deallocate(coeff_2, coeff_2_dx1)
    end subroutine evaluate_batch_splines_2d_many_der

    subroutine evaluate_batch_splines_2d_der(spl, x, y_batch, dy_batch)
        type(BatchSplineData2D), intent(in) :: spl
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: y_batch(:)     ! (n_quantities)
        real(dp), intent(out) :: dy_batch(:,:)  ! (2, n_quantities)

        real(dp) :: x_arr(2, 1)
        real(dp) :: y_arr(spl%num_quantities, 1)
        real(dp) :: dy_arr(2, spl%num_quantities, 1)

        x_arr(:, 1) = x
        call evaluate_batch_splines_2d_many_der(spl, x_arr, y_arr, dy_arr)
        y_batch(1:spl%num_quantities) = y_arr(:, 1)
        dy_batch(1:2, 1:spl%num_quantities) = dy_arr(:, :, 1)
    end subroutine evaluate_batch_splines_2d_der
    
end module batch_interpolate_2d
