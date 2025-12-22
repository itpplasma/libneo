module batch_interpolate_3d
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use batch_interpolate_types, only: BatchSplineData3D
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
    public :: construct_batch_splines_3d
    public :: construct_batch_splines_3d_lines
    public :: construct_batch_splines_3d_resident
    public :: construct_batch_splines_3d_resident_device
    public :: destroy_batch_splines_3d

#ifdef LIBNEO_ENABLE_SPLINE_ORACLE
    public :: construct_batch_splines_3d_legacy
#endif
    
    ! Export batch spline evaluation routines
    public :: evaluate_batch_splines_3d
    public :: evaluate_batch_splines_3d_der
    public :: evaluate_batch_splines_3d_der2
    public :: evaluate_batch_splines_3d_many
    public :: evaluate_batch_splines_3d_many_resident
    
contains
    
    subroutine construct_batch_splines_3d(x_min, x_max, y_batch, order, periodic, spl)
        real(dp), intent(in) :: x_min(:), x_max(:)
        real(dp), intent(in) :: y_batch(:,:,:,:)  ! (n1, n2, n3, n_quantities)
        integer, intent(in) :: order(3)
        logical, intent(in) :: periodic(3)
        type(BatchSplineData3D), intent(out) :: spl

        call construct_batch_splines_3d_legacy(x_min, x_max, y_batch, order, periodic, spl)

#ifdef _OPENACC
        !$acc enter data copyin(spl%coeff(1:spl%num_quantities, 0:spl%order(1), 0:spl%order(2), &
        !$acc&                         0:spl%order(3), 1:spl%num_points(1), 1:spl%num_points(2), &
        !$acc&                         1:spl%num_points(3)))
#endif
    end subroutine construct_batch_splines_3d

    subroutine construct_batch_splines_3d_legacy(x_min, x_max, y_batch, order, periodic, spl)
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
    end subroutine construct_batch_splines_3d_legacy

    subroutine construct_batch_splines_3d_lines(x_min, x_max, y_batch, order, periodic, spl)
        real(dp), intent(in) :: x_min(3), x_max(3)
        real(dp), intent(in) :: y_batch(:, :, :, :)  ! (n1, n2, n3, n_quantities)
        integer, intent(in) :: order(3)
        logical, intent(in) :: periodic(3)
        type(BatchSplineData3D), intent(out) :: spl

        integer :: n1, n2, n3, n_quantities
        integer :: N1_order, N2_order, N3_order
        integer :: istat
        integer :: i1, i2, i3, iq, k1, k2, k3
        integer :: line, line2, line3
        real(dp), allocatable :: work3(:, :, :)
        real(dp), allocatable :: work2(:, :, :)
        real(dp), allocatable :: work1(:, :, :)
        real(dp) :: h1, h2, h3

        n1 = size(y_batch, 1)
        n2 = size(y_batch, 2)
        n3 = size(y_batch, 3)
        n_quantities = size(y_batch, 4)
        N1_order = order(1)
        N2_order = order(2)
        N3_order = order(3)

        if (n1 < 2 .or. n2 < 2 .or. n3 < 2) then
            error stop "construct_batch_splines_3d_lines: Need at least 2 points"
        end if
        if (n_quantities < 1) then
            error stop "construct_batch_splines_3d_lines: Need at least 1 quantity"
        end if
        if (any(order < 3) .or. any(order > 5)) then
            error stop "construct_batch_splines_3d_lines: Order must be between 3 and 5"
        end if
        if (N1_order == 3) then
            if (periodic(1) .and. n1 < 3) then
                error stop "construct_batch_splines_3d_lines: Need at least 3 points " // &
                           "for periodic order(1)=3"
            end if
        else if (N1_order == 4) then
            if (n1 < 5) then
                error stop "construct_batch_splines_3d_lines: Need at least 5 points for order(1)=4"
            end if
        else
            if (n1 < 6) then
                error stop "construct_batch_splines_3d_lines: Need at least 6 points for order(1)=5"
            end if
        end if
        if (N2_order == 3) then
            if (periodic(2) .and. n2 < 3) then
                error stop "construct_batch_splines_3d_lines: Need at least 3 points " // &
                           "for periodic order(2)=3"
            end if
        else if (N2_order == 4) then
            if (n2 < 5) then
                error stop "construct_batch_splines_3d_lines: Need at least 5 points for order(2)=4"
            end if
        else
            if (n2 < 6) then
                error stop "construct_batch_splines_3d_lines: Need at least 6 points for order(2)=5"
            end if
        end if
        if (N3_order == 3) then
            if (periodic(3) .and. n3 < 3) then
                error stop "construct_batch_splines_3d_lines: Need at least 3 points " // &
                           "for periodic order(3)=3"
            end if
        else if (N3_order == 4) then
            if (n3 < 5) then
                error stop "construct_batch_splines_3d_lines: Need at least 5 points for order(3)=4"
            end if
        else
            if (n3 < 6) then
                error stop "construct_batch_splines_3d_lines: Need at least 6 points for order(3)=5"
            end if
        end if

        spl%order = order
        spl%num_points = [n1, n2, n3]
        spl%periodic = periodic
        spl%h_step = [(x_max(1) - x_min(1))/dble(n1 - 1), &
                      (x_max(2) - x_min(2))/dble(n2 - 1), &
                      (x_max(3) - x_min(3))/dble(n3 - 1)]
        spl%x_min = x_min
        spl%num_quantities = n_quantities

        allocate(spl%coeff(n_quantities, 0:N1_order, 0:N2_order, 0:N3_order, &
                           n1, n2, n3), stat=istat)
        if (istat /= 0) then
            error stop "construct_batch_splines_3d_lines: Allocation failed for coeff"
        end if

        allocate(work3(n3, n1*n2*n_quantities, 0:N3_order), stat=istat)
        if (istat /= 0) then
            error stop "construct_batch_splines_3d_lines: Allocation failed for work3"
        end if
        allocate(work2(n2, n1*n3*n_quantities, 0:N2_order), stat=istat)
        if (istat /= 0) then
            error stop "construct_batch_splines_3d_lines: Allocation failed for work2"
        end if
        allocate(work1(n1, n2*n3*n_quantities, 0:N1_order), stat=istat)
        if (istat /= 0) then
            error stop "construct_batch_splines_3d_lines: Allocation failed for work1"
        end if

        do iq = 1, n_quantities
            do i2 = 1, n2
                do i1 = 1, n1
                    line = i1 + (i2 - 1)*n1 + (iq - 1)*n1*n2
                    work3(:, line, 0) = y_batch(i1, i2, :, iq)
                end do
            end do
        end do

#ifdef _OPENMP
        !$omp parallel do default(none) shared(work3, n1, n2, n3, n_quantities, N3_order, periodic, &
        !$omp& spl) private(line)
#endif
        do line = 1, n1*n2*n_quantities
            if (N3_order == 3) then
                if (periodic(3)) then
                    call spl_three_per_line(n3, spl%h_step(3), work3(:, line, 0), &
                                            work3(:, line, 1), work3(:, line, 2), &
                                            work3(:, line, 3))
                else
                    call spl_three_reg_line(n3, spl%h_step(3), work3(:, line, 0), &
                                            work3(:, line, 1), work3(:, line, 2), &
                                            work3(:, line, 3))
                end if
            else if (N3_order == 4) then
                if (periodic(3)) then
                    call spl_four_per_line(n3, spl%h_step(3), work3(:, line, 0), &
                                           work3(:, line, 1), work3(:, line, 2), &
                                           work3(:, line, 3), work3(:, line, 4))
                else
                    call spl_four_reg_line(n3, spl%h_step(3), work3(:, line, 0), &
                                           work3(:, line, 1), work3(:, line, 2), &
                                           work3(:, line, 3), work3(:, line, 4))
                end if
            else
                if (periodic(3)) then
                    call spl_five_per_line(n3, spl%h_step(3), work3(:, line, 0), &
                                           work3(:, line, 1), work3(:, line, 2), &
                                           work3(:, line, 3), work3(:, line, 4), &
                                           work3(:, line, 5))
                else
                    call spl_five_reg_line(n3, spl%h_step(3), work3(:, line, 0), &
                                           work3(:, line, 1), work3(:, line, 2), &
                                           work3(:, line, 3), work3(:, line, 4), &
                                           work3(:, line, 5))
                end if
            end if
        end do
#ifdef _OPENMP
        !$omp end parallel do
#endif

        do k3 = 0, N3_order
            do iq = 1, n_quantities
                do i3 = 1, n3
                    do i1 = 1, n1
                        line2 = i1 + (i3 - 1)*n1 + (iq - 1)*n1*n3
                        do i2 = 1, n2
                            line3 = i1 + (i2 - 1)*n1 + (iq - 1)*n1*n2
                            work2(i2, line2, 0) = work3(i3, line3, k3)
                        end do
                    end do
                end do
            end do

#ifdef _OPENMP
            !$omp parallel do default(none) shared(work2, n1, n2, n3, n_quantities, N2_order, periodic, &
            !$omp& spl) private(line)
#endif
            do line = 1, n1*n3*n_quantities
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
                    do i3 = 1, n3
                        do i2 = 1, n2
                            line = i2 + (i3 - 1)*n2 + (iq - 1)*n2*n3
                            do i1 = 1, n1
                                line2 = i1 + (i3 - 1)*n1 + (iq - 1)*n1*n3
                                work1(i1, line, 0) = work2(i2, line2, k2)
                            end do
                        end do
                    end do
                end do

#ifdef _OPENMP
                !$omp parallel do default(none) shared(work1, n1, n2, n3, n_quantities, N1_order, periodic, &
                !$omp& spl) private(line)
#endif
                do line = 1, n2*n3*n_quantities
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
                    do i3 = 1, n3
                        do i2 = 1, n2
                            line = i2 + (i3 - 1)*n2 + (iq - 1)*n2*n3
                            do i1 = 1, n1
                                do k1 = 0, N1_order
                                    spl%coeff(iq, k1, k2, k3, i1, i2, i3) = &
                                        work1(i1, line, k1)
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do

        deallocate(work1)
        deallocate(work2)
        deallocate(work3)
    end subroutine construct_batch_splines_3d_lines

    subroutine construct_batch_splines_3d_resident(x_min, x_max, y_batch, order, &
                                                   periodic, spl)
        real(dp), intent(in) :: x_min(:), x_max(:)
        real(dp), intent(in) :: y_batch(:, :, :, :)  ! (n1, n2, n3, n_quantities)
        integer, intent(in) :: order(:)
        logical, intent(in) :: periodic(:)
        type(BatchSplineData3D), intent(out) :: spl

        call construct_batch_splines_3d(x_min, x_max, y_batch, order, periodic, spl)
    end subroutine construct_batch_splines_3d_resident

    subroutine construct_batch_splines_3d_resident_device(x_min, x_max, y_batch, &
                                                          order, periodic, spl, &
                                                          update_host, &
                                                          assume_y_present)
        real(dp), intent(in) :: x_min(:), x_max(:)
        real(dp), intent(in) :: y_batch(:, :, :, :)  ! (n1, n2, n3, n_quantities)
        integer, intent(in) :: order(3)
        logical, intent(in) :: periodic(3)
        type(BatchSplineData3D), intent(out) :: spl
        logical, intent(in), optional :: update_host
        logical, intent(in), optional :: assume_y_present

        integer :: n1, n2, n3, n_quantities
        logical :: do_update, do_assume_present

        do_update = .true.
        if (present(update_host)) do_update = update_host
        do_assume_present = .false.
        if (present(assume_y_present)) do_assume_present = assume_y_present

        n1 = size(y_batch, 1)
        n2 = size(y_batch, 2)
        n3 = size(y_batch, 3)
        n_quantities = size(y_batch, 4)

#ifdef _OPENACC
        call construct_batch_splines_3d_resident_device_impl(x_min, x_max, n1, n2, &
                                                            n3, n_quantities, y_batch, &
                                                            order, periodic, spl, &
                                                            do_update, do_assume_present)
#else
        call construct_batch_splines_3d_resident(x_min, x_max, y_batch, order, &
                                                 periodic, spl)
#endif
    end subroutine construct_batch_splines_3d_resident_device

    subroutine construct_batch_splines_3d_resident_device_impl(x_min, x_max, n1, n2, &
                                                               n3, n_quantities, &
                                                               y_batch, &
                                                               order, periodic, &
                                                               spl, &
                                                               do_update, &
                                                               do_assume_present)
        real(dp), intent(in) :: x_min(3), x_max(3)
        integer, intent(in) :: n1, n2, n3, n_quantities
        real(dp), intent(in) :: y_batch(n1, n2, n3, n_quantities)
        integer, intent(in) :: order(3)
        logical, intent(in) :: periodic(3)
        type(BatchSplineData3D), intent(out) :: spl
        logical, intent(in) :: do_update
        logical, intent(in) :: do_assume_present

        integer :: istat
        integer :: N1_order, N2_order, N3_order
        integer :: order1, order2, order3
        integer :: i1, i2, i3, iq, k1, k2, k3
        integer :: line, line2, line3
        real(dp), allocatable :: work3(:, :, :)
        real(dp), allocatable :: work2(:, :, :)
        real(dp), allocatable :: work1(:, :, :)
        real(dp) :: h1, h2, h3
        logical :: periodic1, periodic2, periodic3

        N1_order = order(1)
        N2_order = order(2)
        N3_order = order(3)
        order1 = N1_order
        order2 = N2_order
        order3 = N3_order
        periodic1 = periodic(1)
        periodic2 = periodic(2)
        periodic3 = periodic(3)

        if (n1 < 2 .or. n2 < 2 .or. n3 < 2) then
            error stop "construct_batch_splines_3d_resident_device:" // &
                " Need at least 2 points"
        end if
        if (n_quantities < 1) then
            error stop "construct_batch_splines_3d_resident_device:" // &
                " Need at least 1 quantity"
        end if
        if (any(order < 3) .or. any(order > 5)) then
            error stop "construct_batch_splines_3d_resident_device:" // &
                " Order must be between 3 and 5"
        end if
        if (N1_order == 3) then
            if (periodic(1) .and. n1 < 3) then
                error stop "construct_batch_splines_3d_resident_device:" // &
                    " Need at least 3 points for periodic order(1)=3"
            end if
        else if (N1_order == 4) then
            if (n1 < 5) then
                error stop "construct_batch_splines_3d_resident_device:" // &
                    " Need at least 5 points for order(1)=4"
            end if
        else
            if (n1 < 6) then
                error stop "construct_batch_splines_3d_resident_device:" // &
                    " Need at least 6 points for order(1)=5"
            end if
        end if
        if (N2_order == 3) then
            if (periodic(2) .and. n2 < 3) then
                error stop "construct_batch_splines_3d_resident_device:" // &
                    " Need at least 3 points for periodic order(2)=3"
            end if
        else if (N2_order == 4) then
            if (n2 < 5) then
                error stop "construct_batch_splines_3d_resident_device:" // &
                    " Need at least 5 points for order(2)=4"
            end if
        else
            if (n2 < 6) then
                error stop "construct_batch_splines_3d_resident_device:" // &
                    " Need at least 6 points for order(2)=5"
            end if
        end if
        if (N3_order == 3) then
            if (periodic(3) .and. n3 < 3) then
                error stop "construct_batch_splines_3d_resident_device:" // &
                    " Need at least 3 points for periodic order(3)=3"
            end if
        else if (N3_order == 4) then
            if (n3 < 5) then
                error stop "construct_batch_splines_3d_resident_device:" // &
                    " Need at least 5 points for order(3)=4"
            end if
        else
            if (n3 < 6) then
                error stop "construct_batch_splines_3d_resident_device:" // &
                    " Need at least 6 points for order(3)=5"
            end if
        end if

        spl%order = order
        spl%num_points = [n1, n2, n3]
        spl%periodic = periodic
        spl%h_step = [(x_max(1) - x_min(1))/dble(n1 - 1), &
                      (x_max(2) - x_min(2))/dble(n2 - 1), &
                      (x_max(3) - x_min(3))/dble(n3 - 1)]
        spl%x_min = x_min
        spl%num_quantities = n_quantities

        allocate(spl%coeff(n_quantities, 0:N1_order, 0:N2_order, 0:N3_order, &
            n1, n2, n3), stat=istat)
        if (istat /= 0) then
            error stop "construct_batch_splines_3d_resident_device:" // &
                " Allocation failed for coeff"
        end if

        !$acc enter data create(spl%coeff(1:n_quantities, 0:N1_order, 0:N2_order, &
        !$acc&                          0:N3_order, 1:n1, 1:n2, 1:n3))

        allocate(work3(n3, n1*n2*n_quantities, 0:N3_order), stat=istat)
        if (istat /= 0) then
            error stop "construct_batch_splines_3d_resident_device:" // &
                " Allocation failed for work3"
        end if
        allocate(work2(n2, n1*n3*n_quantities, 0:N2_order), stat=istat)
        if (istat /= 0) then
            error stop "construct_batch_splines_3d_resident_device:" // &
                " Allocation failed for work2"
        end if
        allocate(work1(n1, n2*n3*n_quantities, 0:N1_order), stat=istat)
        if (istat /= 0) then
            error stop "construct_batch_splines_3d_resident_device:" // &
                " Allocation failed for work1"
        end if

        h1 = spl%h_step(1)
        h2 = spl%h_step(2)
        h3 = spl%h_step(3)

#ifdef _OPENACC
        block
            logical :: y_was_present
            logical :: did_copyin_y

            y_was_present = acc_is_present(y_batch)
            if (do_assume_present .and. .not. y_was_present) then
                error stop "construct_batch_splines_3d_resident_device:" // &
                    " assume_y_present=T but y_batch is not present"
            end if
            did_copyin_y = .false.
            if (.not. y_was_present) then
                !$acc enter data copyin(y_batch(1:n1, 1:n2, 1:n3, 1:n_quantities))
                did_copyin_y = .true.
            end if

            !$acc enter data create(work3, work2, work1)

            !$acc parallel loop collapse(4) gang present(y_batch, work3)
            do iq = 1, n_quantities
                do i2 = 1, n2
                    do i1 = 1, n1
                        do i3 = 1, n3
                            line = i1 + (i2 - 1)*n1 + (iq - 1)*n1*n2
                            work3(i3, line, 0) = y_batch(i1, i2, i3, iq)
                        end do
                    end do
                end do
            end do

            !$acc parallel loop gang present(work3)
            do line = 1, n1*n2*n_quantities
                call spl_build_line_inplace(order3, periodic3, n3, h3, work3, line)
            end do

            do k3 = 0, N3_order
                !$acc parallel loop collapse(3) gang present(work2, work3)
                do iq = 1, n_quantities
                    do i3 = 1, n3
                        do i1 = 1, n1
                            line2 = i1 + (i3 - 1)*n1 + (iq - 1)*n1*n3
                            !$acc loop vector
                            do i2 = 1, n2
                                line3 = i1 + (i2 - 1)*n1 + (iq - 1)*n1*n2
                                work2(i2, line2, 0) = work3(i3, line3, k3)
                            end do
                        end do
                    end do
                end do

                !$acc parallel loop gang present(work2)
                do line = 1, n1*n3*n_quantities
                    call spl_build_line_inplace(order2, periodic2, n2, h2, work2, &
                                                line)
                end do

                do k2 = 0, N2_order
                    !$acc parallel loop collapse(3) gang present(work1, work2)
                    do iq = 1, n_quantities
                        do i3 = 1, n3
                            do i2 = 1, n2
                                line = i2 + (i3 - 1)*n2 + (iq - 1)*n2*n3
                                !$acc loop vector
                                do i1 = 1, n1
                                    line2 = i1 + (i3 - 1)*n1 + (iq - 1)*n1*n3
                                    work1(i1, line, 0) = work2(i2, line2, k2)
                                end do
                            end do
                        end do
                    end do

                    !$acc parallel loop gang present(work1)
                    do line = 1, n2*n3*n_quantities
                        call spl_build_line_inplace(order1, periodic1, n1, h1, &
                                                    work1, line)
                    end do

                    associate (coeff => spl%coeff)
                        !$acc parallel loop collapse(4) gang present(work1, coeff)
                        do iq = 1, n_quantities
                            do i3 = 1, n3
                                do i2 = 1, n2
                                    do i1 = 1, n1
                                        line = i2 + (i3 - 1)*n2 + (iq - 1)*n2*n3
                                        !$acc loop vector
                                        do k1 = 0, N1_order
                                            coeff(iq, k1, k2, k3, i1, i2, i3) = &
                                                work1(i1, line, k1)
                                        end do
                                    end do
                                end do
                            end do
                        end do
                    end associate
                end do
            end do

            !$acc exit data delete(work3, work2, work1)
            if (did_copyin_y) then
                !$acc exit data delete(y_batch(1:n1, 1:n2, 1:n3, 1:n_quantities))
            end if
        end block
#endif

        deallocate(work1)
        deallocate(work2)
        deallocate(work3)
        if (do_update) then
            !$acc update self(spl%coeff(1:n_quantities, 0:N1_order, 0:N2_order, &
            !$acc&                        0:N3_order, 1:n1, 1:n2, 1:n3))
        end if
    end subroutine construct_batch_splines_3d_resident_device_impl
    
    
    subroutine destroy_batch_splines_3d(spl)
        type(BatchSplineData3D), intent(inout) :: spl

#ifdef _OPENACC
        if (allocated(spl%coeff)) then
            if (acc_is_present(spl%coeff(1:spl%num_quantities, 0:spl%order(1), &
                                         0:spl%order(2), 0:spl%order(3), &
                                         1:spl%num_points(1), 1:spl%num_points(2), &
                                         1:spl%num_points(3)))) then
                !$acc exit data delete(spl%coeff(1:spl%num_quantities, 0:spl%order(1), &
                !$acc&                         0:spl%order(2), 0:spl%order(3), &
                !$acc&                         1:spl%num_points(1), 1:spl%num_points(2), &
                !$acc&                         1:spl%num_points(3)))
            end if
        end if
#endif
        if (allocated(spl%coeff)) deallocate(spl%coeff)
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
        do k3 = 0, N3
            do k2 = 0, N2
                !$omp simd
                do iq = 1, spl%num_quantities
                    coeff_23(iq, k2, k3) = spl%coeff(iq, N1, k2, k3, &
                        interval_index(1)+1, interval_index(2)+1, interval_index(3)+1)
                end do
            end do
        end do
        
        ! Apply backward Horner method (standard polynomial evaluation)
        do k1 = N1-1, 0, -1
            do k3 = 0, N3
                do k2 = 0, N2
                    !$omp simd
                    do iq = 1, spl%num_quantities
                        coeff_23(iq, k2, k3) = spl%coeff(iq, k1, k2, k3, &
                            interval_index(1)+1, interval_index(2)+1, &
                            interval_index(3)+1) + x_local(1)*coeff_23(iq, k2, k3)
                    end do
                end do
            end do
        end do
        
        ! Second reduction: interpolation over x2 using backward Horner
        ! Initialize with highest order
        do k3 = 0, N3
            !$omp simd
            do iq = 1, spl%num_quantities
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

    subroutine evaluate_batch_splines_3d_many(spl, x, y_batch)
        type(BatchSplineData3D), intent(in) :: spl
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: y_batch(:, :)  ! (n_quantities, n_points)

        integer :: ipt, iq, k1, k2, k3, i1, i2, i3, k_wrap
        integer :: nq, order1, order2, order3
        integer :: num_points(3)
        real(dp) :: xj1, xj2, xj3
        real(dp) :: x_norm1, x_norm2, x_norm3
        real(dp) :: x_local1, x_local2, x_local3
        real(dp) :: x_min(3), h_step(3), period(3)
        real(dp) :: t, w, v, w2, yq

        if (size(x, 1) /= 3) then
            error stop "evaluate_batch_splines_3d_many: First dimension of x must be 3"
        end if
        if (size(y_batch, 1) < spl%num_quantities) then
            error stop "evaluate_batch_splines_3d_many: First dimension too small"
        end if
        if (size(y_batch, 2) /= size(x, 2)) then
            error stop "evaluate_batch_splines_3d_many: y_batch second dim mismatch"
        end if

#ifdef _OPENACC
        if (acc_is_present(spl%coeff, size(spl%coeff)*8)) then
            !$acc data copyin(x) copy(y_batch)
            call evaluate_batch_splines_3d_many_resident(spl, x, y_batch)
            !$acc end data
            return
        end if
#endif

        nq = spl%num_quantities
        num_points = spl%num_points
        order1 = spl%order(1)
        order2 = spl%order(2)
        order3 = spl%order(3)
        x_min = spl%x_min
        h_step = spl%h_step
        period(1) = h_step(1)*real(num_points(1) - 1, dp)
        period(2) = h_step(2)*real(num_points(2) - 1, dp)
        period(3) = h_step(3)*real(num_points(3) - 1, dp)

        do ipt = 1, size(x, 2)
            include "spline3d_many_point_body.inc"
        end do
    end subroutine evaluate_batch_splines_3d_many

    subroutine evaluate_batch_splines_3d_many_resident(spl, x, y_batch)
        type(BatchSplineData3D), intent(in) :: spl
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(inout) :: y_batch(:, :)

        integer :: ipt, iq, k1, k2, k3, i1, i2, i3, k_wrap
        integer :: nq, order1, order2, order3
        integer :: num_points(3)
        real(dp) :: xj1, xj2, xj3
        real(dp) :: x_norm1, x_norm2, x_norm3
        real(dp) :: x_local1, x_local2, x_local3
        real(dp) :: x_min(3), h_step(3), period(3)
        real(dp) :: t, w, v, w2, yq

        if (size(x, 1) /= 3) then
            error stop "evaluate_batch_splines_3d_many_resident: x first dim must be 3"
        end if
        if (size(y_batch, 1) < spl%num_quantities) then
            error stop "evaluate_batch_splines_3d_many_resident: y_batch dim1 too small"
        end if
        if (size(y_batch, 2) /= size(x, 2)) then
            error stop "evaluate_batch_splines_3d_many_resident: y_batch dim2 mismatch"
        end if

        nq = spl%num_quantities
        num_points = spl%num_points
        order1 = spl%order(1)
        order2 = spl%order(2)
        order3 = spl%order(3)
        x_min = spl%x_min
        h_step = spl%h_step
        period(1) = h_step(1)*real(num_points(1) - 1, dp)
        period(2) = h_step(2)*real(num_points(2) - 1, dp)
        period(3) = h_step(3)*real(num_points(3) - 1, dp)

        !$acc parallel loop present(spl%coeff, x, y_batch) &
        !$acc& private(ipt, iq, k1, k2, k3, i1, i2, i3, k_wrap) &
        !$acc& private(xj1, xj2, xj3, x_norm1, x_norm2, x_norm3) &
        !$acc& private(x_local1, x_local2, x_local3, t, w, v, w2, yq)
        do ipt = 1, size(x, 2)
            include "spline3d_many_point_body.inc"
        end do
        !$acc end parallel loop
    end subroutine evaluate_batch_splines_3d_many_resident
    
    
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
        do k3 = 0, N3
            do k2 = 0, N2
                !$omp simd
                do iq = 1, spl%num_quantities
                    coeff_23(iq, k2, k3) = spl%coeff(iq, N1, k2, k3, &
                        interval_index(1)+1, interval_index(2)+1, &
                        interval_index(3)+1)
                end do
            end do
        end do
        
        ! Apply backward Horner for value
        do k1 = N1-1, 0, -1
            do k3 = 0, N3
                do k2 = 0, N2
                    !$omp simd
                    do iq = 1, spl%num_quantities
                        coeff_23(iq, k2, k3) = spl%coeff(iq, k1, k2, k3, &
                            interval_index(1)+1, interval_index(2)+1, &
                            interval_index(3)+1) + x_local(1)*coeff_23(iq, k2, k3)
                    end do
                end do
            end do
        end do
        
        ! First derivative over x1 using backward Horner for derivative
        do k3 = 0, N3
            do k2 = 0, N2
                !$omp simd
                do iq = 1, spl%num_quantities
                    coeff_23_dx1(iq, k2, k3) = N1 * spl%coeff(iq, N1, k2, k3, &
                        interval_index(1)+1, interval_index(2)+1, &
                        interval_index(3)+1)
                end do
            end do
        end do
        
        do k1 = N1-1, 1, -1
            do k3 = 0, N3
                do k2 = 0, N2
                    !$omp simd
                    do iq = 1, spl%num_quantities
                        coeff_23_dx1(iq, k2, k3) = k1 * spl%coeff(iq, k1, k2, k3, &
                            interval_index(1)+1, interval_index(2)+1, &
                            interval_index(3)+1) + x_local(1)*coeff_23_dx1(iq, k2, k3)
                    end do
                end do
            end do
        end do
        
        ! Second reduction: interpolation over x2
        ! Initialize with highest order
        do k3 = 0, N3
            !$omp simd
            do iq = 1, spl%num_quantities
                coeff_3(iq, k3) = coeff_23(iq, N2, k3)
                coeff_3_dx1(iq, k3) = coeff_23_dx1(iq, N2, k3)
                coeff_3_dx2(iq, k3) = N2 * coeff_23(iq, N2, k3)
            end do
        end do
        
        ! Apply Horner for value and dx1
        do k2 = N2-1, 0, -1
            do k3 = 0, N3
                !$omp simd
                do iq = 1, spl%num_quantities
                    coeff_3(iq, k3) = coeff_23(iq, k2, k3) + &
                        x_local(2)*coeff_3(iq, k3)
                    coeff_3_dx1(iq, k3) = coeff_23_dx1(iq, k2, k3) + &
                        x_local(2)*coeff_3_dx1(iq, k3)
                end do
            end do
        end do
        
        ! Derivative over x2
        do k2 = N2-1, 1, -1
            do k3 = 0, N3
                !$omp simd
                do iq = 1, spl%num_quantities
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
        do k3 = 0, N3
            do k2 = 0, N2
                !$omp simd
                do iq = 1, spl%num_quantities
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
            do k3 = 0, N3
                do k2 = 0, N2
                    !$omp simd
                    do iq = 1, spl%num_quantities
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
            do k3 = 0, N3
                do k2 = 0, N2
                    !$omp simd
                    do iq = 1, spl%num_quantities
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
        do k3 = 0, N3
            do k2 = 0, N2
                !$omp simd
                do iq = 1, spl%num_quantities
                    c = spl%coeff(iq, k1, k2, k3, interval_index(1)+1, &
                        interval_index(2)+1, interval_index(3)+1)
                    coeff_23(iq, k2, k3) = c + x_local(1)*coeff_23(iq, k2, k3)
                end do
            end do
        end do
        
        ! Second reduction: x2 interpolation with fused derivative computation
        ! Initialize all polynomials with highest order coefficient (k2 = N2)
        do k3 = 0, N3
            !$omp simd
            do iq = 1, spl%num_quantities
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
            do k3 = 0, N3
                !$omp simd
                do iq = 1, spl%num_quantities
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
            do k3 = 0, N3
                !$omp simd
                do iq = 1, spl%num_quantities
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
        do k3 = 0, N3
            !$omp simd
            do iq = 1, spl%num_quantities
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
