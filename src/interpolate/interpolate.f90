module interpolate
    use spl_three_to_five_sub
    use batch_interpolate
    use cgls_dense, only : cgls_dense_solve
    
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    
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

    interface disp
        module procedure disp_3d
    end interface disp

contains

    subroutine build_design_matrix_1d(x_min, x_max, order, periodic, num_points, &
        x_data, phi)
        real(dp), intent(in) :: x_min, x_max
        integer, intent(in) :: order, num_points
        logical, intent(in) :: periodic
        real(dp), intent(in) :: x_data(:)
        real(dp), intent(out) :: phi(:,:)

        integer :: n_data, i, j
        real(dp), allocatable :: y_basis(:)
        type(SplineData1D) :: spl

        n_data = size(x_data)
        if (size(phi, 1) /= n_data .or. size(phi, 2) /= num_points) then
            error stop "build_design_matrix_1d: phi has wrong shape"
        end if

        allocate(y_basis(num_points))

        !$omp parallel do default(shared) private(j, y_basis, spl, i)
        do j = 1, num_points
            y_basis = 0.0_dp
            y_basis(j) = 1.0_dp
            call construct_splines_1d(x_min, x_max, y_basis, order, periodic, spl)
            do i = 1, n_data
                call evaluate_splines_1d(spl, x_data(i), phi(i, j))
            end do
            call destroy_splines_1d(spl)
        end do
        !$omp end parallel do

        deallocate(y_basis)
    end subroutine build_design_matrix_1d


    subroutine build_design_matrix_2d(x_min, x_max, order, periodic, num_points, &
        x_data, y_data, phi)
        real(dp), intent(in) :: x_min(2), x_max(2)
        integer, intent(in) :: order(2), num_points(2)
        logical, intent(in) :: periodic(2)
        real(dp), intent(in) :: x_data(:), y_data(:)
        real(dp), intent(out) :: phi(:,:)

        integer :: n_data, i, i1, i2, col
        real(dp), allocatable :: y_basis(:,:)
        type(SplineData2D) :: spl

        n_data = size(x_data)
        if (size(y_data) /= n_data) then
            error stop "build_design_matrix_2d: data size mismatch"
        end if
        if (size(phi, 1) /= n_data .or. size(phi, 2) /= num_points(1)*num_points(2)) then
            error stop "build_design_matrix_2d: phi wrong shape"
        end if

        allocate(y_basis(num_points(1), num_points(2)))

        !$omp parallel do default(shared) private(i1, i2, col, y_basis, spl, i)
        do col = 1, num_points(1)*num_points(2)
            i2 = (col - 1)/num_points(1) + 1
            i1 = col - (i2 - 1)*num_points(1)
            y_basis = 0.0_dp
            y_basis(i1, i2) = 1.0_dp
            call construct_splines_2d(x_min, x_max, y_basis, order, periodic, spl)
            do i = 1, n_data
                call evaluate_splines_2d(spl, [x_data(i), y_data(i)], phi(i, col))
            end do
            call destroy_splines_2d(spl)
        end do
        !$omp end parallel do

        deallocate(y_basis)
    end subroutine build_design_matrix_2d


    subroutine build_design_matrix_3d(x_min, x_max, order, periodic, num_points, &
        x_data, y_data, z_data, phi)
        real(dp), intent(in) :: x_min(3), x_max(3)
        integer, intent(in) :: order(3), num_points(3)
        logical, intent(in) :: periodic(3)
        real(dp), intent(in) :: x_data(:), y_data(:), z_data(:)
        real(dp), intent(out) :: phi(:,:)

        integer :: n_data, i, i1, i2, i3, col
        real(dp), allocatable :: y_basis(:,:,:)
        type(SplineData3D) :: spl

        n_data = size(x_data)
        if (size(y_data) /= n_data .or. size(z_data) /= n_data) then
            error stop "build_design_matrix_3d: data size mismatch"
        end if
        if (size(phi, 1) /= n_data .or. size(phi, 2) /= &
            num_points(1)*num_points(2)*num_points(3)) then
            error stop "build_design_matrix_3d: phi wrong shape"
        end if

        allocate(y_basis(num_points(1), num_points(2), num_points(3)))

        !$omp parallel do default(shared) private(i1, i2, i3, col, y_basis, spl, i)
        do col = 1, num_points(1)*num_points(2)*num_points(3)
            i3 = (col - 1)/(num_points(1)*num_points(2)) + 1
            i2 = (col - (i3 - 1)*num_points(1)*num_points(2) - 1)/ &
                num_points(1) + 1
            i1 = col - (i3 - 1)*num_points(1)*num_points(2) - (i2 - 1)* &
                num_points(1)
            y_basis = 0.0_dp
            y_basis(i1, i2, i3) = 1.0_dp
            call construct_splines_3d(x_min, x_max, y_basis, order, periodic, spl)
            do i = 1, n_data
                call evaluate_splines_3d(spl, [x_data(i), y_data(i), z_data(i)], &
                    phi(i, col))
            end do
            call destroy_splines_3d(spl)
        end do
        !$omp end parallel do

        deallocate(y_basis)
    end subroutine build_design_matrix_3d



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

    subroutine construct_splines_1d_lsq(x_min, x_max, order, periodic, &
        num_points, x_data, f_data, spl, weights)
        real(dp), intent(in) :: x_min, x_max
        integer, intent(in) :: order
        logical, intent(in) :: periodic
        integer, intent(in) :: num_points
        real(dp), intent(in) :: x_data(:)
        real(dp), intent(in) :: f_data(:)
        real(dp), intent(in), optional :: weights(:)
        type(SplineData1D), intent(out) :: spl

        integer :: n_data, n_keep, i
        real(dp), allocatable :: x_used(:), f_used(:), w_used(:)
        real(dp), allocatable :: phi(:,:), y(:)
        real(dp) :: h_ref, tol_x
        logical :: use_direct

        n_data = size(x_data)
        if (n_data /= size(f_data)) then
            error stop "construct_splines_1d_lsq: x_data and f_data size mismatch"
        end if
        if (num_points < 2) then
            error stop "construct_splines_1d_lsq: need at least 2 grid points"
        end if

        use_direct = .false.
        tol_x = 1.0d-12
        if (.not. periodic) then
            if (n_data == num_points) then
                h_ref = (x_max - x_min)/dble(num_points - 1)
                use_direct = .true.
                do i = 1, num_points
                    if (abs(x_data(i) - (x_min + h_ref*dble(i - 1))) > tol_x) then
                        use_direct = .false.
                        exit
                    end if
                end do
            end if
        end if

        if (use_direct) then
            call construct_splines_1d(x_min, x_max, f_data, order, periodic, spl)
            return
        end if

        n_keep = 0
        do i = 1, n_data
            if (.not. periodic) then
                if (x_data(i) < x_min .or. x_data(i) > x_max) cycle
            end if
            if (present(weights)) then
                if (weights(i) == 0.0_dp) cycle
            end if
            n_keep = n_keep + 1
        end do

        if (n_keep == 0) then
            error stop "construct_splines_1d_lsq: no usable data"
        end if

        allocate(x_used(n_keep), f_used(n_keep))
        if (present(weights)) allocate(w_used(n_keep))

        n_keep = 0
        do i = 1, n_data
            if (.not. periodic) then
                if (x_data(i) < x_min .or. x_data(i) > x_max) cycle
            end if
            if (present(weights)) then
                if (weights(i) == 0.0_dp) cycle
                w_used(n_keep + 1) = weights(i)
            end if
            n_keep = n_keep + 1
            x_used(n_keep) = x_data(i)
            f_used(n_keep) = f_data(i)
        end do

        allocate(phi(n_keep, num_points))
        call build_design_matrix_1d(x_min, x_max, order, periodic, num_points, &
            x_used, phi)

        allocate(y(num_points))
        if (present(weights)) then
            call cgls_dense_solve(phi, f_used, y, w_used, &
                max_iter = 4*num_points, tol = 1.0d-14)
        else
            call cgls_dense_solve(phi, f_used, y, &
                max_iter = 4*num_points, tol = 1.0d-14)
        end if

        call construct_splines_1d(x_min, x_max, y, order, periodic, spl)

        deallocate(phi, y, x_used, f_used)
        if (present(weights)) deallocate(w_used)
    end subroutine construct_splines_1d_lsq


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

    subroutine construct_splines_2d_lsq(x_min, x_max, order, periodic, &
        num_points, x_data, y_data, f_data, spl, weights)
        real(dp), intent(in) :: x_min(2), x_max(2)
        integer, intent(in) :: order(2)
        logical, intent(in) :: periodic(2)
        integer, intent(in) :: num_points(2)
        real(dp), intent(in) :: x_data(:)
        real(dp), intent(in) :: y_data(:)
        real(dp), intent(in) :: f_data(:)
        real(dp), intent(in), optional :: weights(:)
        type(SplineData2D), intent(out) :: spl

        integer :: n_data, n_keep, i
        real(dp), allocatable :: x_used(:), y_used(:), f_used(:), w_used(:)
        real(dp), allocatable :: phi(:,:), y_vec(:), y_grid(:,:)
        integer :: n_cols

        n_data = size(x_data)
        if (n_data /= size(y_data) .or. n_data /= size(f_data)) then
            error stop "construct_splines_2d_lsq: data size mismatch"
        end if
        if (num_points(1) < 2 .or. num_points(2) < 2) then
            error stop "construct_splines_2d_lsq: need at least 2 grid points"
        end if

        n_keep = 0
        do i = 1, n_data
            if (.not. periodic(1)) then
                if (x_data(i) < x_min(1) .or. x_data(i) > x_max(1)) cycle
            end if
            if (.not. periodic(2)) then
                if (y_data(i) < x_min(2) .or. y_data(i) > x_max(2)) cycle
            end if
            if (present(weights)) then
                if (weights(i) == 0.0_dp) cycle
            end if
            n_keep = n_keep + 1
        end do

        if (n_keep == 0) then
            error stop "construct_splines_2d_lsq: no usable data"
        end if

        allocate(x_used(n_keep), y_used(n_keep), f_used(n_keep))
        if (present(weights)) allocate(w_used(n_keep))

        n_keep = 0
        do i = 1, n_data
            if (.not. periodic(1)) then
                if (x_data(i) < x_min(1) .or. x_data(i) > x_max(1)) cycle
            end if
            if (.not. periodic(2)) then
                if (y_data(i) < x_min(2) .or. y_data(i) > x_max(2)) cycle
            end if
            if (present(weights)) then
                if (weights(i) == 0.0_dp) cycle
                w_used(n_keep + 1) = weights(i)
            end if
            n_keep = n_keep + 1
            x_used(n_keep) = x_data(i)
            y_used(n_keep) = y_data(i)
            f_used(n_keep) = f_data(i)
        end do

        n_cols = num_points(1)*num_points(2)
        allocate(phi(n_keep, n_cols))
        call build_design_matrix_2d(x_min, x_max, order, periodic, num_points, &
            x_used, y_used, phi)

        allocate(y_vec(n_cols))
        if (present(weights)) then
            call cgls_dense_solve(phi, f_used, y_vec, w_used, &
                max_iter = 4*n_cols, tol = 1.0d-14)
        else
            call cgls_dense_solve(phi, f_used, y_vec, &
                max_iter = 4*n_cols, tol = 1.0d-14)
        end if

        allocate(y_grid(num_points(1), num_points(2)))
        y_grid = reshape(y_vec, shape=[num_points(1), num_points(2)])
        call construct_splines_2d(x_min, x_max, y_grid, order, periodic, spl)

        deallocate(phi, y_vec, y_grid, x_used, y_used, f_used)
        if (present(weights)) deallocate(w_used)
    end subroutine construct_splines_2d_lsq


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


    subroutine construct_splines_3d_lsq(x_min, x_max, order, periodic, &
        num_points, x_data, y_data, z_data, f_data, spl, weights)
        real(dp), intent(in) :: x_min(3), x_max(3)
        integer, intent(in) :: order(3)
        logical, intent(in) :: periodic(3)
        integer, intent(in) :: num_points(3)
        real(dp), intent(in) :: x_data(:)
        real(dp), intent(in) :: y_data(:)
        real(dp), intent(in) :: z_data(:)
        real(dp), intent(in) :: f_data(:)
        real(dp), intent(in), optional :: weights(:)
        type(SplineData3D), intent(out) :: spl

        integer :: n_data, n_keep, i
        real(dp), allocatable :: x_used(:), y_used(:), z_used(:), f_used(:)
        real(dp), allocatable :: w_used(:)
        real(dp), allocatable :: phi(:,:), y_vec(:), y_grid(:,:,:)
        integer :: n_cols

        n_data = size(x_data)
        if (n_data /= size(y_data) .or. n_data /= size(z_data) .or. &
            n_data /= size(f_data)) then
            error stop "construct_splines_3d_lsq: data size mismatch"
        end if
        if (num_points(1) < 2 .or. num_points(2) < 2 .or. num_points(3) < 2) then
            error stop "construct_splines_3d_lsq: need at least 2 grid points"
        end if

        n_keep = 0
        do i = 1, n_data
            if (.not. periodic(1)) then
                if (x_data(i) < x_min(1) .or. x_data(i) > x_max(1)) cycle
            end if
            if (.not. periodic(2)) then
                if (y_data(i) < x_min(2) .or. y_data(i) > x_max(2)) cycle
            end if
            if (.not. periodic(3)) then
                if (z_data(i) < x_min(3) .or. z_data(i) > x_max(3)) cycle
            end if
            if (present(weights)) then
                if (weights(i) == 0.0_dp) cycle
            end if
            n_keep = n_keep + 1
        end do

        if (n_keep == 0) then
            error stop "construct_splines_3d_lsq: no usable data"
        end if

        allocate(x_used(n_keep), y_used(n_keep), z_used(n_keep), f_used(n_keep))
        if (present(weights)) allocate(w_used(n_keep))

        n_keep = 0
        do i = 1, n_data
            if (.not. periodic(1)) then
                if (x_data(i) < x_min(1) .or. x_data(i) > x_max(1)) cycle
            end if
            if (.not. periodic(2)) then
                if (y_data(i) < x_min(2) .or. y_data(i) > x_max(2)) cycle
            end if
            if (.not. periodic(3)) then
                if (z_data(i) < x_min(3) .or. z_data(i) > x_max(3)) cycle
            end if
            if (present(weights)) then
                if (weights(i) == 0.0_dp) cycle
                w_used(n_keep + 1) = weights(i)
            end if
            n_keep = n_keep + 1
            x_used(n_keep) = x_data(i)
            y_used(n_keep) = y_data(i)
            z_used(n_keep) = z_data(i)
            f_used(n_keep) = f_data(i)
        end do

        n_cols = num_points(1)*num_points(2)*num_points(3)
        allocate(phi(n_keep, n_cols))
        call build_design_matrix_3d(x_min, x_max, order, periodic, num_points, &
            x_used, y_used, z_used, phi)

        allocate(y_vec(n_cols))
        if (present(weights)) then
            call cgls_dense_solve(phi, f_used, y_vec, w_used, &
                max_iter = 4*n_cols, tol = 1.0d-14)
        else
            call cgls_dense_solve(phi, f_used, y_vec, &
                max_iter = 4*n_cols, tol = 1.0d-14)
        end if

        allocate(y_grid(num_points(1), num_points(2), num_points(3)))
        y_grid = reshape(y_vec, shape=[num_points(1), num_points(2), num_points(3)])
        call construct_splines_3d(x_min, x_max, y_grid, order, periodic, spl)

        deallocate(phi, y_vec, y_grid, x_used, y_used, z_used, f_used)
        if (present(weights)) deallocate(w_used)

    end subroutine construct_splines_3d_lsq



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

end module interpolate
