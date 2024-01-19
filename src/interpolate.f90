module interpolate
    implicit none

    type :: SplineData1D
        integer :: order
        integer :: num_points
        logical :: periodic
        real(8) :: h_step
        real(8), dimension(:,:), allocatable :: coeff
    end type SplineData1D

    type :: SplineData2D
        integer :: order
        integer :: num_points
        logical :: periodic(2)
        real(8) :: h_step
        real(8), dimension(:,:,:,:), allocatable :: coeff
    end type SplineData2D

    type :: SplineData3D
        integer :: order
        integer :: num_points
        logical :: periodic(3)
        real(8) :: h_step
        real(8), dimension(:,:,:,:,:,:,:), allocatable :: coeff
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
        call spl_reg(spl%order, spl%num_points, spl%h_step, spl%coeff)
    end subroutine construct_splines_1d


    subroutine destroy_splines_1d(spl)
        type(SplineData1D), intent(inout) :: spl

        deallocate(spl%coeff)
    end subroutine destroy_splines_1d


    subroutine evaluate_splines_1d(x, spl, y)
        real(8), intent(in) :: x
        type(SplineData1D), intent(in) :: spl
        real(8), intent(out) :: y

        real(8) :: x_norm, x_local, local_coeff(0:spl%order)
        integer :: interval_index, k_power

        x_norm = x/spl%h_step
        interval_index = max(0, min(spl%num_points-1, int(x_norm)))
        x_local = (x_norm - dble(interval_index))*spl%h_step  ! Distance to grid point

        local_coeff = spl%coeff(:, interval_index+1)

        ! Start with largest power and then multiply recursively
        y = local_coeff(spl%order+1)
        do k_power = spl%order, 0, -1
            y = local_coeff(k_power) + x_local*y
        enddo
    end subroutine evaluate_splines_1d


    subroutine construct_splines_3d(spl_data, spl_order, n_r, n_z, n_phi, h_r, h_z, h_phi)
        integer, intent(in) :: spl_order
        integer, intent(in) :: n_r, n_z, n_phi
        real(8), intent(in) :: h_r, h_z, h_phi
        real(8), intent(inout) :: spl_data(:,:,:,:,:,:,:)

        real(8), dimension(:,:), allocatable  :: splcoe

        integer :: n_data
        integer :: i_r, i_z, i_phi            ! Loop counters for grid
        integer :: i_r_z, i_r_phi, i_data, k  ! Loop counters for splines

        n_data = size(spl_data, 1)

        ! Spline over $\varphi$
        allocate(splcoe(0:spl_order,n_phi))
        do i_r=1,n_r
        do i_z=1,n_z
            do i_data=1,n_data
                splcoe(0,:)=spl_data(i_data,1,1,1,i_r,i_z,:)
                call spl_per(spl_order,n_phi,h_phi,splcoe)
                do k=1,spl_order
                    spl_data(i_data,1,1,k+1,i_r,i_z,:)=splcoe(k,:)
                enddo
            enddo
        enddo
        enddo
        deallocate(splcoe)

        ! Spline over $Z$
        allocate(splcoe(0:spl_order,n_z))
        do i_r=1,n_r
        do i_phi=1,n_phi
            do i_r_phi=1,spl_order+1
                do i_data=1,n_data
                    splcoe(0,:)=spl_data(i_data,1,1,i_r_phi,i_r,:,i_phi)
                    call spl_reg(spl_order,n_z,h_z,splcoe)
                    do k=1,spl_order
                        spl_data(i_data,1,k+1,i_r_phi,i_r,:,i_phi)=splcoe(k,:)
                    enddo
                enddo
            enddo
        enddo
        enddo
        deallocate(splcoe)

        ! Spline over $R$
        allocate(splcoe(0:spl_order,n_r))
        do i_z=1,n_z
        do i_phi=1,n_phi
            do i_r_z=1,spl_order+1
            do i_r_phi=1,spl_order+1
                do i_data=1,n_data
                    splcoe(0,:)=spl_data(i_data,1,i_r_z,i_r_phi,:,i_z,i_phi)
                    call spl_reg(spl_order,n_r,h_r,splcoe)
                    do k=1,spl_order
                        spl_data(i_data,k+1,i_r_z,i_r_phi,:,i_z,i_phi)=splcoe(k,:)
                    enddo
                enddo
            enddo
            enddo
        enddo
        enddo
        deallocate(splcoe)

    end subroutine construct_splines_3d

end module interpolate
