module interpolate
    implicit none

contains

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
