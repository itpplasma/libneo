program test_chartmap_full_torus
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_coordinates, only: coordinate_system_t, &
        make_chartmap_coordinate_system
    use netcdf
    implicit none

    class(coordinate_system_t), allocatable :: coords
    real(dp), parameter :: pi = acos(-1.0_dp), period = 2.0_dp*pi/3.0_dp
    real(dp) :: u(3), x(3), base(3), basis(3, 3), base_basis(3, 3), rot(3, 3)
    integer :: k

    call write_torus('cyl_full_torus.nc')
    call make_chartmap_coordinate_system(coords, 'cyl_full_torus.nc')
    u = [0.5_dp, pi/2.0_dp, 0.3_dp]
    call coords%evaluate_cart(u, base)
    call coords%covariant_basis(u, base_basis)
    do k = -2, 4
        rot = 0.0_dp
        rot(1, 1) = cos(k*period)
        rot(1, 2) = -sin(k*period)
        rot(2, 1) = sin(k*period)
        rot(2, 2) = cos(k*period)
        rot(3, 3) = 1.0_dp
        u(3) = 0.3_dp + k*period
        call coords%evaluate_cart(u, x)
        call coords%covariant_basis(u, basis)
        if (norm2(x - matmul(rot, base)) > 1.0e-8_dp) &
            error stop 'Cartesian field-period rotation failed'
        if (maxval(abs(basis - matmul(rot, base_basis))) > 1.0e-8_dp) &
            error stop 'Covariant basis field-period rotation failed'
    end do
    ! A small physical step across a seam must remain a small Cartesian segment.
    u(3) = period - 1.0e-6_dp
    call coords%evaluate_cart(u, base)
    u(3) = period + 1.0e-6_dp
    call coords%evaluate_cart(u, x)
    if (norm2(x - base) > 1.0e-3_dp) error stop 'Cartesian seam discontinuity'
    print *, 'PASS: cylindrical chartmap Cartesian positions and bases on full torus'

contains

    subroutine check(status)
        integer, intent(in) :: status
        if (status /= nf90_noerr) then
            print *, nf90_strerror(status)
            error stop 1
        end if
    end subroutine check

    subroutine write_torus(path)
        character(*), intent(in) :: path
        integer, parameter :: nr = 5, nt = 16, nz = 9
        real(dp) :: rho(nr), theta(nt), zeta(nz), xyz(nr, nt, nz, 3), radius
        integer :: nc, dims(3), vars(3), pos(3), nfp, i, j, l
        character(len=5), parameter :: names(3) = ['rho  ', 'theta', 'zeta ']
        character(len=1), parameter :: axes(3) = ['x', 'y', 'z']

        rho = [(real(i - 1, dp)/(nr - 1), i=1, nr)]
        theta = [(2.0_dp*pi*real(i - 1, dp)/nt, i=1, nt)]
        zeta = [(period*real(i - 1, dp)/nz, i=1, nz)]
        do l = 1, nz
            do j = 1, nt
                do i = 1, nr
                    radius = 170.0_dp + 35.0_dp*rho(i)*cos(theta(j))
                    xyz(i, j, l, :) = [radius*cos(zeta(l)), radius*sin(zeta(l)), &
                        35.0_dp*rho(i)*sin(theta(j))]
                end do
            end do
        end do
        call check(nf90_create(path, nf90_clobber, nc))
        call check(nf90_put_att(nc, nf90_global, 'zeta_convention', 'cyl'))
        call check(nf90_def_dim(nc, 'rho', nr, dims(1)))
        call check(nf90_def_dim(nc, 'theta', nt, dims(2)))
        call check(nf90_def_dim(nc, 'zeta', nz, dims(3)))
        do i = 1, 3
            call check(nf90_def_var(nc, trim(names(i)), nf90_double, &
                [dims(i)], vars(i)))
            call check(nf90_def_var(nc, axes(i), nf90_double, dims, pos(i)))
            call check(nf90_put_att(nc, pos(i), 'units', 'cm'))
        end do
        call check(nf90_def_var(nc, 'num_field_periods', nf90_int, nfp))
        call check(nf90_enddef(nc))
        call check(nf90_put_var(nc, vars(1), rho))
        call check(nf90_put_var(nc, vars(2), theta))
        call check(nf90_put_var(nc, vars(3), zeta))
        call check(nf90_put_var(nc, nfp, 3))
        do i = 1, 3
            call check(nf90_put_var(nc, pos(i), xyz(:, :, :, i)))
        end do
        call check(nf90_close(nc))
    end subroutine write_torus
end program test_chartmap_full_torus
