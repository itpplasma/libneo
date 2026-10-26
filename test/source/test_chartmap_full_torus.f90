program test_chartmap_full_torus
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_coordinates, only: coordinate_system_t, chartmap_coordinate_system_t, &
        make_chartmap_coordinate_system
    use netcdf
    implicit none

    class(coordinate_system_t), allocatable :: coords
    real(dp), parameter :: pi = acos(-1.0_dp), period = 2.0_dp*pi/3.0_dp
    real(dp) :: u(3), x(3), base(3), basis(3, 3), base_basis(3, 3), rot(3, 3)
    real(dp) :: recovered(3), roundtrip(3), cyl(3)
    integer :: k, ierr

    call write_torus('cyl_full_torus.nc')
    call make_chartmap_coordinate_system(coords, 'cyl_full_torus.nc')
    select type (coords)
    type is (chartmap_coordinate_system_t)
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
            call coords%from_cart(x, recovered, ierr)
            if (ierr /= 0) error stop 'Full-torus Cartesian inverse failed'
            call coords%evaluate_cart(recovered, roundtrip)
            if (norm2(x - roundtrip) > 1.0e-7_dp) &
                error stop 'Full-torus Cartesian round trip failed'
            call coords%evaluate_cyl(u, cyl)
            roundtrip = [cyl(1)*cos(cyl(2)), cyl(1)*sin(cyl(2)), cyl(3)]
            if (norm2(x - roundtrip) > 1.0e-8_dp) &
                error stop 'Cylindrical and Cartesian positions disagree'
            call coords%from_cyl(cyl, recovered, ierr)
            if (ierr /= 0) error stop 'Full-torus cylindrical inverse failed'
            call coords%evaluate_cart(recovered, roundtrip)
            if (norm2(x - roundtrip) > 1.0e-7_dp) &
                error stop 'Full-torus cylindrical round trip failed'
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
    class default
        error stop 'Unexpected chartmap coordinate type'
    end select

    call check_cold_cartesian_inverse()

contains

    subroutine check_cold_cartesian_inverse()
        class(coordinate_system_t), allocatable :: cold_coords

        call write_torus('cyl_cold_cartesian.nc', 1000.0_dp, 300.0_dp, 2)
        call make_chartmap_coordinate_system(cold_coords, 'cyl_cold_cartesian.nc')
        select type (cold_coords)
        type is (chartmap_coordinate_system_t)
            call check_cold_points(cold_coords)
        class default
            error stop 'Unexpected cold chartmap coordinate type'
        end select
    end subroutine check_cold_cartesian_inverse

    subroutine check_cold_points(cold_coords)
        type(chartmap_coordinate_system_t), intent(in) :: cold_coords
        real(dp), parameter :: rho_values(3) = [0.2_dp, 0.7_dp, 0.93_dp]
        real(dp), parameter :: phi_values(8) = [-0.1_dp, 1.0e-6_dp, &
            pi - 1.0e-6_dp, pi + 1.0e-6_dp, 1.2_dp*pi, 1.6_dp*pi, &
            2.0_dp*pi - 1.0e-6_dp, 2.0_dp*pi + 0.1_dp]
        real(dp), parameter :: outside_phi(2) = [-0.1_dp, pi + 0.1_dp]
        real(dp), parameter :: tol_cm = 1.0e-7_dp, cm_to_m = 0.01_dp
        real(dp), parameter :: physical_tol_m = 1.0e-7_dp
        real(dp) :: point(3), target(3), target_cyl(3), cart_u(3), cyl_u(3)
        real(dp) :: x_back(3), cart_res_cm, cyl_res_cm, max_cart_cm, max_cyl_cm
        integer :: ir, it, ip, cart_status, cyl_status, failures, rows, output

        open(newunit=output, file='cyl_cold_cartesian_results.dat', status='replace')
        write(output, '(a)') '# kind rho theta phi cart_status cyl_status ' // &
            'cart_res_cm cyl_res_cm cart_res_m cyl_res_m ' // &
            'cart_0p1um_residual_pass cyl_0p1um_residual_pass'
        failures = 0
        rows = 0
        max_cart_cm = 0.0_dp
        max_cyl_cm = 0.0_dp
        do ip = 1, size(phi_values)
            do it = 1, 8
                do ir = 1, size(rho_values)
                    point = [rho_values(ir), 0.17_dp + (it - 1)*pi/4.0_dp, &
                        phi_values(ip)]
                    call cold_coords%evaluate_cart(point, target)
                    target_cyl = [sqrt(target(1)**2 + target(2)**2), &
                        atan2(target(2), target(1)), target(3)]
                    call cold_coords%from_cart(target, cart_u, cart_status)
                    call cold_coords%from_cyl(target_cyl, cyl_u, cyl_status)
                    cart_res_cm = huge(1.0_dp)
                    cyl_res_cm = huge(1.0_dp)
                    if (cart_status == 0) then
                        call cold_coords%evaluate_cart(cart_u, x_back)
                        cart_res_cm = norm2(target - x_back)
                    end if
                    if (cyl_status == 0) then
                        call cold_coords%evaluate_cart(cyl_u, x_back)
                        cyl_res_cm = norm2(target - x_back)
                    end if
                    write(output, *) 'interior', point, cart_status, cyl_status, &
                        cart_res_cm, cyl_res_cm, cm_to_m*cart_res_cm, &
                        cm_to_m*cyl_res_cm, &
                        cart_status == 0 .and. cm_to_m*cart_res_cm <= physical_tol_m, &
                        cyl_status == 0 .and. cm_to_m*cyl_res_cm <= physical_tol_m
                    rows = rows + 1
                    max_cart_cm = max(max_cart_cm, cart_res_cm)
                    max_cyl_cm = max(max_cyl_cm, cyl_res_cm)
                    if (cart_status /= 0 .or. .not. (cart_res_cm <= tol_cm)) &
                        failures = failures + 1
                    if (cyl_status /= 0 .or. .not. (cyl_res_cm <= tol_cm)) &
                        failures = failures + 1
                end do
            end do
        end do
        if (rows /= 192) error stop 'Unexpected cold Cartesian case count'

        ! These physical targets lie beyond the analytic circular last surface.
        ! Construct them directly because evaluate_cart clamps rho to the chart.
        do ip = 1, size(outside_phi)
            target = [1315.0_dp*cos(outside_phi(ip)), &
                1315.0_dp*sin(outside_phi(ip)), 0.0_dp]
            target_cyl = [1315.0_dp, outside_phi(ip), 0.0_dp]
            call cold_coords%from_cart(target, cart_u, cart_status)
            call cold_coords%from_cyl(target_cyl, cyl_u, cyl_status)
            cart_res_cm = huge(1.0_dp)
            cyl_res_cm = huge(1.0_dp)
            if (cart_status == 0) then
                call cold_coords%evaluate_cart(cart_u, x_back)
                cart_res_cm = norm2(target - x_back)
                failures = failures + 1
            end if
            if (cyl_status == 0) then
                call cold_coords%evaluate_cart(cyl_u, x_back)
                cyl_res_cm = norm2(target - x_back)
                failures = failures + 1
            end if
            write(output, *) 'outside', 1.05_dp, 0.0_dp, outside_phi(ip), &
                cart_status, cyl_status, cart_res_cm, cyl_res_cm, &
                cm_to_m*cart_res_cm, cm_to_m*cyl_res_cm, &
                cart_status == 0 .and. cm_to_m*cart_res_cm <= physical_tol_m, &
                cyl_status == 0 .and. cm_to_m*cyl_res_cm <= physical_tol_m
        end do
        close(output)
        print *, 'Cold inverse interior cases:', rows, ' failures:', failures
        print *, 'Maximum Cartesian/cylindrical roundtrip residuals [cm]:', &
            max_cart_cm, max_cyl_cm
        if (failures /= 0) error stop 'Cold Cartesian/cylindrical inverse check failed'
    end subroutine check_cold_points

    subroutine check(status)
        integer, intent(in) :: status
        if (status /= nf90_noerr) then
            print *, nf90_strerror(status)
            error stop 1
        end if
    end subroutine check

    subroutine write_torus(path, major_radius, minor_radius, field_periods)
        character(*), intent(in) :: path
        real(dp), intent(in), optional :: major_radius, minor_radius
        integer, intent(in), optional :: field_periods
        integer, parameter :: nr = 5, nt = 16, nz = 9
        real(dp) :: rho(nr), theta(nt), zeta(nz), xyz(nr, nt, nz, 3), radius
        real(dp) :: r_major, r_minor, torus_period
        integer :: nperiods
        integer :: nc, dims(3), vars(3), pos(3), nfp, i, j, l
        character(len=5), parameter :: names(3) = ['rho  ', 'theta', 'zeta ']
        character(len=1), parameter :: axes(3) = ['x', 'y', 'z']

        r_major = 170.0_dp
        r_minor = 35.0_dp
        nperiods = 3
        if (present(major_radius)) r_major = major_radius
        if (present(minor_radius)) r_minor = minor_radius
        if (present(field_periods)) nperiods = field_periods
        torus_period = 2.0_dp*pi/real(nperiods, dp)
        rho = [(real(i - 1, dp)/(nr - 1), i=1, nr)]
        theta = [(2.0_dp*pi*real(i - 1, dp)/nt, i=1, nt)]
        zeta = [(torus_period*real(i - 1, dp)/nz, i=1, nz)]
        do l = 1, nz
            do j = 1, nt
                do i = 1, nr
                    radius = r_major + r_minor*rho(i)*cos(theta(j))
                    xyz(i, j, l, :) = [radius*cos(zeta(l)), radius*sin(zeta(l)), &
                        r_minor*rho(i)*sin(theta(j))]
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
        call check(nf90_put_var(nc, nfp, nperiods))
        do i = 1, 3
            call check(nf90_put_var(nc, pos(i), xyz(:, :, :, i)))
        end do
        call check(nf90_close(nc))
    end subroutine write_torus
end program test_chartmap_full_torus
