program test_coordinate_basis_jets
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_coordinates, only: chartmap_coordinate_system_t
    use libneo_coordinate_jets, only: cylindrical_field_direction_der
    use interpolate, only: construct_batch_splines_3d
    implicit none
    type(chartmap_coordinate_system_t) :: chart
    integer, parameter :: nr = 8, nt = 129, np = 9
    real(dp), parameter :: pi = acos(-1.0_dp), period = 2*pi/3
    real(dp) :: values(nr, nt, np, 3), u(3), b(3, 3), db(3, 3, 3)
    real(dp) :: bp(3, 3), bm(3, 3), shifted(3), fd(3, 3), rho, theta, step
    real(dp) :: expected(3, 3, 3), c, s, cp, sp, err
    integer :: i, j, k, kind, axis

    call check_direction_derivative()
    do kind = 1, 2
        do k = 1, np
            do j = 1, nt
                theta = 2*pi*real(j - 1, dp)/(nt - 1)
                do i = 1, nr
                    rho = real(i - 1, dp)/(nr - 1)
                    values(i, j, k, 1) = 170.0_dp + 35*rho*cos(theta)
                    values(i, j, k, 2) = 35*rho*sin(theta)
                    values(i, j, k, 3) = 0.0_dp
                    if (kind == 2) then
                        values(i, j, k, 3) = values(i, j, k, 2)
                        values(i, j, k, 2) = 0.0_dp
                    end if
                end do
            end do
        end do
        chart%num_field_periods = 3
        chart%has_spl_rz = kind == 1
        if (kind == 1) then
            call construct_batch_splines_3d([0.0_dp, 0.0_dp, 0.0_dp], &
                              [1.0_dp, 2*pi, period], values(:, :, :, 1:2), [5, 5, 5], &
                                            [.false., .true., .true.], chart%spl_rz)
        else
            call construct_batch_splines_3d([0.0_dp, 0.0_dp, 0.0_dp], &
                                            [1.0_dp, 2*pi, period], values, [5, 5, 5], &
                                            [.false., .true., .true.], chart%spl_cart)
        end if
        do i = 1, 3
            u = [10.0_dp**(-i - 1), 0.37_dp, 0.42_dp + period]
            call chart%covariant_basis_with_der(u, b, db)
            c = cos(u(2)); s = sin(u(2)); cp = cos(u(3)); sp = sin(u(3))
            expected = 0.0_dp
            expected(:, 1, 2) = [-35*s*cp, -35*s*sp, 35*c]
            expected(:, 1, 3) = [-35*c*sp, 35*c*cp, 0.0_dp]
            expected(:, 2, 2) = [-35*u(1)*c*cp, -35*u(1)*c*sp, -35*u(1)*s]
            expected(:, 2, 3) = [35*u(1)*s*sp, -35*u(1)*s*cp, 0.0_dp]
            expected(:, 3, 3) = [-(170 + 35*u(1)*c)*cp, &
                                 -(170 + 35*u(1)*c)*sp, 0.0_dp]
            expected(:, 2, 1) = expected(:, 1, 2)
            expected(:, 3, 1) = expected(:, 1, 3)
            expected(:, 3, 2) = expected(:, 2, 3)
            err = maxval(abs(db - expected))
            print *, 'kind', kind, 'rho', u(1), 'analytic Hessian error', err
            if (err > 2.0e-6_dp) error stop 'analytic torus basis derivative'
            do axis = 1, 3
                step = min(1.0e-5_dp, 0.1_dp*u(1))
                shifted = u; shifted(axis) = u(axis) + step
                call chart%covariant_basis(shifted, bp)
                shifted(axis) = u(axis) - step
                call chart%covariant_basis(shifted, bm)
                fd = (bp - bm)/(2*step)
                if (maxval(abs(db(:, :, axis) - fd)) > 2.0e-5_dp) &
                    error stop 'basis derivative disagrees with native basis'
            end do
        end do
    end do
contains
    subroutine check_direction_derivative()
        real(dp) :: position(3), dposition(3, 3), field(3), dfield(3, 3)
        real(dp) :: direction(3), derivative(3, 3), plus(3), minus(3), reference(3)
        real(dp), parameter :: delta = 1.0e-5_dp
        integer :: j
        position = [147.0_dp, -39.0_dp, 12.0_dp]
        dposition = reshape([2.0_dp, -3.0_dp, 0.4_dp, -0.7_dp, 1.3_dp, 2.0_dp, &
                             0.2_dp, 0.5_dp, -0.1_dp], [3, 3])
        field = [0.3_dp, 0.014_dp, 1.2_dp]
        dfield = reshape([0.02_dp, -0.003_dp, 0.03_dp, 0.01_dp, 0.001_dp, &
                          -0.02_dp, -0.03_dp, 0.002_dp, 0.01_dp], [3, 3])
        call cylindrical_field_direction_der(position, dposition, field, dfield, &
                                              direction, derivative)
        reference = explicit_direction(position, field)
        if (maxval(abs(direction - reference)) > 1.0e-13_dp) &
            error stop 'physical cylindrical direction disagrees with Cartesian field'
        if (abs(sum(direction**2) - 1.0_dp) > 1.0e-13_dp) &
            error stop 'physical cylindrical direction is not normalized'
        do j = 1, 3
            plus = explicit_direction(position + delta*dposition(:, j), &
                                      field + delta*dfield(:, j))
            minus = explicit_direction(position - delta*dposition(:, j), &
                                       field - delta*dfield(:, j))
            if (maxval(abs(derivative(:, j) - (plus - minus)/(2*delta))) > 1.0e-8_dp) &
                error stop 'physical cylindrical direction derivative'
            if (abs(dot_product(direction, derivative(:, j))) > 1.0e-13_dp) &
                error stop 'unit-direction derivative is not tangent'
        end do
    end subroutine

    function explicit_direction(position, field) result(direction)
        real(dp), intent(in) :: position(3), field(3)
        real(dp) :: direction(3), radius
        radius = sqrt(sum(position(1:2)**2))
        direction = [position(1)*field(1)/radius - position(2)*field(2), &
                     position(2)*field(1)/radius + position(1)*field(2), field(3)]
        direction = direction/sqrt(sum(direction**2))
    end function
end program test_coordinate_basis_jets
