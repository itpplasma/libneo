program test_spectre_inverse_status
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use spectre_reader, only: spectre_data_t
    use libneo_coordinates, only: coordinate_system_t, &
                                  spectre_coordinate_system_t, &
                                  make_spectre_coordinate_system, SPECTRE_SUCCESS, &
                                  SPECTRE_OUTSIDE_DOMAIN, SPECTRE_NOT_CONVERGED, &
                                  SPECTRE_DOMAIN_UNDETERMINED, SPECTRE_INVALID_VOLUME
    implicit none

    type(spectre_data_t) :: data
    class(coordinate_system_t), allocatable :: cs
    real(dp) :: u(3), saved(3), point(3), angle, radius, phi, center(2)
    integer :: exterior_status, singular_status, status, i, j, k

    ! Synthetic analytic domain: (R - 10)**2 + Z**2 <= 1.
    data%Mvol = 1
    data%mn = 2
    allocate (data%im(2), data%in(2))
    data%im = [0, 1]
    data%in = [0, 0]
    allocate (data%Rbc(2, 0:1), data%Zbs(2, 0:1))
    allocate (data%Rbs(2, 0:1), data%Zbc(2, 0:1))
    data%Rbs = 0.0_dp
    data%Zbc = 0.0_dp
    data%Zbs = 0.0_dp
    data%Rbc(:, 0) = [10.0_dp, 0.0_dp]
    data%Rbc(:, 1) = [10.0_dp, 1.0_dp]
    data%Zbs(:, 1) = [0.0_dp, -1.0_dp]
    call make_spectre_coordinate_system(cs, data)

    call cs%from_cyl([12.0_dp, 0.0_dp, 0.0_dp], u, exterior_status)
    u = [0.0_dp, 0.0_dp, 0.0_dp]
    select type (cs)
    type is (spectre_coordinate_system_t)
        call cs%from_cyl_warm([10.5_dp, 0.0_dp, 0.0_dp], u, 1, singular_status)
    class default
        error stop 'unexpected chart type'
    end select
    print *, 'exterior and interior singular-seed statuses:', &
        exterior_status, singular_status
    call check_status('cold exterior', exterior_status, SPECTRE_OUTSIDE_DOMAIN)
    call check_status('singular interior seed', singular_status, SPECTRE_NOT_CONVERGED)

    select type (cs)
    type is (spectre_coordinate_system_t)
        saved = 0.0_dp
        u = saved
        call cs%from_cyl_warm([12.0_dp, 0.0_dp, 0.0_dp], u, 1, status)
        call check_status('warm exterior', status, SPECTRE_OUTSIDE_DOMAIN)
        if (any(u /= saved)) error stop 'failed warm inverse changed seed'
        call cs%from_cyl_warm([10.5_dp, 0.0_dp, 0.0_dp], u, 0, status)
        call check_status('invalid volume', status, SPECTRE_INVALID_VOLUME)
        call cs%from_cyl_warm([10.5_dp, 0.0_dp, 0.0_dp], u, 2, status)
        call check_status('volume above Mvol', status, SPECTRE_INVALID_VOLUME)

        call cs%from_cyl([10.5_dp, 0.4_dp, 0.0_dp], u, status)
        call check_status('cold interior', status, SPECTRE_SUCCESS)
        if (abs(u(1) - 0.5_dp) > 1.0e-9_dp) error stop 'wrong inverse radius'
        call cs%from_cyl([11.0_dp, 0.4_dp, 0.0_dp], u, status)
        call check_status('cold boundary', status, SPECTRE_SUCCESS)
        u = saved
        call cs%from_cyl_warm([11.0_dp, 0.0_dp, 0.0_dp], u, 1, status)
        call check_status('failed boundary seed', status, SPECTRE_DOMAIN_UNDETERMINED)
        status = cs%domain_status([11.0_dp + 1.0e-12_dp, 0.0_dp, 0.0_dp])
        call check_status('unresolved near boundary', status, &
                          SPECTRE_DOMAIN_UNDETERMINED)

        ! Elliptic cross-sections translated with phi, including Rbs/Zbc.
        cs%mn = 3
        deallocate (cs%im, cs%in, cs%Rbc, cs%Rbs, cs%Zbc, cs%Zbs)
        allocate (cs%im(3), cs%in(3))
        allocate (cs%Rbc(3, 0:1), cs%Rbs(3, 0:1))
        allocate (cs%Zbc(3, 0:1), cs%Zbs(3, 0:1))
        cs%im = [0, 1, 0]
        cs%in = [0, 0, 1]
        cs%Rbc = 0.0_dp
        cs%Rbs = 0.0_dp
        cs%Zbc = 0.0_dp
        cs%Zbs = 0.0_dp
        cs%Rbc(1, :) = 10.0_dp
        cs%Rbc(3, :) = 0.4_dp
        cs%Zbs(3, :) = -0.3_dp
        cs%Rbc(2, 1) = 1.5_dp
        cs%Rbs(2, 1) = 0.2_dp
        cs%Zbc(2, 1) = 0.1_dp
        cs%Zbs(2, 1) = -0.8_dp
        do k = 0, 7
            phi = real(k, dp)*acos(-1.0_dp)/4.0_dp
            center(1) = 10.0_dp + 0.4_dp*cos(phi)
            center(2) = 0.3_dp*sin(phi)
            do j = 1, 2
                radius = 0.5_dp
                status = SPECTRE_SUCCESS
                if (j == 2) then
                    radius = 1.5_dp
                    status = SPECTRE_OUTSIDE_DOMAIN
                end if
                do i = 0, 15
                    angle = real(i, dp)*acos(-1.0_dp)/8.0_dp
                    point(1) = center(1) + radius*(1.5_dp*cos(angle) &
                                                   + 0.2_dp*sin(angle))
                    point(2) = phi
                    point(3) = center(2) + radius*(0.1_dp*cos(angle) &
                                                   - 0.8_dp*sin(angle))
                    call check_status('shaped toroidal domain', &
                                      cs%domain_status(point), status)
                end do
            end do
        end do
    end select

    ! Two concentric volumes: failure in volume 1 is not global exterior.
    data%Mvol = 2
    deallocate (data%Rbc, data%Rbs, data%Zbc, data%Zbs)
    allocate (data%Rbc(2, 0:2), data%Zbs(2, 0:2))
    allocate (data%Rbs(2, 0:2), data%Zbc(2, 0:2))
    data%Rbc(1, :) = 10.0_dp
    data%Rbc(2, :) = [0.0_dp, 1.0_dp, 2.0_dp]
    data%Zbs = 0.0_dp
    data%Zbs(2, :) = [0.0_dp, -1.0_dp, -2.0_dp]
    data%Rbs = 0.0_dp
    data%Zbc = 0.0_dp
    call make_spectre_coordinate_system(cs, data)
    select type (cs)
    type is (spectre_coordinate_system_t)
        u = 0.0_dp
        call cs%from_cyl_warm([11.5_dp, 0.0_dp, 0.0_dp], u, 1, status)
        call check_status('wrong warm volume', status, SPECTRE_NOT_CONVERGED)
        call cs%from_cyl([11.5_dp, 0.0_dp, 0.0_dp], u, status)
        call check_status('cold finds other volume', status, SPECTRE_SUCCESS)
        call check_status('outermost volume domain', &
                          cs%domain_status([11.5_dp, 0.0_dp, 0.0_dp]), SPECTRE_SUCCESS)
    end select

    ! A star-shaped contour r(theta) = 1 + 0.1*cos(64*theta). A fixed
    ! 32-point polygon aliases the corrugation; the arc bound must refine it.
    data%Mvol = 1
    data%mn = 4
    deallocate (data%im, data%in, data%Rbc, data%Rbs, data%Zbc, data%Zbs)
    allocate (data%im(4), data%in(4))
    allocate (data%Rbc(4, 0:1), data%Rbs(4, 0:1))
    allocate (data%Zbc(4, 0:1), data%Zbs(4, 0:1))
    data%im = [0, 1, 63, 65]
    data%in = 0
    data%Rbc = 0.0_dp
    data%Rbs = 0.0_dp
    data%Zbc = 0.0_dp
    data%Zbs = 0.0_dp
    data%Rbc(1, :) = 10.0_dp
    data%Rbc(2:4, 1) = [1.0_dp, 0.05_dp, 0.05_dp]
    data%Zbs(2:4, 1) = [-1.0_dp, 0.05_dp, -0.05_dp]
    call make_spectre_coordinate_system(cs, data)
    select type (cs)
    type is (spectre_coordinate_system_t)
        call check_status('resolved corrugation inside', &
                          cs%domain_status([11.05_dp, 0.0_dp, 0.0_dp]), SPECTRE_SUCCESS)
        angle = acos(-1.0_dp)/64.0_dp
        point(1) = 10.0_dp + 0.95_dp*cos(angle)
        point(2) = 0.0_dp
        point(3) = -0.95_dp*sin(angle)
        call check_status('resolved corrugation outside', &
                          cs%domain_status(point), SPECTRE_OUTSIDE_DOMAIN)
    end select
contains
    subroutine check_status(label, actual, expected)
        character(len=*), intent(in) :: label
        integer, intent(in) :: actual, expected

        if (actual == expected) return
        print *, label, ': actual =', actual, ' expected =', expected
        error stop 'incorrect SPECTRE status'
    end subroutine check_status
end program test_spectre_inverse_status
