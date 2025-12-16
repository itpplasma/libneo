program test_chartmap_basis_debug
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_coordinates, only: coordinate_system_t, &
                                  make_chartmap_coordinate_system, &
                                  chartmap_coordinate_system_t
    implicit none

    class(coordinate_system_t), allocatable :: cs
    real(dp) :: u(3), e_cov(3, 3), e_cov_fd(3, 3)
    real(dp) :: u_plus(3), u_minus(3), x_plus(3), x_minus(3)
    real(dp) :: rel_err, ratio
    real(dp), parameter :: h = 1.0e-6_dp
    integer :: k, i, j
    character(len=256) :: filename

    call get_command_argument(1, filename)
    if (len_trim(filename) == 0) then
        print *, "Usage: test_chartmap_basis_debug <chartmap.nc>"
        error stop 1
    end if

    print *, "Loading chartmap: ", trim(filename)
    call make_chartmap_coordinate_system(cs, trim(filename))

    u = [0.5_dp, 1.0_dp, 0.5_dp]  ! Use phi=0.5 to avoid periodicity edge
    print *, "Test point u = ", u

    select type (ccs => cs)
    type is (chartmap_coordinate_system_t)
        call ccs%evaluate_cart(u, x_plus)
        print *, "Cartesian position x = ", x_plus
        print *, "R = ", sqrt(x_plus(1)**2 + x_plus(2)**2)

        call ccs%covariant_basis(u, e_cov)

        do k = 1, 3
            u_plus = u
            u_minus = u
            u_plus(k) = u(k) + h
            u_minus(k) = u(k) - h
            call ccs%evaluate_cart(u_plus, x_plus)
            call ccs%evaluate_cart(u_minus, x_minus)
            e_cov_fd(:, k) = (x_plus - x_minus) / (2.0_dp * h)
            if (k == 3) then
                print *, "FD for phi derivative:"
                print *, "  u_plus = ", u_plus
                print *, "  u_minus = ", u_minus
                print *, "  x_plus = ", x_plus
                print *, "  x_minus = ", x_minus
                print *, "  x_diff = ", x_plus - x_minus
                print *, "  2h = ", 2.0_dp * h
            end if
        end do

        print *, ""
        print *, "Covariant basis from covariant_basis():"
        print *, "  e_cov(:, 1) = dx/drho   = ", e_cov(:, 1)
        print *, "  e_cov(:, 2) = dx/dtheta = ", e_cov(:, 2)
        print *, "  e_cov(:, 3) = dx/dphi   = ", e_cov(:, 3)

        print *, ""
        print *, "Covariant basis from finite-difference:"
        print *, "  e_cov_fd(:, 1) = dx/drho   = ", e_cov_fd(:, 1)
        print *, "  e_cov_fd(:, 2) = dx/dtheta = ", e_cov_fd(:, 2)
        print *, "  e_cov_fd(:, 3) = dx/dphi   = ", e_cov_fd(:, 3)

        print *, ""
        print *, "Difference (e_cov - e_cov_fd):"
        print *, "  diff(:, 1) = ", e_cov(:, 1) - e_cov_fd(:, 1)
        print *, "  diff(:, 2) = ", e_cov(:, 2) - e_cov_fd(:, 2)
        print *, "  diff(:, 3) = ", e_cov(:, 3) - e_cov_fd(:, 3)

        print *, ""
        print *, "Ratio e_cov / e_cov_fd (element-wise):"
        do j = 1, 3
            do i = 1, 3
                if (abs(e_cov_fd(i, j)) > 1.0e-10_dp) then
                    ratio = e_cov(i, j) / e_cov_fd(i, j)
                    print *, "  e_cov(", i, ",", j, ") / e_cov_fd = ", ratio
                end if
            end do
        end do

        rel_err = maxval(abs(e_cov - e_cov_fd)) / &
                  max(maxval(abs(e_cov_fd)), 1.0e-10_dp)
        print *, ""
        print *, "Relative error: ", rel_err

        if (rel_err > 1.0e-3_dp) then
            print *, "FAIL: relative error too large"
            error stop 1
        else
            print *, "PASS: relative error acceptable"
        end if

    class default
        print *, "ERROR: not a chartmap coordinate system"
        error stop 1
    end select

end program test_chartmap_basis_debug
