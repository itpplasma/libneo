program test_vector_conversion
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_coordinates, only: coordinate_system_t, &
                                  make_chartmap_coordinate_system, &
                                  chartmap_coordinate_system_t
    use math_constants, only: TWOPI
    implicit none

    integer :: nerrors
    class(coordinate_system_t), allocatable :: cs
    character(len=*), parameter :: chartmap_file = "chartmap.nc"

    nerrors = 0

    call make_chartmap_coordinate_system(cs, chartmap_file)
    if (.not. allocated(cs)) error stop "Failed to load chartmap"

    select type (ccs => cs)
    type is (chartmap_coordinate_system_t)
        call test_ctr_to_cart_basis_vectors(ccs, nerrors)
        call test_ctr_to_cart_linearity(ccs, nerrors)
        call test_cov_to_cart_norm_preservation(ccs, nerrors)
        call test_cov_ctr_cart_consistency(ccs, nerrors)
        call test_cart_ctr_cart_roundtrip(ccs, nerrors)
        call test_cart_cov_cart_roundtrip(ccs, nerrors)
    class default
        error stop "Unexpected coordinate system type"
    end select

    if (nerrors > 0) then
        print *, "FAILED:", nerrors, "error(s) in vector conversion tests"
        error stop 1
    end if
    print *, "All vector conversion tests passed!"

contains

    subroutine test_ctr_to_cart_basis_vectors(cs, nerrors)
        class(coordinate_system_t), intent(in) :: cs
        integer, intent(inout) :: nerrors

        real(dp) :: u(3), e_cov(3, 3)
        real(dp) :: v_ctr(3), v_cart(3)
        real(dp) :: tol, err
        integer :: k

        tol = 1.0e-12_dp
        u = [0.4_dp, 1.2_dp, 0.8_dp]

        call cs%covariant_basis(u, e_cov)

        do k = 1, 3
            v_ctr = 0.0_dp
            v_ctr(k) = 1.0_dp
            call cs%ctr_to_cart(u, v_ctr, v_cart)
            err = maxval(abs(v_cart - e_cov(:, k)))
            if (err > tol) then
                print *, "  FAIL: ctr_to_cart unit vector", k, "err=", err
                nerrors = nerrors + 1
            end if
        end do

        if (nerrors == 0) then
            print *, "  PASS: ctr_to_cart maps unit vectors to basis vectors"
        end if
    end subroutine test_ctr_to_cart_basis_vectors

    subroutine test_ctr_to_cart_linearity(cs, nerrors)
        class(coordinate_system_t), intent(in) :: cs
        integer, intent(inout) :: nerrors

        real(dp) :: u(3)
        real(dp) :: v1_ctr(3), v2_ctr(3), v_combined_ctr(3)
        real(dp) :: v1_cart(3), v2_cart(3), v_combined_cart(3)
        real(dp) :: v_expected(3)
        real(dp) :: a, b, tol, err

        tol = 1.0e-12_dp
        u = [0.5_dp, 2.5_dp, 1.0_dp]
        v1_ctr = [1.3_dp, -0.7_dp, 2.1_dp]
        v2_ctr = [-0.5_dp, 1.8_dp, 0.3_dp]
        a = 2.5_dp
        b = -1.3_dp

        call cs%ctr_to_cart(u, v1_ctr, v1_cart)
        call cs%ctr_to_cart(u, v2_ctr, v2_cart)

        v_combined_ctr = a * v1_ctr + b * v2_ctr
        call cs%ctr_to_cart(u, v_combined_ctr, v_combined_cart)

        v_expected = a * v1_cart + b * v2_cart
        err = maxval(abs(v_combined_cart - v_expected))

        if (err > tol) then
            print *, "  FAIL: ctr_to_cart linearity err=", err
            nerrors = nerrors + 1
        else
            print *, "  PASS: ctr_to_cart is linear"
        end if
    end subroutine test_ctr_to_cart_linearity

    subroutine test_cov_to_cart_norm_preservation(cs, nerrors)
        class(coordinate_system_t), intent(in) :: cs
        integer, intent(inout) :: nerrors

        real(dp) :: u(3)
        real(dp) :: g(3, 3), ginv(3, 3), sqrtg
        real(dp) :: v_cov(3), v_cart(3)
        real(dp) :: norm_sq_curv, norm_sq_cart
        real(dp) :: tol, rel_err
        integer :: i, j, n_tests, n_passed
        real(dp) :: test_points(3, 4)

        tol = 1.0e-10_dp
        n_tests = 0
        n_passed = 0

        test_points(:, 1) = [0.2_dp, 0.5_dp, 0.3_dp]
        test_points(:, 2) = [0.7_dp, 3.0_dp, 1.5_dp]
        test_points(:, 3) = [0.5_dp, 5.0_dp, 0.1_dp]
        test_points(:, 4) = [0.9_dp, 1.0_dp, 2.0_dp]

        do i = 1, 4
            u = test_points(:, i)
            call cs%metric_tensor(u, g, ginv, sqrtg)

            v_cov = [1.5_dp, -0.8_dp, 2.3_dp]

            norm_sq_curv = 0.0_dp
            do j = 1, 3
                norm_sq_curv = norm_sq_curv + v_cov(j) * sum(ginv(j, :) * v_cov)
            end do

            call cs%cov_to_cart(u, v_cov, v_cart)
            norm_sq_cart = sum(v_cart**2)

            n_tests = n_tests + 1
            rel_err = abs(norm_sq_cart - norm_sq_curv) / max(norm_sq_curv, 1.0e-15_dp)

            if (rel_err < tol) then
                n_passed = n_passed + 1
            else
                print *, "  FAIL: norm mismatch at u=", u
                print *, "    |v|²_curv=", norm_sq_curv, " |v|²_cart=", norm_sq_cart
            end if
        end do

        if (n_passed == n_tests) then
            print *, "  PASS: cov_to_cart preserves vector norm (", n_tests, " points)"
        else
            print *, "  FAIL: norm preservation failed", n_tests - n_passed, "/", n_tests
            nerrors = nerrors + 1
        end if
    end subroutine test_cov_to_cart_norm_preservation

    subroutine test_cov_ctr_cart_consistency(cs, nerrors)
        class(coordinate_system_t), intent(in) :: cs
        integer, intent(inout) :: nerrors

        real(dp) :: u(3)
        real(dp) :: g(3, 3), ginv(3, 3), sqrtg
        real(dp) :: v_ctr(3), v_cov(3)
        real(dp) :: v_cart_from_ctr(3), v_cart_from_cov(3)
        real(dp) :: tol, err
        integer :: i, n_tests, n_passed
        real(dp) :: test_points(3, 3)

        tol = 1.0e-11_dp
        n_tests = 0
        n_passed = 0

        test_points(:, 1) = [0.3_dp, 1.0_dp, 0.5_dp]
        test_points(:, 2) = [0.6_dp, 4.0_dp, 1.2_dp]
        test_points(:, 3) = [0.8_dp, 2.5_dp, 0.8_dp]

        do i = 1, 3
            u = test_points(:, i)

            v_ctr = [0.7_dp, -1.2_dp, 0.9_dp]

            call cs%ctr_to_cart(u, v_ctr, v_cart_from_ctr)

            call cs%metric_tensor(u, g, ginv, sqrtg)
            v_cov = matmul(g, v_ctr)

            call cs%cov_to_cart(u, v_cov, v_cart_from_cov)

            n_tests = n_tests + 1
            err = maxval(abs(v_cart_from_ctr - v_cart_from_cov))

            if (err < tol) then
                n_passed = n_passed + 1
            else
                print *, "  FAIL: cov/ctr inconsistency at u=", u, " err=", err
            end if
        end do

        if (n_passed == n_tests) then
            print *, "  PASS: cov_to_cart and ctr_to_cart are consistent (", n_tests, " points)"
        else
            print *, "  FAIL: cov/ctr consistency failed", n_tests - n_passed, "/", n_tests
            nerrors = nerrors + 1
        end if
    end subroutine test_cov_ctr_cart_consistency

    subroutine test_cart_ctr_cart_roundtrip(cs, nerrors)
        class(coordinate_system_t), intent(in) :: cs
        integer, intent(inout) :: nerrors

        real(dp) :: u(3), e_cov(3, 3)
        real(dp) :: v_cart(3), v_ctr(3), v_cart_back(3)
        real(dp) :: tol, err
        integer :: i, n_passed
        real(dp) :: test_points(3, 4)
        real(dp) :: test_vectors(3, 3)
        integer :: iv

        tol = 1.0e-10_dp
        n_passed = 0

        test_points(:, 1) = [0.3_dp, 0.8_dp, 0.4_dp]
        test_points(:, 2) = [0.6_dp, 2.5_dp, 1.0_dp]
        test_points(:, 3) = [0.5_dp, 4.0_dp, 0.2_dp]
        test_points(:, 4) = [0.8_dp, 1.5_dp, 1.8_dp]

        test_vectors(:, 1) = [1.0_dp, 0.0_dp, 0.0_dp]
        test_vectors(:, 2) = [0.5_dp, -0.3_dp, 0.8_dp]
        test_vectors(:, 3) = [-0.2_dp, 1.5_dp, 0.7_dp]

        do i = 1, 4
            u = test_points(:, i)
            call cs%covariant_basis(u, e_cov)

            do iv = 1, 3
                v_cart = test_vectors(:, iv)
                call solve_3x3(e_cov, v_cart, v_ctr)
                call cs%ctr_to_cart(u, v_ctr, v_cart_back)

                err = maxval(abs(v_cart_back - v_cart))
                if (err < tol) then
                    n_passed = n_passed + 1
                else
                    print *, "  FAIL: cart->ctr->cart roundtrip at u=", u, " err=", err
                end if
            end do
        end do

        if (n_passed == 12) then
            print *, "  PASS: cart->ctr->cart roundtrip (12 tests)"
        else
            print *, "  FAIL: cart->ctr->cart roundtrip", 12 - n_passed, "/12"
            nerrors = nerrors + 1
        end if
    end subroutine test_cart_ctr_cart_roundtrip

    subroutine test_cart_cov_cart_roundtrip(cs, nerrors)
        class(coordinate_system_t), intent(in) :: cs
        integer, intent(inout) :: nerrors

        real(dp) :: u(3), e_cov(3, 3)
        real(dp) :: g(3, 3), ginv(3, 3), sqrtg
        real(dp) :: v_cart(3), v_ctr(3), v_cov(3), v_cart_back(3)
        real(dp) :: tol, err
        integer :: i, n_passed
        real(dp) :: test_points(3, 4)
        real(dp) :: test_vectors(3, 3)
        integer :: iv

        tol = 1.0e-10_dp
        n_passed = 0

        test_points(:, 1) = [0.25_dp, 1.0_dp, 0.5_dp]
        test_points(:, 2) = [0.55_dp, 3.0_dp, 0.8_dp]
        test_points(:, 3) = [0.45_dp, 5.5_dp, 0.3_dp]
        test_points(:, 4) = [0.75_dp, 2.0_dp, 1.5_dp]

        test_vectors(:, 1) = [0.0_dp, 1.0_dp, 0.0_dp]
        test_vectors(:, 2) = [0.8_dp, 0.2_dp, -0.5_dp]
        test_vectors(:, 3) = [-0.4_dp, 0.9_dp, 1.2_dp]

        do i = 1, 4
            u = test_points(:, i)
            call cs%covariant_basis(u, e_cov)
            call cs%metric_tensor(u, g, ginv, sqrtg)

            do iv = 1, 3
                v_cart = test_vectors(:, iv)
                call solve_3x3(e_cov, v_cart, v_ctr)
                v_cov = matmul(g, v_ctr)
                call cs%cov_to_cart(u, v_cov, v_cart_back)

                err = maxval(abs(v_cart_back - v_cart))
                if (err < tol) then
                    n_passed = n_passed + 1
                else
                    print *, "  FAIL: cart->cov->cart roundtrip at u=", u, " err=", err
                end if
            end do
        end do

        if (n_passed == 12) then
            print *, "  PASS: cart->cov->cart roundtrip (12 tests)"
        else
            print *, "  FAIL: cart->cov->cart roundtrip", 12 - n_passed, "/12"
            nerrors = nerrors + 1
        end if
    end subroutine test_cart_cov_cart_roundtrip

    subroutine solve_3x3(A, b, x)
        real(dp), intent(in) :: A(3, 3), b(3)
        real(dp), intent(out) :: x(3)

        real(dp) :: det, Ainv(3, 3)

        det = A(1,1)*(A(2,2)*A(3,3) - A(2,3)*A(3,2)) &
            - A(1,2)*(A(2,1)*A(3,3) - A(2,3)*A(3,1)) &
            + A(1,3)*(A(2,1)*A(3,2) - A(2,2)*A(3,1))

        Ainv(1,1) = (A(2,2)*A(3,3) - A(2,3)*A(3,2))/det
        Ainv(1,2) = (A(1,3)*A(3,2) - A(1,2)*A(3,3))/det
        Ainv(1,3) = (A(1,2)*A(2,3) - A(1,3)*A(2,2))/det
        Ainv(2,1) = (A(2,3)*A(3,1) - A(2,1)*A(3,3))/det
        Ainv(2,2) = (A(1,1)*A(3,3) - A(1,3)*A(3,1))/det
        Ainv(2,3) = (A(1,3)*A(2,1) - A(1,1)*A(2,3))/det
        Ainv(3,1) = (A(2,1)*A(3,2) - A(2,2)*A(3,1))/det
        Ainv(3,2) = (A(1,2)*A(3,1) - A(1,1)*A(3,2))/det
        Ainv(3,3) = (A(1,1)*A(2,2) - A(1,2)*A(2,1))/det

        x = matmul(Ainv, b)
    end subroutine solve_3x3

end program test_vector_conversion
