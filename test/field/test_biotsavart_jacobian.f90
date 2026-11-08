program test_biotsavart_jacobian
    use neo_biotsavart, only: coils_t, coils_init, compute_vector_potential, &
                              compute_magnetic_field, compute_field_jacobian
    implicit none

    integer, parameter :: dp = kind(1.0d0)
    real(dp), parameter :: pi = 3.141592653589793238462643383279502884197d0
    integer, parameter :: n_points = 40
    type(coils_t) :: coils
    real(dp) :: x(3), A(3), B(3), dA(3, 3), dB(3, 3), points(3, n_points)
    real(dp) :: fd(3, 3), scale_b, scale_db, curl_a(3)
    real(dp) :: err_value, err_fd, err_div, err_curl_b, err_curl_a
    integer :: k

    call tilted_polygon_loop(coils)
    call random_seed(put=[(20261108 + k, k=1, 64)])
    call random_number(points)
    points = 300.0_dp*(points - 0.5_dp)
    err_value = 0
    err_fd = 0
    err_div = 0
    err_curl_b = 0
    err_curl_a = 0
    do k = 1, n_points
        x = points(:, k)
        call compute_field_jacobian(coils, x, A, dA, B, dB)
        scale_b = norm2(B)
        scale_db = maxval(abs(dB))
        err_value = max(err_value, &
                        norm2(A - compute_vector_potential(coils, x))/norm2(A), &
                        norm2(B - compute_magnetic_field(coils, x))/scale_b)
        call central_difference(coils, x, fd)
        err_fd = max(err_fd, maxval(abs(dB - fd))/scale_db)
        err_div = max(err_div, abs(dB(1, 1) + dB(2, 2) + dB(3, 3))/scale_db)
        err_curl_b = max(err_curl_b, maxval(abs(curl(dB)))/scale_db)
        curl_a = curl(dA)
        err_curl_a = max(err_curl_a, norm2(curl_a - B)/scale_b)
    end do
    call check('A and B equal the existing routines', err_value, 1.0e-13_dp)
    call check('dB equals central differences of B', err_fd, 1.0e-8_dp)
    call check('div B vanishes', err_div, 1.0e-12_dp)
    call check('curl B vanishes off the wire', err_curl_b, 1.0e-12_dp)
    call check('curl A equals B', err_curl_a, 1.0e-12_dp)

contains

    subroutine tilted_polygon_loop(coils)
        type(coils_t), intent(out) :: coils
        integer, parameter :: n = 64
        real(dp) :: t(n + 1), xs(n + 1), ys(n + 1), zs(n + 1), cur(n + 1)
        integer :: i

        t = [(2*pi*i/n, i=0, n)]
        xs = 100.0_dp*cos(t) + 20.0_dp
        ys = 80.0_dp*sin(t)
        zs = 30.0_dp*cos(t) + 15.0_dp*sin(2*t)
        cur = 1.0e6_dp
        cur(n + 1) = 0
        call coils_init(xs, ys, zs, cur, coils)
    end subroutine tilted_polygon_loop

    subroutine central_difference(coils, x, fd)
        type(coils_t), intent(in) :: coils
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: fd(3, 3)
        real(dp), parameter :: h = 1.0e-4_dp
        real(dp) :: e(3)
        integer :: j

        do j = 1, 3
            e = 0
            e(j) = h
            fd(:, j) = (compute_magnetic_field(coils, x + e) &
                        - compute_magnetic_field(coils, x - e))/(2*h)
        end do
    end subroutine central_difference

    pure function curl(d) result(c)
        real(dp), intent(in) :: d(3, 3)
        real(dp) :: c(3)

        c = [d(3, 2) - d(2, 3), d(1, 3) - d(3, 1), d(2, 1) - d(1, 2)]
    end function curl

    subroutine check(label, err, tol)
        character(*), intent(in) :: label
        real(dp), intent(in) :: err, tol

        print '(a, es10.2, a, es9.2)', label//': ', err, ' <= ', tol
        if (.not. err <= tol) error stop 'test_biotsavart_jacobian failed'
    end subroutine check

end program test_biotsavart_jacobian
