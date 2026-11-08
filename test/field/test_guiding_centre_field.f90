program test_guiding_centre_field
    ! Guiding-centre quantities of a coil field in the analytic toroidal chart
    ! x = ((R0 + rho cos t) cos z, (R0 + rho cos t) sin z, rho sin t): sqrt(g)
    ! against its closed form, d_j log|B| and (curl h)^j against central
    ! differences of log|B| and h_j along the chart coordinates.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neo_biotsavart, only: coils_t, coils_init, compute_field_jacobian
    use neo_guiding_centre_field, only: guiding_centre_quantities, GC_BMOD, &
                                   GC_SQRTG, GC_BDER, GC_HCOV, GC_HCTR, GC_HCURL, GC_COUNT
    implicit none

    real(dp), parameter :: pi = 3.141592653589793238462643383279502884_dp
    real(dp), parameter :: r0 = 300.0_dp, h = 2.0e-4_dp
    integer, parameter :: npts = 30
    type(coils_t) :: coils
    real(dp) :: u(3, npts), q(GC_COUNT), qp(GC_COUNT), qm(GC_COUNT), qpp(GC_COUNT), &
                qmm(GC_COUNT), dh(3, 3), &
                dlog(3), curl(3), err_g, err_bder, err_curl, err_dual
    integer :: k, j

    call helical_coil(coils)
    call random_seed(put=[(20261108 + k, k=1, 64)])
    call random_number(u)
    u(1, :) = 20.0_dp + 60.0_dp*u(1, :)
    u(2, :) = 2*pi*u(2, :)
    u(3, :) = 2*pi*u(3, :)
    err_g = 0
    err_bder = 0
    err_curl = 0
    err_dual = 0
    do k = 1, npts
        call quantities(u(:, k), q)
        err_g = max(err_g, abs(abs(q(GC_SQRTG)) &
                               - u(1, k)*(r0 + u(1, k)*cos(u(2, k))))/abs(q(GC_SQRTG)))
        err_dual = max(err_dual, abs(dot_product(q(GC_HCOV:GC_HCOV + 2), &
                                                 q(GC_HCTR:GC_HCTR + 2)) - 1))
        do j = 1, 3
            call quantities(u(:, k) + h*unit(j), qp)
            call quantities(u(:, k) - h*unit(j), qm)
            call quantities(u(:, k) + 2*h*unit(j), qpp)
            call quantities(u(:, k) - 2*h*unit(j), qmm)
            ! fourth-order central differences
            dlog(j) = (8*(log(qp(GC_BMOD)) - log(qm(GC_BMOD))) &
                       - (log(qpp(GC_BMOD)) - log(qmm(GC_BMOD))))/(12*h)
            dh(j, :) = (8*(qp(GC_HCOV:GC_HCOV + 2) - qm(GC_HCOV:GC_HCOV + 2)) &
                        - (qpp(GC_HCOV:GC_HCOV + 2) - qmm(GC_HCOV:GC_HCOV + 2)))/(12*h)
        end do
        curl = [dh(2, 3) - dh(3, 2), dh(3, 1) - dh(1, 3), dh(1, 2) - dh(2, 1)] &
               /q(GC_SQRTG)
        err_bder = max(err_bder, maxval(abs(q(GC_BDER:GC_BDER + 2) - dlog)) &
                       /maxval(abs(dlog)))
        err_curl = max(err_curl, maxval(abs(q(GC_HCURL:GC_HCURL + 2) - curl)) &
                       /maxval(abs(curl)))
    end do
    call check('sqrt(g) equals rho (R0 + rho cos t)', err_g, 1.0e-14_dp)
    call check('h_j h^j = 1', err_dual, 1.0e-14_dp)
    call check('d_j log|B| equals central differences', err_bder, 1.0e-7_dp)
    call check('(curl h)^j equals central differences', err_curl, 1.0e-7_dp)

contains

    subroutine quantities(u, q)
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: q(GC_COUNT)
        real(dp) :: x(3), e(3, 3), a(3), da(3, 3), b(3), db(3, 3), rr

        rr = r0 + u(1)*cos(u(2))
        x = [rr*cos(u(3)), rr*sin(u(3)), u(1)*sin(u(2))]
        e(:, 1) = [cos(u(2))*cos(u(3)), cos(u(2))*sin(u(3)), sin(u(2))]
        e(:, 2) = [-u(1)*sin(u(2))*cos(u(3)), -u(1)*sin(u(2))*sin(u(3)), &
                   u(1)*cos(u(2))]
        e(:, 3) = [-rr*sin(u(3)), rr*cos(u(3)), 0.0_dp]
        call compute_field_jacobian(coils, x, a, da, b, db)
        call guiding_centre_quantities(e, b, db, q)
    end subroutine quantities

    pure function unit(j) result(v)
        integer, intent(in) :: j
        real(dp) :: v(3)

        v = 0
        v(j) = 1
    end function unit

    subroutine helical_coil(coils)
        ! A closed helical winding around the torus plus a vertical return
        ! wire on the axis, so that B has toroidal and poloidal parts.
        type(coils_t), intent(out) :: coils
        integer, parameter :: n = 2000
        real(dp) :: t(n + 1), xs(n + 3), ys(n + 3), zs(n + 3), cur(n + 3)
        integer :: i

        t = [(2*pi*i/n, i=0, n)]
        xs(1:n + 1) = (r0 + 120.0_dp*cos(5*t))*cos(t)
        ys(1:n + 1) = (r0 + 120.0_dp*cos(5*t))*sin(t)
        zs(1:n + 1) = 120.0_dp*sin(5*t)
        cur(1:n + 1) = 1.0e7_dp
        cur(n + 1) = 0
        xs(n + 2:n + 3) = 0
        ys(n + 2:n + 3) = 0
        zs(n + 2:n + 3) = [-3.0e3_dp, 3.0e3_dp]
        cur(n + 2:n + 3) = [5.0e7_dp, 0.0_dp]
        call coils_init(xs, ys, zs, cur, coils)
    end subroutine helical_coil

    subroutine check(label, err, tol)
        character(*), intent(in) :: label
        real(dp), intent(in) :: err, tol

        print '(a, es10.2, a, es9.2)', label//': ', err, ' <= ', tol
        if (.not. err <= tol) error stop 'test_guiding_centre_field failed'
    end subroutine check

end program test_guiding_centre_field
