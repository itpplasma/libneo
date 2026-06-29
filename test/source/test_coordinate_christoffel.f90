module coordinate_christoffel_charts
    !> Synthetic coordinate_system_t charts with known closed-form Christoffel
    !> symbols, used to exercise the christoffel binding and the shared
    !> finite-difference helper christoffel_fd:
    !>   - cart_chart_t: flat Cartesian chart, identity metric (all symbols 0),
    !>   - cyl_chart_t: cylindrical chart u = (R, phi, Z) with g = diag(1,R^2,1),
    !>     where only Gamma^R_{phiphi} = -R and Gamma^phi_{R phi} = 1/R differ
    !>     from zero.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_coordinates_base, only: coordinate_system_t
    implicit none
    private
    public :: cart_chart_t, cyl_chart_t
    public :: test_flat_cartesian, test_cylindrical, test_symmetry_and_compat

    type, extends(coordinate_system_t) :: cart_chart_t
    contains
        procedure :: evaluate_cart => cart_eval
        procedure :: evaluate_cyl => cart_eval
        procedure :: covariant_basis => cart_basis
        procedure :: metric_tensor => cart_metric
        procedure :: christoffel => cart_christoffel
        procedure :: from_cyl => cart_from_cyl
    end type cart_chart_t

    type, extends(coordinate_system_t) :: cyl_chart_t
    contains
        procedure :: evaluate_cart => cyl_eval
        procedure :: evaluate_cyl => cyl_eval
        procedure :: covariant_basis => cyl_basis
        procedure :: metric_tensor => cyl_metric
        procedure :: christoffel => cyl_christoffel
        procedure :: from_cyl => cyl_from_cyl
    end type cyl_chart_t

contains

    subroutine test_flat_cartesian(cs)
        type(cart_chart_t), intent(in) :: cs
        real(dp) :: u(3), Gamma(3, 3, 3), err
        integer :: i, m, n

        u = [0.3_dp, 1.1_dp, -0.4_dp]
        call cs%christoffel(u, Gamma)
        err = 0.0_dp
        do i = 1, 3
            do m = 1, 3
                do n = 1, 3
                    err = max(err, abs(Gamma(i, m, n)))
                end do
            end do
        end do
        print '(a,es12.4)', '  flat Cartesian max |Gamma| = ', err
        if (err > 1.0e-8_dp) error stop 'flat Cartesian Christoffel not zero'
    end subroutine test_flat_cartesian

    subroutine test_cylindrical(cs)
        type(cyl_chart_t), intent(in) :: cs
        real(dp) :: u(3), Gamma(3, 3, 3), R, err
        integer :: i, m, n

        u = [1.7_dp, 0.6_dp, 0.9_dp] ! (R, phi, Z)
        R = u(1)
        call cs%christoffel(u, Gamma)

        err = abs(Gamma(1, 2, 2) - (-R))
        print '(a,es12.4)', '  cyl Gamma^R_phiphi err = ', err
        if (err > 1.0e-6_dp) error stop 'cyl Gamma^R_phiphi wrong'

        err = max(abs(Gamma(2, 1, 2) - 1.0_dp/R), abs(Gamma(2, 2, 1) - 1.0_dp/R))
        print '(a,es12.4)', '  cyl Gamma^phi_Rphi err = ', err
        if (err > 1.0e-6_dp) error stop 'cyl Gamma^phi_Rphi wrong'

        ! All other components must be (near) zero.
        do i = 1, 3
            do m = 1, 3
                do n = 1, 3
                    if (i == 1 .and. m == 2 .and. n == 2) cycle
                    if (i == 2 .and. m == 1 .and. n == 2) cycle
                    if (i == 2 .and. m == 2 .and. n == 1) cycle
                    if (abs(Gamma(i, m, n)) > 1.0e-6_dp) then
                        print '(a,3i2,es12.4)', '  spurious Gamma ', i, m, n, &
                            Gamma(i, m, n)
                        error stop 'cyl Christoffel has spurious nonzero entry'
                    end if
                end do
            end do
        end do
    end subroutine test_cylindrical

    subroutine test_symmetry_and_compat(cs)
        type(cyl_chart_t), intent(in) :: cs
        real(dp) :: u(3), Gamma(3, 3, 3)
        real(dp) :: g(3, 3), ginv(3, 3), sqrtg
        real(dp) :: gp(3, 3), gm(3, 3), dummy(3, 3), ddet
        real(dp) :: dg(3, 3, 3), up(3), um(3), h, nabla, err
        integer :: i, j, k, l, m, n

        u = [1.3_dp, 0.4_dp, 0.7_dp]
        h = 1.0e-5_dp
        call cs%christoffel(u, Gamma)

        err = 0.0_dp
        do i = 1, 3
            do m = 1, 3
                do n = 1, 3
                    err = max(err, abs(Gamma(i, m, n) - Gamma(i, n, m)))
                end do
            end do
        end do
        print '(a,es12.4)', '  symmetry max |Gamma^i_mn - Gamma^i_nm| = ', err
        if (err > 1.0e-12_dp) error stop 'coordinate christoffel not symmetric'

        call cs%metric_tensor(u, g, ginv, sqrtg)
        do k = 1, 3
            up = u
            um = u
            up(k) = u(k) + h
            um(k) = u(k) - h
            call cs%metric_tensor(up, gp, dummy, ddet)
            call cs%metric_tensor(um, gm, dummy, ddet)
            dg(:, :, k) = (gp - gm)/(2.0_dp*h)
        end do

        err = 0.0_dp
        do k = 1, 3
            do i = 1, 3
                do j = 1, 3
                    nabla = dg(i, j, k)
                    do l = 1, 3
                        nabla = nabla - Gamma(l, k, i)*g(l, j) - Gamma(l, k, j)*g(i, l)
                    end do
                    err = max(err, abs(nabla))
                end do
            end do
        end do
        print '(a,es12.4)', '  metric compatibility max |nabla_k g_ij| = ', err
        if (err > 1.0e-5_dp) error stop 'coordinate christoffel metric compatibility failed'
    end subroutine test_symmetry_and_compat

    ! ----- flat Cartesian chart: identity metric -----

    subroutine cart_eval(self, u, x)
        class(cart_chart_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: x(3)

        associate (dummy => self)
        end associate

        x = u
    end subroutine cart_eval

    subroutine cart_basis(self, u, e_cov)
        class(cart_chart_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: e_cov(3, 3)
        integer :: i

        associate (dummy => self)
        end associate
        associate (dummy => u)
        end associate

        e_cov = 0.0_dp
        do i = 1, 3
            e_cov(i, i) = 1.0_dp
        end do
    end subroutine cart_basis

    subroutine cart_metric(self, u, g, ginv, sqrtg)
        class(cart_chart_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: g(3, 3), ginv(3, 3), sqrtg
        integer :: i

        associate (dummy => self)
        end associate
        associate (dummy => u)
        end associate

        g = 0.0_dp
        ginv = 0.0_dp
        do i = 1, 3
            g(i, i) = 1.0_dp
            ginv(i, i) = 1.0_dp
        end do
        sqrtg = 1.0_dp
    end subroutine cart_metric

    subroutine cart_christoffel(self, u, Gamma)
        class(cart_chart_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: Gamma(3, 3, 3)

        call self%christoffel_fd(u, Gamma)
    end subroutine cart_christoffel

    subroutine cart_from_cyl(self, xcyl, u, ierr)
        class(cart_chart_t), intent(in) :: self
        real(dp), intent(in) :: xcyl(3)
        real(dp), intent(out) :: u(3)
        integer, intent(out) :: ierr

        associate (dummy => self)
        end associate

        u = xcyl
        ierr = 0
    end subroutine cart_from_cyl

    ! ----- cylindrical chart u = (R, phi, Z): g = diag(1, R^2, 1) -----

    subroutine cyl_eval(self, u, x)
        class(cyl_chart_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: x(3)

        associate (dummy => self)
        end associate

        x(1) = u(1)*cos(u(2))
        x(2) = u(1)*sin(u(2))
        x(3) = u(3)
    end subroutine cyl_eval

    subroutine cyl_basis(self, u, e_cov)
        class(cyl_chart_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: e_cov(3, 3)
        real(dp) :: R, phi

        associate (dummy => self)
        end associate

        R = u(1)
        phi = u(2)
        e_cov = 0.0_dp
        e_cov(1, 1) = cos(phi)
        e_cov(2, 1) = sin(phi)
        e_cov(1, 2) = -R*sin(phi)
        e_cov(2, 2) = R*cos(phi)
        e_cov(3, 3) = 1.0_dp
    end subroutine cyl_basis

    subroutine cyl_metric(self, u, g, ginv, sqrtg)
        class(cyl_chart_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: g(3, 3), ginv(3, 3), sqrtg
        real(dp) :: R

        associate (dummy => self)
        end associate

        R = u(1)
        g = 0.0_dp
        ginv = 0.0_dp
        g(1, 1) = 1.0_dp
        g(2, 2) = R*R
        g(3, 3) = 1.0_dp
        ginv(1, 1) = 1.0_dp
        ginv(2, 2) = 1.0_dp/(R*R)
        ginv(3, 3) = 1.0_dp
        sqrtg = R
    end subroutine cyl_metric

    subroutine cyl_christoffel(self, u, Gamma)
        class(cyl_chart_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: Gamma(3, 3, 3)

        call self%christoffel_fd(u, Gamma)
    end subroutine cyl_christoffel

    subroutine cyl_from_cyl(self, xcyl, u, ierr)
        class(cyl_chart_t), intent(in) :: self
        real(dp), intent(in) :: xcyl(3)
        real(dp), intent(out) :: u(3)
        integer, intent(out) :: ierr

        associate (dummy => self)
        end associate

        u = xcyl
        ierr = 0
    end subroutine cyl_from_cyl

end module coordinate_christoffel_charts

program test_coordinate_christoffel
    use coordinate_christoffel_charts, only: cart_chart_t, cyl_chart_t, &
        test_flat_cartesian, test_cylindrical, test_symmetry_and_compat
    implicit none

    type(cart_chart_t) :: cart
    type(cyl_chart_t) :: cyl

    call test_flat_cartesian(cart)
    call test_cylindrical(cyl)
    call test_symmetry_and_compat(cyl)

    print *, 'coordinate_system christoffel: all checks passed'
end program test_coordinate_christoffel
