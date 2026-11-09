module chartmap_field_test_field
    ! Smooth field with nfp = 3 dependence in its cylindrical components.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neo_field_base, only: field_t
    implicit none
    type, extends(field_t) :: wavy_field_t
    contains
        procedure :: compute_afield, compute_bfield, compute_abfield
    end type wavy_field_t
contains
    subroutine cylindrical(x, a, b)
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: a(3), b(3)
        real(dp) :: r, p, z
        r = hypot(x(1), x(2)); p = atan2(x(2), x(1)); z = x(3)
        a = [z*cos(3*p), 0.5_dp*r, 1.0e-2_dp*r*z*sin(3*p)]
        b = [1.0e-2_dp*z*cos(3*p), 2.0e4_dp + 10*r, 1.0e3_dp*sin(3*p + z/50)]
    end subroutine cylindrical
    subroutine compute_abfield(self, x, A, B)
        class(wavy_field_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: A(3), B(3)
        real(dp) :: a_c(3), b_c(3), c, s
        call cylindrical(x, a_c, b_c)
        c = x(1)/hypot(x(1), x(2)); s = x(2)/hypot(x(1), x(2))
        A = [c*a_c(1) - s*a_c(2), s*a_c(1) + c*a_c(2), a_c(3)]
        B = [c*b_c(1) - s*b_c(2), s*b_c(1) + c*b_c(2), b_c(3)]
    end subroutine compute_abfield
    subroutine compute_afield(self, x, A)
        class(wavy_field_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: A(3)
        real(dp) :: B(3)
        call self%compute_abfield(x, A, B)
    end subroutine compute_afield
    subroutine compute_bfield(self, x, B)
        class(wavy_field_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: B(3)
        real(dp) :: A(3)
        call self%compute_abfield(x, A, B)
    end subroutine compute_bfield
end module chartmap_field_test_field

program test_chartmap_field
    ! Field on a synthetic nfp = 3 stellarator chartmap (cylindrical angle):
    ! nodes hold the exact cylindrical components, the written file reads back
    ! identically and is detected by content, a plain chartmap is not, and the
    ! spline error at cell centres falls faster than h^3 (the chart map itself
    ! is a cubic spline, which limits the composite smoothness).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_chartmap_vmec_generator, only: write_chartmap_file
    use chartmap_field_test_field, only: wavy_field_t, cylindrical
    use libneo_chartmap_field, only: chartmap_field_t, chartmap_has_field, &
                                     read_chartmap_field, sample_chartmap_field, &
                                     write_chartmap_field
    implicit none
    type(wavy_field_t) :: field
    type(chartmap_field_t) :: f, g, fine
    real(dp) :: a(3), b(3), worst, x(3)
    integer :: i, j, k

    call stellarator_chart('chart_field_coarse.nc', 9, 16, 12)
    call stellarator_chart('chart_field_fine.nc', 17, 32, 24)
    if (chartmap_has_field('chart_field_coarse.nc')) error stop 'plain chart has field'

    call sample_chartmap_field(field, 'chart_field_coarse.nc', 'wavy', f)
    call write_chartmap_field('chart_field_coarse.nc', f, 'chart_field_out.nc')
    if (.not. chartmap_has_field('chart_field_out.nc')) error stop 'field not detected'
    call read_chartmap_field('chart_field_out.nc', g)
    if (any(g%a /= f%a) .or. any(g%b /= f%b)) error stop 'round trip'
    if (g%spline_error_b /= f%spline_error_b .or. g%source /= 'wavy') &
        error stop 'attributes'

    worst = 0
    do k = 1, size(f%zeta)
        do j = 1, size(f%theta)
            do i = 1, size(f%rho)
                x = position(f%rho(i), f%theta(j), f%zeta(k))
                call cylindrical(x, a, b)
                worst = max(worst, maxval(abs(f%a(:, i, j, k) - a)), &
                            maxval(abs(f%b(:, i, j, k) - b))/2.0e4_dp)
            end do
        end do
    end do
    if (worst > 1.0e-12_dp) error stop 'node values'

    call sample_chartmap_field(field, 'chart_field_fine.nc', 'wavy', fine)
    print '(a,2es10.3,a,2es10.3)', 'spline error B coarse, fine', f%spline_error_b, &
        fine%spline_error_b, '   A', f%spline_error_a, fine%spline_error_a
    if (.not. f%spline_error_b > 10*fine%spline_error_b) error stop 'B convergence'
    if (.not. f%spline_error_a > 10*fine%spline_error_a) error stop 'A convergence'
contains
    function position(rho, theta, zeta) result(x)
        real(dp), intent(in) :: rho, theta, zeta
        real(dp) :: x(3), r
        r = 180 + 45*rho*cos(theta) + 4*rho*cos(theta - 3*zeta)
        x = [r*cos(zeta), r*sin(zeta), 45*rho*sin(theta) + 3*rho*sin(3*zeta)]
    end function position

    subroutine stellarator_chart(path, nr, nt, nz)
        character(len=*), intent(in) :: path
        integer, intent(in) :: nr, nt, nz
        real(dp), parameter :: pi = acos(-1.0_dp)
        real(dp) :: rho(nr), theta(nt), zeta(nz), p(3)
        real(dp), allocatable :: x(:, :, :), y(:, :, :), z(:, :, :)
        integer :: i, j, k, ierr
        character(len=256) :: message
        rho = [(real(i - 1, dp)/(nr - 1), i=1, nr)]
        theta = [(2*pi*(j - 1)/nt, j=1, nt)]
        zeta = [(2*pi/3*(k - 1)/nz, k=1, nz)]
        allocate (x(nr, nt, nz), y(nr, nt, nz), z(nr, nt, nz))
        do k = 1, nz
            do j = 1, nt
                do i = 1, nr
                    p = position(rho(i), theta(j), zeta(k))
                    x(i, j, k) = p(1); y(i, j, k) = p(2); z(i, j, k) = p(3)
                end do
            end do
        end do
        call write_chartmap_file(path, rho, theta, zeta, x, y, z, 3, ierr, message)
        if (ierr /= 0) error stop 'write chart'
    end subroutine stellarator_chart
end program test_chartmap_field
