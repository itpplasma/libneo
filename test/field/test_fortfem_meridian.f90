program test_fortfem_meridian
    ! FortFEM builds and solves against libneo's Fortnum: the field of a
    ! localized toroidal source of one mode has the source's sign and decays
    ! away from it.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_meridian_bspline_n, only: meridian_bspline_t, &
                                meridian_bspline_build, meridian_bspline_exterior_build, &
                meridian_bspline_solve, meridian_bspline_load_value, meridian_bspline_eval
    use fortfem_meridian_magnetic_coupling, only: meridian_exterior_t
    implicit none

    type(meridian_bspline_t) :: sp
    type(meridian_exterior_t) :: ex
    complex(dp), allocatable :: f(:), u(:), p(:)
    complex(dp) :: g(3, 1), a(3, 2), b(3, 2)
    real(dp) :: x(2, 1), w(1), xe(2, 2)
    integer :: status

    call meridian_bspline_build(1.0_dp, 3.0_dp, -1.0_dp, 1.0_dp, 12, 12, 2, sp, status)
    call check(status == 0, "build")
    call meridian_bspline_exterior_build(sp, 0, ex, status)
    call check(status == 0, "exterior")
    allocate (f(sp%n1))
    f = (0.0_dp, 0.0_dp)
    x(:, 1) = [2.0_dp, 0.0_dp]
    w = 1.0_dp
    g(:, 1) = [(0.0_dp, 0.0_dp), (1.0_dp, 0.0_dp), (0.0_dp, 0.0_dp)]
    call meridian_bspline_load_value(sp, 0, x, w, g, f, status)
    call check(status == 0, "load")
    call meridian_bspline_solve(sp, ex, f, u, p, status)
    call check(status == 0, "solve")
    xe(:, 1) = [2.0_dp, 0.0_dp]
    xe(:, 2) = [2.0_dp, 0.9_dp]
    call meridian_bspline_eval(sp, 0, u, xe, a, b, status)
    call check(status == 0, "eval")
    call check(real(a(2, 1), dp) > 0.0_dp, "A_phi has the sign of the source")
    call check(real(a(2, 2), dp) < real(a(2, 1), dp), "A_phi decays")
    print *, "PASS"

contains

    subroutine check(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what
        if (.not. ok) then
            print *, "FAIL: ", what
            error stop 1
        end if
    end subroutine check

end program test_fortfem_meridian
