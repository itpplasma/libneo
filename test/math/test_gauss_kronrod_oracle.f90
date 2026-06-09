module gk_oracle_battery
! Shared integrand battery for the GSL oracle comparison; the selector and
! evaluation counter travel through the GSL params pointer, so no module
! state is needed.

    use libneo_kinds, only: dp
    use iso_c_binding, only: c_double, c_int, c_ptr, c_f_pointer

    implicit none

    type, bind(c) :: battery_state_t
        integer(c_int) :: sel
        integer(c_int) :: neval
    end type battery_state_t

contains

    function battery_eval(sel, x) result(fx)
        integer, intent(in) :: sel
        real(dp), intent(in) :: x
        real(dp) :: fx

        select case (sel)
        case (1)
            fx = exp(x)
        case (2)
            fx = sin(50.0d0*x)
        case (3)
            fx = sqrt(x)
        case (4)
            fx = 1.0d0/(x*x + 1.0d-8)
        case (5)
            fx = exp(-400.0d0*(x - 0.3d0)**2)
        case (6)
            fx = 1.0d0/sqrt(x)
        case default
            fx = 0.0d0
        end select
    end function battery_eval

    function f_exp_bench(x) result(fx)
        real(dp), intent(in) :: x
        real(dp) :: fx

        fx = exp(x)
    end function f_exp_bench

    function f_sin50_bench(x) result(fx)
        real(dp), intent(in) :: x
        real(dp) :: fx

        fx = sin(50.0d0*x)
    end function f_sin50_bench

    function gsl_callback(x, params) bind(c) result(fx)
        real(c_double), value :: x
        type(c_ptr), value :: params
        real(c_double) :: fx

        type(battery_state_t), pointer :: state

        call c_f_pointer(params, state)
        state%neval = state%neval + 1
        fx = battery_eval(int(state%sel), x)
    end function gsl_callback

end module gk_oracle_battery

program test_gauss_kronrod_oracle

    use libneo_kinds, only: dp
    use iso_c_binding, only: c_double, c_int, c_size_t, c_ptr, c_funptr, &
                             c_loc, c_funloc
    use neo_gauss_kronrod, only: integrate_gk
    use gk_oracle_battery, only: battery_state_t, battery_eval, gsl_callback, &
                                 f_exp_bench, f_sin50_bench
    use util_for_test, only: print_test, print_ok, print_fail

    implicit none

    type, bind(c) :: gsl_function_t
        type(c_funptr) :: function
        type(c_ptr) :: params
    end type gsl_function_t

    interface
        function gsl_integration_workspace_alloc(n) &
            bind(c, name="gsl_integration_workspace_alloc") result(w)
            import :: c_size_t, c_ptr
            integer(c_size_t), value :: n
            type(c_ptr) :: w
        end function gsl_integration_workspace_alloc

        subroutine gsl_integration_workspace_free(w) &
            bind(c, name="gsl_integration_workspace_free")
            import :: c_ptr
            type(c_ptr), value :: w
        end subroutine gsl_integration_workspace_free

        function gsl_integration_qag(f, a, b, epsabs, epsrel, limit, key, &
                                     workspace, result, abserr) &
            bind(c, name="gsl_integration_qag") result(status)
            import :: c_double, c_int, c_size_t, c_ptr, gsl_function_t
            type(gsl_function_t), intent(in) :: f
            real(c_double), value :: a, b, epsabs, epsrel
            integer(c_size_t), value :: limit
            integer(c_int), value :: key
            type(c_ptr), value :: workspace
            real(c_double), intent(out) :: result, abserr
            integer(c_int) :: status
        end function gsl_integration_qag

        function gsl_set_error_handler_off() &
            bind(c, name="gsl_set_error_handler_off") result(old)
            import :: c_funptr
            type(c_funptr) :: old
        end function gsl_set_error_handler_off
    end interface

    integer, parameter :: ncase = 6
    real(dp), parameter :: a_case(ncase) = [0.0d0, 0.0d0, 0.0d0, -1.0d0, &
                                            0.0d0, 0.0d0]
    real(dp), parameter :: b_case(ncase) = [1.0d0, 1.0d0, 1.0d0, 1.0d0, &
                                            1.0d0, 1.0d0]
    ! references from mpmath.quad at 40 digits
    real(dp), parameter :: ref_case(ncase) = [ &
        1.71828182845904523536028747135266d0, &
        7.0067943015773451862085882197966d-4, &
        0.666666666666666666666666666666667d0, &
        31413.9265359045990512531004997474d0, &
        0.0886226925452758004113398690035196d0, &
        2.0d0]
    character(len=8), parameter :: name_case(ncase) = [ &
        "exp     ", "sin50x  ", "sqrt    ", "lorentz ", "gpeak   ", "invsqrt "]
    integer, parameter :: keys(4) = [15, 21, 31, 61]
    integer, parameter :: gsl_keys(4) = [1, 2, 3, 6]
    real(dp), parameter :: epsabs = 0.0d0, epsrel = 1.0d-10
    integer, parameter :: limit = 200

    type(battery_state_t), target :: state
    type(gsl_function_t) :: gsl_f
    type(c_ptr) :: workspace
    type(c_funptr) :: old_handler
    real(dp) :: max_rel_dev

    old_handler = gsl_set_error_handler_off()
    workspace = gsl_integration_workspace_alloc(int(limit, c_size_t))
    gsl_f%function = c_funloc(gsl_callback)
    gsl_f%params = c_loc(state)
    max_rel_dev = 0.0d0

    call compare_battery
    print '(a, es10.2)', " max |ours - gsl| / |gsl| over battery: ", max_rel_dev
    call benchmark

    call gsl_integration_workspace_free(workspace)

contains

    subroutine compare_battery
        real(dp) :: res, abserr, res_gsl, abserr_gsl, dev, tol
        integer :: ierr, icase, k, neval_ours, status

        call print_test("test_gauss_kronrod_vs_gsl_qag")
        do icase = 1, ncase
            do k = 1, size(keys)
                state%sel = icase
                state%neval = 0
                call integrate_gk(f_current, a_case(icase), b_case(icase), &
                                  epsabs, epsrel, res, abserr, ierr, &
                                  key=keys(k), limit=limit)
                neval_ours = state%neval
                state%neval = 0
                status = gsl_integration_qag(gsl_f, a_case(icase), &
                                             b_case(icase), epsabs, epsrel, &
                                             int(limit, c_size_t), gsl_keys(k), &
                                             workspace, res_gsl, abserr_gsl)
                if (ierr /= 0 .or. status /= 0) then
                    call print_fail
                    print *, name_case(icase), " key ", keys(k), &
                        " ierr ", ierr, " gsl status ", status
                    stop "integration failed"
                end if
                dev = abs(res - res_gsl)
                tol = max(epsabs, epsrel*abs(res_gsl))
                if (dev > tol) then
                    call print_fail
                    print *, name_case(icase), " key ", keys(k), &
                        " ours ", res, " gsl ", res_gsl
                    stop "deviation from GSL exceeds tolerance"
                end if
                if (abs(res_gsl) > 0.0d0) then
                    max_rel_dev = max(max_rel_dev, dev/abs(res_gsl))
                end if
                print '(2a, i3, a, i6, a, i6, a, es10.2)', "    ", &
                    name_case(icase), keys(k), "  neval ours ", neval_ours, &
                    "  gsl ", state%neval, "  |diff| ", dev
            end do
        end do
        call print_ok
    end subroutine compare_battery

    subroutine benchmark
        integer, parameter :: nrep_smooth = 200000, nrep_osc = 20000
        real(dp) :: res, abserr, sink
        integer :: ierr, i, status

        call print_test("benchmark_vs_gsl (ns per integrate call)")
        sink = 0.0d0

        state%sel = 1
        call bench_ours(nrep_smooth, f_exp_bench, "exp    ours: ", sink)
        call bench_gsl(nrep_smooth, "exp    gsl:  ", sink)
        state%sel = 2
        call bench_ours(nrep_osc, f_sin50_bench, "sin50x ours: ", sink)
        call bench_gsl(nrep_osc, "sin50x gsl:  ", sink)

        if (sink /= sink) then
            call print_fail
            stop "benchmark produced NaN"
        end if
        call print_ok
    end subroutine benchmark

    subroutine bench_ours(nrep, fb, label, sink)
        use neo_gauss_kronrod, only: gk_integrand
        integer, intent(in) :: nrep
        procedure(gk_integrand) :: fb
        character(*), intent(in) :: label
        real(dp), intent(inout) :: sink

        real(dp) :: res, abserr, t
        integer :: ierr, i
        integer(8) :: t0, t1, rate

        call system_clock(t0, rate)
        do i = 1, nrep
            call integrate_gk(fb, a_case(state%sel), b_case(state%sel), &
                              epsabs, epsrel, res, abserr, ierr, key=21, &
                              limit=limit)
            sink = sink + res
        end do
        call system_clock(t1)
        t = real(t1 - t0, dp)/real(rate, dp)
        print '(2a, f12.1)', "    ", label, 1.0d9*t/nrep
    end subroutine bench_ours

    subroutine bench_gsl(nrep, label, sink)
        integer, intent(in) :: nrep
        character(*), intent(in) :: label
        real(dp), intent(inout) :: sink

        real(dp) :: res, abserr, t
        integer :: i, status
        integer(8) :: t0, t1, rate

        call system_clock(t0, rate)
        do i = 1, nrep
            status = gsl_integration_qag(gsl_f, a_case(state%sel), &
                                         b_case(state%sel), epsabs, epsrel, &
                                         int(limit, c_size_t), 2, workspace, &
                                         res, abserr)
            sink = sink + res
        end do
        call system_clock(t1)
        t = real(t1 - t0, dp)/real(rate, dp)
        print '(2a, f12.1)', "    ", label, 1.0d9*t/nrep
    end subroutine bench_gsl

    function f_current(x) result(fx)
        real(dp), intent(in) :: x
        real(dp) :: fx

        state%neval = state%neval + 1
        fx = battery_eval(int(state%sel), x)
    end function f_current

end program test_gauss_kronrod_oracle
