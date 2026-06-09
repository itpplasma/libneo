program test_fft_oracle
    ! Compares neo_fft against FFTW3 over a battery of lengths (including the
    ! primes 97, 1009 and typical coil_tools nphi values) and benchmarks
    ! ns/transform with plan reuse on both sides. Build standalone:
    ! gfortran -O2 libneo_kinds.f90 src/math/fft.f90 util_for_test.f90 \
    !     test/math/test_fft_oracle.f90 -L/opt/homebrew/lib -lfftw3

    use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double, c_double_complex
    use, intrinsic :: iso_fortran_env, only: int64
    use libneo_kinds, only: dp
    use neo_fft, only: fft_r2c, neo_fft_plan_t, neo_fft_plan_init
    use util_for_test, only: print_test, print_ok, print_fail

    implicit none

    interface
        function fftw_plan_dft_r2c_1d(n, in, out, flags) &
            bind(c, name="fftw_plan_dft_r2c_1d")
            import :: c_ptr, c_int, c_double, c_double_complex
            integer(c_int), value :: n
            real(c_double), dimension(*) :: in
            complex(c_double_complex), dimension(*) :: out
            integer(c_int), value :: flags
            type(c_ptr) :: fftw_plan_dft_r2c_1d
        end function fftw_plan_dft_r2c_1d
        subroutine fftw_execute_dft_r2c(plan, in, out) &
            bind(c, name="fftw_execute_dft_r2c")
            import :: c_ptr, c_double, c_double_complex
            type(c_ptr), value :: plan
            real(c_double), dimension(*) :: in
            complex(c_double_complex), dimension(*) :: out
        end subroutine fftw_execute_dft_r2c
        subroutine fftw_destroy_plan(plan) bind(c, name="fftw_destroy_plan")
            import :: c_ptr
            type(c_ptr), value :: plan
        end subroutine fftw_destroy_plan
    end interface

    integer(c_int), parameter :: fftw_measure = 0
    integer, parameter :: acc_sizes(23) = [ &
        16, 24, 30, 32, 36, 60, 64, 96, 97, 100, 128, 180, 256, 360, 480, &
        512, 720, 960, 1009, 1024, 2048, 4093, 4096]
    integer, parameter :: bench_sizes(6) = [64, 128, 256, 512, 1024, 4096]
    real(dp) :: worst_scaled_err
    integer :: i

    worst_scaled_err = 0.0d0
    call print_test("fft_r2c vs FFTW over size battery (random + structured)")
    do i = 1, size(acc_sizes)
        call check_against_fftw(acc_sizes(i), worst_scaled_err)
    end do
    call print_ok
    print '(a, es10.3, a)', " worst |err| / (n max|x|) = ", worst_scaled_err, &
        "  (required < 1e-13)"

    print '(/, a)', " benchmark, plan reuse both sides (FFTW_MEASURE):"
    print '(a)', "      n    neo_fft ns/call    fftw ns/call    ratio"
    do i = 1, size(bench_sizes)
        call benchmark(bench_sizes(i))
    end do

contains

    subroutine check_against_fftw(n, worst)
        integer, intent(in) :: n
        real(dp), intent(inout) :: worst
        real(dp), parameter :: tol_fac = 1.0d-13
        real(dp), parameter :: twopi = 8.0d0*atan(1.0d0)
        real(dp) :: x(n), err, tol
        complex(dp) :: c(n/2 + 1), cref(n/2 + 1)
        integer :: trial, j

        do trial = 1, 2
            if (trial == 1) then
                call fill_random(x, n + trial)
                x = x*10.0d0**(8*x)
            else
                do j = 1, n
                    x(j) = 0.3d0 + cos(twopi*modulo(3*(j - 1), n)/real(n, dp)) &
                           - 2.0d0*sin(twopi*modulo((n/3)*(j - 1), n)/real(n, dp))
                end do
                x(1) = x(1) + 5.0d0
                x(n/2 + 1) = x(n/2 + 1) - 3.0d0
            end if
            call fftw_reference(n, x, cref)
            call fft_r2c(x, c)
            err = maxval(abs(c - cref))
            tol = tol_fac*real(n, dp)*maxval(abs(x))
            worst = max(worst, err/(real(n, dp)*maxval(abs(x))))
            if (err > tol) then
                call print_fail
                print *, "n = ", n, ", trial ", trial, ": err ", err, " > tol ", tol
                stop "fft_r2c deviates from FFTW"
            end if
        end do
    end subroutine check_against_fftw

    subroutine fftw_reference(n, x, cref)
        integer, intent(in) :: n
        real(dp), intent(in) :: x(:)
        complex(dp), intent(out) :: cref(:)
        real(c_double) :: xin(n)
        complex(c_double_complex) :: cout(n/2 + 1)
        type(c_ptr) :: plan

        plan = fftw_plan_dft_r2c_1d(int(n, c_int), xin, cout, fftw_measure)
        xin = x
        call fftw_execute_dft_r2c(plan, xin, cout)
        cref = cout
        call fftw_destroy_plan(plan)
    end subroutine fftw_reference

    subroutine benchmark(n)
        integer, intent(in) :: n
        real(dp) :: x(n)
        real(c_double) :: xin(n)
        complex(dp) :: c(n/2 + 1)
        complex(c_double_complex) :: cout(n/2 + 1)
        type(neo_fft_plan_t) :: plan
        type(c_ptr) :: fplan
        integer(int64) :: t0, t1, rate
        integer :: nrep, rep
        real(dp) :: ns_neo, ns_fftw, sink

        nrep = max(100000000/n, 1000)
        call fill_random(x, n)
        call neo_fft_plan_init(plan, n)
        call fft_r2c(x, c, plan)
        call system_clock(t0, rate)
        do rep = 1, nrep
            call fft_r2c(x, c, plan)
        end do
        call system_clock(t1)
        sink = abs(c(2))
        ns_neo = 1.0d9*real(t1 - t0, dp)/real(rate, dp)/real(nrep, dp)

        fplan = fftw_plan_dft_r2c_1d(int(n, c_int), xin, cout, fftw_measure)
        xin = x
        call fftw_execute_dft_r2c(fplan, xin, cout)
        call system_clock(t0)
        do rep = 1, nrep
            call fftw_execute_dft_r2c(fplan, xin, cout)
        end do
        call system_clock(t1)
        sink = sink + abs(cout(2))
        ns_fftw = 1.0d9*real(t1 - t0, dp)/real(rate, dp)/real(nrep, dp)
        call fftw_destroy_plan(fplan)

        print '(i7, f15.1, f17.1, f11.2)', n, ns_neo, ns_fftw, ns_neo/ns_fftw
        if (sink < 0.0d0) print *, sink
    end subroutine benchmark

    subroutine fill_random(x, seed_val)
        real(dp), intent(out) :: x(:)
        integer, intent(in) :: seed_val
        integer, allocatable :: seed(:)
        integer :: nseed

        call random_seed(size=nseed)
        allocate (seed(nseed))
        seed = seed_val
        call random_seed(put=seed)
        call random_number(x)
        x = 2.0d0*x - 1.0d0
    end subroutine fill_random

end program test_fft_oracle
