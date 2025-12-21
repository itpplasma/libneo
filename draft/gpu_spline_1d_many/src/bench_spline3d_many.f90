program bench_spline3d_many
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: iso_c_binding, only: c_ptr, c_loc, c_f_pointer
    use batch_interpolate_types, only: BatchSplineData3D
    use batch_interpolate_3d, only: construct_batch_splines_3d, destroy_batch_splines_3d
    use draft_batch_splines_many_api, only: evaluate_batch_splines_3d_many
    use spline3d_many_offload, only: spline3d_many_setup, &
                                     spline3d_many_teardown, &
                                     spline3d_many_eval_host, &
                                     spline3d_many_eval_resident
    use util_lcg_rng, only: lcg_state_t, lcg_init, lcg_uniform_0_1
    implicit none

    integer, parameter :: order(3) = [5, 3, 3]
    integer, parameter :: num_points(3) = [48, 32, 32]
    integer, parameter :: num_quantities = 8
    integer, parameter :: npts = 200000
    integer, parameter :: niter = 6
    logical, parameter :: periodic(3) = [.true., .true., .true.]

    real(dp), parameter :: x_min(3) = [1.23d0, -0.7d0, 0.5d0]
    real(dp), parameter :: x_max(3) = [x_min(1) + 40.0d0, x_min(2) + 30.0d0, &
                                       x_min(3) + 20.0d0]

    type(BatchSplineData3D) :: spl
    real(dp), allocatable, target :: y_ref(:)
    real(dp), allocatable :: y_out(:)
    real(dp), pointer :: y_ref2d(:, :)
    type(c_ptr) :: y_ref_ptr
    real(dp), allocatable :: x_eval(:, :)
    real(dp), allocatable :: y_grid(:, :, :, :)

    type(lcg_state_t) :: rng
    integer :: i1, i2, i3, iq, ipt
    real(dp) :: best, diff_max

    call lcg_init(rng, 97531_int64)

    allocate (y_grid(num_points(1), num_points(2), num_points(3), num_quantities))
    do iq = 1, num_quantities
        do i3 = 1, num_points(3)
            do i2 = 1, num_points(2)
                do i1 = 1, num_points(1)
                    y_grid(i1, i2, i3, iq) = &
                        cos(x_min(1) + real(i1 - 1, dp)*(x_max(1) - x_min(1))/ &
                            real(num_points(1) - 1, dp) + real(iq - 1, dp)*0.1d0)* &
                        cos(x_min(2) + real(i2 - 1, dp)*(x_max(2) - x_min(2))/ &
                            real(num_points(2) - 1, dp) + real(iq - 1, dp)*0.2d0)* &
                        cos(x_min(3) + real(i3 - 1, dp)*(x_max(3) - x_min(3))/ &
                            real(num_points(3) - 1, dp) + real(iq - 1, dp)*0.05d0)
                end do
            end do
        end do
    end do

    call construct_batch_splines_3d(x_min, x_max, y_grid, order, periodic, spl)

    allocate (x_eval(3, npts))
    do ipt = 1, npts
        x_eval(1, ipt) = x_min(1) + (x_max(1) - x_min(1))*lcg_uniform_0_1(rng)
        x_eval(2, ipt) = x_min(2) + (x_max(2) - x_min(2))*lcg_uniform_0_1(rng)
        x_eval(3, ipt) = x_min(3) + (x_max(3) - x_min(3))*lcg_uniform_0_1(rng)
    end do

    allocate (y_ref(num_quantities*npts))
    allocate (y_out(num_quantities*npts))

    y_ref_ptr = c_loc(y_ref(1))
    call c_f_pointer(y_ref_ptr, y_ref2d, [spl%num_quantities, npts])
    call evaluate_batch_splines_3d_many(spl, x_eval, y_ref2d)

    print *, "Benchmark 3D parameters"
    print *, "order         ", order
    print *, "num_points    ", num_points
    print *, "num_quantities", num_quantities
    print *, "npts          ", npts
    print *, "niter         ", niter
    print *, "periodic      ", periodic

    call bench_cpu()
#if defined(LIBNEO_ENABLE_OPENACC)
    call bench_openacc()
#elif defined(LIBNEO_ENABLE_OPENMP)
    call bench_openmp()
#endif

    call destroy_batch_splines_3d(spl)

contains

    real(dp) function wall_time() result(t)
        integer :: count, rate, max_count
        call system_clock(count, rate, max_count)
        t = real(count, dp)/real(rate, dp)
    end function wall_time

    subroutine bench_cpu()
        integer :: it
        real(dp) :: t0, t1, dt

        best = huge(1.0d0)
        do it = 1, niter
            t0 = wall_time()
            call spline3d_many_eval_host(spl%order, spl%num_points, spl%num_quantities, &
                                         spl%periodic, spl%x_min, spl%h_step, spl%coeff, &
                                         x_eval, y_out)
            t1 = wall_time()
            dt = t1 - t0
            best = min(best, dt)
        end do
        diff_max = maxval(abs(y_out - y_ref))
        print *, "cpu best_s ", best, " pts_per_s ", real(npts, dp)/best, &
            " max_abs_diff ", &
            diff_max
    end subroutine bench_cpu

    subroutine bench_openacc()
        integer :: it
        real(dp) :: t0, t1, dt

        call spline3d_many_setup(spl%coeff, x_eval, y_out)
        call spline3d_many_eval_resident(spl%order, spl%num_points, spl%num_quantities, &
                                         spl%periodic, spl%x_min, spl%h_step, spl%coeff, &
                                         x_eval, y_out)
        !$acc wait

        best = huge(1.0d0)
        do it = 1, niter
            t0 = wall_time()
            call spline3d_many_eval_resident(spl%order, spl%num_points, spl%num_quantities, &
                                             spl%periodic, spl%x_min, spl%h_step, spl%coeff, &
                                             x_eval, y_out)
            !$acc wait
            t1 = wall_time()
            dt = t1 - t0
            best = min(best, dt)
        end do
        !$acc update self(y_out)
        diff_max = maxval(abs(y_out - y_ref))
        print *, "openacc best_s ", best, " pts_per_s ", real(npts, dp)/best, &
            " max_abs_diff ", diff_max
        call spline3d_many_teardown(spl%coeff, x_eval, y_out)
    end subroutine bench_openacc

    subroutine bench_openmp()
        integer :: it
        real(dp) :: t0, t1, dt

        call spline3d_many_setup(spl%coeff, x_eval, y_out)
        call spline3d_many_eval_resident(spl%order, spl%num_points, spl%num_quantities, &
                                         spl%periodic, spl%x_min, spl%h_step, spl%coeff, &
                                         x_eval, y_out)

        best = huge(1.0d0)
        do it = 1, niter
            t0 = wall_time()
            call spline3d_many_eval_resident(spl%order, spl%num_points, spl%num_quantities, &
                                             spl%periodic, spl%x_min, spl%h_step, spl%coeff, &
                                             x_eval, y_out)
            t1 = wall_time()
            dt = t1 - t0
            best = min(best, dt)
        end do

        !$omp target update from(y_out)
        diff_max = maxval(abs(y_out - y_ref))
        print *, "openmp best_s ", best, " pts_per_s ", real(npts, dp)/best, &
            " max_abs_diff ", diff_max

        call spline3d_many_teardown(spl%coeff, x_eval, y_out)
    end subroutine bench_openmp

end program bench_spline3d_many
