program bench_neo_bspline_vs_interpolate
    use, intrinsic :: iso_fortran_env, only: dp_bench => real64
    use, intrinsic :: iso_c_binding, only: c_int, c_long
    use interpolate
    use neo_bspline

    implicit none

    integer, parameter :: N_CASE = 9
    integer, parameter :: N_TOTAL_LIST(N_CASE) = &
        [10, 30, 100, 300, 1000, 3000, 10000, 30000, 100000]
    real(dp_bench), parameter :: PI_BENCH = acos(-1.0_dp_bench)
    real(dp_bench), parameter :: TWOPI_BENCH = 2.0_dp_bench*PI_BENCH

    type, bind(C) :: timespec
        integer(c_long) :: tv_sec
        integer(c_long) :: tv_nsec
    end type timespec

!> clock id for clock_gettime; macOS uses 6 for CLOCK_MONOTONIC, Linux uses 1.
#ifdef __APPLE__
    integer(c_int), parameter :: CLOCK_MONOTONIC = 6_c_int
#else
    integer(c_int), parameter :: CLOCK_MONOTONIC = 1_c_int
#endif

    interface
        function c_clock_gettime(clock_id, tp) bind(C, name="clock_gettime")
            import :: c_int, timespec
            integer(c_int), value :: clock_id
            type(timespec) :: tp
            integer(c_int) :: c_clock_gettime
        end function c_clock_gettime
    end interface

    call run_bench_1d()
    call run_bench_2d()
    call run_bench_3d()

contains

    subroutine bench_time_now(t)
        real(dp_bench), intent(out) :: t
        type(timespec) :: ts
        integer(c_int) :: ierr

        ierr = c_clock_gettime(CLOCK_MONOTONIC, ts)
        if (ierr == 0_c_int) then
            t = real(ts%tv_sec, dp_bench) &
                + real(ts%tv_nsec, dp_bench)*1.0e-9_dp_bench
        else
            ! Fallback for platforms where CLOCK_MONOTONIC is unavailable.
            call cpu_time(t)
        end if
    end subroutine bench_time_now

    subroutine run_bench_1d()
        type(SplineData1D) :: spl_interp
        type(bspline_1d) :: spl_bs, spl_bs_dir
        integer :: i, k, n, degree_bs
        real(dp_bench) :: t0, t1
        real(dp_bench) :: t_create_interp, t_eval_interp
        real(dp_bench) :: t_create_lsq, t_eval_lsq
        real(dp_bench) :: t_create_dir, t_eval_dir
        real(dp_bench), allocatable :: x(:), f(:), coeff_lsq(:), coeff_dir(:)
        real(dp_bench) :: x_min, x_max, y
        integer :: unit, n_ctrl_lsq

        x_min = 1.23_dp_bench
        x_max = x_min + TWOPI_BENCH

        open(newunit=unit, file="bench_neo_bspline_1d.dat", &
            status="replace", action="write", form="formatted")

        do k = 1, N_CASE
            n = N_TOTAL_LIST(k)
            if (n < 8) cycle

            allocate(x(n), f(n))
            do i = 1, n
                x(i) = x_min + (x_max - x_min)*real(i - 1, dp_bench) &
                    /real(n - 1, dp_bench)
                f(i) = cos(2.0_dp_bench*x(i)) + 0.5_dp_bench*sin(3.0_dp_bench*x(i))
            end do

            ! 1) interpolate module: construct + evaluate
            call bench_time_now(t0)
            call construct_splines_1d(x_min, x_max, f, 5, .false., spl_interp)
            call bench_time_now(t1)
            t_create_interp = t1 - t0

            call bench_time_now(t0)
            do i = 1, n
                call evaluate_splines_1d(spl_interp, x(i), y)
            end do
            call bench_time_now(t1)
            t_eval_interp = t1 - t0

            ! 2) neo_bspline LSQ: construct + evaluate
            degree_bs = 5
            n_ctrl_lsq = min(40, n/2)
            n_ctrl_lsq = max(degree_bs + 1, n_ctrl_lsq)

            allocate(coeff_lsq(n_ctrl_lsq))
            coeff_lsq = 0.0_dp_bench

            call bench_time_now(t0)
            call bspline_1d_init_uniform(spl_bs, degree_bs, n_ctrl_lsq, x_min, &
                x_max)
            call bspline_1d_lsq_cgls(spl_bs, x, f, coeff_lsq, max_iter=400, &
                tol=1.0d-10)
            call bench_time_now(t1)
            t_create_lsq = t1 - t0

            call bench_time_now(t0)
            do i = 1, n
                call bspline_1d_eval(spl_bs, coeff_lsq, x(i), y)
            end do
            call bench_time_now(t1)
            t_eval_lsq = t1 - t0

            ! 3) neo_bspline direct interp: construct + evaluate (skip large N)
            allocate(coeff_dir(n))

            call bench_time_now(t0)
            call bspline_1d_init_uniform(spl_bs_dir, degree_bs, n, x_min, x_max)
            call bspline_1d_interp(spl_bs_dir, x, f, coeff_dir)
            call bench_time_now(t1)
            t_create_dir = t1 - t0

            call bench_time_now(t0)
            do i = 1, n
                call bspline_1d_eval(spl_bs_dir, coeff_dir, x(i), y)
            end do
            call bench_time_now(t1)
            t_eval_dir = t1 - t0

            deallocate(coeff_dir)

            write(unit,'(i10,6es16.8)') n, t_create_interp, t_eval_interp, &
                t_create_lsq, t_eval_lsq, t_create_dir, t_eval_dir

            deallocate(x, f, coeff_lsq)
            call destroy_splines_1d(spl_interp)
        end do

        close(unit)
    end subroutine run_bench_1d


    subroutine run_bench_2d()
        type(SplineData2D) :: spl_interp
        type(bspline_2d) :: spl_bs, spl_bs_dir
        integer :: i1, i2, k
        integer :: n_total, n_side, n1, n2
        integer :: degree(2), n_ctrl_lsq(2)
        real(dp_bench) :: t0, t1
        real(dp_bench) :: t_create_interp, t_eval_interp
        real(dp_bench) :: t_create_lsq, t_eval_lsq
        real(dp_bench) :: t_create_dir, t_eval_dir
        real(dp_bench), allocatable :: x1(:), x2(:)
        real(dp_bench), allocatable :: f_grid(:,:)
        real(dp_bench), allocatable :: coeff_lsq(:,:), coeff_dir(:,:)
        real(dp_bench) :: x_min(2), x_max(2)
        real(dp_bench) :: x(2), y
        integer :: unit, n_data

        x_min = [1.23_dp_bench, 1.23_dp_bench]
        x_max = x_min + [TWOPI_BENCH, TWOPI_BENCH]

        open(newunit=unit, file="bench_neo_bspline_2d.dat", &
            status="replace", action="write", form="formatted")

        degree = [3, 3]

        do k = 1, N_CASE
            n_total = N_TOTAL_LIST(k)
            n_side = int(sqrt(real(n_total, dp_bench)))
            if (n_side < degree(1) + 2) n_side = degree(1) + 2

            n1 = n_side
            n2 = n_side
            n_data = n1*n2

            allocate(x1(n1), x2(n2), f_grid(n1, n2))

            do i1 = 1, n1
                x1(i1) = x_min(1) + (x_max(1) - x_min(1))* &
                    real(i1 - 1, dp_bench)/real(n1 - 1, dp_bench)
            end do
            do i2 = 1, n2
                x2(i2) = x_min(2) + (x_max(2) - x_min(2))* &
                    real(i2 - 1, dp_bench)/real(n2 - 1, dp_bench)
            end do

            do i2 = 1, n2
                do i1 = 1, n1
                    f_grid(i1, i2) = cos(x1(i1))*cos(2.0_dp_bench*x2(i2))
                end do
            end do

            ! 1) interpolate module: construct + evaluate
            call bench_time_now(t0)
            call construct_splines_2d(x_min, x_max, f_grid, [5, 5], &
                [.false., .false.], spl_interp)
            call bench_time_now(t1)
            t_create_interp = t1 - t0

            call bench_time_now(t0)
            do i2 = 1, n2
                do i1 = 1, n1
                    x = [x1(i1), x2(i2)]
                    call evaluate_splines_2d(spl_interp, x, y)
                end do
            end do
            call bench_time_now(t1)
            t_eval_interp = t1 - t0

            ! 2) neo_bspline LSQ: construct + evaluate
            n_ctrl_lsq(1) = min(16, n1/2)
            n_ctrl_lsq(2) = min(16, n2/2)
            n_ctrl_lsq(1) = max(degree(1) + 1, n_ctrl_lsq(1))
            n_ctrl_lsq(2) = max(degree(2) + 1, n_ctrl_lsq(2))

            allocate(coeff_lsq(n_ctrl_lsq(1), n_ctrl_lsq(2)))
            coeff_lsq = 0.0_dp_bench

            call bench_time_now(t0)
            call bspline_2d_init_uniform(spl_bs, degree, n_ctrl_lsq, x_min, x_max)
            call bspline_2d_lsq_cgls(spl_bs, x1, x2, f_grid, coeff_lsq, &
                max_iter=800, tol=1.0d-10)
            call bench_time_now(t1)
            t_create_lsq = t1 - t0

            call bench_time_now(t0)
            do i2 = 1, n2
                do i1 = 1, n1
                    x = [x1(i1), x2(i2)]
                    call bspline_2d_eval(spl_bs, coeff_lsq, x, y)
                end do
            end do
            call bench_time_now(t1)
            t_eval_lsq = t1 - t0

            ! 3) neo_bspline direct interp: construct + evaluate (skip large N)
            allocate(coeff_dir(n1, n2))

            call bench_time_now(t0)
            call bspline_2d_init_uniform(spl_bs_dir, degree, [n1, n2], x_min, &
                x_max)
            call bspline_2d_interp(spl_bs_dir, x1, x2, f_grid, coeff_dir)
            call bench_time_now(t1)
            t_create_dir = t1 - t0

            call bench_time_now(t0)
            do i2 = 1, n2
                do i1 = 1, n1
                    x = [x1(i1), x2(i2)]
                    call bspline_2d_eval(spl_bs_dir, coeff_dir, x, y)
                end do
            end do
            call bench_time_now(t1)
            t_eval_dir = t1 - t0

            deallocate(coeff_dir)

            write(unit,'(i10,6es16.8)') n_data, t_create_interp, t_eval_interp, &
                t_create_lsq, t_eval_lsq, t_create_dir, t_eval_dir

            deallocate(x1, x2, f_grid, coeff_lsq)
            call destroy_splines_2d(spl_interp)
        end do

        close(unit)
    end subroutine run_bench_2d


    subroutine run_bench_3d()
        type(SplineData3D) :: spl_interp
        type(bspline_3d) :: spl_bs, spl_bs_dir
        integer :: i1, i2, i3, k
        integer :: n_total, n_side, n1, n2, n3
        integer :: degree(3), n_ctrl_lsq(3)
        real(dp_bench) :: t0, t1
        real(dp_bench) :: t_create_interp, t_eval_interp
        real(dp_bench) :: t_create_lsq, t_eval_lsq
        real(dp_bench) :: t_create_dir, t_eval_dir
        real(dp_bench), allocatable :: x1(:), x2(:), x3(:)
        real(dp_bench), allocatable :: f3d(:,:,:)
        real(dp_bench), allocatable :: coeff_lsq(:,:,:), coeff_dir(:,:,:)
        real(dp_bench) :: x_min(3), x_max(3)
        real(dp_bench) :: x(3), y
        integer :: unit, n_data

        x_min = [1.23_dp_bench, 1.23_dp_bench, 1.23_dp_bench]
        x_max = x_min + [TWOPI_BENCH, TWOPI_BENCH, TWOPI_BENCH]

        open(newunit=unit, file="bench_neo_bspline_3d.dat", &
            status="replace", action="write", form="formatted")

        degree = [3, 3, 3]

        do k = 1, N_CASE
            n_total = N_TOTAL_LIST(k)
            n_side = int(real(n_total, dp_bench)**(1.0_dp_bench/3.0_dp_bench))
            if (n_side < degree(1) + 2) n_side = degree(1) + 2

            n1 = n_side
            n2 = n_side
            n3 = n_side
            n_data = n1*n2*n3

            allocate(x1(n1), x2(n2), x3(n3))
            allocate(f3d(n1, n2, n3))

            do i1 = 1, n1
                x1(i1) = x_min(1) + (x_max(1) - x_min(1))* &
                    real(i1 - 1, dp_bench)/real(n1 - 1, dp_bench)
            end do
            do i2 = 1, n2
                x2(i2) = x_min(2) + (x_max(2) - x_min(2))* &
                    real(i2 - 1, dp_bench)/real(n2 - 1, dp_bench)
            end do
            do i3 = 1, n3
                x3(i3) = x_min(3) + (x_max(3) - x_min(3))* &
                    real(i3 - 1, dp_bench)/real(n3 - 1, dp_bench)
            end do

            do i3 = 1, n3
                do i2 = 1, n2
                    do i1 = 1, n1
                        f3d(i1, i2, i3) = cos(x1(i1))*cos(2.0_dp_bench*x2(i2)) &
                            *cos(3.0_dp_bench*x3(i3))
                    end do
                end do
            end do

            ! 1) interpolate module: construct + evaluate
            call bench_time_now(t0)
            call construct_splines_3d(x_min, x_max, f3d, [5, 5, 5], &
                [.false., .false., .false.], spl_interp)
            call bench_time_now(t1)
            t_create_interp = t1 - t0

            call bench_time_now(t0)
            do i3 = 1, n3
                do i2 = 1, n2
                    do i1 = 1, n1
                        x = [x1(i1), x2(i2), x3(i3)]
                        call evaluate_splines_3d(spl_interp, x, y)
                    end do
                end do
            end do
            call bench_time_now(t1)
            t_eval_interp = t1 - t0

            ! 2) neo_bspline LSQ: construct + evaluate
            n_ctrl_lsq(1) = min(16, n1/2)
            n_ctrl_lsq(2) = min(16, n2/2)
            n_ctrl_lsq(3) = min(16, n3/2)
            n_ctrl_lsq(1) = max(degree(1) + 1, n_ctrl_lsq(1))
            n_ctrl_lsq(2) = max(degree(2) + 1, n_ctrl_lsq(2))
            n_ctrl_lsq(3) = max(degree(3) + 1, n_ctrl_lsq(3))

            allocate(coeff_lsq(n_ctrl_lsq(1), n_ctrl_lsq(2), n_ctrl_lsq(3)))
            coeff_lsq = 0.0_dp_bench

            call bench_time_now(t0)
            call bspline_3d_init_uniform(spl_bs, degree, n_ctrl_lsq, x_min, x_max)
            call bspline_3d_lsq_cgls(spl_bs, x1, x2, x3, f3d, coeff_lsq, &
                max_iter=800, tol=1.0d-10)
            call bench_time_now(t1)
            t_create_lsq = t1 - t0

            call bench_time_now(t0)
            do i3 = 1, n3
                do i2 = 1, n2
                    do i1 = 1, n1
                        x = [x1(i1), x2(i2), x3(i3)]
                        call bspline_3d_eval(spl_bs, coeff_lsq, x, y)
                    end do
                end do
            end do
            call bench_time_now(t1)
            t_eval_lsq = t1 - t0

            ! 3) neo_bspline direct interp: construct + evaluate (skip large N)
            allocate(coeff_dir(n1, n2, n3))

            call bench_time_now(t0)
            call bspline_3d_init_uniform(spl_bs_dir, degree, [n1, n2, n3], &
                x_min, x_max)
            call bspline_3d_interp(spl_bs_dir, x1, x2, x3, f3d, coeff_dir)
            call bench_time_now(t1)
            t_create_dir = t1 - t0

            call bench_time_now(t0)
            do i3 = 1, n3
                do i2 = 1, n2
                    do i1 = 1, n1
                        x = [x1(i1), x2(i2), x3(i3)]
                        call bspline_3d_eval(spl_bs_dir, coeff_dir, x, y)
                    end do
                end do
            end do
            call bench_time_now(t1)
            t_eval_dir = t1 - t0

            deallocate(coeff_dir)

            write(unit,'(i10,6es16.8)') n_data, t_create_interp, t_eval_interp, &
                t_create_lsq, t_eval_lsq, t_create_dir, t_eval_dir

            deallocate(x1, x2, x3, f3d, coeff_lsq)
            if (n_data <= 10000) call destroy_splines_3d(spl_interp)
        end do

        close(unit)
    end subroutine run_bench_3d

end program bench_neo_bspline_vs_interpolate
