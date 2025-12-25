program benchmark_spline_performance
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_util, only: linspace
    use interpolate, only: SplineData1D, SplineData2D, SplineData3D, &
                           construct_splines_1d, construct_splines_2d, &
                           construct_splines_3d, &
                           destroy_splines_1d, destroy_splines_2d, &
                           destroy_splines_3d, &
                           evaluate_splines_1d, evaluate_splines_1d_der, &
                           evaluate_splines_2d, evaluate_splines_2d_der, &
                           evaluate_splines_3d, evaluate_splines_3d_der

    implicit none

    call run_1d_benchmark()
    call run_2d_benchmark()
    call run_3d_benchmark()

contains

    subroutine run_1d_benchmark()
        integer, parameter :: n_grid = 100, n_eval = 10000, n_warmup = 2, n_repeat = 5
        real(dp), parameter :: x_min = 0.0d0, x_max = 6.28d0
        real(dp) :: x_grid(n_grid), y_grid(n_grid)
        real(dp) :: t_start, t_end, t_construct, t_eval
        real(dp) :: x_eval, y, dy, dummy
        type(SplineData1D) :: spl
        integer :: i, ie, rep

        call linspace(x_min, x_max, n_grid, x_grid)
        do i = 1, n_grid
            y_grid(i) = sin(x_grid(i))
        end do

        do rep = 1, n_warmup
            call construct_splines_1d(x_min, x_max, y_grid, 5, .false., spl)
            call destroy_splines_1d(spl)
        end do

        call cpu_time(t_start)
        do rep = 1, n_repeat
            call construct_splines_1d(x_min, x_max, y_grid, 5, .false., spl)
            call destroy_splines_1d(spl)
        end do
        call cpu_time(t_end)
        t_construct = (t_end - t_start) / n_repeat

        call construct_splines_1d(x_min, x_max, y_grid, 5, .false., spl)

        dummy = 0.0d0
        do rep = 1, n_warmup
            do ie = 1, n_eval
                x_eval = x_min + (x_max - x_min) * dble(ie) / dble(n_eval + 1)
                call evaluate_splines_1d_der(spl, x_eval, y, dy)
                dummy = dummy + y
            end do
        end do

        call cpu_time(t_start)
        do rep = 1, n_repeat
            do ie = 1, n_eval
                x_eval = x_min + (x_max - x_min) * dble(ie) / dble(n_eval + 1)
                call evaluate_splines_1d_der(spl, x_eval, y, dy)
                dummy = dummy + y
            end do
        end do
        call cpu_time(t_end)
        t_eval = (t_end - t_start) / n_repeat

        call destroy_splines_1d(spl)

        print '(A)', "1D SPLINE BENCHMARK"
        print '(A,I6,A,I7,A)', "  Grid: ", n_grid, " points, Eval: ", n_eval, " points"
        print '(A,ES12.4,A)', "  Construction: ", t_construct * 1.0d6, " us"
        print '(A,ES12.4,A)', "  Evaluation:   ", t_eval * 1.0d6, " us"
        print '(A,ES12.4,A)', "  Per-point:    ", t_eval / n_eval * 1.0d9, " ns"
        if (dummy == 0.0d0) print *, ""
    end subroutine run_1d_benchmark

    subroutine run_2d_benchmark()
        integer, parameter :: n1 = 50, n2 = 50, n_eval = 5000, n_warmup = 2, n_repeat = 5
        real(dp), parameter :: x_min(2) = [0.0d0, 0.0d0], x_max(2) = [6.28d0, 6.28d0]
        real(dp) :: x1(n1), x2(n2), y_grid(n1, n2)
        real(dp) :: t_start, t_end, t_construct, t_eval
        real(dp) :: x_eval(2), y, dy(2), dummy
        type(SplineData2D) :: spl
        integer :: i1, i2, ie, rep
        integer, parameter :: order(2) = [5, 5]
        logical, parameter :: periodic(2) = [.false., .false.]

        call linspace(x_min(1), x_max(1), n1, x1)
        call linspace(x_min(2), x_max(2), n2, x2)
        do i2 = 1, n2
            do i1 = 1, n1
                y_grid(i1, i2) = sin(x1(i1)) * cos(x2(i2))
            end do
        end do

        do rep = 1, n_warmup
            call construct_splines_2d(x_min, x_max, y_grid, order, periodic, spl)
            call destroy_splines_2d(spl)
        end do

        call cpu_time(t_start)
        do rep = 1, n_repeat
            call construct_splines_2d(x_min, x_max, y_grid, order, periodic, spl)
            call destroy_splines_2d(spl)
        end do
        call cpu_time(t_end)
        t_construct = (t_end - t_start) / n_repeat

        call construct_splines_2d(x_min, x_max, y_grid, order, periodic, spl)

        dummy = 0.0d0
        do rep = 1, n_warmup
            do ie = 1, n_eval
                x_eval(1) = x_min(1) + (x_max(1) - x_min(1)) * dble(ie) / dble(n_eval + 1)
                x_eval(2) = x_min(2) + (x_max(2) - x_min(2)) * dble(mod(ie*7, n_eval)) / dble(n_eval + 1)
                call evaluate_splines_2d_der(spl, x_eval, y, dy)
                dummy = dummy + y
            end do
        end do

        call cpu_time(t_start)
        do rep = 1, n_repeat
            do ie = 1, n_eval
                x_eval(1) = x_min(1) + (x_max(1) - x_min(1)) * dble(ie) / dble(n_eval + 1)
                x_eval(2) = x_min(2) + (x_max(2) - x_min(2)) * dble(mod(ie*7, n_eval)) / dble(n_eval + 1)
                call evaluate_splines_2d_der(spl, x_eval, y, dy)
                dummy = dummy + y
            end do
        end do
        call cpu_time(t_end)
        t_eval = (t_end - t_start) / n_repeat

        call destroy_splines_2d(spl)

        print '(A)', "2D SPLINE BENCHMARK"
        print '(A,I4,A,I4,A,I6,A)', "  Grid: ", n1, "x", n2, ", Eval: ", n_eval, " points"
        print '(A,ES12.4,A)', "  Construction: ", t_construct * 1.0d6, " us"
        print '(A,ES12.4,A)', "  Evaluation:   ", t_eval * 1.0d6, " us"
        print '(A,ES12.4,A)', "  Per-point:    ", t_eval / n_eval * 1.0d9, " ns"
        if (dummy == 0.0d0) print *, ""
    end subroutine run_2d_benchmark

    subroutine run_3d_benchmark()
        integer, parameter :: n1 = 20, n2 = 20, n3 = 20, n_eval = 2000
        integer, parameter :: n_warmup = 2, n_repeat = 5
        real(dp), parameter :: x_min(3) = [0.0d0, 0.0d0, 0.0d0]
        real(dp), parameter :: x_max(3) = [6.28d0, 6.28d0, 6.28d0]
        real(dp) :: x1(n1), x2(n2), x3(n3), y_grid(n1, n2, n3)
        real(dp) :: t_start, t_end, t_construct, t_eval
        real(dp) :: x_eval(3), y, dy(3), dummy
        type(SplineData3D) :: spl
        integer :: i1, i2, i3, ie, rep
        integer, parameter :: order(3) = [5, 5, 5]
        logical, parameter :: periodic(3) = [.false., .false., .false.]

        call linspace(x_min(1), x_max(1), n1, x1)
        call linspace(x_min(2), x_max(2), n2, x2)
        call linspace(x_min(3), x_max(3), n3, x3)
        do i3 = 1, n3
            do i2 = 1, n2
                do i1 = 1, n1
                    y_grid(i1, i2, i3) = sin(x1(i1)) * cos(x2(i2)) * sin(x3(i3))
                end do
            end do
        end do

        do rep = 1, n_warmup
            call construct_splines_3d(x_min, x_max, y_grid, order, periodic, spl)
            call destroy_splines_3d(spl)
        end do

        call cpu_time(t_start)
        do rep = 1, n_repeat
            call construct_splines_3d(x_min, x_max, y_grid, order, periodic, spl)
            call destroy_splines_3d(spl)
        end do
        call cpu_time(t_end)
        t_construct = (t_end - t_start) / n_repeat

        call construct_splines_3d(x_min, x_max, y_grid, order, periodic, spl)

        dummy = 0.0d0
        do rep = 1, n_warmup
            do ie = 1, n_eval
                x_eval(1) = x_min(1) + (x_max(1) - x_min(1)) * dble(ie) / dble(n_eval + 1)
                x_eval(2) = x_min(2) + (x_max(2) - x_min(2)) * dble(mod(ie*7, n_eval)) / dble(n_eval + 1)
                x_eval(3) = x_min(3) + (x_max(3) - x_min(3)) * dble(mod(ie*13, n_eval)) / dble(n_eval + 1)
                call evaluate_splines_3d_der(spl, x_eval, y, dy)
                dummy = dummy + y
            end do
        end do

        call cpu_time(t_start)
        do rep = 1, n_repeat
            do ie = 1, n_eval
                x_eval(1) = x_min(1) + (x_max(1) - x_min(1)) * dble(ie) / dble(n_eval + 1)
                x_eval(2) = x_min(2) + (x_max(2) - x_min(2)) * dble(mod(ie*7, n_eval)) / dble(n_eval + 1)
                x_eval(3) = x_min(3) + (x_max(3) - x_min(3)) * dble(mod(ie*13, n_eval)) / dble(n_eval + 1)
                call evaluate_splines_3d_der(spl, x_eval, y, dy)
                dummy = dummy + y
            end do
        end do
        call cpu_time(t_end)
        t_eval = (t_end - t_start) / n_repeat

        call destroy_splines_3d(spl)

        print '(A)', "3D SPLINE BENCHMARK"
        print '(A,I3,A,I3,A,I3,A,I5,A)', "  Grid: ", n1, "x", n2, "x", n3, ", Eval: ", n_eval, " points"
        print '(A,ES12.4,A)', "  Construction: ", t_construct * 1.0d6, " us"
        print '(A,ES12.4,A)', "  Evaluation:   ", t_eval * 1.0d6, " us"
        print '(A,ES12.4,A)', "  Per-point:    ", t_eval / n_eval * 1.0d9, " ns"
        if (dummy == 0.0d0) print *, ""
    end subroutine run_3d_benchmark

end program benchmark_spline_performance
