program test_spline_performance
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_util, only: linspace
    use interpolate, only: SplineData1D, SplineData2D, SplineData3D, &
                           construct_splines_1d, construct_splines_2d, &
                           construct_splines_3d, &
                           destroy_splines_1d, destroy_splines_2d, &
                           destroy_splines_3d, &
                           evaluate_splines_1d_der, &
                           evaluate_splines_2d_der, &
                           evaluate_splines_3d_der

    implicit none

    logical :: all_passed

    all_passed = .true.
    call test_1d_performance(all_passed)
    call test_2d_performance(all_passed)
    call test_3d_performance(all_passed)

    if (all_passed) then
        print *, "All spline performance tests PASSED!"
    else
        error stop "Spline performance regression detected"
    end if

contains

    subroutine test_1d_performance(all_passed)
        logical, intent(inout) :: all_passed
        integer, parameter :: n_grid = 100, n_eval = 5000, n_repeat = 3
        real(dp), parameter :: max_ns_per_eval = 500.0d0
        real(dp), parameter :: x_min = 0.0d0, x_max = 6.28d0
        real(dp) :: x_grid(n_grid), y_grid(n_grid)
        real(dp) :: t_start, t_end, t_eval, ns_per_eval
        real(dp) :: x_eval, y, dy, dummy
        type(SplineData1D) :: spl
        integer :: i, ie, rep

        call linspace(x_min, x_max, n_grid, x_grid)
        do i = 1, n_grid
            y_grid(i) = sin(x_grid(i))
        end do

        call construct_splines_1d(x_min, x_max, y_grid, 5, .false., spl)

        dummy = 0.0d0
        do rep = 1, 2
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

        ns_per_eval = t_eval / n_eval * 1.0d9
        print '(A,F8.1,A,F8.1,A)', "  1D: ", ns_per_eval, " ns/eval (limit: ", &
            max_ns_per_eval, " ns)"

        if (ns_per_eval > max_ns_per_eval) then
            print *, "  FAILED: 1D evaluation too slow"
            all_passed = .false.
        else
            print *, "  PASSED"
        end if
        if (dummy == 0.0d0) print *, ""
    end subroutine test_1d_performance

    subroutine test_2d_performance(all_passed)
        logical, intent(inout) :: all_passed
        integer, parameter :: n1 = 40, n2 = 40, n_eval = 2000, n_repeat = 3
        real(dp), parameter :: max_ns_per_eval = 2000.0d0
        real(dp), parameter :: x_min(2) = [0.0d0, 0.0d0], x_max(2) = [6.28d0, 6.28d0]
        real(dp) :: x1(n1), x2(n2), y_grid(n1, n2)
        real(dp) :: t_start, t_end, t_eval, ns_per_eval
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

        call construct_splines_2d(x_min, x_max, y_grid, order, periodic, spl)

        dummy = 0.0d0
        do rep = 1, 2
            do ie = 1, n_eval
                x_eval(1) = x_min(1) + (x_max(1) - x_min(1)) * dble(ie) / dble(n_eval + 1)
                x_eval(2) = x_min(2) + (x_max(2) - x_min(2)) * dble(mod(ie*7, n_eval)) / &
                    dble(n_eval + 1)
                call evaluate_splines_2d_der(spl, x_eval, y, dy)
                dummy = dummy + y
            end do
        end do

        call cpu_time(t_start)
        do rep = 1, n_repeat
            do ie = 1, n_eval
                x_eval(1) = x_min(1) + (x_max(1) - x_min(1)) * dble(ie) / dble(n_eval + 1)
                x_eval(2) = x_min(2) + (x_max(2) - x_min(2)) * dble(mod(ie*7, n_eval)) / &
                    dble(n_eval + 1)
                call evaluate_splines_2d_der(spl, x_eval, y, dy)
                dummy = dummy + y
            end do
        end do
        call cpu_time(t_end)
        t_eval = (t_end - t_start) / n_repeat

        call destroy_splines_2d(spl)

        ns_per_eval = t_eval / n_eval * 1.0d9
        print '(A,F8.1,A,F8.1,A)', "  2D: ", ns_per_eval, " ns/eval (limit: ", &
            max_ns_per_eval, " ns)"

        if (ns_per_eval > max_ns_per_eval) then
            print *, "  FAILED: 2D evaluation too slow"
            all_passed = .false.
        else
            print *, "  PASSED"
        end if
        if (dummy == 0.0d0) print *, ""
    end subroutine test_2d_performance

    subroutine test_3d_performance(all_passed)
        logical, intent(inout) :: all_passed
        integer, parameter :: n1 = 15, n2 = 15, n3 = 15, n_eval = 1000, n_repeat = 3
        real(dp), parameter :: max_ns_per_eval = 5000.0d0
        real(dp), parameter :: x_min(3) = [0.0d0, 0.0d0, 0.0d0]
        real(dp), parameter :: x_max(3) = [6.28d0, 6.28d0, 6.28d0]
        real(dp) :: x1(n1), x2(n2), x3(n3), y_grid(n1, n2, n3)
        real(dp) :: t_start, t_end, t_eval, ns_per_eval
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

        call construct_splines_3d(x_min, x_max, y_grid, order, periodic, spl)

        dummy = 0.0d0
        do rep = 1, 2
            do ie = 1, n_eval
                x_eval(1) = x_min(1) + (x_max(1) - x_min(1)) * dble(ie) / dble(n_eval + 1)
                x_eval(2) = x_min(2) + (x_max(2) - x_min(2)) * dble(mod(ie*7, n_eval)) / &
                    dble(n_eval + 1)
                x_eval(3) = x_min(3) + (x_max(3) - x_min(3)) * dble(mod(ie*13, n_eval)) / &
                    dble(n_eval + 1)
                call evaluate_splines_3d_der(spl, x_eval, y, dy)
                dummy = dummy + y
            end do
        end do

        call cpu_time(t_start)
        do rep = 1, n_repeat
            do ie = 1, n_eval
                x_eval(1) = x_min(1) + (x_max(1) - x_min(1)) * dble(ie) / dble(n_eval + 1)
                x_eval(2) = x_min(2) + (x_max(2) - x_min(2)) * dble(mod(ie*7, n_eval)) / &
                    dble(n_eval + 1)
                x_eval(3) = x_min(3) + (x_max(3) - x_min(3)) * dble(mod(ie*13, n_eval)) / &
                    dble(n_eval + 1)
                call evaluate_splines_3d_der(spl, x_eval, y, dy)
                dummy = dummy + y
            end do
        end do
        call cpu_time(t_end)
        t_eval = (t_end - t_start) / n_repeat

        call destroy_splines_3d(spl)

        ns_per_eval = t_eval / n_eval * 1.0d9
        print '(A,F8.1,A,F8.1,A)', "  3D: ", ns_per_eval, " ns/eval (limit: ", &
            max_ns_per_eval, " ns)"

        if (ns_per_eval > max_ns_per_eval) then
            print *, "  FAILED: 3D evaluation too slow"
            all_passed = .false.
        else
            print *, "  PASSED"
        end if
        if (dummy == 0.0d0) print *, ""
    end subroutine test_3d_performance

end program test_spline_performance
