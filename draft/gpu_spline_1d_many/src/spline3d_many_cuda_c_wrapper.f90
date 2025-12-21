module spline3d_many_cuda_c_wrapper
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_null_ptr, c_loc
    implicit none
    private

    public :: cuda3d_c_state_t
    public :: spline3d_many_cuda_c_init
    public :: spline3d_many_cuda_c_free
    public :: spline3d_many_cuda_c_set_x
    public :: spline3d_many_cuda_c_eval_device
    public :: spline3d_many_cuda_c_get_y

    type :: cuda3d_c_state_t
        type(c_ptr) :: handle = c_null_ptr
    end type cuda3d_c_state_t

    interface
        subroutine spline3d_many_cuda_c_init_c(order1, order2, order3, n1, n2, n3, nq, &
                                               periodic1, periodic2, periodic3, &
                                               x_min1, x_min2, &
                                               x_min3, h1, h2, h3, coeff, handle) &
            bind(C, name="spline3d_many_cuda_c_init")
            import :: c_int, c_double, c_ptr
            integer(c_int), value :: order1, order2, order3, n1, n2, n3, nq
            integer(c_int), value :: periodic1, periodic2, periodic3
            real(c_double), value :: x_min1, x_min2, x_min3, h1, h2, h3
            type(c_ptr), value :: coeff
            type(c_ptr) :: handle
        end subroutine spline3d_many_cuda_c_init_c

        subroutine spline3d_many_cuda_c_free_c(handle) &
            bind(C, name="spline3d_many_cuda_c_free")
            import :: c_ptr
            type(c_ptr), value :: handle
        end subroutine spline3d_many_cuda_c_free_c

        subroutine spline3d_many_cuda_c_set_x_c(handle, x, npts) &
            bind(C, name="spline3d_many_cuda_c_set_x")
            import :: c_ptr, c_int
            type(c_ptr), value :: handle
            type(c_ptr), value :: x
            integer(c_int), value :: npts
        end subroutine spline3d_many_cuda_c_set_x_c

        subroutine spline3d_many_cuda_c_eval_device_c(handle, npts) &
            bind(C, name="spline3d_many_cuda_c_eval_device")
            import :: c_ptr, c_int
            type(c_ptr), value :: handle
            integer(c_int), value :: npts
        end subroutine spline3d_many_cuda_c_eval_device_c

        subroutine spline3d_many_cuda_c_get_y_c(handle, y, npts) &
            bind(C, name="spline3d_many_cuda_c_get_y")
            import :: c_ptr, c_int
            type(c_ptr), value :: handle
            type(c_ptr), value :: y
            integer(c_int), value :: npts
        end subroutine spline3d_many_cuda_c_get_y_c
    end interface

contains

    subroutine spline3d_many_cuda_c_init(st, order, num_points, num_quantities, &
                                         periodic, x_min, &
                                         h_step, coeff)
        type(cuda3d_c_state_t), intent(inout) :: st
        integer, intent(in) :: order(3)
        integer, intent(in) :: num_points(3)
        integer, intent(in) :: num_quantities
        logical, intent(in) :: periodic(3)
        real(dp), intent(in) :: x_min(3)
        real(dp), intent(in) :: h_step(3)
        real(dp), intent(in), target :: coeff(num_quantities, 0:order(1), 0:order(2), &
                                              0:order(3), &
                                              num_points(1), num_points(2), &
                                              num_points(3))

        integer(c_int) :: p1, p2, p3

        p1 = 0_c_int
        p2 = 0_c_int
        p3 = 0_c_int
        if (periodic(1)) p1 = 1_c_int
        if (periodic(2)) p2 = 1_c_int
        if (periodic(3)) p3 = 1_c_int

        call spline3d_many_cuda_c_init_c( &
            int(order(1), c_int), int(order(2), c_int), int(order(3), c_int), &
            int(num_points(1), c_int), &
            int(num_points(2), c_int), int(num_points(3), c_int), int(num_quantities, &
                                                                      c_int), &
            p1, p2, &
            p3, &
            real(x_min(1), c_double), real(x_min(2), c_double), real(x_min(3), &
                                                                     c_double), &
            real(h_step(1), c_double), real(h_step(2), c_double), real(h_step(3), &
                                                                       c_double), &
            c_loc(coeff(1, 0, 0, 0, 1, 1, 1)), st%handle &
            )
    end subroutine spline3d_many_cuda_c_init

    subroutine spline3d_many_cuda_c_free(st)
        type(cuda3d_c_state_t), intent(inout) :: st

        if (st%handle /= c_null_ptr) then
            call spline3d_many_cuda_c_free_c(st%handle)
            st%handle = c_null_ptr
        end if
    end subroutine spline3d_many_cuda_c_free

    subroutine spline3d_many_cuda_c_set_x(st, x)
        type(cuda3d_c_state_t), intent(in) :: st
        real(dp), intent(in), target :: x(:, :)
        integer(c_int) :: npts_c

        npts_c = int(size(x, 2), c_int)
        call spline3d_many_cuda_c_set_x_c(st%handle, c_loc(x(1, 1)), npts_c)
    end subroutine spline3d_many_cuda_c_set_x

    subroutine spline3d_many_cuda_c_eval_device(st, npts)
        type(cuda3d_c_state_t), intent(in) :: st
        integer, intent(in) :: npts

        call spline3d_many_cuda_c_eval_device_c(st%handle, int(npts, c_int))
    end subroutine spline3d_many_cuda_c_eval_device

    subroutine spline3d_many_cuda_c_get_y(st, y, npts)
        type(cuda3d_c_state_t), intent(in) :: st
        real(dp), intent(out), target :: y(:)
        integer, intent(in) :: npts

        call spline3d_many_cuda_c_get_y_c(st%handle, c_loc(y(1)), int(npts, c_int))
    end subroutine spline3d_many_cuda_c_get_y

end module spline3d_many_cuda_c_wrapper
