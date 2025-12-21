module spline2d_many_cuda_c_wrapper
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_null_ptr, c_loc
    implicit none
    private

    public :: cuda2d_c_state_t
    public :: spline2d_many_cuda_c_init
    public :: spline2d_many_cuda_c_free
    public :: spline2d_many_cuda_c_set_x
    public :: spline2d_many_cuda_c_eval_device
    public :: spline2d_many_cuda_c_get_y

    type :: cuda2d_c_state_t
        type(c_ptr) :: handle = c_null_ptr
    end type cuda2d_c_state_t

    interface
        subroutine spline2d_many_cuda_c_init_c(order1, order2, n1, n2, nq, periodic1, &
                                               periodic2, &
                                               x_min1, x_min2, h_step1, h_step2, &
                                               coeff, handle) &
            bind(C, name="spline2d_many_cuda_c_init")
            import :: c_int, c_double, c_ptr
            integer(c_int), value :: order1, order2, n1, n2, nq
            integer(c_int), value :: periodic1, periodic2
            real(c_double), value :: x_min1, x_min2, h_step1, h_step2
            type(c_ptr), value :: coeff
            type(c_ptr) :: handle
        end subroutine spline2d_many_cuda_c_init_c

        subroutine spline2d_many_cuda_c_free_c(handle) &
            bind(C, name="spline2d_many_cuda_c_free")
            import :: c_ptr
            type(c_ptr), value :: handle
        end subroutine spline2d_many_cuda_c_free_c

        subroutine spline2d_many_cuda_c_set_x_c(handle, x, npts) &
            bind(C, name="spline2d_many_cuda_c_set_x")
            import :: c_ptr, c_int
            type(c_ptr), value :: handle
            type(c_ptr), value :: x
            integer(c_int), value :: npts
        end subroutine spline2d_many_cuda_c_set_x_c

        subroutine spline2d_many_cuda_c_eval_device_c(handle, npts) &
            bind(C, name="spline2d_many_cuda_c_eval_device")
            import :: c_ptr, c_int
            type(c_ptr), value :: handle
            integer(c_int), value :: npts
        end subroutine spline2d_many_cuda_c_eval_device_c

        subroutine spline2d_many_cuda_c_get_y_c(handle, y, npts) &
            bind(C, name="spline2d_many_cuda_c_get_y")
            import :: c_ptr, c_int
            type(c_ptr), value :: handle
            type(c_ptr), value :: y
            integer(c_int), value :: npts
        end subroutine spline2d_many_cuda_c_get_y_c
    end interface

contains

    subroutine spline2d_many_cuda_c_init(st, order, num_points, num_quantities, &
                                         periodic, x_min, &
                                         h_step, coeff)
        type(cuda2d_c_state_t), intent(inout) :: st
        integer, intent(in) :: order(2)
        integer, intent(in) :: num_points(2)
        integer, intent(in) :: num_quantities
        logical, intent(in) :: periodic(2)
        real(dp), intent(in) :: x_min(2)
        real(dp), intent(in) :: h_step(2)
        real(dp), intent(in), target :: coeff(num_quantities, 0:order(1), 0:order(2), &
                                              num_points(1), num_points(2))

        integer(c_int) :: p1, p2

        if (periodic(1)) then
            p1 = 1_c_int
        else
            p1 = 0_c_int
        end if
        if (periodic(2)) then
            p2 = 1_c_int
        else
            p2 = 0_c_int
        end if

        call spline2d_many_cuda_c_init_c( &
            int(order(1), c_int), int(order(2), c_int), int(num_points(1), c_int), &
            int(num_points(2), c_int), int(num_quantities, c_int), p1, p2, &
            real(x_min(1), c_double), real(x_min(2), c_double), real(h_step(1), &
                                                                     c_double), &
            real(h_step(2), c_double), c_loc(coeff(1, 0, 0, 1, 1)), st%handle &
            )
    end subroutine spline2d_many_cuda_c_init

    subroutine spline2d_many_cuda_c_free(st)
        type(cuda2d_c_state_t), intent(inout) :: st

        if (st%handle /= c_null_ptr) then
            call spline2d_many_cuda_c_free_c(st%handle)
            st%handle = c_null_ptr
        end if
    end subroutine spline2d_many_cuda_c_free

    subroutine spline2d_many_cuda_c_set_x(st, x)
        type(cuda2d_c_state_t), intent(in) :: st
        real(dp), intent(in), target :: x(:, :)

        call spline2d_many_cuda_c_set_x_c(st%handle, c_loc(x(1, 1)), &
                                          int(size(x, 2), c_int))
    end subroutine spline2d_many_cuda_c_set_x

    subroutine spline2d_many_cuda_c_eval_device(st, npts)
        type(cuda2d_c_state_t), intent(in) :: st
        integer, intent(in) :: npts

        call spline2d_many_cuda_c_eval_device_c(st%handle, int(npts, c_int))
    end subroutine spline2d_many_cuda_c_eval_device

    subroutine spline2d_many_cuda_c_get_y(st, y, npts)
        type(cuda2d_c_state_t), intent(in) :: st
        real(dp), intent(out), target :: y(:)
        integer, intent(in) :: npts

        call spline2d_many_cuda_c_get_y_c(st%handle, c_loc(y(1)), int(npts, c_int))
    end subroutine spline2d_many_cuda_c_get_y

end module spline2d_many_cuda_c_wrapper
