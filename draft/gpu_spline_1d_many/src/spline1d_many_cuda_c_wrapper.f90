module spline1d_many_cuda_c_wrapper
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_null_ptr, c_loc
    implicit none
    private

    public :: cuda_c_state_t
    public :: spline1d_many_cuda_c_init
    public :: spline1d_many_cuda_c_free
    public :: spline1d_many_cuda_c_eval
    public :: spline1d_many_cuda_c_set_x
    public :: spline1d_many_cuda_c_eval_device
    public :: spline1d_many_cuda_c_get_y
    public :: spline1d_many_cuda_sync

    type :: cuda_c_state_t
        type(c_ptr) :: handle = c_null_ptr
    end type cuda_c_state_t

    interface
        subroutine spline1d_many_cuda_c_init_c(order, num_points, num_quantities, &
                                               periodic, &
                                               x_min, h_step, coeff, handle) &
            bind(C, name="spline1d_many_cuda_c_init")
            import :: c_int, c_double, c_ptr
            integer(c_int), value :: order
            integer(c_int), value :: num_points
            integer(c_int), value :: num_quantities
            integer(c_int), value :: periodic
            real(c_double), value :: x_min
            real(c_double), value :: h_step
            type(c_ptr), value :: coeff
            type(c_ptr) :: handle
        end subroutine spline1d_many_cuda_c_init_c

        subroutine spline1d_many_cuda_c_free_c(handle) &
            bind(C, name="spline1d_many_cuda_c_free")
            import :: c_ptr
            type(c_ptr), value :: handle
        end subroutine spline1d_many_cuda_c_free_c

        subroutine spline1d_many_cuda_c_eval_c(handle, x, y, npts) &
            bind(C, name="spline1d_many_cuda_c_eval")
            import :: c_ptr, c_int
            type(c_ptr), value :: handle
            type(c_ptr), value :: x
            type(c_ptr), value :: y
            integer(c_int), value :: npts
        end subroutine spline1d_many_cuda_c_eval_c

        subroutine spline1d_many_cuda_c_set_x_c(handle, x, npts) &
            bind(C, name="spline1d_many_cuda_c_set_x")
            import :: c_ptr, c_int
            type(c_ptr), value :: handle
            type(c_ptr), value :: x
            integer(c_int), value :: npts
        end subroutine spline1d_many_cuda_c_set_x_c

        subroutine spline1d_many_cuda_c_eval_device_c(handle, npts) &
            bind(C, name="spline1d_many_cuda_c_eval_device")
            import :: c_ptr, c_int
            type(c_ptr), value :: handle
            integer(c_int), value :: npts
        end subroutine spline1d_many_cuda_c_eval_device_c

        subroutine spline1d_many_cuda_c_get_y_c(handle, y, npts) &
            bind(C, name="spline1d_many_cuda_c_get_y")
            import :: c_ptr, c_int
            type(c_ptr), value :: handle
            type(c_ptr), value :: y
            integer(c_int), value :: npts
        end subroutine spline1d_many_cuda_c_get_y_c

        subroutine spline1d_many_cuda_sync_c() bind(C, name="spline1d_many_cuda_sync")
        end subroutine spline1d_many_cuda_sync_c
    end interface

contains

    subroutine spline1d_many_cuda_c_init(st, order, num_points, num_quantities, &
                                         periodic, &
                                         x_min, h_step, coeff)
        type(cuda_c_state_t), intent(inout) :: st
        integer, intent(in) :: order
        integer, intent(in) :: num_points
        integer, intent(in) :: num_quantities
        logical, intent(in) :: periodic
        real(dp), intent(in) :: x_min
        real(dp), intent(in) :: h_step
        real(dp), intent(in), target :: coeff(num_quantities, 0:order, num_points)

        integer(c_int) :: periodic_int

        if (periodic) then
            periodic_int = 1_c_int
        else
            periodic_int = 0_c_int
        end if

        call spline1d_many_cuda_c_init_c( &
            int(order, c_int), int(num_points, c_int), int(num_quantities, c_int), &
            periodic_int, real(x_min, c_double), real(h_step, c_double), &
            c_loc(coeff(1, 0, 1)), st%handle &
            )
    end subroutine spline1d_many_cuda_c_init

    subroutine spline1d_many_cuda_c_free(st)
        type(cuda_c_state_t), intent(inout) :: st

        if (st%handle /= c_null_ptr) then
            call spline1d_many_cuda_c_free_c(st%handle)
            st%handle = c_null_ptr
        end if
    end subroutine spline1d_many_cuda_c_free

    subroutine spline1d_many_cuda_c_eval(st, x, y)
        type(cuda_c_state_t), intent(in) :: st
        real(dp), intent(in), target :: x(:)
        real(dp), intent(out), target :: y(:)

        call spline1d_many_cuda_c_eval_c(st%handle, c_loc(x(1)), c_loc(y(1)), &
                                         int(size(x), c_int))
    end subroutine spline1d_many_cuda_c_eval

    subroutine spline1d_many_cuda_c_set_x(st, x)
        type(cuda_c_state_t), intent(in) :: st
        real(dp), intent(in), target :: x(:)

        call spline1d_many_cuda_c_set_x_c(st%handle, c_loc(x(1)), int(size(x), c_int))
    end subroutine spline1d_many_cuda_c_set_x

    subroutine spline1d_many_cuda_c_eval_device(st, npts)
        type(cuda_c_state_t), intent(in) :: st
        integer, intent(in) :: npts

        call spline1d_many_cuda_c_eval_device_c(st%handle, int(npts, c_int))
    end subroutine spline1d_many_cuda_c_eval_device

    subroutine spline1d_many_cuda_c_get_y(st, y, npts)
        type(cuda_c_state_t), intent(in) :: st
        real(dp), intent(out), target :: y(:)
        integer, intent(in) :: npts

        call spline1d_many_cuda_c_get_y_c(st%handle, c_loc(y(1)), int(npts, c_int))
    end subroutine spline1d_many_cuda_c_get_y

    subroutine spline1d_many_cuda_sync()
        call spline1d_many_cuda_sync_c()
    end subroutine spline1d_many_cuda_sync

end module spline1d_many_cuda_c_wrapper
