module draft_batch_splines_many_api
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: iso_c_binding, only: c_ptr, c_loc, c_f_pointer
    use batch_interpolate_types, only: BatchSplineData1D, BatchSplineData2D, &
                                       BatchSplineData3D
    use spline1d_many_openacc, only: spline1d_many_openacc_eval_host
    use spline2d_many_openacc, only: spline2d_many_openacc_eval_host
    use spline3d_many_openacc, only: spline3d_many_openacc_eval_host
    implicit none
    private

    public :: evaluate_batch_splines_1d_many
    public :: evaluate_batch_splines_2d_many
    public :: evaluate_batch_splines_3d_many

contains

    subroutine evaluate_batch_splines_1d_many(spl, x, y)
        type(BatchSplineData1D), intent(in) :: spl
        real(dp), intent(in) :: x(:)
        real(dp), intent(out), target :: y(:, :)

        type(c_ptr) :: y_ptr
        real(dp), pointer :: y_flat(:)
        integer :: npts, nq

        npts = size(x)
        nq = spl%num_quantities
        if (size(y, 1) /= nq) error stop
        if (size(y, 2) /= npts) error stop

        y_ptr = c_loc(y(1, 1))
        call c_f_pointer(y_ptr, y_flat, [nq*npts])

        call spline1d_many_openacc_eval_host(spl%order, spl%num_points, &
                                             spl%num_quantities, &
                                             spl%periodic, spl%x_min, spl%h_step, &
                                             spl%coeff, x, &
                                             y_flat)
    end subroutine evaluate_batch_splines_1d_many

    subroutine evaluate_batch_splines_2d_many(spl, x, y)
        type(BatchSplineData2D), intent(in) :: spl
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out), target :: y(:, :)

        type(c_ptr) :: y_ptr
        real(dp), pointer :: y_flat(:)
        integer :: npts, nq

        if (size(x, 1) /= 2) error stop
        npts = size(x, 2)
        nq = spl%num_quantities
        if (size(y, 1) /= nq) error stop
        if (size(y, 2) /= npts) error stop

        y_ptr = c_loc(y(1, 1))
        call c_f_pointer(y_ptr, y_flat, [nq*npts])

        call spline2d_many_openacc_eval_host(spl%order, spl%num_points, &
                                             spl%num_quantities, &
                                             spl%periodic, spl%x_min, spl%h_step, &
                                             spl%coeff, x, &
                                             y_flat)
    end subroutine evaluate_batch_splines_2d_many

    subroutine evaluate_batch_splines_3d_many(spl, x, y)
        type(BatchSplineData3D), intent(in) :: spl
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out), target :: y(:, :)

        type(c_ptr) :: y_ptr
        real(dp), pointer :: y_flat(:)
        integer :: npts, nq

        if (size(x, 1) /= 3) error stop
        npts = size(x, 2)
        nq = spl%num_quantities
        if (size(y, 1) /= nq) error stop
        if (size(y, 2) /= npts) error stop

        y_ptr = c_loc(y(1, 1))
        call c_f_pointer(y_ptr, y_flat, [nq*npts])

        call spline3d_many_openacc_eval_host(spl%order, spl%num_points, &
                                             spl%num_quantities, &
                                             spl%periodic, spl%x_min, spl%h_step, &
                                             spl%coeff, x, &
                                             y_flat)
    end subroutine evaluate_batch_splines_3d_many

end module draft_batch_splines_many_api
