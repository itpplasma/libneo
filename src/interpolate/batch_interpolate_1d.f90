module batch_interpolate_1d
   use, intrinsic :: iso_fortran_env, only: dp => real64
   use batch_interpolate_types, only: BatchSplineData1D
   use spline_build_lines, only: spl_three_per_line, spl_three_reg_line, &
                                 spl_four_per_line, spl_four_reg_line, &
                                 spl_five_per_line, spl_five_reg_line
   use spl_three_to_five_sub, only: spl_per, spl_reg
#ifdef _OPENACC
   use openacc, only: acc_is_present
#endif

   implicit none
   private

   ! Export batch spline construction/destruction routines
   public :: construct_batch_splines_1d
   public :: construct_batch_splines_1d_lines
   public :: construct_batch_splines_1d_resident
   public :: construct_batch_splines_1d_resident_device
   public :: destroy_batch_splines_1d

#ifdef LIBNEO_ENABLE_SPLINE_ORACLE
   public :: construct_batch_splines_1d_legacy
#endif

   ! Export batch spline evaluation routines
   public :: evaluate_batch_splines_1d
   public :: evaluate_batch_splines_1d_single
   public :: evaluate_batch_splines_1d_many
   public :: evaluate_batch_splines_1d_many_resident
   public :: evaluate_batch_splines_1d_der
   public :: evaluate_batch_splines_1d_der2
   public :: evaluate_batch_splines_1d_der3
   public :: evaluate_batch_splines_1d_many_der
   public :: evaluate_batch_splines_1d_many_der2
   public :: evaluate_batch_splines_1d_many_der3

contains

   subroutine construct_batch_splines_1d(x_min, x_max, y_batch, order, periodic, spl)
      real(dp), intent(in) :: x_min, x_max
      real(dp), intent(in) :: y_batch(:, :)  ! (n_points, n_quantities)
      integer, intent(in) :: order
      logical, intent(in) :: periodic
      type(BatchSplineData1D), intent(out) :: spl

      call construct_batch_splines_1d_legacy(x_min, x_max, y_batch, order, periodic, spl)

#ifdef _OPENACC
      ! Map only the allocatable component, not the whole derived type
      ! This is compatible with both gfortran and nvfortran
      !$acc enter data copyin(spl%coeff)
#endif
   end subroutine construct_batch_splines_1d

   subroutine construct_batch_splines_1d_legacy(x_min, x_max, y_batch, order, periodic, spl)
      real(dp), intent(in) :: x_min, x_max
      real(dp), intent(in) :: y_batch(:, :)  ! (n_points, n_quantities)
      integer, intent(in) :: order
      logical, intent(in) :: periodic
      type(BatchSplineData1D), intent(out) :: spl

      integer :: iq, n_points, n_quantities, istat
      real(dp), dimension(:, :), allocatable :: splcoe_temp

      n_points = size(y_batch, 1)
      n_quantities = size(y_batch, 2)

      ! Validate input
      if (n_points < 2) then
         error stop "construct_batch_splines_1d: Need at least 2 points"
      end if
      if (n_quantities < 1) then
         error stop "construct_batch_splines_1d: Need at least 1 quantity"
      end if
      if (order < 3 .or. order > 5) then
         error stop "construct_batch_splines_1d: Order must be between 3 and 5"
      end if

      ! Store metadata
      spl%order = order
      spl%num_points = n_points
      spl%periodic = periodic
      spl%x_min = x_min
      spl%h_step = (x_max - x_min)/dble(n_points - 1)
      spl%inv_h_step = 1.0_dp/spl%h_step
      spl%period = spl%h_step*real(n_points - 1, dp)
      spl%num_quantities = n_quantities

      ! Allocate batch coefficients with cache-friendly layout
      allocate (spl%coeff(n_quantities, 0:order, n_points), stat=istat)
      if (istat /= 0) then
         error stop "construct_batch_splines_1d: Allocation failed for coeff"
      end if

      ! Temporary array for single spline construction
      allocate (splcoe_temp(0:order, n_points), stat=istat)
      if (istat /= 0) then
         error stop "construct_batch_splines_1d: Allocation failed for splcoe_temp"
      end if

      ! Process each quantity
      do iq = 1, n_quantities
         ! Set the data values
         splcoe_temp(0, :) = y_batch(:, iq)

         ! Construct spline coefficients
         if (periodic) then
            call spl_per(order, n_points, spl%h_step, splcoe_temp)
         else
            call spl_reg(order, n_points, spl%h_step, splcoe_temp)
         end if

         ! Store coefficients in batch array with quantities first
         spl%coeff(iq, :, :) = splcoe_temp(:, :)
      end do

      deallocate (splcoe_temp)
   end subroutine construct_batch_splines_1d_legacy

   subroutine construct_batch_splines_1d_lines(x_min, x_max, y_batch, order, periodic, spl)
      real(dp), intent(in) :: x_min, x_max
      real(dp), intent(in) :: y_batch(:, :)  ! (n_points, n_quantities)
      integer, intent(in) :: order
      logical, intent(in) :: periodic
      type(BatchSplineData1D), intent(out) :: spl

      integer :: iq, n_points, n_quantities, istat
      real(dp), allocatable :: work(:, :, :)

      n_points = size(y_batch, 1)
      n_quantities = size(y_batch, 2)

      if (n_points < 2) then
         error stop "construct_batch_splines_1d_lines: Need at least 2 points"
      end if
      if (n_quantities < 1) then
         error stop "construct_batch_splines_1d_lines: Need at least 1 quantity"
      end if
      if (order < 3 .or. order > 5) then
         error stop "construct_batch_splines_1d_lines: Order must be between 3 and 5"
      end if
      if (order == 3) then
         if (periodic .and. n_points < 3) then
            error stop "construct_batch_splines_1d_lines: Need at least 3 points "// &
               "for periodic order=3"
         end if
      else if (order == 4) then
         if (n_points < 5) then
            error stop "construct_batch_splines_1d_lines: Need at least 5 points for order=4"
         end if
      else
         if (n_points < 6) then
            error stop "construct_batch_splines_1d_lines: Need at least 6 points for order=5"
         end if
      end if

      spl%order = order
      spl%num_points = n_points
      spl%periodic = periodic
      spl%x_min = x_min
      spl%h_step = (x_max - x_min)/dble(n_points - 1)
      spl%inv_h_step = 1.0_dp/spl%h_step
      spl%period = spl%h_step*real(n_points - 1, dp)
      spl%num_quantities = n_quantities

      allocate (spl%coeff(n_quantities, 0:order, n_points), stat=istat)
      if (istat /= 0) then
         error stop "construct_batch_splines_1d_lines: Allocation failed for coeff"
      end if

      allocate (work(n_points, n_quantities, 0:order), stat=istat)
      if (istat /= 0) then
         error stop "construct_batch_splines_1d_lines: Allocation failed for work"
      end if

      do iq = 1, n_quantities
         work(:, iq, 0) = y_batch(:, iq)
      end do

#ifdef _OPENMP
!$omp parallel do default(none) shared(work, n_points, n_quantities, order, periodic, &
!$omp& spl) private(iq)
#endif
      do iq = 1, n_quantities
         if (order == 3) then
            if (periodic) then
               call spl_three_per_line(n_points, spl%h_step, work(:, iq, 0), &
                                       work(:, iq, 1), work(:, iq, 2), &
                                       work(:, iq, 3))
            else
               call spl_three_reg_line(n_points, spl%h_step, work(:, iq, 0), &
                                       work(:, iq, 1), work(:, iq, 2), &
                                       work(:, iq, 3))
            end if
         else if (order == 4) then
            if (periodic) then
               call spl_four_per_line(n_points, spl%h_step, work(:, iq, 0), &
                                      work(:, iq, 1), work(:, iq, 2), &
                                      work(:, iq, 3), work(:, iq, 4))
            else
               call spl_four_reg_line(n_points, spl%h_step, work(:, iq, 0), &
                                      work(:, iq, 1), work(:, iq, 2), &
                                      work(:, iq, 3), work(:, iq, 4))
            end if
         else
            if (periodic) then
               call spl_five_per_line(n_points, spl%h_step, work(:, iq, 0), &
                                      work(:, iq, 1), work(:, iq, 2), &
                                      work(:, iq, 3), work(:, iq, 4), &
                                      work(:, iq, 5))
            else
               call spl_five_reg_line(n_points, spl%h_step, work(:, iq, 0), &
                                      work(:, iq, 1), work(:, iq, 2), &
                                      work(:, iq, 3), work(:, iq, 4), &
                                      work(:, iq, 5))
            end if
         end if
      end do
#ifdef _OPENMP
!$omp end parallel do
#endif

      do iq = 1, n_quantities
         spl%coeff(iq, :, :) = transpose(work(:, iq, :))
      end do

      deallocate (work)
   end subroutine construct_batch_splines_1d_lines

   subroutine construct_batch_splines_1d_resident(x_min, x_max, y_batch, order, &
                                                  periodic, spl)
      real(dp), intent(in) :: x_min, x_max
      real(dp), intent(in) :: y_batch(:, :)  ! (n_points, n_quantities)
      integer, intent(in) :: order
      logical, intent(in) :: periodic
      type(BatchSplineData1D), intent(out) :: spl

      call construct_batch_splines_1d(x_min, x_max, y_batch, order, periodic, spl)
   end subroutine construct_batch_splines_1d_resident

   subroutine construct_batch_splines_1d_resident_device(x_min, x_max, y_batch, &
                                                         order, periodic, spl, &
                                                         update_host, &
                                                         assume_y_present)
      real(dp), intent(in) :: x_min, x_max
      real(dp), intent(in) :: y_batch(:, :)  ! (n_points, n_quantities)
      integer, intent(in) :: order
      logical, intent(in) :: periodic
      type(BatchSplineData1D), intent(out) :: spl
      logical, intent(in), optional :: update_host
      logical, intent(in), optional :: assume_y_present

      integer :: n_points, n_quantities
      logical :: do_update, do_assume_present

      do_update = .true.
      if (present(update_host)) do_update = update_host
      do_assume_present = .false.
      if (present(assume_y_present)) do_assume_present = assume_y_present

      n_points = size(y_batch, 1)
      n_quantities = size(y_batch, 2)

#ifdef _OPENACC
      call construct_batch_splines_1d_resident_device_impl(x_min, x_max, n_points, &
                                                           n_quantities, y_batch, &
                                                           order, periodic, spl, &
                                                           do_update, do_assume_present)
#else
      call construct_batch_splines_1d_resident(x_min, x_max, y_batch, order, &
                                               periodic, spl)
#endif
   end subroutine construct_batch_splines_1d_resident_device

   subroutine construct_batch_splines_1d_resident_device_impl(x_min, x_max, n_points, &
                                                              n_quantities, y_batch, &
                                                              order, periodic, spl, &
                                                              do_update, &
                                                              do_assume_present)
      real(dp), intent(in) :: x_min, x_max
      integer, intent(in) :: n_points, n_quantities
      real(dp), intent(in) :: y_batch(n_points, n_quantities)
      integer, intent(in) :: order
      logical, intent(in) :: periodic
      type(BatchSplineData1D), intent(out) :: spl
      logical, intent(in) :: do_update
      logical, intent(in) :: do_assume_present

      integer :: istat
      integer :: ip, iq, k
      real(dp), allocatable :: work(:, :, :)
      real(dp) :: h_step

      if (n_points < 2) then
         error stop "construct_batch_splines_1d_resident_device:"// &
            " Need at least 2 points"
      end if
      if (n_quantities < 1) then
         error stop "construct_batch_splines_1d_resident_device:"// &
            " Need at least 1 quantity"
      end if
      if (order < 3 .or. order > 5) then
         error stop "construct_batch_splines_1d_resident_device:"// &
            " Order must be between 3 and 5"
      end if
      if (order == 3) then
         if (periodic .and. n_points < 3) then
            error stop "construct_batch_splines_1d_resident_device:"// &
               " Need at least 3 points for periodic order=3"
         end if
      else if (order == 4) then
         if (n_points < 5) then
            error stop "construct_batch_splines_1d_resident_device:"// &
               " Need at least 5 points for order=4"
         end if
      else
         if (n_points < 6) then
            error stop "construct_batch_splines_1d_resident_device:"// &
               " Need at least 6 points for order=5"
         end if
      end if

      spl%order = order
      spl%num_points = n_points
      spl%periodic = periodic
      spl%x_min = x_min
      spl%h_step = (x_max - x_min)/dble(n_points - 1)
      spl%inv_h_step = 1.0_dp/spl%h_step
      spl%period = spl%h_step*real(n_points - 1, dp)
      spl%num_quantities = n_quantities

      allocate (spl%coeff(n_quantities, 0:order, n_points), stat=istat)
      if (istat /= 0) then
         error stop "construct_batch_splines_1d_resident_device:"// &
            " Allocation failed for coeff"
      end if

      !$acc enter data create(spl%coeff)

      allocate (work(n_points, n_quantities, 0:order), stat=istat)
      if (istat /= 0) then
         error stop "construct_batch_splines_1d_resident_device:"// &
            " Allocation failed for work"
      end if

      h_step = spl%h_step

#ifdef _OPENACC
      if (do_assume_present) then
         if (.not. acc_is_present(y_batch)) then
            error stop "construct_batch_splines_1d_resident_device:"// &
               " assume_y_present=T but y_batch is not present"
         end if
      end if
#endif
      !$acc data present_or_copyin(y_batch(1:n_points, 1:n_quantities)) create(work)
      !$acc parallel loop collapse(2) gang
      do iq = 1, n_quantities
         do ip = 1, n_points
            work(ip, iq, 0) = y_batch(ip, iq)
         end do
      end do

      !$acc parallel loop gang
      do iq = 1, n_quantities
         if (order == 3) then
            if (periodic) then
               call spl_three_per_line(n_points, h_step, work(:, iq, 0), &
                                       work(:, iq, 1), work(:, iq, 2), &
                                       work(:, iq, 3))
            else
               call spl_three_reg_line(n_points, h_step, work(:, iq, 0), &
                                       work(:, iq, 1), work(:, iq, 2), &
                                       work(:, iq, 3))
            end if
         else if (order == 4) then
            if (periodic) then
               call spl_four_per_line(n_points, h_step, work(:, iq, 0), &
                                      work(:, iq, 1), work(:, iq, 2), &
                                      work(:, iq, 3), work(:, iq, 4))
            else
               call spl_four_reg_line(n_points, h_step, work(:, iq, 0), &
                                      work(:, iq, 1), work(:, iq, 2), &
                                      work(:, iq, 3), work(:, iq, 4))
            end if
         else
            if (periodic) then
               call spl_five_per_line(n_points, h_step, work(:, iq, 0), &
                                      work(:, iq, 1), work(:, iq, 2), &
                                      work(:, iq, 3), work(:, iq, 4), &
                                      work(:, iq, 5))
            else
               call spl_five_reg_line(n_points, h_step, work(:, iq, 0), &
                                      work(:, iq, 1), work(:, iq, 2), &
                                      work(:, iq, 3), work(:, iq, 4), &
                                      work(:, iq, 5))
            end if
         end if
      end do

      call store_coeff_1d(n_points, n_quantities, order, work, spl%coeff)
      !$acc end data

      deallocate (work)
      if (do_update) then
         !$acc update self(spl%coeff(1:n_quantities, 0:order, 1:n_points))
      end if
   end subroutine construct_batch_splines_1d_resident_device_impl

   subroutine store_coeff_1d(n_points, n_quantities, order, work, coeff)
      integer, intent(in) :: n_points, n_quantities, order
      real(dp), intent(in) :: work(n_points, n_quantities, 0:order)
      real(dp), intent(out) :: coeff(n_quantities, 0:order, n_points)

      integer :: ip, iq, k

#ifdef _OPENACC
      !$acc parallel loop collapse(2) gang present(work, coeff)
#endif
      do ip = 1, n_points
         do k = 0, order
#ifdef _OPENACC
            !$acc loop vector
#endif
            do iq = 1, n_quantities
               coeff(iq, k, ip) = work(ip, iq, k)
            end do
         end do
      end do
   end subroutine store_coeff_1d

   subroutine destroy_batch_splines_1d(spl)
      type(BatchSplineData1D), intent(inout) :: spl

	#ifdef _OPENACC
	      if (allocated(spl%coeff)) then
	         if (acc_is_present(spl%coeff)) then
	            !$acc exit data delete(spl%coeff)
	            !$acc wait
	         end if
	      end if
	#endif
      if (allocated(spl%coeff)) deallocate (spl%coeff)
   end subroutine destroy_batch_splines_1d

   subroutine evaluate_batch_splines_1d(spl, x, y_batch)
      type(BatchSplineData1D), intent(in) :: spl
      real(dp), intent(in) :: x
      real(dp), intent(out) :: y_batch(:)  ! (n_quantities)

      real(dp) :: x_arr(1)
      real(dp) :: y_arr(spl%num_quantities, 1)

      x_arr(1) = x
      call evaluate_batch_splines_1d_many(spl, x_arr, y_arr)
      y_batch(1:spl%num_quantities) = y_arr(:, 1)
   end subroutine evaluate_batch_splines_1d

   subroutine evaluate_batch_splines_1d_single(spl, x, iq, y)
      type(BatchSplineData1D), intent(in) :: spl
      real(dp), intent(in) :: x
      integer, intent(in) :: iq  ! quantity index
      real(dp), intent(out) :: y

      real(dp) :: y_all(spl%num_quantities)

      if (iq < 1 .or. iq > spl%num_quantities) then
         error stop "evaluate_batch_splines_1d_single: Invalid quantity index"
      end if

      call evaluate_batch_splines_1d(spl, x, y_all)
      y = y_all(iq)
   end subroutine evaluate_batch_splines_1d_single

   subroutine evaluate_batch_splines_1d_many(spl, x, y_batch)
      type(BatchSplineData1D), intent(in) :: spl
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: y_batch(:, :)  ! (n_quantities, n_points)

      integer :: ipt, iq, k_power, idx, k_wrap
      integer :: num_points, nq, order
      logical :: periodic
      real(dp) :: xj, x_norm, x_local, x_min, h_step, period, t, w

      if (size(y_batch, 1) < spl%num_quantities) then
         error stop "evaluate_batch_splines_1d_many: First dimension too small"
      end if
      if (size(y_batch, 2) /= size(x)) then
         error stop "evaluate_batch_splines_1d_many: y_batch second dim mismatch"
      end if

#ifdef _OPENACC
      if (acc_is_present(spl%coeff)) then
         !$acc data copyin(x) copy(y_batch)
         call evaluate_batch_splines_1d_many_resident(spl, x, y_batch)
         !$acc end data
         return
      end if
#endif

      order = spl%order
      num_points = spl%num_points
      nq = spl%num_quantities
      periodic = spl%periodic
      x_min = spl%x_min
      h_step = spl%h_step
      period = h_step*real(num_points - 1, dp)

      do ipt = 1, size(x)
         include "spline1d_many_point_body.inc"
      end do
   end subroutine evaluate_batch_splines_1d_many

   subroutine evaluate_batch_splines_1d_many_resident(spl, x, y_batch)
      type(BatchSplineData1D), intent(in) :: spl
      real(dp), intent(in) :: x(:)
      real(dp), intent(inout) :: y_batch(:, :)

      integer :: ipt, iq, k_power, idx, k_wrap
      integer :: num_points, nq, order
      logical :: periodic
      real(dp) :: xj, x_norm, x_local, x_min, h_step, period, t, w

      if (size(y_batch, 1) < spl%num_quantities) then
         error stop "evaluate_batch_splines_1d_many_resident: y_batch dim1 too small"
      end if
      if (size(y_batch, 2) /= size(x)) then
         error stop "evaluate_batch_splines_1d_many_resident: y_batch dim2 mismatch"
      end if

      order = spl%order
      num_points = spl%num_points
      nq = spl%num_quantities
      periodic = spl%periodic
      x_min = spl%x_min
      h_step = spl%h_step
      period = h_step*real(num_points - 1, dp)

      !$acc parallel loop present(spl%coeff, x, y_batch) &
      !$acc& private(ipt, iq, k_power, idx, k_wrap, xj, x_norm, x_local, t, w)
      do ipt = 1, size(x)
         include "spline1d_many_point_body.inc"
      end do
      !$acc end parallel loop
   end subroutine evaluate_batch_splines_1d_many_resident

   subroutine evaluate_batch_splines_1d_many_der(spl, x, y_batch, dy_batch)
      type(BatchSplineData1D), intent(in) :: spl
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: y_batch(:, :)   ! (nq, npts)
      real(dp), intent(out) :: dy_batch(:, :)  ! (nq, npts)

      integer :: ipt, iq, k, idx, npts, nq, N
      real(dp) :: xj, x_norm, x_local, x_min, h_step, period

      npts = size(x)
      nq = spl%num_quantities
      N = spl%order
      x_min = spl%x_min
      h_step = spl%h_step
      period = h_step*real(spl%num_points - 1, dp)

      do ipt = 1, npts
         if (spl%periodic) then
            xj = modulo(x(ipt) - x_min, period) + x_min
         else
            xj = x(ipt)
         end if
         x_norm = (xj - x_min)/h_step
         idx = max(0, min(spl%num_points - 2, int(x_norm))) + 1
         x_local = (x_norm - real(idx - 1, dp))*h_step

         do iq = 1, nq
            y_batch(iq, ipt) = spl%coeff(iq, N, idx)
            dy_batch(iq, ipt) = N*spl%coeff(iq, N, idx)
         end do

         do k = N - 1, 1, -1
            do iq = 1, nq
               dy_batch(iq, ipt) = k*spl%coeff(iq, k, idx) + &
                                   x_local*dy_batch(iq, ipt)
            end do
         end do

         do k = N - 1, 0, -1
            do iq = 1, nq
               y_batch(iq, ipt) = spl%coeff(iq, k, idx) + &
                                  x_local*y_batch(iq, ipt)
            end do
         end do
      end do
   end subroutine evaluate_batch_splines_1d_many_der

   subroutine evaluate_batch_splines_1d_many_der2(spl, x, y_batch, dy_batch, d2y_batch)
      type(BatchSplineData1D), intent(in) :: spl
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: y_batch(:, :), dy_batch(:, :), d2y_batch(:, :)

      integer :: ipt, iq, k, idx, npts, nq, N
      real(dp) :: xj, x_norm, x_local, x_min, h_step, period

      npts = size(x)
      nq = spl%num_quantities
      N = spl%order
      x_min = spl%x_min
      h_step = spl%h_step
      period = h_step*real(spl%num_points - 1, dp)

      do ipt = 1, npts
         if (spl%periodic) then
            xj = modulo(x(ipt) - x_min, period) + x_min
         else
            xj = x(ipt)
         end if
         x_norm = (xj - x_min)/h_step
         idx = max(0, min(spl%num_points - 2, int(x_norm))) + 1
         x_local = (x_norm - real(idx - 1, dp))*h_step

         do iq = 1, nq
            y_batch(iq, ipt) = spl%coeff(iq, N, idx)
            dy_batch(iq, ipt) = N*spl%coeff(iq, N, idx)
            d2y_batch(iq, ipt) = N*(N - 1)*spl%coeff(iq, N, idx)
         end do

         do k = N - 1, 2, -1
            do iq = 1, nq
               d2y_batch(iq, ipt) = k*(k - 1)*spl%coeff(iq, k, idx) + &
                                    x_local*d2y_batch(iq, ipt)
            end do
         end do

         do k = N - 1, 1, -1
            do iq = 1, nq
               dy_batch(iq, ipt) = k*spl%coeff(iq, k, idx) + &
                                   x_local*dy_batch(iq, ipt)
            end do
         end do

         do k = N - 1, 0, -1
            do iq = 1, nq
               y_batch(iq, ipt) = spl%coeff(iq, k, idx) + &
                                  x_local*y_batch(iq, ipt)
            end do
         end do
      end do
   end subroutine evaluate_batch_splines_1d_many_der2

   subroutine evaluate_batch_splines_1d_many_der3(spl, x, y_batch, dy_batch, d2y_batch, &
                                                  d3y_batch)
      type(BatchSplineData1D), intent(in) :: spl
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: y_batch(:, :), dy_batch(:, :), d2y_batch(:, :), d3y_batch(:, :)

      integer :: ipt, iq, k, idx, npts, nq, N
      real(dp) :: xj, x_norm, x_local, x_min, h_step, period

      npts = size(x)
      nq = spl%num_quantities
      N = spl%order
      x_min = spl%x_min
      h_step = spl%h_step
      period = h_step*real(spl%num_points - 1, dp)

      do ipt = 1, npts
         if (spl%periodic) then
            xj = modulo(x(ipt) - x_min, period) + x_min
         else
            xj = x(ipt)
         end if
         x_norm = (xj - x_min)/h_step
         idx = max(0, min(spl%num_points - 2, int(x_norm))) + 1
         x_local = (x_norm - real(idx - 1, dp))*h_step

         do iq = 1, nq
            y_batch(iq, ipt) = spl%coeff(iq, N, idx)
            dy_batch(iq, ipt) = N*spl%coeff(iq, N, idx)
            d2y_batch(iq, ipt) = N*(N - 1)*spl%coeff(iq, N, idx)
            d3y_batch(iq, ipt) = N*(N - 1)*(N - 2)*spl%coeff(iq, N, idx)
         end do

         do k = N - 1, 3, -1
            do iq = 1, nq
               d3y_batch(iq, ipt) = k*(k - 1)*(k - 2)*spl%coeff(iq, k, idx) + &
                                    x_local*d3y_batch(iq, ipt)
            end do
         end do

         do k = N - 1, 2, -1
            do iq = 1, nq
               d2y_batch(iq, ipt) = k*(k - 1)*spl%coeff(iq, k, idx) + &
                                    x_local*d2y_batch(iq, ipt)
            end do
         end do

         do k = N - 1, 1, -1
            do iq = 1, nq
               dy_batch(iq, ipt) = k*spl%coeff(iq, k, idx) + &
                                   x_local*dy_batch(iq, ipt)
            end do
         end do

         do k = N - 1, 0, -1
            do iq = 1, nq
               y_batch(iq, ipt) = spl%coeff(iq, k, idx) + &
                                  x_local*y_batch(iq, ipt)
            end do
         end do
      end do
   end subroutine evaluate_batch_splines_1d_many_der3

   subroutine evaluate_batch_splines_1d_der(spl, x, y_batch, dy_batch)
      type(BatchSplineData1D), intent(in) :: spl
      real(dp), intent(in) :: x
      real(dp), intent(out) :: y_batch(:)   ! (n_quantities)
      real(dp), intent(out) :: dy_batch(:)  ! (n_quantities)

      real(dp) :: x_arr(1)
      real(dp) :: y_arr(spl%num_quantities, 1), dy_arr(spl%num_quantities, 1)

      x_arr(1) = x
      call evaluate_batch_splines_1d_many_der(spl, x_arr, y_arr, dy_arr)
      y_batch(1:spl%num_quantities) = y_arr(:, 1)
      dy_batch(1:spl%num_quantities) = dy_arr(:, 1)
   end subroutine evaluate_batch_splines_1d_der

   subroutine evaluate_batch_splines_1d_der2(spl, x, y_batch, dy_batch, d2y_batch)
      type(BatchSplineData1D), intent(in) :: spl
      real(dp), intent(in) :: x
      real(dp), intent(out) :: y_batch(:), dy_batch(:), d2y_batch(:)

      integer :: iq, k, idx, nq, N
      real(dp) :: xj, x_norm, x_local, x_min, h_step, period

      nq = spl%num_quantities
      N = spl%order
      x_min = spl%x_min
      h_step = spl%h_step
      period = h_step*real(spl%num_points - 1, dp)

      if (spl%periodic) then
         xj = x
         if (xj < x_min .or. xj >= x_min + period) then
            xj = modulo(xj - x_min, period) + x_min
         end if
      else
         xj = x
      end if

      x_norm = (xj - x_min)/h_step
      idx = max(0, min(spl%num_points - 2, int(x_norm))) + 1
      x_local = (x_norm - real(idx - 1, dp))*h_step

      do iq = 1, nq
         y_batch(iq) = spl%coeff(iq, N, idx)
         dy_batch(iq) = N*spl%coeff(iq, N, idx)
         d2y_batch(iq) = N*(N - 1)*spl%coeff(iq, N, idx)
      end do

      do k = N - 1, 2, -1
         do iq = 1, nq
            d2y_batch(iq) = k*(k - 1)*spl%coeff(iq, k, idx) + &
                            x_local*d2y_batch(iq)
         end do
      end do

      do k = N - 1, 1, -1
         do iq = 1, nq
            dy_batch(iq) = k*spl%coeff(iq, k, idx) + &
                           x_local*dy_batch(iq)
         end do
      end do

      do k = N - 1, 0, -1
         do iq = 1, nq
            y_batch(iq) = spl%coeff(iq, k, idx) + &
                          x_local*y_batch(iq)
         end do
      end do
   end subroutine evaluate_batch_splines_1d_der2

   subroutine evaluate_batch_splines_1d_der3(spl, x, y_batch, dy_batch, d2y_batch, &
                                             d3y_batch)
      type(BatchSplineData1D), intent(in) :: spl
      real(dp), intent(in) :: x
      real(dp), intent(out) :: y_batch(:), dy_batch(:), d2y_batch(:), d3y_batch(:)

      integer :: iq, k, idx, nq, N
      real(dp) :: xj, x_norm, x_local, x_min, h_step, period

      nq = spl%num_quantities
      N = spl%order
      x_min = spl%x_min
      h_step = spl%h_step
      period = h_step*real(spl%num_points - 1, dp)

      if (spl%periodic) then
         xj = x
         if (xj < x_min .or. xj >= x_min + period) then
            xj = modulo(xj - x_min, period) + x_min
         end if
      else
         xj = x
      end if

      x_norm = (xj - x_min)/h_step
      idx = max(0, min(spl%num_points - 2, int(x_norm))) + 1
      x_local = (x_norm - real(idx - 1, dp))*h_step

      do iq = 1, nq
         y_batch(iq) = spl%coeff(iq, N, idx)
         dy_batch(iq) = N*spl%coeff(iq, N, idx)
         d2y_batch(iq) = N*(N - 1)*spl%coeff(iq, N, idx)
         d3y_batch(iq) = N*(N - 1)*(N - 2)*spl%coeff(iq, N, idx)
      end do

      do k = N - 1, 3, -1
         do iq = 1, nq
            d3y_batch(iq) = k*(k - 1)*(k - 2)*spl%coeff(iq, k, idx) + &
                            x_local*d3y_batch(iq)
         end do
      end do

      do k = N - 1, 2, -1
         do iq = 1, nq
            d2y_batch(iq) = k*(k - 1)*spl%coeff(iq, k, idx) + &
                            x_local*d2y_batch(iq)
         end do
      end do

      do k = N - 1, 1, -1
         do iq = 1, nq
            dy_batch(iq) = k*spl%coeff(iq, k, idx) + &
                           x_local*dy_batch(iq)
         end do
      end do

      do k = N - 1, 0, -1
         do iq = 1, nq
            y_batch(iq) = spl%coeff(iq, k, idx) + &
                          x_local*y_batch(iq)
         end do
      end do
   end subroutine evaluate_batch_splines_1d_der3

end module batch_interpolate_1d
