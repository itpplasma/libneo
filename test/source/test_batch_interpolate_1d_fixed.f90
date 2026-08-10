program test_batch_interpolate_1d_fixed
   use, intrinsic :: iso_fortran_env, only: dp => real64
   use batch_interpolate_types, only: BatchSplineData1D
   use batch_interpolate_1d, only: evaluate_batch_spline_1d_scalar_cubic_der, &
                                   evaluate_batch_spline_1d_pair_cubic_der, &
                                   evaluate_batch_spline_1d_scalar_quintic_der, &
                                   evaluate_batch_spline_1d_pair_quintic_der
   implicit none

   real(dp), parameter :: TOL = 2.0e-14_dp
   type(BatchSplineData1D) :: scalar_spl, pair_spl
   real(dp) :: y, dy, y1, dy1, y2, dy2
   real(dp) :: expected_y, expected_dy
   integer :: iq, k

   call test_cubic

   call make_polynomial_spline(scalar_spl, 1, 5, .false.)
   call evaluate_batch_spline_1d_scalar_quintic_der( &
      scalar_spl, scalar_spl%x_min + 0.37_dp*scalar_spl%h_step, y, dy)
   call polynomial_oracle(scalar_spl, 1, 0.37_dp*scalar_spl%h_step, &
      expected_y, expected_dy)
   call assert_close(y, expected_y, "scalar value")
   call assert_close(dy, expected_dy, "scalar derivative")

   call make_polynomial_spline(pair_spl, 2, 5, .true.)
   call evaluate_batch_spline_1d_pair_quintic_der( &
      pair_spl, pair_spl%x_min + pair_spl%period + &
      0.61_dp*pair_spl%h_step, y1, dy1, y2, dy2)
   call polynomial_oracle(pair_spl, 1, 0.61_dp*pair_spl%h_step, &
      expected_y, expected_dy)
   call assert_close(y1, expected_y, "pair first value")
   call assert_close(dy1, expected_dy, "pair first derivative")
   call polynomial_oracle(pair_spl, 2, 0.61_dp*pair_spl%h_step, &
      expected_y, expected_dy)
   call assert_close(y2, expected_y, "pair second value")
   call assert_close(dy2, expected_dy, "pair second derivative")

   deallocate (scalar_spl%coeff, pair_spl%coeff)
   print *, "PASSED: fixed cubic/quintic 1D evaluators satisfy polynomial oracle"

contains

   subroutine test_cubic
      type(BatchSplineData1D) :: scalar, pair
      real(dp) :: got_y, got_dy, got_y1, got_dy1, got_y2, got_dy2
      real(dp) :: ref_y, ref_dy

      call make_polynomial_spline(scalar, 1, 3, .false.)
      call evaluate_batch_spline_1d_scalar_cubic_der( &
         scalar, scalar%x_min + 0.29_dp*scalar%h_step, got_y, got_dy)
      call polynomial_oracle(scalar, 1, 0.29_dp*scalar%h_step, ref_y, ref_dy)
      call assert_close(got_y, ref_y, "cubic scalar value")
      call assert_close(got_dy, ref_dy, "cubic scalar derivative")

      call make_polynomial_spline(pair, 2, 3, .true.)
      call evaluate_batch_spline_1d_pair_cubic_der( &
         pair, pair%x_min - pair%period + 0.73_dp*pair%h_step, &
         got_y1, got_dy1, got_y2, got_dy2)
      call polynomial_oracle(pair, 1, 0.73_dp*pair%h_step, ref_y, ref_dy)
      call assert_close(got_y1, ref_y, "cubic pair first value")
      call assert_close(got_dy1, ref_dy, "cubic pair first derivative")
      call polynomial_oracle(pair, 2, 0.73_dp*pair%h_step, ref_y, ref_dy)
      call assert_close(got_y2, ref_y, "cubic pair second value")
      call assert_close(got_dy2, ref_dy, "cubic pair second derivative")
      deallocate (scalar%coeff, pair%coeff)
   end subroutine test_cubic

   subroutine make_polynomial_spline(spl, nq, order, periodic)
      type(BatchSplineData1D), intent(out) :: spl
      integer, intent(in) :: nq, order
      logical, intent(in) :: periodic

      spl%order = order
      spl%num_points = 2
      spl%periodic = periodic
      spl%x_min = -0.4_dp
      spl%h_step = 1.7_dp
      spl%inv_h_step = 1.0_dp/spl%h_step
      spl%period = spl%h_step
      spl%num_quantities = nq
      allocate (spl%coeff(nq, 0:order, 2))
      spl%coeff = 0.0_dp
      do iq = 1, nq
         do k = 0, order
            spl%coeff(iq, k, 1) = 0.0625_dp*real(3*iq + 2*k - 5, dp)
         end do
      end do
   end subroutine make_polynomial_spline

   subroutine polynomial_oracle(spl, quantity, x_local, value, derivative)
      type(BatchSplineData1D), intent(in) :: spl
      integer, intent(in) :: quantity
      real(dp), intent(in) :: x_local
      real(dp), intent(out) :: value, derivative
      integer :: power

      value = 0.0_dp
      derivative = 0.0_dp
      do power = 0, spl%order
         value = value + spl%coeff(quantity, power, 1)*x_local**power
      end do
      do power = 1, spl%order
         derivative = derivative + real(power, dp)* &
            spl%coeff(quantity, power, 1)*x_local**(power - 1)
      end do
   end subroutine polynomial_oracle

   subroutine assert_close(got, expected, label)
      real(dp), intent(in) :: got, expected
      character(*), intent(in) :: label
      if (abs(got - expected) > TOL*max(1.0_dp, abs(expected))) then
         print *, trim(label), got, expected
         error stop "fixed quintic 1D evaluator failed polynomial oracle"
      end if
   end subroutine assert_close

end program test_batch_interpolate_1d_fixed
