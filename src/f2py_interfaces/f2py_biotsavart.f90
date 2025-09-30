! f2py interface for biotsavart module
! Provides Python-friendly wrappers for direct Biot-Savart field computation
! This is a standalone wrapper that does not expose the module internals to f2py

subroutine compute_field_direct_biotsavart(x_coil, y_coil, z_coil, current_coil, n_segments, &
                                           x_eval, y_eval, z_eval, n_points, Bx, By, Bz)
  implicit none

  integer, parameter :: dp = kind(1.0d0)

  ! Input: coil geometry and currents
  integer, intent(in) :: n_segments
  real(dp), intent(in) :: x_coil(n_segments), y_coil(n_segments), z_coil(n_segments)
  real(dp), intent(in) :: current_coil(n_segments)

  ! Input: evaluation points
  integer, intent(in) :: n_points
  real(dp), intent(in) :: x_eval(n_points), y_eval(n_points), z_eval(n_points)

  ! Output: magnetic field at evaluation points
  real(dp), intent(out) :: Bx(n_points), By(n_points), Bz(n_points)

  ! Local variables
  real(dp) :: pos(3), B_vec(3)
  integer :: i

  ! Call internal implementation
  call compute_field_direct_biotsavart_impl(x_coil, y_coil, z_coil, current_coil, n_segments, &
                                             x_eval, y_eval, z_eval, n_points, Bx, By, Bz)

end subroutine compute_field_direct_biotsavart


! Internal implementation that uses the bio module
subroutine compute_field_direct_biotsavart_impl(x_coil, y_coil, z_coil, current_coil, n_segments, &
                                                 x_eval, y_eval, z_eval, n_points, Bx, By, Bz)
  use neo_biotsavart, only: coils_t, coils_init, coils_deinit, compute_magnetic_field
  implicit none

  integer, parameter :: dp = kind(1.0d0)

  integer, intent(in) :: n_segments, n_points
  real(dp), intent(in) :: x_coil(n_segments), y_coil(n_segments), z_coil(n_segments)
  real(dp), intent(in) :: current_coil(n_segments)
  real(dp), intent(in) :: x_eval(n_points), y_eval(n_points), z_eval(n_points)
  real(dp), intent(out) :: Bx(n_points), By(n_points), Bz(n_points)

  type(coils_t) :: coils
  real(dp) :: pos(3), B_vec(3)
  integer :: i

  call coils_init(x_coil, y_coil, z_coil, current_coil, coils)

  do i = 1, n_points
    pos = [x_eval(i), y_eval(i), z_eval(i)]
    B_vec = compute_magnetic_field(coils, pos)
    Bx(i) = B_vec(1)
    By(i) = B_vec(2)
    Bz(i) = B_vec(3)
  end do

  call coils_deinit(coils)

end subroutine compute_field_direct_biotsavart_impl
