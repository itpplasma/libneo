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


subroutine compute_field_from_gpec_file(filename, coil_currents, n_coils, &
                                         x_eval, y_eval, z_eval, n_points, Bx, By, Bz)
  implicit none

  integer, parameter :: dp = kind(1.0d0)

  character(*), intent(in) :: filename
  integer, intent(in) :: n_coils, n_points
  real(dp), intent(in) :: coil_currents(n_coils)
  real(dp), intent(in) :: x_eval(n_points), y_eval(n_points), z_eval(n_points)
  real(dp), intent(out) :: Bx(n_points), By(n_points), Bz(n_points)

  ! Call internal implementation
  call compute_field_from_gpec_file_impl(filename, coil_currents, n_coils, &
                                         x_eval, y_eval, z_eval, n_points, Bx, By, Bz)

end subroutine compute_field_from_gpec_file


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


subroutine compute_field_from_gpec_file_impl(filename, coil_currents, n_coils, &
                                              x_eval, y_eval, z_eval, n_points, Bx, By, Bz)
  use neo_biotsavart, only: coils_t, load_coils_from_gpec_file, coils_deinit, compute_magnetic_field
  implicit none

  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: clight = 2.99792458d10

  character(*), intent(in) :: filename
  integer, intent(in) :: n_coils, n_points
  real(dp), intent(in) :: coil_currents(n_coils)
  real(dp), intent(in) :: x_eval(n_points), y_eval(n_points), z_eval(n_points)
  real(dp), intent(out) :: Bx(n_points), By(n_points), Bz(n_points)

  type(coils_t) :: coils
  real(dp) :: pos(3), B_vec(3)
  real(dp) :: ampere_to_statampere
  integer :: i, ic, seg_start, seg_end, n_segments_per_coil

  ! Convert current from Amperes to statamperes (CGS-Gaussian units)
  ampere_to_statampere = clight / 10.0d0

  ! Load coils from GPEC file
  call load_coils_from_gpec_file(filename, coils)

  ! Figure out points per coil
  n_segments_per_coil = size(coils%current) / n_coils

  ! Multiply ALL currents by coil-specific value and convert units
  ! Points with current=0 (closing points) will remain 0
  do ic = 1, n_coils
    seg_start = (ic - 1) * n_segments_per_coil + 1
    seg_end = ic * n_segments_per_coil
    coils%current(seg_start:seg_end) = coils%current(seg_start:seg_end) * coil_currents(ic) * ampere_to_statampere
  end do

  ! Compute field at each evaluation point
  do i = 1, n_points
    pos = [x_eval(i), y_eval(i), z_eval(i)]
    B_vec = compute_magnetic_field(coils, pos)
    Bx(i) = B_vec(1)
    By(i) = B_vec(2)
    Bz(i) = B_vec(3)
  end do

  call coils_deinit(coils)

end subroutine compute_field_from_gpec_file_impl
