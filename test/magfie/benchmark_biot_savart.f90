program benchmark_biot_savart
  use iso_fortran_env, only: dp => real64
  use math_constants, only: current_si_to_cgs, C
  use coil_tools, only: coil_t, coils_read_GPEC, coils_append, coil_deinit, grid_from_bounding_box, &
       biot_savart_sum_coils
  use neo_biotsavart, only: coils_t, load_coils_from_gpec_file, coils_deinit, compute_magnetic_field, clight

  implicit none

  character(len=*), parameter :: coil_files(2) = [character(len=32) :: 'aug_bu.dat', 'aug_bl.dat']
  character(len=*), parameter :: currents_name = 'aug_currents.txt'
  character(len=*), parameter :: grid_name = 'vacfield_AUG_lowres.in'

  character(len=512) :: data_dir
  integer :: argc

  type(coil_t), allocatable :: coils_ct(:), tmp_coils(:)
  type(coils_t) :: coils_neo, tmp_neo

  real(dp), allocatable :: currents(:)
  real(dp), allocatable :: currents_ct(:)
  real(dp), allocatable :: B_neo(:, :, :, :), B_ct(:, :, :, :)
  real(dp), allocatable :: R(:), Z(:), phi(:)
  real(dp), allocatable :: cosphi(:), sinphi(:)
  real(dp), dimension(:, :, :, :), allocatable :: Bvac_ct

  integer, allocatable :: seg_per_coil(:)
  real(dp), allocatable :: windings(:)

  integer :: nc_per_file(2), nseg_per_file(2)
  real(dp) :: wind_per_file(2)

  integer :: total_coils, ic, idx, nphi, nR, nZ
  real(dp) :: Rmin, Rmax, Zmin, Zmax
  integer :: nphi_in, nR_in, nZ_in, nmax_dummy
  real(dp) :: prefactor, ampere_to_statampere
  integer :: seg_start, seg_end

  real(dp) :: t_start, t_end, time_neo, time_ct
  real(dp) :: max_abs_diff, max_rel_diff

  namelist /coil_field/ Rmin, Rmax, Zmin, Zmax, nR_in, nZ_in, nphi_in, nmax_dummy

  character(len=512) :: path_coil, path_currents, path_grid

  argc = command_argument_count()
  if (argc >= 1) then
    call get_command_argument(1, data_dir)
  else
    data_dir = 'test/magfie/test_data'
  end if

  call ensure_trailing_slash(data_dir)

  path_currents = trim(data_dir) // currents_name
  path_grid = trim(data_dir) // grid_name

  call read_gpec_headers(trim(data_dir), coil_files, nc_per_file, nseg_per_file, wind_per_file)

  total_coils = sum(nc_per_file)
  allocate(currents(total_coils))
  allocate(currents_ct(total_coils))
  allocate(seg_per_coil(total_coils))
  allocate(windings(total_coils))

  idx = 0
  do ic = 1, size(coil_files)
    seg_per_coil(idx+1:idx+nc_per_file(ic)) = nseg_per_file(ic) + 1
    windings(idx+1:idx+nc_per_file(ic)) = wind_per_file(ic)
    idx = idx + nc_per_file(ic)
  end do

  call read_currents_file(path_currents, total_coils, currents)

  call read_grid_file(path_grid, Rmin, Rmax, Zmin, Zmax, nR_in, nZ_in, nphi_in)

  nR = nR_in
  nZ = nZ_in
  nphi = nphi_in

  allocate(R(nR), Z(nZ), phi(nphi))
  call grid_from_bounding_box(Rmin, Rmax, nR, R, Zmin, Zmax, nZ, Z, nphi, phi)
  allocate(cosphi(nphi), sinphi(nphi))
  cosphi = cos(phi)
  sinphi = sin(phi)

  call load_coils_for_coiltools(trim(data_dir), coil_files, coils_ct)
  call load_coils_for_neo(trim(data_dir), coil_files, coils_neo)

  ampere_to_statampere = clight / 10.0_dp
  idx = 0
  do ic = 1, total_coils
    seg_start = idx + 1
    seg_end = idx + seg_per_coil(ic)
    coils_neo%current(seg_start:seg_end) = coils_neo%current(seg_start:seg_end) * &
      currents(ic) * windings(ic) * ampere_to_statampere
    idx = seg_end
  end do

  prefactor = current_si_to_cgs / C
  do ic = 1, total_coils
    currents_ct(ic) = prefactor * currents(ic) * windings(ic)
  end do

  allocate(B_neo(3, nZ, nphi, nR))
  allocate(B_ct(3, nZ, nphi, nR))

  call cpu_time(t_start)
  call evaluate_neo(coils_neo, R, Z, cosphi, sinphi, B_neo)
  call cpu_time(t_end)
  time_neo = t_end - t_start

  call cpu_time(t_start)
  call biot_savart_sum_coils(coils_ct, currents_ct, Rmin, Rmax, Zmin, Zmax, nR, nphi, nZ, Bvac_ct)
  call cpu_time(t_end)
  time_ct = t_end - t_start

  B_ct = Bvac_ct

  max_abs_diff = maxval(abs(B_neo - B_ct))
  max_rel_diff = maxval(abs(B_neo - B_ct) / max(abs(B_neo), 1.0e-20_dp))

  if (max_abs_diff > 1.0e-9_dp .or. max_rel_diff > 1.0e-9_dp) then
    write(*, '(A,1pe12.3)') 'Max abs diff: ', max_abs_diff
    write(*, '(A,1pe12.3)') 'Max rel diff: ', max_rel_diff
    error stop 'Biot-Savart implementations disagree beyond tolerance'
  end if

  write(*, '(A,F8.4)') 'neo_biotsavart runtime (s):   ', time_neo
  write(*, '(A,F8.4)') 'coil_tools runtime (s):        ', time_ct
  write(*, '(A,1pe12.3)') 'Max abs difference (Gauss):   ', max_abs_diff
  write(*, '(A,1pe12.3)') 'Max relative difference:      ', max_rel_diff

  call cleanup_coils(coils_ct)
  call coils_deinit(coils_neo)
  if (allocated(Bvac_ct)) deallocate(Bvac_ct)

contains

  subroutine ensure_trailing_slash(path)
    character(len=*), intent(inout) :: path
    integer :: l
    l = len_trim(path)
    if (path(l:l) /= '/' .and. path(l:l) /= '\\') then
      path = trim(path) // '/'
    end if
  end subroutine ensure_trailing_slash

  subroutine read_gpec_headers(base, files, nc_list, nseg_list, wind_list)
    character(len=*), intent(in) :: base
    character(len=*), intent(in) :: files(:)
    integer, intent(out) :: nc_list(:), nseg_list(:)
    real(dp), intent(out) :: wind_list(:)
    integer :: i, unit
    character(len=512) :: fullpath
    integer :: idum
    real(dp) :: ddum

    do i = 1, size(files)
      fullpath = trim(base) // trim(files(i))
      open(newunit=unit, file=fullpath, status='old', action='read')
      read(unit, *) nc_list(i), idum, nseg_list(i), ddum
      wind_list(i) = ddum
      close(unit)
    end do
  end subroutine read_gpec_headers

  subroutine read_currents_file(filename, n, values)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: n
    real(dp), intent(out) :: values(n)
    integer :: unit

    open(newunit=unit, file=filename, status='old', action='read')
    read(unit, *) values
    close(unit)
  end subroutine read_currents_file

  subroutine read_grid_file(filename, Rmin, Rmax, Zmin, Zmax, nR, nZ, nphi)
    character(len=*), intent(in) :: filename
    real(dp), intent(out) :: Rmin, Rmax, Zmin, Zmax
    integer, intent(out) :: nR, nZ, nphi
    integer :: nmax_dummy, unit
    namelist /coil_field/ Rmin, Rmax, Zmin, Zmax, nR, nZ, nphi, nmax_dummy

    open(newunit=unit, file=filename, status='old', action='read')
    read(unit, nml=coil_field)
    close(unit)
  end subroutine read_grid_file

  subroutine load_coils_for_coiltools(base, files, coils)
    character(len=*), intent(in) :: base
    character(len=*), intent(in) :: files(:)
    type(coil_t), allocatable, intent(out) :: coils(:)
    type(coil_t), allocatable :: temp(:)
    integer :: i

    do i = 1, size(files)
      if (i == 1) then
        call coils_read_GPEC(trim(base) // trim(files(i)), coils)
      else
        call coils_read_GPEC(trim(base) // trim(files(i)), temp)
        call coils_append(coils, temp)
      end if
    end do
  end subroutine load_coils_for_coiltools

  subroutine load_coils_for_neo(base, files, coils)
    character(len=*), intent(in) :: base
    character(len=*), intent(in) :: files(:)
    type(coils_t), intent(out) :: coils
    type(coils_t) :: temp
    integer :: i

    if (allocated(coils%x)) call coils_deinit(coils)
    do i = 1, size(files)
      call load_coils_from_gpec_file(trim(base) // trim(files(i)), temp)
      call append_neo_coils(coils, temp)
      call coils_deinit(temp)
    end do
  end subroutine load_coils_for_neo

  subroutine append_neo_coils(dst, src)
    type(coils_t), intent(inout) :: dst
    type(coils_t), intent(inout) :: src
    integer :: n_old, n_new
    real(dp), allocatable :: tmp(:)

    if (.not. allocated(dst%x)) then
      call move_alloc(src%x, dst%x)
      call move_alloc(src%y, dst%y)
      call move_alloc(src%z, dst%z)
      call move_alloc(src%current, dst%current)
      return
    end if

    n_old = size(dst%x)
    n_new = n_old + size(src%x)

    allocate(tmp(n_new))
    tmp(1:n_old) = dst%x
    tmp(n_old+1:n_new) = src%x
    call move_alloc(tmp, dst%x)

    allocate(tmp(n_new))
    tmp(1:n_old) = dst%y
    tmp(n_old+1:n_new) = src%y
    call move_alloc(tmp, dst%y)

    allocate(tmp(n_new))
    tmp(1:n_old) = dst%z
    tmp(n_old+1:n_new) = src%z
    call move_alloc(tmp, dst%z)

    allocate(tmp(n_new))
    tmp(1:n_old) = dst%current
    tmp(n_old+1:n_new) = src%current
    call move_alloc(tmp, dst%current)
  end subroutine append_neo_coils

  subroutine evaluate_neo(coils, R, Z, cosphi, sinphi, Bcyl)
    type(coils_t), intent(in) :: coils
    real(dp), intent(in), dimension(:) :: R, Z, cosphi, sinphi
    real(dp), intent(out), dimension(3, size(Z), size(cosphi), size(R)) :: Bcyl
    integer :: kR, kZ, kphi
    real(dp) :: eval_point(3), Bxyz(3)

    do kZ = 1, size(Z)
      eval_point(3) = Z(kZ)
      do kphi = 1, size(cosphi)
        do kR = 1, size(R)
          eval_point(1) = R(kR) * cosphi(kphi)
          eval_point(2) = R(kR) * sinphi(kphi)
          Bxyz = compute_magnetic_field(coils, eval_point)
          Bcyl(1, kZ, kphi, kR) = Bxyz(1) * cosphi(kphi) + Bxyz(2) * sinphi(kphi)
          Bcyl(2, kZ, kphi, kR) = -Bxyz(1) * sinphi(kphi) + Bxyz(2) * cosphi(kphi)
          Bcyl(3, kZ, kphi, kR) = Bxyz(3)
        end do
      end do
    end do
  end subroutine evaluate_neo

  subroutine cleanup_coils(coils)
    type(coil_t), allocatable, intent(inout) :: coils(:)
    integer :: i
    if (.not. allocated(coils)) return
    do i = 1, size(coils)
      call coil_deinit(coils(i))
    end do
    deallocate(coils)
  end subroutine cleanup_coils

end program benchmark_biot_savart
