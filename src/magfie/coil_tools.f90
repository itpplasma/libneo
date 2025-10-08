module coil_tools

  use iso_fortran_env, only : error_unit
  use libneo_kinds, only : dp

  implicit none

  private

  public :: coil_t, coil_init, coil_deinit, coils_append, &
    process_fixed_number_of_args, check_number_of_args, &
    coils_write_AUG, coils_read_AUG, &
    coils_write_nemov, coils_read_nemov, &
    coils_write_GPEC, coils_read_GPEC, &
    read_currents, grid_from_bounding_box, &
    biot_savart_sum_coils, biot_savart_fourier, &
    write_Bvac_nemov, write_Bnvac_fourier, read_Bnvac_fourier, &
    vector_potential_biot_savart_fourier, write_Anvac_fourier, read_Anvac_fourier, &
    sum_coils_gauge_single_mode_Anvac, gauged_Anvac_from_Bnvac, &
    segment_vector_potential_contribution, segment_eccentricity, &
    segment_gradient_kernel, segment_gradient_contribution

  type :: coil_t
    integer :: nseg = 0
    integer :: nwind = 0
    real(dp), allocatable :: XYZ(:, :)
  end type coil_t

  character(len = *), parameter :: arg_size_fmt = &
    '("Argument size mismatch in ", a, ": ", a, " = ", i0, ", ", a, " = ", i0)'

contains

  function linspace(lo, hi, cnt, excl_lo, excl_hi)
    real(dp), intent(in) :: lo, hi
    integer, intent(in) :: cnt
    integer, intent(in), optional :: excl_lo, excl_hi
    real(dp) :: linspace(cnt)
    real(dp) :: step
    integer :: k, omit_lo, omit_hi

    omit_lo = 0
    if (present(excl_lo)) omit_lo = excl_lo
    omit_hi = 0
    if (present(excl_hi)) omit_hi = excl_hi
    step = (hi - lo) / dble(cnt - 1 + omit_lo + omit_hi)
    ! gfortran 14.2.1 thinks that the equivalent
    ! implied do loop may be uninitialized here
    do k = omit_lo, cnt - 1 + omit_lo
      linspace(k + 1) = lo + k * step
    end do
    if (omit_hi == 0) linspace(cnt) = hi
  end function linspace

  subroutine check_number_of_args(expected)
    integer, intent(in) :: expected
    integer :: argc
    character(len = *), parameter :: argc_err_fmt = &
      '("Expected at least ", i0, " command line arguments, but received only ", i0, ".")'

    argc = command_argument_count()
    if (argc < expected) then
      write (error_unit, argc_err_fmt) expected, argc
      error stop
    end if
  end subroutine check_number_of_args

  subroutine process_fixed_number_of_args(offset, number, args)
    integer, intent(in) :: offset
    integer, intent(out) :: number
    character(len = 1024), dimension(:), allocatable, intent(out) :: args
    character(len = 1024) :: decimal_number, err_msg
    integer :: k, status

    call check_number_of_args(offset)
    call get_command_argument(offset, decimal_number)
    read (decimal_number, *, iostat = status, iomsg = err_msg) number
    if (status /= 0) then
      write (error_unit, '(a)') trim(err_msg)
      error stop
    end if
    if (number < 1) then
      write (error_unit, '("Number of further command line arguments must be positive.")')
      error stop
    end if
    call check_number_of_args(offset + number)
    allocate(character(len = 1024) :: args(number))
    do k = 1, number
      call get_command_argument(offset + k, args(k))
    end do
  end subroutine process_fixed_number_of_args

  subroutine coil_init(coil, nseg, nwind)
    type(coil_t), intent(inout) :: coil
    integer, intent(in) :: nseg
    integer, intent(in) :: nwind

    call coil_deinit(coil)
    allocate(coil%XYZ(3, nseg))
    coil%XYZ(:, :) = 0d0
    coil%nseg = nseg
    coil%nwind = nwind
  end subroutine coil_init

  subroutine coil_deinit(coil)
    type(coil_t), intent(inout) :: coil

    if (allocated(coil%XYZ)) deallocate(coil%XYZ)
    coil%nseg = 0
    coil%nwind = 0
  end subroutine coil_deinit

  subroutine coils_append(coils_head, coils_tail)
    type(coil_t), dimension(:), allocatable, intent(inout) :: coils_head
    type(coil_t), dimension(:), allocatable, intent(inout) :: coils_tail
    type(coil_t), dimension(:), allocatable :: coils_temp

    if (.not. allocated(coils_tail)) then
      return
    end if
    if (.not. allocated(coils_head)) then
      call move_alloc(coils_tail, coils_head)
      return
    end if
    if (0 == size(coils_tail)) then
      deallocate(coils_tail)
      return
    end if
    if (0 == size(coils_head)) then
      deallocate(coils_head)
      call move_alloc(coils_tail, coils_head)
      return
    end if
    allocate(coils_temp(size(coils_head) + size(coils_tail)))
    coils_temp(:size(coils_head)) = coils_head
    coils_temp(size(coils_head) + 1:) = coils_tail
    deallocate(coils_head, coils_tail)
    call move_alloc(coils_temp, coils_head)
  end subroutine coils_append

  subroutine coils_write_AUG(filename, coil)
    use math_constants, only: length_si_to_cgs
    character(len = *), intent(in) :: filename
    type(coil_t), intent(in) :: coil
    integer :: fid, ks

    open(newunit = fid, file = filename, status = 'replace', &
      action = 'write', form = 'formatted')
    do ks = 1, coil%nseg
      write (fid, '(3(es24.16e3, :, 1x))') &
        hypot(coil%XYZ(2, ks), coil%XYZ(1, ks)) / length_si_to_cgs, &
        coil%XYZ(3, ks) / length_si_to_cgs, &
        atan2(coil%XYZ(2, ks), coil%XYZ(1, ks))
    end do
    close(fid)
  end subroutine coils_write_AUG

  subroutine coils_read_AUG(filename, coil)
    use math_constants, only: length_si_to_cgs
    character(len = *), intent(in) :: filename
    type(coil_t), intent(inout) :: coil
    integer :: fid, status, nseg, ks
    real(dp), dimension(:, :), allocatable :: R_Z_phi

    ! count number of lines = number of coil segments
    open(newunit = fid, file = filename, status = 'old', &
      action = 'read', form = 'formatted')
    nseg = 0
    do
      read(fid, *, iostat = status)
      if (status /= 0) exit
      nseg = nseg + 1
    end do
    call coil_init(coil, nseg, 1)
    allocate(R_Z_phi(nseg, 3))
    rewind fid
    do ks = 1, nseg
      read (fid, *) R_Z_phi(ks, :)
    end do
    close(fid)
    coil%XYZ(1, :) = length_si_to_cgs * R_Z_phi(:, 1) * cos(R_Z_phi(:, 3))
    coil%XYZ(2, :) = length_si_to_cgs * R_Z_phi(:, 1) * sin(R_Z_phi(:, 3))
    coil%XYZ(3, :) = length_si_to_cgs * R_Z_phi(:, 2)
    deallocate(R_Z_phi)
  end subroutine coils_read_AUG

  subroutine coils_write_nemov(filename, coils)
    character(len = *), intent(in) :: filename
    type(coil_t), intent(in), dimension(:) :: coils
    integer :: fid, kc, ks

    open(newunit = fid, file = filename, status = 'replace', &
      action = 'write', form = 'formatted')
    ! first point is repeated for each coil
    write (fid, '(i0)') sum(coils(:)%nseg + 1)
    do kc = 1, size(coils)
      do ks = 1, coils(kc)%nseg
        write (fid, '(3(es24.16e3, 1x), f3.1, 1x, i0)') &
          coils(kc)%XYZ(:, ks), 1.0, kc
      end do
      write (fid, '(3(es24.16e3, 1x), f3.1, 1x, i0)') &
        coils(kc)%XYZ(:, 1), 0.0, kc
    end do
    close(fid)
  end subroutine coils_write_nemov

  subroutine coils_read_nemov(filename, coils)
    character(len = *), intent(in) :: filename
    type(coil_t), intent(out), allocatable, dimension(:) :: coils
    integer :: fid, total, k, ncoil, kc, nseg, ks, idum
    real(dp) :: X, Y, Z, cur

    open(newunit = fid, file = filename, status = 'old', &
      action = 'read', form = 'formatted')
    read (fid, *) total
    ncoil = 0
    do k = 1, total !figure out the number of coils in this loop
      read (fid, *) X, Y, Z, cur, kc
      ncoil = max(ncoil, kc)
    end do
    allocate(coils(ncoil))
    rewind fid
    read (fid, *) total
    nseg = 0
    do k = 1, total !figure out number of segments per coil and initialise the vertex coordinates in this loop
      read (fid, *) X, Y, Z, cur, kc
      if (abs(cur) > 0d0) then
        nseg = nseg + 1
      else
        call coil_init(coils(kc), nseg, 1)
        nseg = 0
      end if
    end do
    rewind fid
    read (fid, *) total
    k = 1
    do kc = 1, ncoil !save the vertex coordinates in coils in this loop
      do ks = 1, coils(kc)%nseg
        read (fid, *) coils(kc)%XYZ(:, ks), cur, idum
        k = k + 1
        if (idum /= kc) then
          write (error_unit, '("coils_read_nemov: expected coil index ", ' // &
            'i0, " in ", a, " at line ", i0, ", but got ", i0)') kc, filename, k, idum
          error stop
        end if
      end do
      read (fid, *) X, Y, Z, cur, idum
    end do
    close(fid)
  end subroutine coils_read_nemov

  subroutine coils_write_GPEC(filename, coils)
    use math_constants, only: length_si_to_cgs
    character(len = *), intent(in) :: filename
    type(coil_t), intent(in), dimension(:) :: coils
    integer :: fid, ncoil, kc, ks

    if (any(coils(:)%nseg /= coils(1)%nseg)) then
      write (error_unit, '("coils_write_GPEC: not all coils have the same number of segments.")')
      error stop
    end if
    if (any(coils(:)%nwind /= coils(1)%nwind)) then
      write (error_unit, '("coils_write_GPEC: not all coils have the same winding number.")')
      error stop
    end if
    ncoil = size(coils)
    open(newunit = fid, file = filename, status = 'replace', &
      action = 'write', form = 'formatted')
    write (fid, '(3(i0, 1x), es24.16e3)') ncoil, 1, coils(1)%nseg + 1, dble(coils(1)%nwind)
    do kc = 1, ncoil
      do ks = 1, coils(kc)%nseg
        write (fid, '(3(es24.16e3, :, 1x))') coils(kc)%XYZ(:, ks) / length_si_to_cgs
      end do
      write (fid, '(3(es24.16e3, :, 1x))') coils(kc)%XYZ(:, 1) / length_si_to_cgs
    end do
    close(fid)
  end subroutine coils_write_GPEC

  subroutine coils_read_GPEC(filename, coils)
    use math_constants, only: length_si_to_cgs
    character(len = *), intent(in) :: filename
    type(coil_t), intent(out), allocatable, dimension(:) :: coils
    integer :: fid, ncoil, nseg, nwind, kc, ks, idum
    real(dp) :: ddum

    open(newunit = fid, file = filename, status = 'old', &
      action = 'read', form = 'formatted')
    read (fid, *) ncoil, idum, nseg, ddum
    nseg = nseg - 1
    nwind = int(ddum)
    allocate(coils(ncoil))
    do kc = 1, ncoil
      call coil_init(coils(kc), nseg, nwind)
      do ks = 1, nseg
        read (fid, *) coils(kc)%XYZ(:, ks)
      end do
      read (fid, *)
      coils(kc)%XYZ(:, :) = length_si_to_cgs * coils(kc)%XYZ
    end do
    close(fid)
  end subroutine coils_read_GPEC

  !> Reads all floating-point numbers from one line.
  subroutine read_currents(filename, Ic)
    character(len = *), intent(in) :: filename
    real(dp), intent(out), allocatable :: Ic(:)
    integer :: fid, status, ncoil

    open(newunit = fid, file = filename, status = 'old', action = 'read', form = 'formatted')
    ncoil = 1
    status = 0
    do while (status == 0)
      allocate(Ic(ncoil))
      read (fid, *, iostat = status) Ic(:)
      rewind fid
      deallocate(Ic)
      if (status == 0) then
        ncoil = ncoil + 1
      else
        ncoil = ncoil - 1
      end if
    end do
    allocate(Ic(ncoil))
    read (fid, *) Ic
    close(fid)
  end subroutine read_currents

  subroutine grid_from_bounding_box(Rmin, Rmax, nR, R, Zmin, Zmax, nZ, Z, nphi, phi)
    use math_constants, only: pi
    real(dp), intent(in) :: Rmin, Rmax
    integer, intent(in) :: nR
    real(dp), dimension(:), intent(inout) :: R
    real(dp), intent(in) :: Zmin, Zmax
    integer, intent(in) :: nZ
    real(dp), dimension(:), intent(inout) :: Z
    integer, intent(in), optional :: nphi
    real(dp), dimension(:), intent(inout), optional :: phi

    R(:) = linspace(Rmin, Rmax, nR, 0, 0)
    Z(:) = linspace(Zmin, Zmax, nZ, 0, 0)
    if (present(nphi) .and. present(phi)) then
      ! half-open interval: do not repeat phi = 0 at phi = 2 pi
      phi(:) = linspace(0d0, 2d0 * pi, nphi, 0, 1)
    end if
  end subroutine grid_from_bounding_box

  subroutine biot_savart_sum_coils(coils, Ic, &
    Rmin, Rmax, Zmin, Zmax, nR, nphi, nZ, Bvac)
    type(coil_t), intent(in), dimension(:) :: coils
    real(dp), intent(in), dimension(:) :: Ic
    real(dp), intent(in) :: Rmin, Rmax, Zmin, Zmax
    integer, intent(in) :: nR, nphi, nZ
    real(dp), intent(out), dimension(:, :, :, :), allocatable :: Bvac
    integer :: ncoil, kc, ks, kR, kphi, kZ
    real(dp), dimension(nphi) :: phi, cosphi, sinphi
    real(dp) :: R(nR), Z(nZ), XYZ_r(3), XYZ_i(3), XYZ_f(3), dist_i, dist_f, BXYZ_c(3), BXYZ(3), denom

    if (size(coils) /= size(Ic)) then
      write (error_unit, arg_size_fmt) 'biot_savart_sum_coils', &
        'size(coils)', size(coils), 'size(Ic)', size(Ic)
      error stop
    end if
    ncoil = size(coils)
    call grid_from_bounding_box(Rmin, Rmax, nR, R, Zmin, Zmax, nZ, Z, nphi, phi)
    cosphi(:) = cos(phi)
    sinphi(:) = sin(phi)
    allocate(Bvac(3, nZ, nphi, nR))
    Bvac(:, :, :, :) = 0d0
    !$omp parallel do schedule(static) collapse(3) default(none) &
    !$omp private(kr, kphi, kZ, kc, ks, XYZ_r, XYZ_i, XYZ_f, dist_i, dist_f, BXYZ, BXYZ_c, denom) &
    !$omp shared(nR, nphi, nZ, ncoil, R, Z, cosphi, sinphi, coils, Ic, Bvac)
    do kZ = 1, nZ
      do kphi = 1, nphi
        do kR = 1, nR
          XYZ_r(:) = [R(kR) * cosphi(kphi), R(kR) * sinphi(kphi), Z(kZ)]
          ! Biot-Savart integral over coil segments
          BXYZ(:) = 0d0
          do kc = 1, ncoil
            BXYZ_c(:) = 0d0
            XYZ_f(:) = coils(kc)%XYZ(:, coils(kc)%nseg) - XYZ_r
            dist_f = sqrt(sum(XYZ_f * XYZ_f))
            do ks = 1, coils(kc)%nseg
              XYZ_i(:) = XYZ_f
              dist_i = max(1.0d-12, sqrt(sum(XYZ_i * XYZ_i)))
              XYZ_f(:) = coils(kc)%XYZ(:, ks) - XYZ_r
              dist_f = max(1.0d-12, sqrt(sum(XYZ_f * XYZ_f)))
              denom = dist_i * dist_f * (dist_i * dist_f + sum(XYZ_i * XYZ_f))
              if (denom == 0.0d0) denom = 1.0d-12
              BXYZ_c(:) = BXYZ_c + &
                (XYZ_i([2, 3, 1]) * XYZ_f([3, 1, 2]) - XYZ_i([3, 1, 2]) * XYZ_f([2, 3, 1])) * &
                (dist_i + dist_f) / denom
            end do
            BXYZ(:) = BXYZ + Ic(kc) * BXYZ_c
          end do
          Bvac(1, kZ, kphi, kR) = BXYZ(1) * cosphi(kphi) + BXYZ(2) * sinphi(kphi)
          Bvac(2, kZ, kphi, kR) = BXYZ(2) * cosphi(kphi) - BXYZ(1) * sinphi(kphi)
          Bvac(3, kZ, kphi, kR) = BXYZ(3)
        end do
      end do
    end do
    !$omp end parallel do
  end subroutine biot_savart_sum_coils

  subroutine biot_savart_fourier(coils, nmax, &
    Rmin, Rmax, Zmin, Zmax, nR, nphi, nZ, Bn)
    use iso_c_binding, only: c_ptr, c_double, c_double_complex, c_size_t, c_f_pointer, c_null_ptr
    !$ use omp_lib, only: omp_get_max_threads
    use FFTW3, only: fftw_init_threads, fftw_plan_with_nthreads, fftw_cleanup_threads, &
      fftw_alloc_real, fftw_alloc_complex, fftw_plan_dft_r2c_1d, FFTW_PATIENT, &
      FFTW_DESTROY_INPUT, fftw_execute_dft_r2c, fftw_destroy_plan, fftw_free
    type(coil_t), intent(in), dimension(:) :: coils
    integer, intent(in) :: nmax
    real(dp), intent(in) :: Rmin, Rmax, Zmin, Zmax
    integer, intent(in) :: nR, nphi, nZ
    complex(dp), intent(out), dimension(:, :, :, :, :), allocatable :: Bn
    integer :: nfft, ncoil, kc, ks, kR, kphi, kZ
    real(dp), dimension(nphi) :: phi, cosphi, sinphi
    real(dp) :: R(nR), Z(nZ), XYZ_r(3), XYZ_i(3), XYZ_f(3), dist_i, dist_f, BXYZ(3), denom
    type(c_ptr) :: plan_nphi, p_BR, p_Bphi, p_BZ, p_BnR, p_Bnphi, p_BnZ
    real(c_double), dimension(:), pointer :: BR, Bphi, BZ
    complex(c_double_complex), dimension(:), pointer :: BnR, Bnphi, BnZ

    if (nmax > nphi / 4) then
      write (error_unit, '("biot_savart_fourier: requested nmax = ", ' // &
        'i0, ", but only ", i0, " modes available.")') nmax, nphi / 4
      error stop
    end if
    nfft = nphi / 2 + 1
    call grid_from_bounding_box(Rmin, Rmax, nR, R, Zmin, Zmax, nZ, Z, nphi, phi)
    cosphi(:) = cos(phi)
    sinphi(:) = sin(phi)
    ncoil = size(coils)
    allocate(Bn(0:nmax, 3, nR, nZ, ncoil))
    ! prepare FFTW
    !$ if (fftw_init_threads() == 0) error stop 'OpenMP support in FFTW could not be initialized'
    !$ call fftw_plan_with_nthreads(omp_get_max_threads())
    p_BR = fftw_alloc_real(int(nphi, c_size_t))
    call c_f_pointer(p_BR, BR, [nphi])
    p_Bphi = fftw_alloc_real(int(nphi, c_size_t))
    call c_f_pointer(p_Bphi, Bphi, [nphi])
    p_BZ = fftw_alloc_real(int(nphi, c_size_t))
    call c_f_pointer(p_BZ, BZ, [nphi])
    p_BnR = fftw_alloc_complex(int(nfft, c_size_t))
    call c_f_pointer(p_BnR, BnR, [nfft])
    p_Bnphi = fftw_alloc_complex(int(nfft, c_size_t))
    call c_f_pointer(p_Bnphi, Bnphi, [nfft])
    p_BnZ = fftw_alloc_complex(int(nfft, c_size_t))
    call c_f_pointer(p_BnZ, BnZ, [nfft])
    plan_nphi = fftw_plan_dft_r2c_1d(nphi, BR, BnR, ior(FFTW_PATIENT, FFTW_DESTROY_INPUT))
    do kc = 1, ncoil
      do kZ = 1, nZ
        do kR = 1, nR
          !$omp parallel do schedule(static) default(none) &
          !$omp private(kphi, ks, XYZ_r, XYZ_i, XYZ_f, dist_i, dist_f, BXYZ, denom) &
          !$omp shared(nphi, kc, coils, R, kR, Z, kZ, cosphi, sinphi, BR, Bphi, BZ)
          do kphi = 1, nphi
            XYZ_r(:) = [R(kR) * cosphi(kphi), R(kR) * sinphi(kphi), Z(kZ)]
            ! Biot-Savart integral over coil segments
            BXYZ(:) = 0d0
            XYZ_f(:) = coils(kc)%XYZ(:, coils(kc)%nseg) - XYZ_r
            dist_f = sqrt(sum(XYZ_f * XYZ_f))
            do ks = 1, coils(kc)%nseg
              XYZ_i(:) = XYZ_f
              dist_i = max(1.0d-12, sqrt(sum(XYZ_i * XYZ_i)))
              XYZ_f(:) = coils(kc)%XYZ(:, ks) - XYZ_r
              dist_f = max(1.0d-12, sqrt(sum(XYZ_f * XYZ_f)))
              denom = dist_i * dist_f * (dist_i * dist_f + sum(XYZ_i * XYZ_f))
              if (denom == 0.0d0) denom = 1.0d-12
              BXYZ(:) = BXYZ + &
                (XYZ_i([2, 3, 1]) * XYZ_f([3, 1, 2]) - XYZ_i([3, 1, 2]) * XYZ_f([2, 3, 1])) * &
                (dist_i + dist_f) / denom
            end do
            BR(kphi) = BXYZ(1) * cosphi(kphi) + BXYZ(2) * sinphi(kphi)
            Bphi(kphi) = BXYZ(2) * cosphi(kphi) - BXYZ(1) * sinphi(kphi)
            BZ(kphi) = BXYZ(3)
          end do
          !$omp end parallel do
          call fftw_execute_dft_r2c(plan_nphi, BR, BnR)
          call fftw_execute_dft_r2c(plan_nphi, Bphi, Bnphi)
          call fftw_execute_dft_r2c(plan_nphi, BZ, BnZ)
          Bn(0:nmax, 1, kR, kZ, kc) = BnR(1:nmax+1) / dble(nphi)
          Bn(0:nmax, 2, kR, kZ, kc) = Bnphi(1:nmax+1) / dble(nphi)
          Bn(0:nmax, 3, kR, kZ, kc) = BnZ(1:nmax+1) / dble(nphi)
        end do
      end do
    end do
    call fftw_destroy_plan(plan_nphi)
    call fftw_free(p_BR)
    call fftw_free(p_Bphi)
    call fftw_free(p_BZ)
    call fftw_free(p_BnR)
    call fftw_free(p_Bnphi)
    call fftw_free(p_BnZ)
    !$ call fftw_cleanup_threads()
    ! nullify pointers past this point
  end subroutine biot_savart_fourier

  subroutine vector_potential_biot_savart_fourier(coils, nmax, min_distance, max_eccentricity, use_convex_wall, &
    Rmin, Rmax, Zmin, Zmax, nR, nphi, nZ, AnR, Anphi, AnZ, dAnphi_dR, &
    dAnphi_dZ, AnX_raw, AnY_raw, AnZ_raw)
    use iso_c_binding, only: c_ptr, c_double, c_double_complex, c_size_t, c_f_pointer, c_null_ptr
    !$ use omp_lib, only: omp_get_max_threads
    use FFTW3, only: fftw_init_threads, fftw_plan_with_nthreads, fftw_cleanup_threads, &
      fftw_alloc_real, fftw_alloc_complex, fftw_plan_dft_r2c_1d, FFTW_PATIENT, &
      FFTW_DESTROY_INPUT, fftw_execute_dft_r2c, fftw_destroy_plan, fftw_free
    use field_sub, only: read_field_input, stretch_coords

    type(coil_t), intent(in), dimension(:) :: coils
    integer, intent(in) :: nmax
    real(dp), intent(in) :: min_distance, max_eccentricity
    logical, intent(in) :: use_convex_wall
    real(dp), intent(in) :: Rmin, Rmax, Zmin, Zmax
    integer, intent(in) :: nR, nphi, nZ
    complex(dp), intent(out), dimension(:, :, :, :), allocatable :: &
      AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ
    complex(dp), intent(out), dimension(:, :, :, :), allocatable, optional :: &
      AnX_raw, AnY_raw, AnZ_raw
    integer :: nfft, ncoil, kc, ks, ks_prev, kR, kphi, kZ
    real(dp), dimension(nphi) :: phi, cosphi, sinphi
    real(dp) :: R(nR), Z(nZ), actual_R, actual_Z, XYZ_r(3), XYZ_i(3), XYZ_f(3), XYZ_if(3), dist_i, dist_f, dist_if, &
      AXYZ(3), grad_AX(3), grad_AY(3), grad_AZ(3), eccentricity, common_gradient_term(3)
    type(c_ptr) :: plan_nphi, p_AR, p_Aphi, p_AZ, p_AX, p_AY, p_fft_output
    real(c_double), dimension(:), pointer :: AR, Aphi, AZ, AX, AY
    complex(c_double_complex), dimension(:), pointer :: fft_output
    real(dp), dimension(:, :, :, :), allocatable :: debug_grad_AX, debug_grad_AY, debug_grad_AZ
    logical :: export_raw

    export_raw = present(AnX_raw) .and. present(AnY_raw) .and. present(AnZ_raw)
    p_AX = c_null_ptr
    p_AY = c_null_ptr

    if (nmax > nphi / 4) then
      write (error_unit, '("biot_savart_fourier: requested nmax = ", ' // &
        'i0, ", but only ", i0, " modes available.")') nmax, nphi / 4
      error stop
    end if
    nfft = nphi / 2 + 1
    ncoil = size(coils)
    call grid_from_bounding_box(Rmin, Rmax, nR, R, Zmin, Zmax, nZ, Z, nphi, phi)
    cosphi(:) = cos(phi)
    sinphi(:) = sin(phi)
    allocate(AnR(0:nmax, nR, nZ, ncoil))
    allocate(Anphi(0:nmax, nR, nZ, ncoil))
    allocate(AnZ(0:nmax, nR, nZ, ncoil))
    if (export_raw) then
      allocate(AnX_raw(0:nmax, nR, nZ, ncoil))
      allocate(AnY_raw(0:nmax, nR, nZ, ncoil))
      allocate(AnZ_raw(0:nmax, nR, nZ, ncoil))
    end if
    allocate(dAnphi_dR(0:nmax, nR, nZ, ncoil))
    allocate(dAnphi_dZ(0:nmax, nR, nZ, ncoil))
    allocate(debug_grad_AX(3, nR, nZ, ncoil))
    allocate(debug_grad_AY(3, nR, nZ, ncoil))
    allocate(debug_grad_AZ(3, nR, nZ, ncoil))
    ! prepare FFTW
    !$ if (fftw_init_threads() == 0) error stop 'OpenMP support in FFTW could not be initialized'
    !$ call fftw_plan_with_nthreads(omp_get_max_threads())
    p_AR = fftw_alloc_real(int(nphi, c_size_t))
    call c_f_pointer(p_AR, AR, [nphi])
    p_Aphi = fftw_alloc_real(int(nphi, c_size_t))
    call c_f_pointer(p_Aphi, Aphi, [nphi])
    if (export_raw) then
      p_AX = fftw_alloc_real(int(nphi, c_size_t))
      call c_f_pointer(p_AX, AX, [nphi])
      p_AY = fftw_alloc_real(int(nphi, c_size_t))
      call c_f_pointer(p_AY, AY, [nphi])
    end if
    p_AZ = fftw_alloc_real(int(nphi, c_size_t))
    call c_f_pointer(p_AZ, AZ, [nphi])
    p_fft_output = fftw_alloc_complex(int(nfft, c_size_t))
    call c_f_pointer(p_fft_output, fft_output, [nfft])
    plan_nphi = fftw_plan_dft_r2c_1d(nphi, AR, AnR, ior(FFTW_PATIENT, FFTW_DESTROY_INPUT))

    if (use_convex_wall) then
      call read_field_input  ! read convex wall
    end if

    do kc = 1, ncoil
      write (*, '("Computig Biot-Savart field for coil ", i0, " of ", i0, "...")') kc, ncoil
      do kZ = 1, nZ
        do kR = 1, nR
          if (use_convex_wall) then
            call stretch_coords(R(kR), Z(kZ), actual_R, actual_Z)
          else
            actual_R = R(kR)
            actual_Z = Z(kZ)
          end if
          !$omp parallel do schedule(static) default(none) &
          !$omp private(kphi, ks, ks_prev, XYZ_r, XYZ_i, XYZ_f, XYZ_if, dist_i, dist_f, dist_if, &
          !$omp AXYZ, grad_AX, grad_AY, grad_AZ, eccentricity, common_gradient_term) &
          !$omp shared(nphi, kc, coils, R, kr, Z, kZ, cosphi, sinphi, AR, Aphi, AX, AY, AZ, &
          !$omp actual_R, actual_Z, min_distance, max_eccentricity, debug_grad_AX, debug_grad_AY, debug_grad_AZ, &
          !$omp export_raw)
          do kphi = 1, nphi
            XYZ_r(:) = [actual_R * cosphi(kphi), actual_R * sinphi(kphi), actual_Z]
            AXYZ(:) = 0d0
            grad_AX(:) = 0d0
            grad_AY(:) = 0d0
            grad_AZ(:) = 0d0
            ks_prev = coils(kc)%nseg
            XYZ_f(:) = coils(kc)%XYZ(:, ks_prev) - XYZ_r
            dist_f = max(min_distance, sqrt(sum(XYZ_f * XYZ_f)))
            do ks = 1, coils(kc)%nseg
              XYZ_i(:) = XYZ_f
              dist_i = dist_f
              XYZ_f(:) = coils(kc)%XYZ(:, ks) - XYZ_r
              dist_f = max(min_distance, sqrt(sum(XYZ_f * XYZ_f)))
              XYZ_if = coils(kc)%XYZ(:, ks) - coils(kc)%XYZ(:, ks_prev)
              dist_if = sqrt(sum(XYZ_if * XYZ_if))
              eccentricity = segment_eccentricity(dist_if, dist_i, dist_f, max_eccentricity)
              AXYZ(:) = AXYZ + segment_vector_potential_contribution(XYZ_if, dist_if, eccentricity)
              common_gradient_term(:) = segment_gradient_kernel(XYZ_i, XYZ_f, dist_i, dist_f)
              call segment_gradient_contribution(XYZ_if, common_gradient_term, grad_AX, grad_AY, grad_AZ)
              ks_prev = ks
            end do
            if (export_raw) then
              AX(kphi) = AXYZ(1)
              AY(kphi) = AXYZ(2)
            end if
            AR(kphi) = AXYZ(1) * cosphi(kphi) + AXYZ(2) * sinphi(kphi)
            Aphi(kphi) = AXYZ(2) * cosphi(kphi) - AXYZ(1) * sinphi(kphi)
            AZ(kphi) = AXYZ(3)
            if (kphi == 1) then
              debug_grad_AX(:, kR, kZ, kc) = grad_AX
              debug_grad_AY(:, kR, kZ, kc) = grad_AY
              debug_grad_AZ(:, kR, kZ, kc) = grad_AZ
            end if
          end do
          !$omp end parallel do
          call fftw_execute_dft_r2c(plan_nphi, AR, fft_output)
          AnR(0:nmax, kR, kZ, kc) = fft_output(1:nmax+1) / dble(nphi)
          call fftw_execute_dft_r2c(plan_nphi, Aphi, fft_output)
          Anphi(0:nmax, kR, kZ, kc) = fft_output(1:nmax+1) / dble(nphi)
          call fftw_execute_dft_r2c(plan_nphi, AZ, fft_output)
          AnZ(0:nmax, kR, kZ, kc) = fft_output(1:nmax+1) / dble(nphi)
          if (export_raw) then
            call fftw_execute_dft_r2c(plan_nphi, AX, fft_output)
            if (present(AnX_raw)) then
              AnX_raw(0:nmax, kR, kZ, kc) = fft_output(1:nmax+1) / dble(nphi)
            end if
            call fftw_execute_dft_r2c(plan_nphi, AY, fft_output)
            if (present(AnY_raw)) then
              AnY_raw(0:nmax, kR, kZ, kc) = fft_output(1:nmax+1) / dble(nphi)
            end if
            if (present(AnZ_raw)) then
              AnZ_raw(0:nmax, kR, kZ, kc) = AnZ(0:nmax, kR, kZ, kc)
            end if
          end if
        end do
      end do
    end do
    call compute_fourier_derivatives(Anphi, dAnphi_dR, dAnphi_dZ, Rmin, Rmax, nR, Zmin, Zmax, nZ)
    call write_debug_gradients(debug_grad_AX, debug_grad_AY, debug_grad_AZ)
    call fftw_destroy_plan(plan_nphi)
    call fftw_free(p_AR)
    call fftw_free(p_Aphi)
    if (export_raw) then
      call fftw_free(p_AX)
      call fftw_free(p_AY)
    end if
    call fftw_free(p_AZ)
    call fftw_free(p_fft_output)
    !$ call fftw_cleanup_threads()
    ! nullify pointers past this point
  end subroutine vector_potential_biot_savart_fourier

  subroutine compute_fourier_derivatives(Anphi, dAnphi_dR, dAnphi_dZ, Rmin, Rmax, nR, Zmin, Zmax, nZ)
    complex(dp), intent(in) :: Anphi(:, :, :, :)
    complex(dp), intent(out) :: dAnphi_dR(:, :, :, :)
    complex(dp), intent(out) :: dAnphi_dZ(:, :, :, :)
    real(dp), intent(in) :: Rmin, Rmax, Zmin, Zmax
    integer, intent(in) :: nR, nZ
    real(dp) :: deltaR, deltaZ
    integer :: mode_lo, mode_hi, nmode, kR, kZ, kcoil

    if (nR > 1) then
      deltaR = (Rmax - Rmin) / real(nR - 1, dp)
    else
      deltaR = 1.0_dp
    end if
    if (nZ > 1) then
      deltaZ = (Zmax - Zmin) / real(nZ - 1, dp)
    else
      deltaZ = 1.0_dp
    end if

    mode_lo = lbound(Anphi, 1)
    mode_hi = ubound(Anphi, 1)

    dAnphi_dR(:, :, :, :) = (0.0_dp, 0.0_dp)
    dAnphi_dZ(:, :, :, :) = (0.0_dp, 0.0_dp)

    if (nR >= 2) then
      do kcoil = 1, size(Anphi, 4)
        do kZ = 1, nZ
          do nmode = mode_lo, mode_hi
            ! Forward difference at lower boundary
            if (nR == 2) then
              dAnphi_dR(nmode, 1, kZ, kcoil) = (Anphi(nmode, 2, kZ, kcoil) - Anphi(nmode, 1, kZ, kcoil)) / deltaR
              dAnphi_dR(nmode, 2, kZ, kcoil) = dAnphi_dR(nmode, 1, kZ, kcoil)
            else
              dAnphi_dR(nmode, 1, kZ, kcoil) = (-3.0_dp * Anphi(nmode, 1, kZ, kcoil) + &
                4.0_dp * Anphi(nmode, 2, kZ, kcoil) - Anphi(nmode, 3, kZ, kcoil)) / (2.0_dp * deltaR)
              do kR = 2, nR - 1
                dAnphi_dR(nmode, kR, kZ, kcoil) = (Anphi(nmode, kR + 1, kZ, kcoil) - &
                  Anphi(nmode, kR - 1, kZ, kcoil)) / (2.0_dp * deltaR)
              end do
              dAnphi_dR(nmode, nR, kZ, kcoil) = (3.0_dp * Anphi(nmode, nR, kZ, kcoil) - &
                4.0_dp * Anphi(nmode, nR - 1, kZ, kcoil) + Anphi(nmode, nR - 2, kZ, kcoil)) / (2.0_dp * deltaR)
            end if
          end do
        end do
      end do
    end if

    if (nZ >= 2) then
      do kcoil = 1, size(Anphi, 4)
        do kR = 1, nR
          do nmode = mode_lo, mode_hi
            if (nZ == 2) then
              dAnphi_dZ(nmode, kR, 1, kcoil) = (Anphi(nmode, kR, 2, kcoil) - Anphi(nmode, kR, 1, kcoil)) / deltaZ
              dAnphi_dZ(nmode, kR, 2, kcoil) = dAnphi_dZ(nmode, kR, 1, kcoil)
            else
              dAnphi_dZ(nmode, kR, 1, kcoil) = (-3.0_dp * Anphi(nmode, kR, 1, kcoil) + &
                4.0_dp * Anphi(nmode, kR, 2, kcoil) - Anphi(nmode, kR, 3, kcoil)) / (2.0_dp * deltaZ)
              do kZ = 2, nZ - 1
                dAnphi_dZ(nmode, kR, kZ, kcoil) = (Anphi(nmode, kR, kZ + 1, kcoil) - &
                  Anphi(nmode, kR, kZ - 1, kcoil)) / (2.0_dp * deltaZ)
              end do
              dAnphi_dZ(nmode, kR, nZ, kcoil) = (3.0_dp * Anphi(nmode, kR, nZ, kcoil) - &
                4.0_dp * Anphi(nmode, kR, nZ - 1, kcoil) + Anphi(nmode, kR, nZ - 2, kcoil)) / (2.0_dp * deltaZ)
            end if
          end do
        end do
      end do
    end if
  end subroutine compute_fourier_derivatives

  subroutine write_debug_gradients(grad_AX, grad_AY, grad_AZ)
    use, intrinsic :: iso_fortran_env, only: output_unit
    real(dp), dimension(:, :, :, :), intent(in) :: grad_AX, grad_AY, grad_AZ
    integer :: fid
    open (newunit=fid, file='debug_cartesian_gradients.dat', status='replace', action='write', form='unformatted')
    write (fid) grad_AX
    write (fid) grad_AY
    write (fid) grad_AZ
    close (fid)
    write (output_unit, '(a)') 'Wrote debug Cartesian gradients to debug_cartesian_gradients.dat'
  end subroutine write_debug_gradients

  subroutine sum_coils_gauge_single_mode_Anvac(AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ, &
    Ic, ntor, Rmin, Rmax, nR, Zmin, Zmax, nZ, gauged_AnR, gauged_AnZ)
    use math_constants, only: imun => IMUN
    complex(dp), dimension(:, :, :, :), intent(in) :: AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ
    real(dp), dimension(:), intent(in) :: Ic
    integer, intent(in) :: ntor
    real(dp), intent(in) :: Rmin, Rmax
    integer, intent(in) :: nR
    real(dp), intent(in) :: Zmin, Zmax
    integer, intent(in) :: nZ
    complex(dp), dimension(:, :), allocatable, intent(inout) :: gauged_AnR, gauged_AnZ
    integer :: ncoil, kcoil, kZ
    real(dp) :: R(nR), Z(nZ)

    ! TODO: check array sizes
    ncoil = size(Ic)
    call grid_from_bounding_box(Rmin, Rmax, nR, R, Zmin, Zmax, nZ, Z)
    if (allocated(gauged_AnR)) deallocate(gauged_AnR)
    if (allocated(gauged_AnZ)) deallocate(gauged_AnZ)
    allocate(gauged_AnR(nR, nZ), gauged_AnZ(nR, nZ))
    gauged_AnR(:, :) = (0d0, 0d0)
    gauged_AnZ(:, :) = (0d0, 0d0)
    ! TODO: OpenMP?
    do kcoil = 1, ncoil
      do kZ = 1, nZ
        if (ntor == 0) then
          gauged_AnR(:, kZ) = gauged_AnR(:, kZ) + Ic(kcoil) * AnR(ntor, :, kZ, kcoil)
          gauged_AnZ(:, kZ) = gauged_AnZ(:, kZ) + Ic(kcoil) * AnZ(ntor, :, kZ, kcoil)
        else
          gauged_AnR(:, kZ) = gauged_AnR(:, kZ) + Ic(kcoil) * &
            (AnR(ntor, :, kZ, kcoil) + imun / ntor * Anphi(ntor, :, kZ, kcoil))
          gauged_AnZ(:, kZ) = gauged_AnZ(:, kZ) + Ic(kcoil) * AnZ(ntor, :, kZ, kcoil)
        end if
      end do
    end do
  end subroutine sum_coils_gauge_single_mode_Anvac

  subroutine gauged_Anvac_from_Bnvac(BnR, BnZ, ntor, Rmin, Rmax, nR, Zmin, Zmax, nZ, &
    gauged_AnR, gauged_AnZ)
    use math_constants, only: imun => IMUN
    complex(dp), dimension(:, :), intent(in) :: BnR, BnZ
    integer, intent(in) :: ntor
    real(dp), intent(in) :: Rmin, Rmax
    integer, intent(in) :: nR
    real(dp), intent(in) :: Zmin, Zmax
    integer, intent(in) :: nZ
    complex(dp), dimension(:, :), allocatable, intent(inout) :: gauged_AnR, gauged_AnZ
    integer :: kZ
    real(dp) :: R(nR), Z(nZ)

    ! TODO: check array sizes
    call grid_from_bounding_box(Rmin, Rmax, nR, R, Zmin, Zmax, nZ, Z)
    if (allocated(gauged_AnR)) deallocate(gauged_AnR)
    if (allocated(gauged_AnZ)) deallocate(gauged_AnZ)
    allocate(gauged_AnR(nR, nZ), gauged_AnZ(nR, nZ))
    gauged_AnR(:, :) = (0d0, 0d0)
    gauged_AnZ(:, :) = (0d0, 0d0)
    do kZ = 1, nZ
      gauged_AnR(:, kZ) = imun / ntor * R * BnZ(:, kZ)
      gauged_AnZ(:, kZ) = -imun / ntor * R * BnR(:, kZ)
    end do
  end subroutine gauged_Anvac_from_Bnvac

  subroutine write_Anvac_fourier(filename, ncoil, nmax, &
    Rmin, Rmax, Zmin, Zmax, nR, nphi, nZ, AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ, &
    AnX_raw, AnY_raw, AnZ_raw)
    use netcdf
    character(len = *), intent(in) :: filename
    integer, intent(in) :: ncoil, nmax
    real(dp), intent(in) :: Rmin, Rmax, Zmin, Zmax
    integer, intent(in) :: nR, nphi, nZ
    complex(dp), intent(in), dimension(:, :, :, :) :: AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ, &
      AnX_raw, AnY_raw, AnZ_raw
    real(dp) :: R(nR), Z(nZ)
    integer :: coil_number(ncoil), ntor(nmax + 1)
    integer :: ncid
    integer :: dimid_R, dimid_Z, dimid_tor, dimid_coil
    integer :: varid_R, varid_Z, varid_ntor, varid_coils, varid_nR, varid_nphi, varid_nZ
    integer :: k

    call grid_from_bounding_box(Rmin, Rmax, nR, R, Zmin, Zmax, nZ, Z)
    coil_number = [(k, k = 1, ncoil)]
    ntor = [(k, k = 0, nmax)]

    call nc_check('create', nf90_create(filename, NF90_NETCDF4, ncid))

    ! define dimensions metadata
    call nc_check('def_dim', nf90_def_dim(ncid, 'R', nR, dimid_R))
    call nc_check('def_dim', nf90_def_dim(ncid, 'Z', nZ, dimid_Z))
    call nc_check('def_dim', nf90_def_dim(ncid, 'ntor', nmax+1, dimid_tor))
    call nc_check('def_dim', nf90_def_dim(ncid, 'coil_number', ncoil, dimid_coil))

    ! define variables metadata
    call nc_check('def_var', nf90_def_var(ncid, 'R', NF90_DOUBLE, [dimid_R], varid_R))
    call nc_check('def_var', nf90_def_var(ncid, 'Z', NF90_DOUBLE, [dimid_Z], varid_Z))
    call nc_check('def_var', nf90_def_var(ncid, 'ntor', NF90_INT, [dimid_tor], varid_ntor))
    call nc_check('def_var', nf90_def_var(ncid, 'coil_number', NF90_INT, [dimid_coil], varid_coils))
    call nc_check('def_var', nf90_def_var(ncid, 'nR', NF90_INT, varid_nR))
    call nc_check('def_var', nf90_def_var(ncid, 'nphi', NF90_INT, varid_nphi))
    call nc_check('def_var', nf90_def_var(ncid, 'nZ', NF90_INT, varid_nZ))

    ! write variables and comments metadata
    call nc_check('put_var', nf90_put_var(ncid, varid_R, R))
    call nc_check('put_att', nf90_put_att(ncid, varid_R, 'comment', 'R components of grid in cm'))
    call nc_check('put_var', nf90_put_var(ncid, varid_Z, Z))
    call nc_check('put_att', nf90_put_att(ncid, varid_Z, 'comment', 'Z components of grid in cm'))
    call nc_check('put_var', nf90_put_var(ncid, varid_ntor, ntor))
    call nc_check('put_att', nf90_put_att(ncid, varid_ntor, 'comment', 'toroidal mode numbers'))
    call nc_check('put_var', nf90_put_var(ncid, varid_coils, coil_number))
    call nc_check('put_att', nf90_put_att(ncid, varid_coils, 'comment', 'coil numbers'))
    call nc_check('put_var', nf90_put_var(ncid, varid_nR, nR))
    call nc_check('put_att', nf90_put_att(ncid, varid_nR, 'comment', 'number of grid points in R direction'))
    call nc_check('put_var', nf90_put_var(ncid, varid_nphi, nphi))
    call nc_check('put_att', nf90_put_att(ncid, varid_nphi, 'comment', 'number of grid points in phi direction'))
    call nc_check('put_var', nf90_put_var(ncid, varid_nZ, nZ))
    call nc_check('put_att', nf90_put_att(ncid, varid_nZ, 'comment', 'number of grid points in Z direction'))

    ! process actual data
    call write_actual_data(AnR, 'AnR', 'R')
    call write_actual_data(Anphi, 'Anphi', 'phi')
    call write_actual_data(AnZ, 'AnZ', 'Z')
    call write_actual_data(AnX_raw, 'AnX_raw', 'ungauged Cartesian X')
    call write_actual_data(AnY_raw, 'AnY_raw', 'ungauged Cartesian Y')
    call write_actual_data(AnZ_raw, 'AnZ_raw', 'ungauged Cartesian Z')
    call write_actual_data(dAnphi_dR, 'dAnphi_dR', 'derivative w.r.t. R of phi')
    call write_actual_data(dAnphi_dZ, 'dAnphi_dZ', 'derivative w.r.t. Z of phi')

    call nc_check('close', nf90_close(ncid))

  contains
    subroutine write_actual_data(var, name, component)
      complex(dp), intent(in), dimension(:, :, :, :) :: var
      character(len = *), intent(in) :: name
      character(len = *), intent(in) :: component
      integer :: varid_actual_data

      call nc_check('def_var', nf90_def_var(ncid, name // '_real', NF90_DOUBLE, &
        [dimid_tor, dimid_R, dimid_Z, dimid_coil], varid_actual_data))
      call nc_check('put_var', nf90_put_var(ncid, varid_actual_data, var%Re))
      call nc_check('put_att', nf90_put_att(ncid, varid_actual_data, 'comment', &
        'real part of toroidal Fourier mode of ' // component // &
        ' component of vector potential'))
      call nc_check('def_var', nf90_def_var(ncid, name // '_imag', NF90_DOUBLE, &
        [dimid_tor, dimid_R, dimid_Z, dimid_coil], varid_actual_data))
      call nc_check('put_var', nf90_put_var(ncid, varid_actual_data, var%Im))
      call nc_check('put_att', nf90_put_att(ncid, varid_actual_data, 'comment', &
        'imaginary part of toroidal Fourier mode of ' // component // &
        ' component of vector potential'))
    end subroutine write_actual_data
  end subroutine write_Anvac_fourier

  subroutine read_Anvac_fourier(filename, ncoil, nmax, &
    Rmin, Rmax, Zmin, Zmax, nR, nphi, nZ, AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ, &
    AnX_raw, AnY_raw, AnZ_raw)
    use netcdf
    character(len = *), intent(in) :: filename
    integer, intent(out) :: ncoil, nmax
    real(dp), intent(out) :: Rmin, Rmax, Zmin, Zmax
    integer, intent(out) :: nR, nphi, nZ
    complex(dp), dimension(:, :, :, :), allocatable, intent(inout) :: &
      AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ, AnX_raw, AnY_raw, AnZ_raw
    integer :: ncid
    integer :: dimid_R, dimid_Z, dimid_tor, dimid_coil
    integer :: varid_R, varid_Z, varid_nphi
    real(dp), dimension(:), allocatable :: R, Z

    call nc_check('open', nf90_open(filename, NF90_NOWRITE, ncid))

    call nc_check('inq_dimid', nf90_inq_dimid(ncid, 'R', dimid_R))
    call nc_check('inq_dimid', nf90_inq_dimid(ncid, 'Z', dimid_Z))
    call nc_check('inq_dimid', nf90_inq_dimid(ncid, 'ntor', dimid_tor))
    call nc_check('inq_dimid', nf90_inq_dimid(ncid, 'coil_number', dimid_coil))

    call nc_check('inquire_dimension', &
      nf90_inquire_dimension(ncid, dimid_R, len = nR))
    call nc_check('inquire_dimension', &
      nf90_inquire_dimension(ncid, dimid_Z, len = nZ))
    call nc_check('inquire_dimension', &
      nf90_inquire_dimension(ncid, dimid_tor, len = nmax))
    nmax = nmax - 1
    call nc_check('inquire_dimension', &
      nf90_inquire_dimension(ncid, dimid_coil, len = ncoil))

    call nc_check('inq_varid', nf90_inq_varid(ncid, 'nphi', varid_nphi))
    call nc_check('get_var', nf90_get_var(ncid, varid_nphi, nphi))

    call nc_check('inq_varid', nf90_inq_varid(ncid, 'R', varid_R))
    call nc_check('inq_varid', nf90_inq_varid(ncid, 'Z', varid_Z))
    allocate(R(nR), Z(nZ))
    call nc_check('get_var', nf90_get_var(ncid, varid_R, R))
    call nc_check('get_var', nf90_get_var(ncid, varid_Z, Z))
    Rmin = R(1)
    Rmax = R(nR)
    Zmin = Z(1)
    Zmax = Z(nZ)

    ! process actual data
    call read_actual_data(AnR, 'AnR')
    call read_actual_data(Anphi, 'Anphi')
    call read_actual_data(AnZ, 'AnZ')
    call read_optional_data(AnX_raw, 'AnX_raw')
    call read_optional_data(AnY_raw, 'AnY_raw')
    call read_optional_data(AnZ_raw, 'AnZ_raw', AnZ)
    call read_actual_data(dAnphi_dR, 'dAnphi_dR')
    call read_actual_data(dAnphi_dZ, 'dAnphi_dZ')

    call nc_check('close', nf90_close(ncid))
    deallocate(R, Z)

  contains
    subroutine read_optional_data(var, name, fallback)
      complex(dp), dimension(:, :, :, :), allocatable, intent(inout) :: var
      character(len = *), intent(in) :: name
      complex(dp), dimension(:, :, :, :), intent(in), optional :: fallback
      integer :: status_real, status_imag, varid_actual_data

      allocate(var(0:nmax, nR, nZ, ncoil))
      status_real = nf90_inq_varid(ncid, name // '_real', varid_actual_data)
      if (status_real /= NF90_NOERR) then
        if (present(fallback)) then
          var(:, :, :, :) = fallback
        else
          var(:, :, :, :) = (0d0, 0d0)
        end if
        return
      end if

      call nc_check('get_var', nf90_get_var(ncid, varid_actual_data, var%Re))
      status_imag = nf90_inq_varid(ncid, name // '_imag', varid_actual_data)
      if (status_imag == NF90_NOERR) then
        call nc_check('get_var', nf90_get_var(ncid, varid_actual_data, var%Im))
      else
        var%Im = 0d0
      end if
    end subroutine read_optional_data

    subroutine read_actual_data(var, name)
      complex(dp), dimension(:, :, :, :), allocatable, intent(inout) :: var
      character(len = *), intent(in) :: name
      integer :: varid_actual_data

      allocate(var(0:nmax, nR, nZ, ncoil))
      call nc_check('inq_varid', nf90_inq_varid(ncid, name // '_real', varid_actual_data))
      call nc_check('get_var', nf90_get_var(ncid, varid_actual_data, var%Re))
      call nc_check('inq_varid', nf90_inq_varid(ncid, name // '_imag', varid_actual_data))
      call nc_check('get_var', nf90_get_var(ncid, varid_actual_data, var%Im))
    end subroutine read_actual_data
  end subroutine read_Anvac_fourier

  subroutine nc_check(operation, status)
    use netcdf, only: NF90_NOERR, nf90_strerror
    character(len=*), intent(in) :: operation
    integer, intent(in) :: status

    if (status == NF90_NOERR) return
    write (*, '("Error encountered during nf90_", a, ": ", a)') operation, nf90_strerror(status)
    error stop
  end subroutine nc_check

  subroutine write_Bvac_nemov(filename, Rmin, Rmax, Zmin, Zmax, Bvac)
    use math_constants, only: pi
    character(len = *), intent(in) :: filename
    real(dp), intent(in) :: Rmin, Rmax, Zmin, Zmax
    real(dp), intent(in), dimension(:, :, :, :) :: Bvac
    integer :: fid, nR, nphi, nZ, kR, kphi, kZ

    if (size(Bvac, 1) /= 3) then
      write (error_unit, arg_size_fmt) 'write_Bnvac_nemov', 'size(Bvac, 1)', size(Bvac, 1), '3', 3
      error stop
    end if
    nZ = size(Bvac, 2)
    nphi = size(Bvac, 3)
    nR = size(Bvac, 4)
    open(newunit = fid, file = filename, status = 'replace', action = 'write', form = 'formatted')
    ! kisslinger_asdex.f90 writes points at phi = 0 and at phi = 2 pi
    write (fid, '(4(i0, :, 1x))') nR, nphi + 1, nZ, 1
    write (fid, '(2(es24.16e3, :, 1x))') Rmin, Rmax
    write (fid, '(2(es24.16e3, :, 1x))') 0d0, 2d0 * pi
    write (fid, '(2(es24.16e3, :, 1x))') Zmin, Zmax
    do kR = 1, nR
      do kphi = 1, nphi
        do kZ = 1, nZ
          write (fid, '(3(es24.16e3, :, 1x))') Bvac(:, kZ, kphi, kR)
        end do
      end do
      ! kisslinger_asdex.f90 writes points at phi = 0 and at phi = 2 pi
      kphi = 1
      do kZ = 1, nZ
        write (fid, '(3(es24.16e3, :, 1x))') Bvac(:, kZ, kphi, kR)
      end do
    end do
    close(fid)
  end subroutine write_Bvac_nemov

  subroutine write_Bnvac_fourier(filename, Bn, Rmin, Rmax, Zmin, Zmax)
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    character(len = *), intent(in) :: filename
    complex(dp), dimension(0:, :, :, :, :), intent(in) :: Bn
    real(dp), intent(in) :: Rmin, Rmax, Zmin, Zmax
    integer :: nR, nZ, nmax, ncoil, ntor, kc
    integer(HID_T) :: h5id_root
    character(len = 7) :: modename, coilname

    if (size(Bn, 2) /= 3) then
      write (error_unit, arg_size_fmt) 'write_Bnvac_fourier', &
        'size(Bn, 2)', size(Bn, 2), '3', 3
      error stop
    end if
    nmax = ubound(Bn, 1)
    nR = size(Bn, 3)
    nZ = size(Bn, 4)
    ncoil = size(Bn, 5)
    call h5_open_rw(filename, h5id_root)
    do ntor = 0, nmax
      write (modename, '("ntor_", i2.2)') ntor
      call h5_create_parent_groups(h5id_root, modename // '/')
      call h5_add(h5id_root, modename // '/R_min', Rmin, &
        unit = 'cm', comment = 'minimal R coordinate of computational grid')
      call h5_add(h5id_root, modename // '/R_max', Rmax, &
        unit = 'cm', comment = 'maximal R coordinate of computational grid')
      call h5_add(h5id_root, modename // '/Z_min', Zmin, &
        unit = 'cm', comment = 'minimal Z coordinate of computational grid')
      call h5_add(h5id_root, modename // '/Z_max', Zmax, &
        unit = 'cm', comment = 'maximal Z coordinate of computational grid')
      call h5_add(h5id_root, modename // '/nR', nR, 'number of grid points in R direction')
      call h5_add(h5id_root, modename // '/nZ', nZ, 'number of grid points in Z direction')
      call h5_add(h5id_root, modename // '/ncoil', ncoil, 'number of coils')
      do kc = 1, ncoil
        write (coilname, '("coil_", i2.2)') kc
        call h5_create_parent_groups(h5id_root, modename // '/' // coilname // '/')
        call h5_add(h5id_root, modename // '/' // coilname // '/Bn_R', &
          Bn(ntor, 1, :, :, kc), [1, 1], [nR, nZ], &
          unit = 'G', comment = 'R component of coil field for I_c = 0.1 A')
        call h5_add(h5id_root, modename // '/' // coilname // '/Bn_phi', &
          Bn(ntor, 2, :, :, kc), [1, 1], [nR, nZ], &
          unit = 'G', comment = 'physical phi component of coil field for I_c = 0.1 A')
        call h5_add(h5id_root, modename // '/' // coilname // '/Bn_Z', &
          Bn(ntor, 3, :, :, kc), [1, 1], [nR, nZ], &
          unit = 'G', comment = 'Z component of coil field for I_c = 0.1 A')
      end do
    end do
    call h5_close(h5id_root)
  end subroutine write_Bnvac_fourier

  subroutine read_Bnvac_fourier(filename, ntor, Ic, nR, nZ, Rmin, Rmax, Zmin, Zmax, &
    Bnvac_R, Bnvac_Z)
    use hdf5_tools, only: HID_T, h5_open, h5_get, h5_close
    character(len = *), intent(in) :: filename
    integer, intent(in) :: ntor
    real(dp), intent(in), dimension(:) :: Ic
    integer, intent(out) :: nR, nZ
    real(dp), intent(out) :: Rmin, Rmax, Zmin, Zmax
    complex(dp), intent(out), dimension(:, :), allocatable :: Bnvac_R, Bnvac_Z
    logical :: file_exists
    integer :: kc, ncoil
    integer(HID_T) :: h5id_root
    character(len = 7) :: modename, coilname
    complex(dp), dimension(:, :), allocatable :: Bn

    inquire(file = filename, exist = file_exists)
    if (.not. file_exists) then
      write (error_unit, '("read_Bnvac_fourier: File ", a, " not found.")') trim(filename)
      error stop
    end if
    call h5_open(filename, h5id_root)
    write (modename, '("ntor_", i2.2)') ntor
    call h5_get(h5id_root, modename // '/R_min', Rmin)
    call h5_get(h5id_root, modename // '/R_max', Rmax)
    call h5_get(h5id_root, modename // '/Z_min', Zmin)
    call h5_get(h5id_root, modename // '/Z_max', Zmax)
    call h5_get(h5id_root, modename // '/nR', nR)
    call h5_get(h5id_root, modename // '/nZ', nZ)
    call h5_get(h5id_root, modename // '/ncoil', ncoil)
    if (ncoil /= size(Ic)) then
      write (error_unit, arg_size_fmt) 'read_Bnvac_fourier', &
        'ncoil', ncoil, 'size(Ic)', size(Ic)
      error stop
    end if
    allocate(Bnvac_R(nR, nZ), Bnvac_Z(nR, nZ), Bn(nR, nZ))
    Bnvac_R(:, :) = (0d0, 0d0)
    Bnvac_Z(:, :) = (0d0, 0d0)
    do kc = 1, ncoil
      write (coilname, '("coil_", i2.2)') kc
      call h5_get(h5id_root, modename // '/' // coilname // '/Bn_R', Bn)
      Bnvac_R(:, :) = Bnvac_R + Ic(kc) * Bn
      call h5_get(h5id_root, modename // '/' // coilname // '/Bn_Z', Bn)
      Bnvac_Z(:, :) = Bnvac_Z + Ic(kc) * Bn
    end do
    call h5_close(h5id_root)
    deallocate(Bn)
  end subroutine read_Bnvac_fourier

  ! Small helper functions for Biot-Savart segment contributions
  ! These are extracted for testing and validation

  pure function segment_vector_potential_contribution(XYZ_segment, dist_segment, eccentricity) result(dA)
    real(dp), intent(in) :: XYZ_segment(3)
    real(dp), intent(in) :: dist_segment
    real(dp), intent(in) :: eccentricity
    real(dp) :: dA(3)
    dA = XYZ_segment / dist_segment * log((1.0_dp + eccentricity) / (1.0_dp - eccentricity))
  end function segment_vector_potential_contribution

  pure function segment_eccentricity(dist_segment, dist_i, dist_f, max_ecc) result(ecc)
    real(dp), intent(in) :: dist_segment, dist_i, dist_f, max_ecc
    real(dp) :: ecc
    ecc = min(max_ecc, dist_segment / (dist_i + dist_f))
  end function segment_eccentricity

  pure function segment_gradient_kernel(XYZ_i, XYZ_f, dist_i, dist_f) result(grad_kernel)
    real(dp), intent(in) :: XYZ_i(3), XYZ_f(3)
    real(dp), intent(in) :: dist_i, dist_f
    real(dp) :: grad_kernel(3)
    grad_kernel = (XYZ_i / dist_i + XYZ_f / dist_f) / (dist_i * dist_f + sum(XYZ_i * XYZ_f))
  end function segment_gradient_kernel

  pure subroutine segment_gradient_contribution(XYZ_segment, grad_kernel, grad_AX, grad_AY, grad_AZ)
    real(dp), intent(in) :: XYZ_segment(3)
    real(dp), intent(in) :: grad_kernel(3)
    real(dp), intent(inout) :: grad_AX(3), grad_AY(3), grad_AZ(3)
    grad_AX = grad_AX + XYZ_segment(1) * grad_kernel
    grad_AY = grad_AY + XYZ_segment(2) * grad_kernel
    grad_AZ = grad_AZ + XYZ_segment(3) * grad_kernel
  end subroutine segment_gradient_contribution

end module coil_tools
