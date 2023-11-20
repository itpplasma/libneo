module coil_tools

  use iso_fortran_env, only: dp => real64, error_unit

  implicit none

  private

  public :: coil_t, coil_init, coil_deinit, coils_append, &
    process_fixed_number_of_args, check_number_of_args, &
    coils_write_AUG, coils_read_AUG, &
    coils_write_Nemov, coils_read_Nemov, &
    coils_write_GPEC, coils_read_GPEC, &
    read_currents, Biot_Savart_sum_coils, Biot_Savart_Fourier, &
    write_Bvac_Nemov, write_Bnvac_Fourier, read_Bnvac_Fourier

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
    linspace = lo + [(k * step, k = omit_lo, cnt - 1 + omit_lo)]
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
    character(len = :), dimension(:), allocatable, intent(out) :: args
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

  subroutine coils_write_Nemov(filename, coils)
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
  end subroutine coils_write_Nemov

  subroutine coils_read_Nemov(filename, coils)
    character(len = *), intent(in) :: filename
    type(coil_t), intent(out), allocatable, dimension(:) :: coils
    integer :: fid, total, k, ncoil, kc, nseg, ks, idum
    real(dp) :: X, Y, Z, cur

    open(newunit = fid, file = filename, status = 'old', &
      action = 'read', form = 'formatted')
    read (fid, *) total
    ncoil = 0
    do k = 1, total
      read (fid, *) X, Y, Z, cur, kc
      ncoil = max(ncoil, kc)
    end do
    allocate(coils(ncoil))
    rewind fid
    read (fid, *) total
    nseg = 0
    do k = 1, total
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
    do kc = 1, ncoil
      do ks = 1, coils(kc)%nseg
        read (fid, *) coils(kc)%XYZ(:, ks), cur, idum
        k = k + 1
        if (idum /= kc) then
          write (error_unit, '("coils_read_Nemov: expected coil index ", ' // &
            'i0, " in ", a, " at line ", i0, ", but got ", i0)') kc, filename, k, idum
          error stop
        end if
      end do
      read (fid, *) X, Y, Z, cur, idum
    end do
    close(fid)
  end subroutine coils_read_Nemov

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

  subroutine Biot_Savart_sum_coils(coils, Ic, &
    Rmin, Rmax, Zmin, Zmax, nR, nphi, nZ, Bvac)
    use math_constants, only: pi
    type(coil_t), intent(in), dimension(:) :: coils
    real(dp), intent(in), dimension(:) :: Ic
    real(dp), intent(in) :: Rmin, Rmax, Zmin, Zmax
    integer, intent(in) :: nR, nphi, nZ
    real(dp), intent(out), dimension(:, :, :, :), allocatable :: Bvac
    integer :: ncoil, kc, ks, kR, kphi, kZ
    real(dp), dimension(nphi) :: phi, cosphi, sinphi
    real(dp) :: R(nR), Z(nZ), XYZ_r(3), XYZ_i(3), XYZ_f(3), dist_i, dist_f, BXYZ_c(3), BXYZ(3)

    if (size(coils) /= size(Ic)) then
      write (error_unit, arg_size_fmt) 'Biot_Savart_sum_coils', &
        'size(coils)', size(coils), 'size(Ic)', size(Ic)
      error stop
    end if
    ncoil = size(coils)
    R(:) = linspace(Rmin, Rmax, nR, 0, 0)
    phi(:) = linspace(0d0, 2d0 * pi, nphi, 0, 1)  ! half-open interval: do not repeat phi = 0 at phi = 2 pi
    cosphi(:) = cos(phi)
    sinphi(:) = sin(phi)
    Z(:) = linspace(Zmin, Zmax, nZ, 0, 0)
    allocate(Bvac(3, nZ, nphi, nR))
    Bvac(:, :, :, :) = 0d0
    !$omp parallel do schedule(static) collapse(3) default(none) &
    !$omp private(kr, kphi, kZ, kc, ks, XYZ_r, XYZ_i, XYZ_f, dist_i, dist_f, BXYZ, BXYZ_c) &
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
              dist_i = dist_f
              XYZ_f(:) = coils(kc)%XYZ(:, ks) - XYZ_r
              dist_f = sqrt(sum(XYZ_f * XYZ_f))
              BXYZ_c(:) = BXYZ_c + &
                (XYZ_i([2, 3, 1]) * XYZ_f([3, 1, 2]) - XYZ_i([3, 1, 2]) * XYZ_f([2, 3, 1])) * &
                (dist_i + dist_f) / (dist_i * dist_f * (dist_i * dist_f + sum(XYZ_i * XYZ_f)))
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
  end subroutine Biot_Savart_sum_coils

  subroutine Biot_Savart_Fourier(coils, nmax, &
    Rmin, Rmax, Zmin, Zmax, nR, nphi, nZ, Bn)
    use iso_c_binding, only: c_ptr, c_double, c_double_complex, c_size_t, c_f_pointer
    !$ use omp_lib, only: omp_get_max_threads
    use FFTW3, only: fftw_init_threads, fftw_plan_with_nthreads, fftw_cleanup_threads, &
      fftw_alloc_real, fftw_alloc_complex, fftw_plan_dft_r2c_1d, FFTW_PATIENT, &
      FFTW_DESTROY_INPUT, fftw_execute_dft_r2c, fftw_destroy_plan, fftw_free
    use math_constants, only: pi
    type(coil_t), intent(in), dimension(:) :: coils
    integer, intent(in) :: nmax
    real(dp), intent(in) :: Rmin, Rmax, Zmin, Zmax
    integer, intent(in) :: nR, nphi, nZ
    complex(dp), intent(out), dimension(:, :, :, :, :), allocatable :: Bn
    integer :: nfft, ncoil, kc, ks, kR, kphi, kZ
    real(dp), dimension(nphi) :: phi, cosphi, sinphi
    real(dp) :: R(nR), Z(nZ), XYZ_r(3), XYZ_i(3), XYZ_f(3), dist_i, dist_f, BXYZ(3)
    type(c_ptr) :: plan_nphi, p_BR, p_Bphi, p_BZ, p_BnR, p_Bnphi, p_BnZ
    real(c_double), dimension(:), pointer :: BR, Bphi, BZ
    complex(c_double_complex), dimension(:), pointer :: BnR, Bnphi, BnZ

    if (nmax > nphi / 4) then
      write (error_unit, '("Biot_Savart_Fourier: requested nmax = ", ' // &
        'i0, ", but only ", i0, " modes available.")') nmax, nphi / 4
      error stop
    end if
    nfft = nphi / 2 + 1
    R(:) = linspace(Rmin, Rmax, nR, 0, 0)
    phi(:) = linspace(0d0, 2d0 * pi, nphi, 0, 1)  ! half-open interval: do not repeat phi = 0 at phi = 2 pi
    cosphi(:) = cos(phi)
    sinphi(:) = sin(phi)
    Z(:) = linspace(Zmin, Zmax, nZ, 0, 0)
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
          !$omp private(kphi, ks, XYZ_r, XYZ_i, XYZ_f, dist_i, dist_f, BXYZ) &
          !$omp shared(nphi, kc, coils, R, kR, Z, kZ, cosphi, sinphi, BR, Bphi, BZ)
          do kphi = 1, nphi
            XYZ_r(:) = [R(kR) * cosphi(kphi), R(kR) * sinphi(kphi), Z(kZ)]
            ! Biot-Savart integral over coil segments
            BXYZ(:) = 0d0
            XYZ_f(:) = coils(kc)%XYZ(:, coils(kc)%nseg) - XYZ_r
            dist_f = sqrt(sum(XYZ_f * XYZ_f))
            do ks = 1, coils(kc)%nseg
              XYZ_i(:) = XYZ_f
              dist_i = dist_f
              XYZ_f(:) = coils(kc)%XYZ(:, ks) - XYZ_r
              dist_f = sqrt(sum(XYZ_f * XYZ_f))
              BXYZ(:) = BXYZ + &
                (XYZ_i([2, 3, 1]) * XYZ_f([3, 1, 2]) - XYZ_i([3, 1, 2]) * XYZ_f([2, 3, 1])) * &
                (dist_i + dist_f) / (dist_i * dist_f * (dist_i * dist_f + sum(XYZ_i * XYZ_f)))
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
  end subroutine Biot_Savart_Fourier

  subroutine write_Bvac_Nemov(filename, Rmin, Rmax, Zmin, Zmax, Bvac)
    use math_constants, only: pi
    character(len = *), intent(in) :: filename
    real(dp), intent(in) :: Rmin, Rmax, Zmin, Zmax
    real(dp), intent(in), dimension(:, :, :, :) :: Bvac
    integer :: fid, nR, nphi, nZ, kR, kphi, kZ

    if (size(Bvac, 1) /= 3) then
      write (error_unit, arg_size_fmt) 'write_Bnvac_Nemov', 'size(Bvac, 1)', size(Bvac, 1), '3', 3
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
  end subroutine write_Bvac_Nemov

  subroutine write_Bnvac_Fourier(filename, Bn, Rmin, Rmax, Zmin, Zmax)
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    character(len = *), intent(in) :: filename
    complex(dp), dimension(0:, :, :, :, :), intent(in) :: Bn
    real(dp), intent(in) :: Rmin, Rmax, Zmin, Zmax
    integer :: nR, nZ, nmax, ncoil, ntor, kc
    integer(HID_T) :: h5id_root
    character(len = 7) :: modename, coilname

    if (size(Bn, 2) /= 3) then
      write (error_unit, arg_size_fmt) 'write_Bnvac_Fourier', &
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
  end subroutine write_Bnvac_Fourier

  subroutine read_Bnvac_Fourier(filename, ntor, Ic, nR, nZ, Rmin, Rmax, Zmin, Zmax, &
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
      write (error_unit, arg_size_fmt) 'read_Bnvac_Fourier', &
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
  end subroutine read_Bnvac_Fourier

end module coil_tools
