module coil_tools

  use iso_fortran_env, only: dp => real64, error_unit

  implicit none

  private

  public :: AUG_coils_read, AUG_coils_write_Nemov, AUG_coils_read_Nemov, &
       AUG_coils_write_GPEC, AUG_coils_read_GPEC, AUG_coils_write_Fourier, &
       read_currents_Nemov, Biot_Savart_sum_coils, write_Bvac_Nemov

  character(len = *), parameter :: arg_size_fmt = &
       '("Argument size mismatch in ", a, ": ", a, " = ", i0, ", ", a, " = ", i0)'

contains

  function linspace(lo, hi, cnt, excl_lo, excl_hi)
    real(dp), intent(in) :: lo, hi
    integer, intent(in) :: cnt, excl_lo, excl_hi
    real(dp) :: linspace(cnt)
    real(dp) :: step
    integer :: k

    step = (hi - lo) / dble(cnt - 1 + excl_lo + excl_hi)
    linspace = lo + [(k * step, k = excl_lo, cnt - 1 + excl_lo)]
    if (excl_hi == 0) linspace(cnt) = hi
  end function linspace

  subroutine AUG_coils_read(directory, ncoil, nseg, nwind, XYZ)
    character(len = *), intent(in) :: directory
    integer, intent(out) :: ncoil, nseg, nwind
    real(dp), intent(out), dimension(:, :, :), allocatable :: XYZ
    character(len = 8) :: filename
    integer :: fid, status, kc, ks
    real(dp), dimension(:, :), allocatable :: R_Z_phi

    ncoil = 8
    nwind = 5
    ! count number of lines = number of coil segments
    filename = 'Bu1n.asc'
    open(newunit = fid, file = directory // '/' // filename, status = 'old', &
         action = 'read', form = 'formatted')
    nseg = 0
    do
       read(fid, *, iostat = status)
       if (status /= 0) exit
       nseg = nseg + 1
    end do
    close(fid)
    ! allocate coordinates and read all coil data, starting with upper B coil set
    allocate(R_Z_phi(nseg, 3))
    allocate(XYZ(3, nseg, 2 * ncoil))
    do kc = 1, ncoil
       write (filename, '("Bu", i1, "n.asc")') kc
       open(newunit = fid, file = directory // '/' // filename, status = 'old', &
            action = 'read', form = 'formatted')
       do ks = 1, nseg
          read (fid, '(f8.6, 1x, f8.5, 1x, f8.5)') R_Z_phi(ks, 1), R_Z_phi(ks, 2), R_Z_phi(ks, 3)
       end do
       close(fid)
       XYZ(1, :, kc) = 1d2 * R_Z_phi(:, 1) * cos(R_Z_phi(:, 3))
       XYZ(2, :, kc) = 1d2 * R_Z_phi(:, 1) * sin(R_Z_phi(:, 3))
       XYZ(3, :, kc) = 1d2 * R_Z_phi(:, 2)
    end do
    do kc = 1, ncoil
       write (filename, '("Bl", i1, "n.asc")') kc
       open(newunit = fid, file = directory // '/' // filename, status = 'old', &
            action = 'read', form = 'formatted')
       do ks = 1, nseg
          read (fid, '(f8.6, 1x, f8.5, 1x, f8.5)') R_Z_phi(ks, 1), R_Z_phi(ks, 2), R_Z_phi(ks, 3)
       end do
       close(fid)
       XYZ(1, :, kc + ncoil) = 1d2 * R_Z_phi(:, 1) * cos(R_Z_phi(:, 3))
       XYZ(2, :, kc + ncoil) = 1d2 * R_Z_phi(:, 1) * sin(R_Z_phi(:, 3))
       XYZ(3, :, kc + ncoil) = 1d2 * R_Z_phi(:, 2)
    end do
    deallocate(R_Z_phi)
  end subroutine AUG_coils_read

  subroutine AUG_coils_write_Nemov(directory, ncoil, nseg, XYZ)
    character(len = *), intent(in) :: directory
    integer, intent(in) :: ncoil, nseg
    real(dp), intent(in), dimension(:, :, :) :: XYZ
    integer :: fid, kc, ks

    if (3 /= size(XYZ, 1)) then
       write (error_unit, arg_size_fmt) 'AUG_coils_write_Nemov', '3', 'size(XYZ, 1)', 3, size(XYZ, 1)
       error stop
    end if
    if (nseg /= size(XYZ, 2)) then
       write (error_unit, arg_size_fmt) 'AUG_coils_write_Nemov', 'nseg', 'size(XYZ, 2)', nseg, size(XYZ, 2)
       error stop
    end if
    if (2 * ncoil /= size(XYZ, 3)) then
       write (error_unit, arg_size_fmt) 'AUG_coils_write_Nemov', '2 * ncoil', 'size(XYZ, 3)', &
            2 * ncoil, size(XYZ, 3)
       error stop
    end if
    open(newunit = fid, file = directory // '/co_asd.dd', status = 'replace', &
         action = 'write', form = 'formatted')
    write (fid, '(1x, i6)') 2 * ncoil * (nseg + 1)
    do kc = 1, 2 * ncoil
       do ks = 1, nseg
          write (fid, '(3(1x, es17.9e2), 1x, es12.4e2, 1x, i3)') &
               XYZ(1, ks, kc), XYZ(2, ks, kc), XYZ(3, ks, kc), 1.0, kc
       end do
       write (fid, '(3(1x, es17.9e2), 1x, es12.4e2, 1x, i3)') &
            XYZ(1, 1, kc), XYZ(2, 1, kc), XYZ(3, 1, kc), 0.0, kc
    end do
    close(fid)
  end subroutine AUG_coils_write_Nemov

  subroutine AUG_coils_read_Nemov(directory, ncoil, nseg, nwind, XYZ)
    character(len = *), intent(in) :: directory
    integer, intent(out) :: ncoil, nseg, nwind
    real(dp), intent(out), allocatable, dimension(:, :, :) :: XYZ
    integer :: fid, kc, ks, temp
    real(dp) :: rdum

    ncoil = 8
    nwind = 5
    open(newunit = fid, file = directory // '/co_asd.dd', status = 'old', &
         action = 'read', form = 'formatted')
    read (fid, *) temp
    nseg = temp / (2 * ncoil) - 1
    allocate(XYZ(3, nseg, 2 * ncoil))
    do kc = 1, 2 * ncoil
       do ks = 1, nseg
          read (fid, *) XYZ(1, ks, kc), XYZ(2, ks, kc), XYZ(3, ks, kc), rdum, temp
          if (temp /= kc) then
             write (error_unit, '("Expected coil index ", i0, " in co_asd.dd at line ", ' // &
                  'i0, ", but got ", i0)') kc, (kc - 1) * (nseg + 1) + ks + 1, temp
             error stop
          end if
       end do
       read (fid, *)
    end do
    close(fid)
  end subroutine AUG_coils_read_Nemov

  subroutine AUG_coils_write_GPEC(directory, ncoil, nseg, nwind, XYZ)
    character(len = *), intent(in) :: directory
    integer, intent(in) :: ncoil, nseg, nwind
    real(dp), intent(in), dimension(:, :, :) :: XYZ
    integer :: fid, kc, ks

    if (3 /= size(XYZ, 1)) then
       write (error_unit, arg_size_fmt) 'AUG_coils_write_GPEC', '3', 'size(XYZ, 1)', 3, size(XYZ, 1)
       error stop
    end if
    if (nseg /= size(XYZ, 2)) then
       write (error_unit, arg_size_fmt) 'AUG_coils_write_GPEC', 'nseg', 'size(XYZ, 2)', nseg, size(XYZ, 2)
       error stop
    end if
    if (2 * ncoil /= size(XYZ, 3)) then
       write (error_unit, arg_size_fmt) 'AUG_coils_write_GPEC', '2 * ncoil', 'size(XYZ, 3)', &
            2 * ncoil, size(XYZ, 3)
       error stop
    end if
    open(newunit = fid, file = directory // '/aug_bu.dat', status = 'replace', &
         action = 'write', form = 'formatted')
    write (fid, '(3(1x, i4), 1x, f7.2)') ncoil, 1, nseg + 1, dble(nwind)
    do kc = 1, ncoil
       do ks = 1, nseg
          write (fid, '(3(1x, es12.4e2))') &
               1d-2 * XYZ(1, ks, kc), 1d-2 * XYZ(2, ks, kc), 1d-2 * XYZ(3, ks, kc)
       end do
       write (fid, '(3(1x, es12.4e2))') &
            1d-2 * XYZ(1, 1, kc), 1d-2 * XYZ(2, 1, kc), 1d-2 * XYZ(3, 1, kc)
    end do
    close(fid)
    open(newunit = fid, file = directory // '/aug_bl.dat', status = 'replace', &
         action = 'write', form = 'formatted')
    write (fid, '(3(1x, i4), 1x, f7.2)') ncoil, 1, nseg + 1, dble(nwind)
    do kc = ncoil + 1, 2 * ncoil
       do ks = 1, nseg
          write (fid, '(3(1x, es12.4e2))') &
               1d-2 * XYZ(1, ks, kc), 1d-2 * XYZ(2, ks, kc), 1d-2 * XYZ(3, ks, kc)
       end do
       write (fid, '(3(1x, es12.4e2))') &
            1d-2 * XYZ(1, 1, kc), 1d-2 * XYZ(2, 1, kc), 1d-2 * XYZ(3, 1, kc)
    end do
    close(fid)
  end subroutine AUG_coils_write_GPEC

  subroutine AUG_coils_read_GPEC(directory, ncoil, nseg, nwind, XYZ)
    character(len = *), intent(in) :: directory
    integer, intent(out) :: ncoil, nseg, nwind
    real(dp), intent(out), allocatable, dimension(:, :, :) :: XYZ
    integer :: fid, kc, ks, idum
    real(dp) :: ddum

    open(newunit = fid, file = directory // '/aug_bu.dat', status = 'old', &
         action = 'read', form = 'formatted')
    read (fid, '(3(1x, i4), 1x, f7.2)') ncoil, idum, nseg, ddum
    nseg = nseg - 1
    nwind = int(ddum)
    allocate(XYZ(3, nseg, 2 * ncoil))
    do kc = 1, ncoil
       do ks = 1, nseg
          read (fid, '(3(1x, es12.4e2))') XYZ(1, ks, kc), XYZ(2, ks, kc), XYZ(3, ks, kc)
       end do
       read (fid, *)
    end do
    close(fid)
    open(newunit = fid, file = directory // '/aug_bl.dat', status = 'old', &
         action = 'read', form = 'formatted')
    read (fid, *)
    do kc = ncoil + 1, 2 * ncoil
       do ks = 1, nseg
          read (fid, '(3(1x, es12.4e2))') XYZ(1, ks, kc), XYZ(2, ks, kc), XYZ(3, ks, kc)
       end do
       read (fid, *)
    end do
    close(fid)
    XYZ(:, :, :) = 1d2 * XYZ
  end subroutine AUG_coils_read_GPEC

  subroutine AUG_coils_write_Fourier(directory, ncoil, nseg, nwind, XYZ, nmax, &
       Rmin, Rmax, Zmin, Zmax, nR, nphi, nZ)
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    character(len = *), intent(in) :: directory
    integer, intent(in) :: ncoil, nseg, nwind
    real(dp), intent(in), dimension(:, :, :) :: XYZ
    integer, intent(in) :: nmax
    real(dp), intent(in) :: Rmin, Rmax, Zmin, Zmax
    integer, intent(in) :: nR, nphi, nZ
    complex(dp), dimension(:, :, :, :, :), allocatable :: Bn
    character(len = 1024) :: filename
    integer :: ntor, kc
    integer(HID_T) :: h5id_root
    character(len = 7) :: modename, coilname

    if (3 /= size(XYZ, 1)) then
       write (error_unit, arg_size_fmt) 'AUG_coils_write_GPEC', '3', 'size(XYZ, 1)', 3, size(XYZ, 1)
       error stop
    end if
    if (nseg /= size(XYZ, 2)) then
       write (error_unit, arg_size_fmt) 'AUG_coils_write_GPEC', 'nseg', 'size(XYZ, 2)', nseg, size(XYZ, 2)
       error stop
    end if
    if (2 * ncoil /= size(XYZ, 3)) then
       write (error_unit, arg_size_fmt) 'AUG_coils_write_GPEC', '2 * ncoil', 'size(XYZ, 3)', &
            2 * ncoil, size(XYZ, 3)
       error stop
    end if
    if (nmax > nphi / 4) then
       write (error_unit, '("Requested nmax = ", i0, ", but only ", i0, " modes available.")') &
            nmax, nphi / 4
       error stop
    end if
    call Biot_Savart_Fourier(ncoil, nseg, nwind, XYZ, nmax, &
         Rmin, Rmax, Zmin, Zmax, nR, nphi, nZ, Bn)
    filename = directory // '/AUG_B_coils.h5'
    call h5_open_rw(trim(filename), h5id_root)
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
       call h5_add(h5id_root, modename // '/ncoil', ncoil, 'number of upper and lower B coils')
       do kc = 1, 2 * ncoil
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
    deallocate(Bn)
  end subroutine AUG_coils_write_Fourier

  subroutine read_currents_Nemov(directory, Ic)
    character(len = *), intent(in) :: directory
    real(dp), intent(out), allocatable :: Ic(:)
    integer, parameter :: nwind = 5, ncoil = 8
    integer :: fid

    allocate(Ic(2 * ncoil))
    open(newunit = fid, file = directory // '/cur_asd.dd', status = 'old', action = 'read', form = 'formatted')
    read (fid, *) Ic(:)
    close(fid)
    Ic(:) = Ic / real(nwind)
  end subroutine read_currents_Nemov

  subroutine Biot_Savart_sum_coils(ncoil, nseg, nwind, XYZ_c, Ic, &
       Rmin, Rmax, Zmin, Zmax, nR, nphi, nZ, Bvac)
    use math_constants, only: pi
    integer, intent(in) :: ncoil, nseg, nwind
    real(dp), intent(in), dimension(:, :, :) :: XYZ_c
    real(dp), intent(in), dimension(:) :: Ic
    real(dp), intent(in) :: Rmin, Rmax, Zmin, Zmax
    integer, intent(in) :: nR, nphi, nZ
    real(dp), intent(out), dimension(:, :, :, :), allocatable :: Bvac
    integer :: kc, ks, kR, kphi, kZ
    real(dp), dimension(nphi) :: phi, cosphi, sinphi
    real(dp) :: R(nR), Z(nZ), XYZ_r(3), XYZ_i(3), XYZ_f(3), dist_i, dist_f, BXYZ_c(3), BXYZ(3)

    if (3 /= size(XYZ_c, 1)) then
       write (error_unit, arg_size_fmt) 'Biot_Savart_sum_coils', &
            '3', 'size(XYZ_c, 1)', 3, size(XYZ_c, 1)
       error stop
    end if
    if (nseg /= size(XYZ_c, 2)) then
       write (error_unit, arg_size_fmt) 'Biot_Savart_sum_coils', &
            'nseg', 'size(XYZ_c, 2)', nseg, size(XYZ_c, 2)
       error stop
    end if
    if (2 * ncoil /= size(XYZ_c, 3)) then
       write (error_unit, arg_size_fmt) 'Biot_Savart_sum_coils', &
            '2 * ncoil', 'size(XYZ_c, 3)', 2 * ncoil, size(XYZ_c, 3)
       error stop
    end if
    if (2 * ncoil /= size(Ic)) then
       write (error_unit, arg_size_fmt) 'Biot_Savart_sum_coils', '2 * ncoil', 'size(Ic)', &
            2 * ncoil, size(Ic)
       error stop
    end if
    R(:) = linspace(Rmin, Rmax, nR, 0, 0)
    phi(:) = linspace(0d0, 2d0 * pi, nphi, 0, 1)  ! half-open interval: do not repeat phi = 0 at phi = 2 pi
    cosphi(:) = cos(phi)
    sinphi(:) = sin(phi)
    Z(:) = linspace(Zmin, Zmax, nZ, 0, 0)
    allocate(Bvac(3, nZ, nphi, nR))
    Bvac(:, :, :, :) = 0d0
    !$omp parallel do schedule(static) collapse(3) default(none) &
    !$omp private(kr, kphi, kZ, kc, ks, XYZ_r, XYZ_i, XYZ_f, dist_i, dist_f, BXYZ) &
    !$omp shared(nR, nphi, nZ, ncoil, nseg, R, Z, cosphi, sinphi, XYZ_c, nwind, Ic, Bvac)
    do kR = 1, nR
       do kphi = 1, nphi
          do kZ = 1, nZ
             XYZ_r(:) = [R(kR) * cosphi(kphi), R(kR) * sinphi(kphi), Z(kZ)]
             ! Biot-Savart integral over coil segments
             BXYZ(:) = 0d0
             do kc = 1, 2 * ncoil
                BXYZ_c(:) = 0d0
                XYZ_f(:) = XYZ_c(:, nseg, kc) - XYZ_r
                dist_f = sqrt(sum(XYZ_f * XYZ_f))
                do ks = 1, nseg
                   XYZ_i(:) = XYZ_f
                   dist_i = dist_f
                   XYZ_f(:) = XYZ_c(:, ks, kc) - XYZ_r
                   dist_f = sqrt(sum(XYZ_f * XYZ_f))
                   BXYZ_c(:) = BXYZ_c + &
                        (XYZ_i([2, 3, 1]) * XYZ_f([3, 1, 2]) - XYZ_i([3, 1, 2]) * XYZ_f([2, 3, 1])) * &
                        (dist_i + dist_f) / (dist_i * dist_f * (dist_i * dist_f + sum(XYZ_i * XYZ_f)))
                end do
                BXYZ(:) = BXYZ + nwind * Ic(kc) * BXYZ_c
             end do
             Bvac(1, kZ, kphi, kR) = BXYZ(1) * cosphi(kphi) + BXYZ(2) * sinphi(kphi)
             Bvac(2, kZ, kphi, kR) = BXYZ(2) * cosphi(kphi) - BXYZ(1) * sinphi(kphi)
             Bvac(3, kZ, kphi, kR) = BXYZ(3)
          end do
       end do
    end do
  end subroutine Biot_Savart_sum_coils

  subroutine Biot_Savart_Fourier(ncoil, nseg, nwind, XYZ_c, nmax, &
       Rmin, Rmax, Zmin, Zmax, nR, nphi, nZ, Bn)
    use iso_c_binding, only: c_ptr, c_double, c_double_complex, c_size_t, c_f_pointer
    !$ use omp_lib, only: omp_get_max_threads
    use FFTW3, only: fftw_init_threads, fftw_plan_with_nthreads, fftw_cleanup_threads, &
         fftw_alloc_real, fftw_alloc_complex, fftw_plan_dft_r2c_1d, FFTW_PATIENT, &
         FFTW_DESTROY_INPUT, fftw_execute_dft_r2c, fftw_destroy_plan, fftw_free
    use math_constants, only: pi
    integer, intent(in) :: ncoil, nseg, nwind
    real(dp), intent(in), dimension(:, :, :) :: XYZ_c
    integer, intent(in) :: nmax
    real(dp), intent(in) :: Rmin, Rmax, Zmin, Zmax
    integer, intent(in) :: nR, nphi, nZ
    complex(dp), intent(out), dimension(:, :, :, :, :), allocatable :: Bn
    integer :: nfft, kc, ks, kR, kphi, kZ
    real(dp), dimension(nphi) :: phi, cosphi, sinphi
    real(dp) :: R(nR), Z(nZ), XYZ_r(3), XYZ_i(3), XYZ_f(3), dist_i, dist_f, BXYZ(3)
    type(c_ptr) :: plan_nphi, p_BR, p_Bphi, p_BZ, p_BnR, p_Bnphi, p_BnZ
    real(c_double), dimension(:), pointer :: BR, Bphi, BZ
    complex(c_double_complex), dimension(:), pointer :: BnR, Bnphi, BnZ

    if (3 /= size(XYZ_c, 1)) then
       write (error_unit, arg_size_fmt) 'Biot_Savart_Fourier', &
            '3', 'size(XYZ_c, 1)', 3, size(XYZ_c, 1)
       error stop
    end if
    if (nseg /= size(XYZ_c, 2)) then
       write (error_unit, arg_size_fmt) 'Biot_Savart_Fourier', &
            'nseg', 'size(XYZ_c, 2)', nseg, size(XYZ_c, 2)
       error stop
    end if
    if (2 * ncoil /= size(XYZ_c, 3)) then
       write (error_unit, arg_size_fmt) 'Biot_Savart_Fourier', &
            '2 * ncoil', 'size(XYZ_c, 3)', 2 * ncoil, size(XYZ_c, 3)
       error stop
    end if
    if (nmax > nphi / 4) then
       write (error_unit, '("Requested nmax = ", i0, ", but only ", i0, " modes available.")') &
            nmax, nphi / 4
       error stop
    end if
    nfft = nphi / 2 + 1
    R(:) = linspace(Rmin, Rmax, nR, 0, 0)
    phi(:) = linspace(0d0, 2d0 * pi, nphi, 0, 1)  ! half-open interval: do not repeat phi = 0 at phi = 2 pi
    cosphi(:) = cos(phi)
    sinphi(:) = sin(phi)
    Z(:) = linspace(Zmin, Zmax, nZ, 0, 0)
    allocate(Bn(0:nmax, 3, nR, nZ, 2 * ncoil))
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
    do kc = 1, 2 * ncoil
       do kZ = 1, nZ
          XYZ_r(3) = Z(kZ)
          do kR = 1, nR
             !$omp parallel do schedule(static) default(none) &
             !$omp private(kphi, ks, XYZ_i, XYZ_f, dist_i, dist_f, BXYZ) firstprivate(XYZ_r) &
             !$omp shared(nphi, nseg, kc, XYZ_c, R, kR, cosphi, sinphi, BR, Bphi, BZ)
             do kphi = 1, nphi
                XYZ_r(1:2) = [R(kR) * cosphi(kphi), R(kR) * sinphi(kphi)]
                ! Biot-Savart integral over coil segments
                XYZ_f(:) = XYZ_c(:, nseg, kc) - XYZ_r
                dist_f = sqrt(sum(XYZ_f * XYZ_f))
                do ks = 1, nseg
                   XYZ_i(:) = XYZ_f
                   dist_i = dist_f
                   XYZ_f(:) = XYZ_c(:, ks, kc) - XYZ_r
                   dist_f = sqrt(sum(XYZ_f * XYZ_f))
                   BXYZ(:) = BXYZ + &
                        (XYZ_i([2, 3, 1]) * XYZ_f([3, 1, 2]) - XYZ_i([3, 1, 2]) * XYZ_f([2, 3, 1])) * &
                        (dist_i + dist_f) / (dist_i * dist_f * (dist_i * dist_f + sum(XYZ_i * XYZ_f)))
                end do
                BR(kphi) = BXYZ(1) * cosphi(kphi) + BXYZ(2) * sinphi(kphi)
                Bphi(kphi) = BXYZ(2) * cosphi(kphi) - BXYZ(1) * sinphi(kphi)
                BZ(kphi) = BXYZ(3)
             end do
             call fftw_execute_dft_r2c(plan_nphi, BR, BnR)
             call fftw_execute_dft_r2c(plan_nphi, Bphi, Bnphi)
             call fftw_execute_dft_r2c(plan_nphi, BZ, BnZ)
             Bn(0:nmax, 1, kR, kZ, kc) = nwind * BnR(1:nmax+1) / dble(nphi)
             Bn(0:nmax, 2, kR, kZ, kc) = nwind * Bnphi(1:nmax+1) / dble(nphi)
             Bn(0:nmax, 3, kR, kZ, kc) = nwind * BnZ(1:nmax+1) / dble(nphi)
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

  subroutine write_Bvac_Nemov(directory, Rmin, Rmax, Zmin, Zmax, Bvac)
    use math_constants, only: pi
    character(len = *), intent(in) :: directory
    real(dp), intent(in) :: Rmin, Rmax, Zmin, Zmax
    real(dp), intent(in), dimension(:, :, :, :) :: Bvac
    integer :: fid, nR, nphi, nZ, kR, kphi, kZ

    if (3 /= size(Bvac, 1)) then
       write (error_unit, arg_size_fmt) 'write_Bnvac_Nemov', '3', 'size(Bvac, 1)', 3, size(Bvac, 1)
       error stop
    end if
    nZ = size(Bvac, 2)
    nphi = size(Bvac, 3)
    nR = size(Bvac, 4)
    open(newunit = fid, file = directory // '/field.dat', status = 'replace', &
         action = 'write', form = 'formatted')
    ! kisslinger_asdex.f90 writes points at phi = 0 and at phi = 2 pi
    write (fid, '(i0, 1x, i0, 1x, i0, 1x, i0)') nR, nphi + 1, nZ, 1
    write (fid, '(es24.16e3, 1x, es24.16e3)') Rmin, Rmax
    write (fid, '(es24.16e3, 1x, es24.16e3)') 0d0, 2d0 * pi
    write (fid, '(es24.16e3, 1x, es24.16e3)') Zmin, Zmax
    do kR = 1, nR
       do kphi = 1, nphi
          do kZ = 1, nZ
             write (fid, '(es24.16e3, 1x, es24.16e3, 1x, es24.16e3)') Bvac(:, kZ, kphi, kR)
          end do
       end do
       ! kisslinger_asdex.f90 writes points at phi = 0 and at phi = 2 pi
       kphi = 1
       do kZ = 1, nZ
          write (fid, '(es24.16e3, 1x, es24.16e3, 1x, es24.16e3)') Bvac(:, kZ, kphi, kR)
       end do
    end do
    close(fid)
  end subroutine write_Bvac_Nemov

  subroutine read_Bnvac_Fourier(directory, ntor, Ic, nR, nZ, Rmin, Rmax, Zmin, Zmax, &
       Bnvac_R, Bnvac_Z)
    use hdf5_tools, only: HID_T, h5_open, h5_get, h5_close
    character(len = *), intent(in) :: directory
    integer, intent(in) :: ntor
    real(dp), intent(in), dimension(:) :: Ic
    integer, intent(out) :: nR, nZ
    real(dp), intent(out) :: Rmin, Rmax, Zmin, Zmax
    complex(dp), intent(out), dimension(:, :), allocatable :: Bnvac_R, Bnvac_Z
    character(len = 1024) :: filename
    logical :: file_exists
    integer :: kc, ncoil
    integer(HID_T) :: h5id_root
    character(len = 7) :: modename, coilname
    complex(dp), dimension(:, :), allocatable :: Bn

    filename = trim(adjustl(directory)) // '/AUG_B_coils.h5'
    inquire(file = filename, exist = file_exists)
    if (.not. file_exists) then
       write (error_unit, '("File ", a, " not found, cannot read vacuum perturbation field.")') &
         trim(filename)
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
    if (2 * ncoil /= size(Ic)) then
       write (error_unit, arg_size_fmt) 'read_Bnvac_Fourier', '2 * ncoil', 'size(Ic)', &
            2 * ncoil, size(Ic)
       error stop
    end if
    allocate(Bnvac_R(nR, nZ), Bnvac_Z(nR, nZ), Bn(nR, nZ))
    Bnvac_R(:, :) = (0d0, 0d0)
    Bnvac_Z(:, :) = (0d0, 0d0)
    do kc = 1, 2 * ncoil
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
