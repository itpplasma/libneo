module geqdsk_tools

  use iso_fortran_env, only : dp => real64
  use math_constants, only : pi, length_si_to_cgs, magneticflux_si_to_cgs, &
       magneticfield_si_to_cgs, clight => C, current_si_to_cgs, &
       pressure_si_to_cgs

  implicit none

  private

  public :: cocos_t, geqdsk_t, geqdsk_deinit, geqdsk_read, geqdsk_write, &
       geqdsk_check_consistency, geqdsk_GS_prefactor, geqdsk_classify, geqdsk_standardise

  type :: cocos_t
     integer :: exp_Bpol, sgn_cyl, sgn_dpsi, sgn_Btor, sgn_Itor, &
          sgn_F, sgn_q, sgn_Bpol, sgn_pol, index
  end type cocos_t

  type :: geqdsk_t
     type(cocos_t) :: cocos
     character(len = 1024) :: fname
     character(len = 48) :: header
     integer :: nw, nh, nbbbs, limitr
     real(dp) :: rdim, zdim, rcentr, rleft, zmid, rmaxis, zmaxis, simag, sibry, bcentr, current
     real(dp), dimension(:), allocatable :: fpol, pres, ffprim, pprime, qpsi, &
          rbbbs, zbbbs, rlim, zlim, psi_eqd, R_eqd, Z_eqd, fprime
     real(dp), dimension(:, :), allocatable :: psirz
  end type geqdsk_t

  character(len = *), parameter :: format_2000 = '(a48, 3i4)'
  character(len = *), parameter :: format_2020 = '(5es16.9)'
  character(len = *), parameter :: format_2022 = '(2i5)'

  character(len = *), parameter :: &
       incons_fmt = '("Signs of ", a, " and ", a, " are inconsistent.")', &
       invert_fmt = '("Inverting sign of ", a, "...")', &
       unscaled_fmt = '("Flux in ", a, " is not normalized by 2 pi.")', &
       rescale_fmt = '("Rescaling ", a, "...")'

contains

  subroutine geqdsk_deinit(geqdsk)
    type(geqdsk_t), intent(inout) :: geqdsk

    geqdsk%fname = ''
    if (allocated(geqdsk%fpol)) deallocate(geqdsk%fpol)
    if (allocated(geqdsk%pres)) deallocate(geqdsk%pres)
    if (allocated(geqdsk%ffprim)) deallocate(geqdsk%ffprim)
    if (allocated(geqdsk%pprime)) deallocate(geqdsk%pprime)
    if (allocated(geqdsk%qpsi)) deallocate(geqdsk%qpsi)
    if (allocated(geqdsk%rbbbs)) deallocate(geqdsk%rbbbs)
    if (allocated(geqdsk%zbbbs)) deallocate(geqdsk%zbbbs)
    if (allocated(geqdsk%rlim)) deallocate(geqdsk%rlim)
    if (allocated(geqdsk%zlim)) deallocate(geqdsk%zlim)
    if (allocated(geqdsk%psi_eqd)) deallocate(geqdsk%psi_eqd)
    if (allocated(geqdsk%R_eqd)) deallocate(geqdsk%R_eqd)
    if (allocated(geqdsk%Z_eqd)) deallocate(geqdsk%Z_eqd)
    if (allocated(geqdsk%fprime)) deallocate(geqdsk%fprime)
    if (allocated(geqdsk%psirz)) deallocate(geqdsk%psirz)
  end subroutine geqdsk_deinit

  subroutine geqdsk_read(geqdsk, fname)
    type(geqdsk_t), intent(inout) :: geqdsk
    character(len = *), intent(in) :: fname
    integer :: fid, kw, kh, idum
    real(dp) :: xdum

    call geqdsk_deinit(geqdsk)
    geqdsk%fname = fname
    open(newunit = fid, file = geqdsk%fname, status = 'old', form = 'formatted', action = 'read')
    read (fid, format_2000) geqdsk%header, idum, geqdsk%nw, geqdsk%nh
    allocate(geqdsk%fpol(geqdsk%nw))
    allocate(geqdsk%pres(geqdsk%nw))
    allocate(geqdsk%ffprim(geqdsk%nw))
    allocate(geqdsk%pprime(geqdsk%nw))
    allocate(geqdsk%qpsi(geqdsk%nw))
    allocate(geqdsk%psirz(geqdsk%nw, geqdsk%nh))
    read (fid, format_2020) geqdsk%rdim, geqdsk%zdim, geqdsk%rcentr, geqdsk%rleft, geqdsk%zmid
    read (fid, format_2020) geqdsk%rmaxis, geqdsk%zmaxis, geqdsk%simag, geqdsk%sibry, geqdsk%bcentr
    read (fid, format_2020) geqdsk%current, (xdum, kw = 1, 4)
    read (fid, format_2020) (xdum, kw = 1, 5)
    read (fid, format_2020) (geqdsk%fpol(kw), kw = 1, geqdsk%nw)
    read (fid, format_2020) (geqdsk%pres(kw), kw = 1, geqdsk%nw)
    read (fid, format_2020) (geqdsk%ffprim(kw), kw = 1, geqdsk%nw)
    read (fid, format_2020) (geqdsk%pprime(kw), kw = 1, geqdsk%nw)
    read (fid, format_2020) ((geqdsk%psirz(kw, kh), kw = 1, geqdsk%nw), kh = 1, geqdsk%nh)
    read (fid, format_2020) (geqdsk%qpsi(kw), kw = 1, geqdsk%nw)
    read (fid, format_2022) geqdsk%nbbbs, geqdsk%limitr
    allocate(geqdsk%rbbbs(geqdsk%nbbbs))
    allocate(geqdsk%zbbbs(geqdsk%nbbbs))
    allocate(geqdsk%rlim(geqdsk%limitr))
    allocate(geqdsk%zlim(geqdsk%limitr))
    read (fid, format_2020) (geqdsk%rbbbs(kw), geqdsk%zbbbs(kw), kw = 1, geqdsk%nbbbs)
    read (fid, format_2020) (geqdsk%rlim(kw), geqdsk%zlim(kw), kw = 1, geqdsk%limitr)
    close(fid)

    geqdsk%rdim = geqdsk%rdim * length_si_to_cgs
    geqdsk%zdim = geqdsk%zdim * length_si_to_cgs
    geqdsk%rcentr = geqdsk%rcentr * length_si_to_cgs
    geqdsk%rleft = geqdsk%rleft * length_si_to_cgs
    geqdsk%zmid = geqdsk%zmid * length_si_to_cgs
    geqdsk%rmaxis = geqdsk%rmaxis * length_si_to_cgs
    geqdsk%zmaxis = geqdsk%zmaxis * length_si_to_cgs

    geqdsk%simag = geqdsk%simag * magneticflux_si_to_cgs
    geqdsk%sibry = geqdsk%sibry * magneticflux_si_to_cgs

    geqdsk%bcentr = geqdsk%bcentr * magneticfield_si_to_cgs

    geqdsk%current = geqdsk%current * current_si_to_cgs

    geqdsk%fpol = geqdsk%fpol * magneticfield_si_to_cgs * length_si_to_cgs

    geqdsk%pres = geqdsk%pres * pressure_si_to_cgs

    geqdsk%ffprim = geqdsk%ffprim * magneticfield_si_to_cgs

    geqdsk%pprime = geqdsk%pprime * pressure_si_to_cgs/magneticflux_si_to_cgs

    geqdsk%psirz = geqdsk%psirz * magneticflux_si_to_cgs

    geqdsk%rbbbs = geqdsk%rbbbs * length_si_to_cgs
    geqdsk%zbbbs = geqdsk%zbbbs * length_si_to_cgs
    geqdsk%rlim = geqdsk%rlim * length_si_to_cgs
    geqdsk%zlim = geqdsk%zlim * length_si_to_cgs
    ! initialize arrays equidistant values
    allocate(geqdsk%psi_eqd(geqdsk%nw))
    geqdsk%psi_eqd(:) = geqdsk%simag + [(kw, kw = 0, geqdsk%nw - 1)] * (geqdsk%sibry - geqdsk%simag) / dble(geqdsk%nw - 1)
    allocate(geqdsk%R_eqd(geqdsk%nw))
    geqdsk%R_eqd(:) = geqdsk%rleft + [(kw, kw = 0, geqdsk%nw - 1)] * geqdsk%rdim / dble(geqdsk%nw - 1)
    allocate(geqdsk%Z_eqd(geqdsk%nh))
    geqdsk%Z_eqd(:) = geqdsk%zmid - 0.5d0 * geqdsk%zdim + [(kh, kh = 0, geqdsk%nh - 1)] * geqdsk%zdim / dble(geqdsk%nh - 1)
    ! cache repeatedly used values
    allocate(geqdsk%fprime(geqdsk%nw))
    geqdsk%fprime(:) = geqdsk%ffprim / geqdsk%fpol
  end subroutine geqdsk_read

  subroutine geqdsk_write(geqdsk, fname)
    type(geqdsk_t), intent(inout) :: geqdsk
    character(len = *), intent(in) :: fname
    integer :: fid, kw, kh, idum
    real(dp) :: xdum

    idum = 0
    xdum = 0.0d0
    geqdsk%fname = fname
    open(newunit = fid, file = geqdsk%fname, status = 'replace', form = 'formatted', action = 'write')
    write (fid, format_2000) geqdsk%header, idum, geqdsk%nw, geqdsk%nh
    write (fid, format_2020) geqdsk%rdim * 1.0d-2, geqdsk%zdim * 1.0d-2, &
         geqdsk%rcentr * 1.0d-2, geqdsk%rleft * 1.0d-2, geqdsk%zmid * 1.0d-2
    write (fid, format_2020) geqdsk%rmaxis * 1.0d-2, geqdsk%zmaxis * 1.0d-2, &
         geqdsk%simag * 1.0d-8, geqdsk%sibry * 1.0d-8, geqdsk%bcentr * 1.0d-4
    write (fid, format_2020) geqdsk%current * 1.0d1 / clight, geqdsk%simag * 1.0d-8, xdum, geqdsk%rmaxis * 1.0d-2, xdum
    write (fid, format_2020) geqdsk%zmaxis * 1.0d-2, xdum, geqdsk%sibry * 1.0d-8, xdum, xdum
    write (fid, format_2020) (geqdsk%fpol(kw) * 1.0d-6, kw = 1, geqdsk%nw)
    write (fid, format_2020) (geqdsk%pres(kw) * 1.0d-1, kw = 1, geqdsk%nw)
    write (fid, format_2020) (geqdsk%ffprim(kw) * 1.0d-4, kw = 1, geqdsk%nw)
    write (fid, format_2020) (geqdsk%pprime(kw) * 1.0d7, kw = 1, geqdsk%nw)
    write (fid, format_2020) ((geqdsk%psirz(kw, kh) * 1.0d-8, kw = 1, geqdsk%nw), kh = 1, geqdsk%nh)
    write (fid, format_2020) (geqdsk%qpsi(kw), kw = 1, geqdsk%nw)
    write (fid, format_2022) geqdsk%nbbbs, geqdsk%limitr
    write (fid, format_2020) (geqdsk%rbbbs(kw) * 1.0d-2, geqdsk%zbbbs(kw) * 1.0d-2, kw = 1, geqdsk%nbbbs)
    write (fid, format_2020) (geqdsk%rlim(kw) * 1.0d-2, geqdsk%zlim(kw) * 1.0d-2, kw = 1, geqdsk%limitr)
    close(fid)
  end subroutine geqdsk_write

  subroutine geqdsk_check_consistency(geqdsk)
    type(geqdsk_t), intent(inout) :: geqdsk
    integer, parameter :: ignore = 4
    real(dp) :: deriv_eqd(geqdsk%nw), factor
    logical :: mask(geqdsk%nw)

    if (geqdsk%cocos%sgn_Btor /= geqdsk%cocos%sgn_F) then
       write (*, incons_fmt) 'FPOL', 'BCENTR'
       ! only modify array if signs are consistent
       if (geqdsk%cocos%sgn_F /= 0) then
          write (*, invert_fmt) 'FPOL'
          geqdsk%fpol = -geqdsk%fpol
          geqdsk%cocos%sgn_F = -geqdsk%cocos%sgn_F
       end if
    end if
    call resample1d(geqdsk%psi_eqd, geqdsk%pres, geqdsk%psi_eqd, deriv_eqd, 3, .true.)
    ! ignore a few possibly unreliable values on both ends of the interval
    mask = .false.
    mask(ignore+1:geqdsk%nw-ignore-1) = .true.
    factor = sum(deriv_eqd / geqdsk%pprime, mask = mask) / dble(count(mask))
    if (abs(factor) >= sqrt(2.0d0 * pi)) then
       write (*, unscaled_fmt) 'PPRIME'
       write (*, rescale_fmt) 'PPRIME'
       geqdsk%pprime = geqdsk%pprime * (2.0d0 * pi)
    end if
    if (factor < 0.0d0) then
       write (*, incons_fmt) 'PPRIME', 'PRES/(SIBRY-SIMAG)'
       write (*, invert_fmt) 'PPRIME'
       geqdsk%pprime = -geqdsk%pprime
    end if
    call resample1d(geqdsk%psi_eqd, geqdsk%fpol, geqdsk%psi_eqd, deriv_eqd, 3, .true.)
    deriv_eqd(:) = deriv_eqd * geqdsk%fpol
    factor = sum(deriv_eqd / geqdsk%ffprim, mask = mask) / dble(count(mask))
    if (abs(factor) >= sqrt(2.0d0 * pi)) then
       write (*, unscaled_fmt) 'FFPRIM'
       write (*, rescale_fmt) 'FFPRIM'
       geqdsk%ffprim = geqdsk%ffprim * (2.0d0 * pi)
    end if
    if (factor < 0.0d0) then
       write (*, incons_fmt) 'FFPRIM', 'FPOL/(SIBRY-SIMAG)'
       write (*, invert_fmt) 'FFPRIM'
       geqdsk%ffprim = -geqdsk%ffprim
    end if
  end subroutine geqdsk_check_consistency

  !> Estimates terms of Grad-Shafranov equation to determine cocos_t::exp_bpol.
  function geqdsk_GS_prefactor(geqdsk)
    type(geqdsk_t), intent(in) :: geqdsk
    real(dp) :: geqdsk_GS_prefactor
    real(dp) :: Delta_R, Delta_Z, psi, gs_rhs, gs_lhs, dp_dpsi, FdF_dpsi
    integer :: kw, kh

    ! find evaluation point to the upper right of the magnetic axis
    call binsearch(geqdsk%R_eqd, lbound(geqdsk%R_eqd, 1), geqdsk%rmaxis, kw)
    call binsearch(geqdsk%Z_eqd, lbound(geqdsk%Z_eqd, 1), geqdsk%zmaxis, kh)
    kw = kw + 2
    kh = kh + 2
    psi = geqdsk%psirz(kw, kh)
    ! evaluate flux functions on RHS of GS equation
    dp_dpsi = interp1d(geqdsk%psi_eqd, geqdsk%pprime, psi, 3)
    FdF_dpsi = interp1d(geqdsk%psi_eqd, geqdsk%ffprim, psi, 3)
    gs_rhs = -4.0d0 * pi * geqdsk%R_eqd(kw) ** 2 * dp_dpsi - FdF_dpsi
    ! approximate differential operator on LHS of GS equation
    Delta_R = geqdsk%rdim / dble(geqdsk%nw - 1)
    Delta_Z = geqdsk%zdim / dble(geqdsk%nh - 1)
    gs_lhs = (geqdsk%psirz(kw, kh + 1) - 2.0d0 * psi + geqdsk%psirz(kw, kh - 1)) / Delta_Z ** 2 &
         + (geqdsk%psirz(kw + 1, kh) - 2.0d0 * psi + geqdsk%psirz(kw - 1, kh)) / Delta_R ** 2 &
         - (geqdsk%psirz(kw + 1, kh) - geqdsk%psirz(kw - 1, kh)) / (2.0d0 * Delta_R * geqdsk%R_eqd(kw))
    geqdsk_GS_prefactor = gs_lhs / gs_rhs
  end function geqdsk_GS_prefactor

  function sign_array(array, name, most)
    real(dp), intent(in), dimension(:) :: array
    character(len = *), intent(in) :: name
    logical, intent(in), optional :: most
    logical :: strict
    integer :: sign_array, pos, neg, all

    strict = .true.
    if (present(most)) then
       strict = .not. most
    end if
    pos = count(array > 0.0d0)
    neg = count(array < 0.0d0)
    all = size(array)
    if ((strict .and. all == pos) .or. (.not. strict .and. pos > neg)) then
       sign_array = +1
    elseif ((strict .and. all == neg) .or. (.not. strict .and. neg > pos)) then
       sign_array = -1
    else
       sign_array = 0
       write (*, '("Sign of ", a, " is inconsistent.")') trim(name)
    end if
  end function sign_array

  subroutine geqdsk_classify(geqdsk)
    type(geqdsk_t), intent(inout) :: geqdsk
    real(dp) :: GS_prefactor

    geqdsk%cocos = cocos_t(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    geqdsk%cocos%sgn_cyl = +1  ! specified by GEQDSK format and assumed in field_divB0.f90
    geqdsk%cocos%sgn_dpsi = sign_array([geqdsk%sibry - geqdsk%simag], 'SIBRY-SIMAG')
    geqdsk%cocos%sgn_Btor = sign_array([geqdsk%bcentr], 'BCENTR')
    geqdsk%cocos%sgn_Itor = sign_array([geqdsk%current], 'CURRENT')
    geqdsk%cocos%sgn_F = sign_array(geqdsk%fpol, 'FPOL')
    ! q is often unreliable in gfiles and its sign is only necessary for exact
    ! classification, not for calculations where the poloidal direction is fixed anyway
    geqdsk%cocos%sgn_q = sign_array(geqdsk%qpsi, 'QPSI', .true.)
    geqdsk%cocos%sgn_Bpol = geqdsk%cocos%sgn_dpsi * geqdsk%cocos%sgn_Itor
    geqdsk%cocos%sgn_pol = geqdsk%cocos%sgn_q * geqdsk%cocos%sgn_Itor * geqdsk%cocos%sgn_Btor
    if (geqdsk%cocos%sgn_Bpol == +1) then
       if (geqdsk%cocos%sgn_pol == +1) then
          geqdsk%cocos%index = 1
       elseif (geqdsk%cocos%sgn_pol == -1) then
          geqdsk%cocos%index = 5
       end if
    elseif (geqdsk%cocos%sgn_Bpol == -1) then
       if (geqdsk%cocos%sgn_pol == +1) then
          geqdsk%cocos%index = 7
       elseif (geqdsk%cocos%sgn_pol == -1) then
          geqdsk%cocos%index = 3
       end if
    end if
    if (geqdsk%cocos%index == 0) then
       write (*, '("COCOS index could not be determined for ", a)') trim(geqdsk%fname)
    end if
    call geqdsk_check_consistency(geqdsk)
    GS_prefactor = geqdsk_GS_prefactor(geqdsk)
    write (*, '("Grad-Shafranov equation prefactor estimate: ", es24.16e3)') GS_prefactor
    if (abs(GS_prefactor) >= 2.0d0 * pi) then
       geqdsk%cocos%index = geqdsk%cocos%index + 10
       geqdsk%cocos%exp_Bpol = 1
    else
       ! specified by GEQDSK format and assumed in field_divB0.f90
       geqdsk%cocos%exp_Bpol = 0
    end if
  end subroutine geqdsk_classify

  subroutine geqdsk_standardise(geqdsk)
    type(geqdsk_t), intent(inout) :: geqdsk

    if (geqdsk%cocos%sgn_Bpol == +1) then
       write (*, invert_fmt) 'SIMAG, SIBRY, PSIRZ, PPRIME, FFPRIM'
       geqdsk%simag = -geqdsk%simag
       geqdsk%sibry = -geqdsk%sibry
       geqdsk%psirz(:, :) = -geqdsk%psirz
       geqdsk%psi_eqd(:) = -geqdsk%psi_eqd
       geqdsk%cocos%sgn_dpsi = -geqdsk%cocos%sgn_dpsi
       geqdsk%pprime(:) = -geqdsk%pprime
       geqdsk%ffprim(:) = -geqdsk%ffprim
       geqdsk%fprime(:) = -geqdsk%fprime
       geqdsk%cocos%sgn_Bpol = -geqdsk%cocos%sgn_Bpol
    end if
    if (geqdsk%cocos%sgn_pol == +1) then
       write (*, invert_fmt) 'QPSI'
       geqdsk%qpsi(:) = -geqdsk%qpsi
       geqdsk%cocos%sgn_q = -geqdsk%cocos%sgn_q
       geqdsk%cocos%sgn_pol = -geqdsk%cocos%sgn_pol
    end if
    if (geqdsk%cocos%exp_Bpol == 1) then
       write (*, unscaled_fmt) 'SIMAG, SIBRY, PSIRZ, PPRIME, FFPRIM'
       write (*, rescale_fmt) 'SIMAG, SIBRY, PSIRZ, PPRIME, FFPRIM'
       geqdsk%simag = geqdsk%simag / (2.0d0 * pi)
       geqdsk%sibry = geqdsk%sibry / (2.0d0 * pi)
       geqdsk%psirz(:, :) = geqdsk%psirz / (2.0d0 * pi)
       geqdsk%psi_eqd(:) = geqdsk%psi_eqd / (2.0d0 * pi)
       geqdsk%pprime(:) = geqdsk%pprime * (2.0d0 * pi)
       geqdsk%ffprim(:) = geqdsk%ffprim * (2.0d0 * pi)
       geqdsk%fprime(:) = geqdsk%fprime * (2.0d0 * pi)
       geqdsk%cocos%exp_Bpol = 0
    end if
    geqdsk%cocos%index = 3
  end subroutine geqdsk_standardise

  !> Binary search to find index \p k of ordered array \p x so that \p xi lies between
  !> x(k-1) and x(k). The lower bound of \p x is given by \p lb.
  !>
  !> When \p x is strictly monotonically increasing, \f$ x_{k-1} < \xi < x_{k} \f$.
  !> When \p x is strictly monotonically decreasing, \f$ x_{k-1} > \xi > x_{k} \f$.
  !> It is not checked whether \p x is ordered; only the first and last value are used
  !> to determine wether the array values are increasing or decrasing.
  subroutine binsearch(x, lb, xi, k)
    integer, intent(in) :: lb
    real(dp), intent(in) :: x(lb:)
    real(dp), intent(in) :: xi
    integer, intent(out) :: k
    integer :: k_min, k_max, iter

    k_min = lbound(x, 1)
    k_max = ubound(x, 1)
    if (x(k_min) < x(k_max)) then
       do iter = 1, size(x) - 1
          k = (k_max - k_min) / 2 + k_min
          if (x(k) > xi) then
             k_max = k
          else
             k_min = k
          end if
          if (k_max == k_min + 1) exit
       end do
    else
       do iter = 1, size(x) - 1
          k = (k_max - k_min) / 2 + k_min
          if (x(k) < xi) then
             k_max = k
          else
             k_min = k
          end if
          if (k_max == k_min + 1) exit
       end do
    end if
    k = k_max
  end subroutine binsearch

  !> 1D piecewise Lagrange polynomial interpolator
  function interp1d(sample_x, sample_y, resampled_x, n_lag, deriv) result(resampled_y)
    real(dp), intent(in) :: sample_x(:)
    real(dp), intent(in) :: sample_y(:)
    real(dp), intent(in) :: resampled_x
    real(dp) :: resampled_y
    integer, intent(in) :: n_lag
    logical, intent(in), optional :: deriv
    real(dp) :: lag_coeff(0:1, n_lag + 1)
    integer :: kder, k, k_lo, k_hi

    if (size(sample_x) /= size(sample_y)) then
       write (*, '("Argument size mismatch in interp1d: size(sample_x) = ", i0, ' // &
            '", size(sample_y) = ", i0)') size(sample_x), size(sample_y)
       error stop
    end if
    if (n_lag >= size(sample_x)) then
       write (*, '("interp1d: Lagrange polynomial order ", i0, ' // &
            '" must be lower than number of sample points ", i0)') &
            n_lag, size(sample_x)
       error stop
    end if
    kder = 0
    if (present(deriv)) then
       if (deriv) then
          kder = 1
       end if
    end if
    call binsearch(sample_x, lbound(sample_x, 1), resampled_x, k)
    k_lo = k - (n_lag + 1) / 2
    k_hi = k_lo + n_lag
    ! ensure that polynomial sample points remain within the bounds of sample_x
    if (k_lo < lbound(sample_x, 1)) then
       k_lo = lbound(sample_x, 1)
       k_hi = k_lo + n_lag
    elseif (k_hi > ubound(sample_x, 1)) then
       k_hi = ubound(sample_x, 1)
       k_lo = k_hi - n_lag
    end if
    call plag_coeff(n_lag + 1, 1, resampled_x, sample_x(k_lo:k_hi), lag_coeff)
    resampled_y = sum(sample_y(k_lo:k_hi) * lag_coeff(kder, :))
  end function interp1d

  !> 1D piecewise Lagrange polynomial resampler
  subroutine resample1d(sample_x, sample_y, resampled_x, resampled_y, n_lag, deriv)
    real(dp), intent(in) :: sample_x(:)
    real(dp), intent(in) :: sample_y(:)
    real(dp), intent(in) :: resampled_x(:)
    real(dp), intent(out) :: resampled_y(:)
    integer, intent(in) :: n_lag
    logical, intent(in), optional :: deriv
    real(dp) :: lag_coeff(0:1, n_lag + 1)
    integer :: kder, k, k_lo, k_hi, k_re

    if (size(sample_x) /= size(sample_y)) then
       write (*, '("Argument size mismatch in resample1d: size(sample_x) = ", i0, ' // &
            '", size(sample_y) = ", i0)') size(sample_x), size(sample_y)
       error stop
    end if
    if (size(resampled_x) /= size(resampled_y)) then
       write (*, '("Argument size mismatch in resample1d: size(resampled_x) = ", i0, ' // &
            '", size(resampled_y) = ", i0)') size(resampled_x), size(resampled_y)
       error stop
    end if
    if (n_lag >= size(sample_x)) then
       write (*, '("resample1d: Lagrange polynomial order ", i0, ' // &
            '" must be lower than number of sample points ", i0)') &
            n_lag, size(sample_x)
       error stop
    end if
    kder = 0
    if (present(deriv)) then
       if (deriv) then
          kder = 1
       end if
    end if
    do k_re = lbound(resampled_x, 1), ubound(resampled_x, 1)
       call binsearch(sample_x, lbound(sample_x, 1), resampled_x(k_re), k)
       k_lo = k - (n_lag + 1) / 2
       k_hi = k_lo + n_lag
       ! ensure that polynomial sample points remain within the bounds of sample_x
       if (k_lo < lbound(sample_x, 1)) then
          k_lo = lbound(sample_x, 1)
          k_hi = k_lo + n_lag
       elseif (k_hi > ubound(sample_x, 1)) then
          k_hi = ubound(sample_x, 1)
          k_lo = k_hi - n_lag
       end if
       call plag_coeff(n_lag + 1, 1, resampled_x(k_re), sample_x(k_lo:k_hi), lag_coeff)
       resampled_y(k_re) = sum(sample_y(k_lo:k_hi) * lag_coeff(kder, :))
    end do
  end subroutine resample1d

end module geqdsk_tools
