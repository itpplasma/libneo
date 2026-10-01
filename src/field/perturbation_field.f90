module neo_perturbation_field
    ! Toroidal-harmonic perturbation field on an (R, Z) grid.
    !
    ! Coordinates: right-handed cylinder (R, phi, Z), phi the geometric toroidal
    ! angle, counterclockwise seen from above (+Z). Vector components are physical
    ! (orthonormal) cylindrical components. For each toroidal mode number n the
    ! container holds complex amplitudes dA_n = (dA_R, dA_phi, dA_Z)(R, Z) and,
    ! optionally, a scalar potential dPhi_n(R, Z). The real field is
    !   dA(R, phi, Z) = Re( sum_n dA_n(R, Z) exp(i n phi) ).
    ! dB = curl dA is evaluated from analytic derivatives of the spline:
    !   dB_R   = (i n/R) dA_Z - d_Z dA_phi
    !   dB_phi = d_Z dA_R - d_R dA_Z
    !   dB_Z   = dA_phi/R + d_R dA_phi - (i n/R) dA_R
    ! Units are carried, not converted: UNITS_SI (m, T m, T) or UNITS_GAUSSIAN
    ! (cm, G cm, G). Evaluation outside the grid returns ierr /= 0 and NaN.
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
    use neo_field_base, only: field_t
    use interpolate, only: BatchSplineData2D, construct_batch_splines_2d, &
        destroy_batch_splines_2d, evaluate_batch_splines_2d_der
    implicit none
    private

    integer, parameter :: dp = kind(1.0d0)

    character(len=*), parameter, public :: UNITS_SI = 'SI'
    character(len=*), parameter, public :: UNITS_GAUSSIAN = 'Gaussian'

    integer, parameter, public :: PERTFIELD_OK = 0
    integer, parameter, public :: PERTFIELD_OUTSIDE_GRID = 1
    integer, parameter, public :: PERTFIELD_BAD_INPUT = 2
    integer, parameter, public :: PERTFIELD_NO_POTENTIAL = 3

    integer, parameter :: NQ_A = 6, NQ_PHI = 2

    type, extends(field_t), public :: perturbation_field_t
        integer :: n_modes = 0
        integer, allocatable :: ntor(:)
        real(dp), allocatable :: R(:), Z(:)
        complex(dp), allocatable :: dA(:, :, :, :)
        complex(dp), allocatable :: dPhi(:, :, :)
        logical :: has_potential = .false.
        character(len=16) :: units = ''
        integer :: spline_order = 5
        type(BatchSplineData2D), private :: spl
    contains
        procedure :: init
        procedure :: init_analytic
        procedure :: eval_modes
        procedure :: eval
        procedure :: dbmod_modes
        procedure :: dbmod
        procedure :: inside
        procedure :: compute_afield
        procedure :: compute_bfield
        procedure :: compute_abfield
    end type perturbation_field_t

    abstract interface
        subroutine mode_amplitude_i(n, R, Z, dA, dPhi)
            import :: dp
            integer, intent(in) :: n
            real(dp), intent(in) :: R, Z
            complex(dp), intent(out) :: dA(3), dPhi
        end subroutine mode_amplitude_i
    end interface
    public :: mode_amplitude_i

contains

    ! R(nR), Z(nZ) equidistant and ascending, R(1) > 0. dA(3, nR, nZ, n_modes),
    ! component order (R, phi, Z). dPhi(nR, nZ, n_modes) is optional.
    subroutine init(self, R, Z, ntor, dA, units, ierr, dPhi, order)
        class(perturbation_field_t), intent(inout) :: self
        real(dp), intent(in) :: R(:), Z(:)
        integer, intent(in) :: ntor(:)
        complex(dp), intent(in) :: dA(:, :, :, :)
        character(len=*), intent(in) :: units
        integer, intent(out) :: ierr
        complex(dp), intent(in), optional :: dPhi(:, :, :)
        integer, intent(in), optional :: order

        call clear(self)
        ierr = validate(R, Z, ntor, dA, units, order)
        if (ierr /= PERTFIELD_OK) return
        if (present(dPhi)) then
            if (any(shape(dPhi) /= [size(R), size(Z), size(ntor)])) then
                ierr = PERTFIELD_BAD_INPUT
                return
            end if
            allocate (self%dPhi, source=dPhi)
            self%has_potential = .true.
        end if
        if (present(order)) self%spline_order = order
        self%n_modes = size(ntor)
        allocate (self%ntor, source=ntor)
        allocate (self%R, source=R)
        allocate (self%Z, source=Z)
        allocate (self%dA, source=dA)
        self%units = units
        call build_splines(self)
    end subroutine init

    ! Samples amplitude(n, R, Z) on an equidistant nR x nZ grid.
    subroutine init_analytic(self, R_min, R_max, nR, Z_min, Z_max, nZ, ntor, &
            amplitude, units, ierr, with_potential, order)
        class(perturbation_field_t), intent(inout) :: self
        real(dp), intent(in) :: R_min, R_max, Z_min, Z_max
        integer, intent(in) :: nR, nZ, ntor(:)
        procedure(mode_amplitude_i) :: amplitude
        character(len=*), intent(in) :: units
        integer, intent(out) :: ierr
        logical, intent(in), optional :: with_potential
        integer, intent(in), optional :: order

        real(dp), allocatable :: R(:), Z(:)
        complex(dp), allocatable :: dA(:, :, :, :), dPhi(:, :, :)
        logical :: want_phi
        integer :: i, j, k

        if (nR < 2 .or. nZ < 2) then
            call clear(self)
            ierr = PERTFIELD_BAD_INPUT
            return
        end if
        want_phi = .false.
        if (present(with_potential)) want_phi = with_potential
        allocate (R(nR), Z(nZ), dA(3, nR, nZ, size(ntor)), dPhi(nR, nZ, size(ntor)))
        do i = 1, nR
            R(i) = R_min + (R_max - R_min)*real(i - 1, dp)/real(nR - 1, dp)
        end do
        do j = 1, nZ
            Z(j) = Z_min + (Z_max - Z_min)*real(j - 1, dp)/real(nZ - 1, dp)
        end do
        do k = 1, size(ntor)
            do j = 1, nZ
                do i = 1, nR
                    call amplitude(ntor(k), R(i), Z(j), dA(:, i, j, k), dPhi(i, j, k))
                end do
            end do
        end do
        if (want_phi) then
            call self%init(R, Z, ntor, dA, units, ierr, dPhi=dPhi, order=order)
        else
            call self%init(R, Z, ntor, dA, units, ierr, order=order)
        end if
    end subroutine init_analytic

    integer function validate(R, Z, ntor, dA, units, order) result(ierr)
        real(dp), intent(in) :: R(:), Z(:)
        integer, intent(in) :: ntor(:)
        complex(dp), intent(in) :: dA(:, :, :, :)
        character(len=*), intent(in) :: units
        integer, intent(in), optional :: order

        integer :: k

        ierr = PERTFIELD_BAD_INPUT
        if (units /= UNITS_SI .and. units /= UNITS_GAUSSIAN) return
        if (size(ntor) < 1) return
        if (any(shape(dA) /= [3, size(R), size(Z), size(ntor)])) return
        if (present(order)) then
            if (order < 3 .or. order > 5) return
        end if
        if (size(R) < 6 .or. size(Z) < 6) return
        if (.not. equidistant(R) .or. .not. equidistant(Z)) return
        if (.not. (R(1) > 0.0_dp)) return
        do k = 2, size(ntor)
            if (any(ntor(:k - 1) == ntor(k))) return
        end do
        ierr = PERTFIELD_OK
    end function validate

    logical function equidistant(x)
        real(dp), intent(in) :: x(:)

        real(dp) :: h

        h = (x(size(x)) - x(1))/real(size(x) - 1, dp)
        equidistant = h > 0.0_dp
        if (.not. equidistant) return
        equidistant = all(abs(x(2:) - x(:size(x) - 1) - h) <= 1.0e-10_dp*h)
    end function equidistant

    subroutine clear(self)
        class(perturbation_field_t), intent(inout) :: self

        call destroy_batch_splines_2d(self%spl)
        if (allocated(self%ntor)) deallocate (self%ntor)
        if (allocated(self%R)) deallocate (self%R)
        if (allocated(self%Z)) deallocate (self%Z)
        if (allocated(self%dA)) deallocate (self%dA)
        if (allocated(self%dPhi)) deallocate (self%dPhi)
        self%n_modes = 0
        self%has_potential = .false.
        self%units = ''
        self%spline_order = 5
    end subroutine clear

    pure integer function quantities_per_mode(self) result(nq)
        class(perturbation_field_t), intent(in) :: self

        nq = NQ_A
        if (self%has_potential) nq = NQ_A + NQ_PHI
    end function quantities_per_mode

    subroutine build_splines(self)
        class(perturbation_field_t), intent(inout) :: self

        real(dp), allocatable :: y(:, :, :)
        integer :: k, c, base, nq, nR, nZ

        nR = size(self%R)
        nZ = size(self%Z)
        nq = quantities_per_mode(self)
        allocate (y(nR, nZ, nq*self%n_modes))
        do k = 1, self%n_modes
            base = (k - 1)*nq
            do c = 1, 3
                y(:, :, base + 2*c - 1) = self%dA(c, :, :, k)%re
                y(:, :, base + 2*c) = self%dA(c, :, :, k)%im
            end do
            if (self%has_potential) then
                y(:, :, base + NQ_A + 1) = self%dPhi(:, :, k)%re
                y(:, :, base + NQ_A + 2) = self%dPhi(:, :, k)%im
            end if
        end do
        call construct_batch_splines_2d([self%R(1), self%Z(1)], &
            [self%R(nR), self%Z(nZ)], y, [self%spline_order, self%spline_order], &
            [.false., .false.], self%spl)
    end subroutine build_splines

    logical function inside(self, R, Z)
        class(perturbation_field_t), intent(in) :: self
        real(dp), intent(in) :: R, Z

        inside = .false.
        if (self%n_modes < 1) return
        if (ieee_is_nan(R) .or. ieee_is_nan(Z)) return
        if (.not. (R >= self%R(1) .and. R <= self%R(size(self%R)))) return
        if (.not. (Z >= self%Z(1) .and. Z <= self%Z(size(self%Z)))) return
        inside = .true.
    end function inside

    real(dp) function nan()
        nan = ieee_value(1.0_dp, ieee_quiet_nan)
    end function nan

    ! Complex mode amplitudes dA_n, dB_n = curl(dA_n exp(i n phi)) exp(-i n phi)
    ! and optionally dPhi_n at (R, Z), for all modes in the order of self%ntor.
    subroutine eval_modes(self, R, Z, dA, dB, ierr, dPhi)
        class(perturbation_field_t), intent(in) :: self
        real(dp), intent(in) :: R, Z
        complex(dp), intent(out) :: dA(:, :), dB(:, :)
        integer, intent(out) :: ierr
        complex(dp), intent(out), optional :: dPhi(:)

        real(dp) :: y(quantities_per_mode(self)*self%n_modes)
        real(dp) :: dy(2, quantities_per_mode(self)*self%n_modes)
        complex(dp) :: a(3), a_R(3), a_Z(3), in_R
        integer :: k, c, base, nq

        dA = cmplx(nan(), nan(), dp)
        dB = dA
        if (present(dPhi)) dPhi = cmplx(nan(), nan(), dp)
        ierr = PERTFIELD_OUTSIDE_GRID
        if (.not. self%inside(R, Z)) return
        ierr = PERTFIELD_BAD_INPUT
        if (size(dA, 2) < self%n_modes .or. size(dB, 2) < self%n_modes) return
        ierr = PERTFIELD_NO_POTENTIAL
        if (present(dPhi) .and. .not. self%has_potential) return
        ierr = PERTFIELD_OK

        nq = quantities_per_mode(self)
        call evaluate_batch_splines_2d_der(self%spl, [R, Z], y, dy)
        do k = 1, self%n_modes
            base = (k - 1)*nq
            do c = 1, 3
                a(c) = cmplx(y(base + 2*c - 1), y(base + 2*c), dp)
                a_R(c) = cmplx(dy(1, base + 2*c - 1), dy(1, base + 2*c), dp)
                a_Z(c) = cmplx(dy(2, base + 2*c - 1), dy(2, base + 2*c), dp)
            end do
            in_R = cmplx(0.0_dp, real(self%ntor(k), dp)/R, dp)
            dA(:, k) = a
            dB(1, k) = in_R*a(3) - a_Z(2)
            dB(2, k) = a_Z(1) - a_R(3)
            dB(3, k) = a(2)/R + a_R(2) - in_R*a(1)
            if (present(dPhi)) then
                dPhi(k) = cmplx(y(base + NQ_A + 1), y(base + NQ_A + 2), dp)
            end if
        end do
    end subroutine eval_modes

    ! Real fields dA, dB at x = (R, phi, Z), physical cylindrical components.
    subroutine eval(self, x, dA, dB, ierr)
        class(perturbation_field_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: dA(3), dB(3)
        integer, intent(out) :: ierr

        complex(dp) :: cA(3, self%n_modes), cB(3, self%n_modes), phase
        integer :: k

        dA = nan()
        dB = nan()
        call self%eval_modes(x(1), x(3), cA, cB, ierr)
        if (ierr /= PERTFIELD_OK) return
        dA = 0.0_dp
        dB = 0.0_dp
        do k = 1, self%n_modes
            phase = exp(cmplx(0.0_dp, real(self%ntor(k), dp)*x(2), dp))
            dA = dA + real(cA(:, k)*phase, dp)
            dB = dB + real(cB(:, k)*phase, dp)
        end do
    end subroutine eval

    ! Eulerian d|B|_n = b0 . dB_n per mode, with b0 = B0/|B0| for the given
    ! axisymmetric equilibrium field B0 = (B0_R, B0_phi, B0_Z) at (R, Z).
    subroutine dbmod_modes(self, R, Z, B0, dbmod, ierr)
        class(perturbation_field_t), intent(in) :: self
        real(dp), intent(in) :: R, Z, B0(3)
        complex(dp), intent(out) :: dbmod(:)
        integer, intent(out) :: ierr

        complex(dp) :: cA(3, self%n_modes), cB(3, self%n_modes)
        real(dp) :: b0_unit(3), b0_norm
        integer :: k

        dbmod = cmplx(nan(), nan(), dp)
        call self%eval_modes(R, Z, cA, cB, ierr)
        if (ierr /= PERTFIELD_OK) return
        b0_norm = norm2(B0)
        ierr = PERTFIELD_BAD_INPUT
        if (.not. (b0_norm > 0.0_dp)) return
        if (size(dbmod) < self%n_modes) return
        ierr = PERTFIELD_OK
        b0_unit = B0/b0_norm
        do k = 1, self%n_modes
            dbmod(k) = sum(b0_unit*cB(:, k))
        end do
    end subroutine dbmod_modes

    ! Real Eulerian d|B| = b0 . dB at x = (R, phi, Z).
    subroutine dbmod(self, x, B0, dbm, ierr)
        class(perturbation_field_t), intent(in) :: self
        real(dp), intent(in) :: x(3), B0(3)
        real(dp), intent(out) :: dbm
        integer, intent(out) :: ierr

        complex(dp) :: modes(self%n_modes)
        integer :: k

        dbm = nan()
        call self%dbmod_modes(x(1), x(3), B0, modes, ierr)
        if (ierr /= PERTFIELD_OK) return
        dbm = 0.0_dp
        do k = 1, self%n_modes
            dbm = dbm + real(modes(k)*exp(cmplx(0.0_dp, self%ntor(k)*x(2), dp)), dp)
        end do
    end subroutine dbmod

    ! field_t interface. It has no status argument, so evaluation outside the grid
    ! stops instead of returning a silent value.
    subroutine compute_abfield(self, x, A, B)
        class(perturbation_field_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: A(3), B(3)

        integer :: ierr

        call self%eval(x, A, B, ierr)
        if (ierr /= PERTFIELD_OK) error stop 'perturbation_field_t: point outside grid'
    end subroutine compute_abfield

    subroutine compute_afield(self, x, A)
        class(perturbation_field_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: A(3)

        real(dp) :: B(3)

        call self%compute_abfield(x, A, B)
    end subroutine compute_afield

    subroutine compute_bfield(self, x, B)
        class(perturbation_field_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: B(3)

        real(dp) :: A(3)

        call self%compute_abfield(x, A, B)
    end subroutine compute_bfield

end module neo_perturbation_field
