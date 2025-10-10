module analytical_tokamak_field
    !> Analytical tokamak field evaluation using Cerfon-Freidberg equilibria
    !>
    !> Implementation following Cerfon & Freidberg
    !> (Physics of Plasmas 17, 032502, 2010) and Verena Eslbauer
    !> (Two analytical solutions to the Grad-Shafranov equation,
    !> Bachelor thesis, TU Graz, 2017).
    !> Provides magnetic field components for circular tokamak equilibria
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use analytical_gs_solver
    use neo_field_base, only: field_t
    implicit none

    private
    public :: analytical_circular_eq_t
    public :: compute_analytical_field_cylindrical

    !> Analytical tokamak equilibrium (supports arbitrary shaping)
    type, extends(field_t) :: analytical_circular_eq_t
        real(dp) :: R0         !< Major radius [m]
        real(dp) :: epsilon    !< Inverse aspect ratio
        real(dp) :: kappa      !< Elongation
        real(dp) :: delta      !< Triangularity
        real(dp) :: A_param    !< Shafranov parameter
        real(dp) :: B0         !< Toroidal field on axis [T]
        real(dp) :: psimult    !< Flux scaling factor
        real(dp) :: coeffs(7) !< Boundary condition coefficients
        integer :: Nripple = 0 !< Number of TF coils (0 = no ripple)
        real(dp) :: a0         !< Minor radius for ripple [m]
        real(dp) :: alpha0     !< Ripple radial dependency exponent
        real(dp) :: delta0     !< Ripple amplitude
        real(dp) :: z0         !< Midplane z coordinate for ripple [m]
        logical :: initialized = .false.
    contains
        procedure :: init => init_circular_equilibrium
        procedure :: eval_psi => evaluate_psi
        procedure :: eval_psi_derivatives
        procedure :: eval_bfield => evaluate_bfield
        procedure :: eval_bfield_ripple
        procedure :: cleanup => destroy_equilibrium
        procedure :: compute_afield
        procedure :: compute_bfield
        procedure :: compute_abfield
    end type analytical_circular_eq_t

contains

    !> Initialize tokamak equilibrium
    !>
    !> Parameters:
    !>   R0_in      - Major radius [m]
    !>   epsilon_in - Inverse aspect ratio (a/R0)
    !>   kappa_in   - Elongation (optional, default=1 for circular)
    !>   delta_in   - Triangularity (optional, default=0 for circular)
    !>   A_in       - Shafranov parameter (controls q-profile)
    !>   B0_in      - Toroidal field on axis [T]
    !>   psimult_in - Flux scaling factor (optional, default=200)
    !>   Nripple_in - Number of TF coils for ripple (optional, 0=no ripple)
    !>   a0_in      - Minor radius for ripple calculation (optional)
    !>   alpha0_in  - Ripple radial dependency exponent (optional, default=2.0)
    !>   delta0_in  - Ripple amplitude (optional, default=0.0)
    !>   z0_in      - Midplane z coordinate for ripple (optional, default=0.0)
    subroutine init_circular_equilibrium(self, R0_in, epsilon_in, kappa_in, delta_in, A_in, B0_in, &
                                         psimult_in, Nripple_in, a0_in, alpha0_in, delta0_in, z0_in)
        class(analytical_circular_eq_t), intent(inout) :: self
        real(dp), intent(in) :: R0_in
        real(dp), intent(in) :: epsilon_in
        real(dp), intent(in), optional :: kappa_in
        real(dp), intent(in), optional :: delta_in
        real(dp), intent(in) :: A_in
        real(dp), intent(in) :: B0_in
        real(dp), intent(in), optional :: psimult_in
        integer, intent(in), optional :: Nripple_in
        real(dp), intent(in), optional :: a0_in
        real(dp), intent(in), optional :: alpha0_in
        real(dp), intent(in), optional :: delta0_in
        real(dp), intent(in), optional :: z0_in

        self%R0 = R0_in
        self%epsilon = epsilon_in
        self%A_param = A_in
        self%B0 = B0_in

        if (present(kappa_in)) then
            self%kappa = kappa_in
        else
            self%kappa = 1.0_dp
        end if

        if (present(delta_in)) then
            self%delta = delta_in
        else
            self%delta = 0.0_dp
        end if

        if (present(psimult_in)) then
            self%psimult = psimult_in
        else
            self%psimult = 200.0_dp
        end if

        if (present(Nripple_in)) then
            self%Nripple = Nripple_in
        else
            self%Nripple = 0
        end if

        if (present(a0_in)) then
            self%a0 = a0_in
        else
            self%a0 = R0_in * epsilon_in
        end if

        if (present(alpha0_in)) then
            self%alpha0 = alpha0_in
        else
            self%alpha0 = 2.0_dp
        end if

        if (present(delta0_in)) then
            self%delta0 = delta0_in
        else
            self%delta0 = 0.0_dp
        end if

        if (present(z0_in)) then
            self%z0 = z0_in
        else
            self%z0 = 0.0_dp
        end if

        call solve_coefficients(self%epsilon, self%kappa, self%delta, self%A_param, self%coeffs)

        self%initialized = .true.
    end subroutine init_circular_equilibrium

    !> Evaluate poloidal flux psi(R, Z)
    !>
    !> Parameters:
    !>   R - Major radius coordinate [m]
    !>   Z - Vertical coordinate [m]
    !>
    !> Returns:
    !>   psi - Poloidal flux [Wb] or [T*m^2]
    function evaluate_psi(self, R, Z) result(psi)
        class(analytical_circular_eq_t), intent(in) :: self
        real(dp), intent(in) :: R, Z
        real(dp) :: psi

        real(dp) :: x, y
        real(dp) :: psi_homogeneous, psi_particular
        integer :: i

        if (.not. self%initialized) then
            error stop "analytical_circular_eq_t not initialized"
        end if

        ! Convert to normalized coordinates
        x = R / self%R0
        y = Z / self%R0

        ! Evaluate homogeneous solution (Cerfon-Freidberg basis functions)
        psi_homogeneous = self%coeffs(1) * psi_1(x, y) &
                        + self%coeffs(2) * psi_2(x, y) &
                        + self%coeffs(3) * psi_3(x, y) &
                        + self%coeffs(4) * psi_4(x, y) &
                        + self%coeffs(5) * psi_5(x, y) &
                        + self%coeffs(6) * psi_6(x, y) &
                        + self%coeffs(7) * psi_7(x, y)

        ! Evaluate particular solution
        psi_particular = psi_p(x, y, self%A_param)

        ! Total normalized flux
        psi = psi_homogeneous + psi_particular

        ! Scale to physical units
        psi = psi * self%psimult
    end function evaluate_psi

    !> Evaluate psi derivatives dpsi/dR and dpsi/dZ
    subroutine eval_psi_derivatives(self, R, Z, dpsi_dR, dpsi_dZ)
        class(analytical_circular_eq_t), intent(in) :: self
        real(dp), intent(in) :: R, Z
        real(dp), intent(out) :: dpsi_dR, dpsi_dZ

        real(dp) :: x, y
        real(dp) :: dpsi_dx, dpsi_dy
        real(dp) :: dpsi_homogeneous_dx, dpsi_homogeneous_dy
        real(dp) :: dpsi_particular_dx, dpsi_particular_dy
        real(dp) :: dx_dR, dy_dZ
        integer :: i

        if (.not. self%initialized) then
            error stop "analytical_circular_eq_t not initialized"
        end if

        ! Convert to normalized coordinates
        x = R / self%R0
        y = Z / self%R0

        ! Chain rule: dx/dR = 1/R0, dy/dZ = 1/R0
        dx_dR = 1.0_dp / self%R0
        dy_dZ = 1.0_dp / self%R0

        ! Evaluate homogeneous derivatives (sum over all 7 basis functions)
        dpsi_homogeneous_dx = 0.0_dp
        dpsi_homogeneous_dy = 0.0_dp

        call dpsi_1_derivatives(x, y, dpsi_dx, dpsi_dy)
        dpsi_homogeneous_dx = dpsi_homogeneous_dx + self%coeffs(1) * dpsi_dx
        dpsi_homogeneous_dy = dpsi_homogeneous_dy + self%coeffs(1) * dpsi_dy

        call dpsi_2_derivatives(x, y, dpsi_dx, dpsi_dy)
        dpsi_homogeneous_dx = dpsi_homogeneous_dx + self%coeffs(2) * dpsi_dx
        dpsi_homogeneous_dy = dpsi_homogeneous_dy + self%coeffs(2) * dpsi_dy

        call dpsi_3_derivatives(x, y, dpsi_dx, dpsi_dy)
        dpsi_homogeneous_dx = dpsi_homogeneous_dx + self%coeffs(3) * dpsi_dx
        dpsi_homogeneous_dy = dpsi_homogeneous_dy + self%coeffs(3) * dpsi_dy

        call dpsi_4_derivatives(x, y, dpsi_dx, dpsi_dy)
        dpsi_homogeneous_dx = dpsi_homogeneous_dx + self%coeffs(4) * dpsi_dx
        dpsi_homogeneous_dy = dpsi_homogeneous_dy + self%coeffs(4) * dpsi_dy

        call dpsi_5_derivatives(x, y, dpsi_dx, dpsi_dy)
        dpsi_homogeneous_dx = dpsi_homogeneous_dx + self%coeffs(5) * dpsi_dx
        dpsi_homogeneous_dy = dpsi_homogeneous_dy + self%coeffs(5) * dpsi_dy

        call dpsi_6_derivatives(x, y, dpsi_dx, dpsi_dy)
        dpsi_homogeneous_dx = dpsi_homogeneous_dx + self%coeffs(6) * dpsi_dx
        dpsi_homogeneous_dy = dpsi_homogeneous_dy + self%coeffs(6) * dpsi_dy

        call dpsi_7_derivatives(x, y, dpsi_dx, dpsi_dy)
        dpsi_homogeneous_dx = dpsi_homogeneous_dx + self%coeffs(7) * dpsi_dx
        dpsi_homogeneous_dy = dpsi_homogeneous_dy + self%coeffs(7) * dpsi_dy

        ! Evaluate particular solution derivatives
        dpsi_particular_dx = dpsi_p_dx(x, y, self%A_param)
        dpsi_particular_dy = dpsi_p_dy(x, y, self%A_param)

        ! Total derivatives in normalized coords
        dpsi_dx = dpsi_homogeneous_dx + dpsi_particular_dx
        dpsi_dy = dpsi_homogeneous_dy + dpsi_particular_dy

        ! Transform to physical coordinates
        dpsi_dR = dpsi_dx * dx_dR * self%psimult
        dpsi_dZ = dpsi_dy * dy_dZ * self%psimult
    end subroutine eval_psi_derivatives

    !> Evaluate magnetic field components
    !>
    !> Parameters:
    !>   R - Major radius coordinate [m]
    !>   Z - Vertical coordinate [m]
    !>
    !> Returns:
    !>   B_R   - Radial field component [T]
    !>   B_Z   - Vertical field component [T]
    !>   B_phi - Toroidal field component [T]
    !>   B_mod - Total field magnitude [T]
    subroutine evaluate_bfield(self, R, Z, B_R, B_Z, B_phi, B_mod)
        class(analytical_circular_eq_t), intent(in) :: self
        real(dp), intent(in) :: R, Z
        real(dp), intent(out) :: B_R, B_Z, B_phi, B_mod

        real(dp) :: dpsi_dR, dpsi_dZ
        real(dp) :: F_psi

        if (.not. self%initialized) then
            error stop "analytical_circular_eq_t not initialized"
        end if

        ! Get psi derivatives
        call self%eval_psi_derivatives(R, Z, dpsi_dR, dpsi_dZ)

        ! Poloidal field components from psi
        B_R = -(1.0_dp / R) * dpsi_dZ
        B_Z = (1.0_dp / R) * dpsi_dR

        ! Toroidal field (F(psi) ≈ B0*R0 for low-beta)
        F_psi = self%B0 * self%R0
        B_phi = F_psi / R

        ! Total field magnitude
        B_mod = sqrt(B_R**2 + B_Z**2 + B_phi**2)
    end subroutine evaluate_bfield

    !> Evaluate magnetic field components with TF ripple
    !>
    !> Parameters:
    !>   R   - Major radius coordinate [m]
    !>   phi - Toroidal angle [rad]
    !>   Z   - Vertical coordinate [m]
    !>
    !> Returns:
    !>   B_R   - Radial field component [T]
    !>   B_Z   - Vertical field component [T]
    !>   B_phi - Toroidal field component with ripple [T]
    !>   B_mod - Total field magnitude [T]
    subroutine eval_bfield_ripple(self, R, phi, Z, B_R, B_Z, B_phi, B_mod)
        class(analytical_circular_eq_t), intent(in) :: self
        real(dp), intent(in) :: R, phi, Z
        real(dp), intent(out) :: B_R, B_Z, B_phi, B_mod

        real(dp) :: dpsi_dR, dpsi_dZ
        real(dp) :: F_psi
        real(dp) :: radius, theta, delta_ripple

        if (.not. self%initialized) then
            error stop "analytical_circular_eq_t not initialized"
        end if

        call self%eval_psi_derivatives(R, Z, dpsi_dR, dpsi_dZ)

        B_R = -(1.0_dp / R) * dpsi_dZ
        B_Z = (1.0_dp / R) * dpsi_dR

        F_psi = self%B0 * self%R0
        B_phi = F_psi / R

        if (self%Nripple > 0) then
            radius = sqrt((R - self%R0)**2 + (Z - self%z0)**2)
            theta = atan2(Z - self%z0, R - self%R0)
            delta_ripple = self%delta0 * exp(-0.5_dp * theta**2) &
                         * (radius / self%a0)**self%alpha0
            B_phi = B_phi * (1.0_dp + delta_ripple * cos(real(self%Nripple, dp) * phi))
        end if

        B_mod = sqrt(B_R**2 + B_Z**2 + B_phi**2)
    end subroutine eval_bfield_ripple

    !> Cleanup
    subroutine destroy_equilibrium(self)
        class(analytical_circular_eq_t), intent(inout) :: self
        self%initialized = .false.
    end subroutine destroy_equilibrium

    !> Standalone interface function for external use (e.g., SIMPLE)
    !>
    !> Evaluates magnetic field in cylindrical coordinates with optional ripple.
    !> This function provides a simple interface matching the pattern used by
    !> other libneo field modules (e.g., neo_biotsavart).
    !>
    !> Parameters:
    !>   eq    - Initialized analytical equilibrium
    !>   R     - Major radius [m]
    !>   phi   - Toroidal angle [rad]
    !>   Z     - Vertical coordinate [m]
    !>
    !> Returns:
    !>   B_R   - Radial field component [T]
    !>   B_Z   - Vertical field component [T]
    !>   B_phi - Toroidal field component (with ripple if enabled) [T]
    !>   B_mod - Total field magnitude [T]
    subroutine compute_analytical_field_cylindrical(eq, R, phi, Z, B_R, B_Z, B_phi, B_mod)
        type(analytical_circular_eq_t), intent(in) :: eq
        real(dp), intent(in) :: R, phi, Z
        real(dp), intent(out) :: B_R, B_Z, B_phi, B_mod

        call eq%eval_bfield_ripple(R, phi, Z, B_R, B_Z, B_phi, B_mod)
    end subroutine compute_analytical_field_cylindrical

    !> field_t interface: compute vector potential
    !> x(3) = (R, phi, Z) in cylindrical coordinates
    subroutine compute_afield(self, x, A)
        class(analytical_circular_eq_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: A(3)

        error stop 'analytical_circular_eq_t: compute_afield not implemented'
    end subroutine compute_afield

    !> field_t interface: compute magnetic field
    !> x(3) = (R, phi, Z) in cylindrical coordinates
    subroutine compute_bfield(self, x, B)
        class(analytical_circular_eq_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: B(3)

        real(dp) :: B_R, B_Z, B_phi, B_mod

        call self%eval_bfield_ripple(x(1), x(2), x(3), B_R, B_Z, B_phi, B_mod)

        B(1) = B_R
        B(2) = B_phi
        B(3) = B_Z
    end subroutine compute_bfield

    !> field_t interface: compute vector potential and field
    subroutine compute_abfield(self, x, A, B)
        class(analytical_circular_eq_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: A(3), B(3)

        call self%compute_afield(x, A)
        call self%compute_bfield(x, B)
    end subroutine compute_abfield

end module analytical_tokamak_field
