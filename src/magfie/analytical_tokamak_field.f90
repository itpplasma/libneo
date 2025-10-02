module analytical_tokamak_field
    !> Analytical tokamak field evaluation using Cerfon-Freidberg equilibria
    !>
    !> Clean-room implementation from Cerfon & Freidberg, Phys. Plasmas 17, 032502 (2010)
    !> Provides magnetic field components for circular tokamak equilibria
    use iso_fortran_env, only: dp => real64
    use analytical_gs_solver
    implicit none

    private
    public :: analytical_circular_eq_t

    !> Analytical tokamak equilibrium (supports arbitrary shaping)
    type :: analytical_circular_eq_t
        real(dp) :: R0         !< Major radius [m]
        real(dp) :: epsilon    !< Inverse aspect ratio
        real(dp) :: kappa      !< Elongation
        real(dp) :: delta      !< Triangularity
        real(dp) :: A_param    !< Shafranov parameter
        real(dp) :: B0         !< Toroidal field on axis [T]
        real(dp) :: psimult    !< Flux scaling factor
        real(dp) :: coeffs(7) !< Boundary condition coefficients
        logical :: initialized = .false.
    contains
        procedure :: init => init_circular_equilibrium
        procedure :: eval_psi => evaluate_psi
        procedure :: eval_psi_derivatives
        procedure :: eval_bfield => evaluate_bfield
        procedure :: cleanup => destroy_equilibrium
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
    subroutine init_circular_equilibrium(self, R0_in, epsilon_in, kappa_in, delta_in, A_in, B0_in, psimult_in)
        class(analytical_circular_eq_t), intent(inout) :: self
        real(dp), intent(in) :: R0_in
        real(dp), intent(in) :: epsilon_in
        real(dp), intent(in), optional :: kappa_in
        real(dp), intent(in), optional :: delta_in
        real(dp), intent(in) :: A_in
        real(dp), intent(in) :: B0_in
        real(dp), intent(in), optional :: psimult_in

        self%R0 = R0_in
        self%epsilon = epsilon_in
        self%A_param = A_in
        self%B0 = B0_in

        if (present(kappa_in)) then
            self%kappa = kappa_in
        else
            self%kappa = 1.0_dp  ! Circular
        end if

        if (present(delta_in)) then
            self%delta = delta_in
        else
            self%delta = 0.0_dp  ! Circular
        end if

        if (present(psimult_in)) then
            self%psimult = psimult_in
        else
            self%psimult = 200.0_dp  ! Standard normalization for ITER
        end if

        ! Solve for coefficients with specified shape
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

    !> Cleanup
    subroutine destroy_equilibrium(self)
        class(analytical_circular_eq_t), intent(inout) :: self
        self%initialized = .false.
    end subroutine destroy_equilibrium

end module analytical_tokamak_field
