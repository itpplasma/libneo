module analytical_gs_solver
    !> Analytical Grad-Shafranov solver for tokamak equilibria
    !>
    !> Implementation based on:
    !>   - Cerfon & Freidberg, One size fits all analytic solutions,
    !>     Physics of Plasmas 17, 032502 (2010), DOI: 10.1063/1.3328818
    !>   - Verena Eslbauer, Two analytical solutions to the Grad-Shafranov
    !>     equation using Solovev pressure and poloidal current profiles,
    !>     Bachelor thesis, TU Graz, November 20, 2017
    !>
    !> Fortran port of the MATLAB implementation from the Eslbauer thesis.
    !> Supports general shaped tokamaks with elongation (κ) and triangularity (δ).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none

    private
    public :: psi_1, psi_2, psi_3, psi_4, psi_5, psi_6, psi_7
    public :: dpsi_1_derivatives, dpsi_2_derivatives, dpsi_3_derivatives
    public :: dpsi_4_derivatives, dpsi_5_derivatives, dpsi_6_derivatives, dpsi_7_derivatives
    public :: d2psi_1_dxx, d2psi_1_dyy, d2psi_2_dxx, d2psi_2_dyy
    public :: d2psi_3_dxx, d2psi_3_dyy, d2psi_4_dxx, d2psi_4_dyy
    public :: d2psi_5_dxx, d2psi_5_dyy, d2psi_6_dxx, d2psi_6_dyy
    public :: d2psi_7_dxx, d2psi_7_dyy
    public :: psi_p, dpsi_p_dx, dpsi_p_dy, d2psi_p_dxx, d2psi_p_dyy
    public :: solve_coefficients

contains

    !> Basis function ψ₁ = 1
    pure function psi_1(x, y) result(psi)
        real(dp), intent(in) :: x, y
        real(dp) :: psi
        psi = 1.0_dp
    end function psi_1

    !> Basis function ψ₂ = x²
    pure function psi_2(x, y) result(psi)
        real(dp), intent(in) :: x, y
        real(dp) :: psi
        psi = x**2
    end function psi_2

    !> Basis function ψ₃ = y² - x²ln(x)
    pure function psi_3(x, y) result(psi)
        real(dp), intent(in) :: x, y
        real(dp) :: psi
        psi = y**2 - x**2 * log(x)
    end function psi_3

    !> Basis function ψ₄ = x⁴ - 4x²y²
    pure function psi_4(x, y) result(psi)
        real(dp), intent(in) :: x, y
        real(dp) :: psi
        psi = x**4 - 4.0_dp * x**2 * y**2
    end function psi_4

    !> Basis function ψ₅ = 2y⁴ - 9x²y² + [3x⁴ - 12x²y²]ln(x)
    pure function psi_5(x, y) result(psi)
        real(dp), intent(in) :: x, y
        real(dp) :: psi
        real(dp) :: logx
        logx = log(x)
        psi = 2.0_dp * y**4 - 9.0_dp * x**2 * y**2 &
            + (3.0_dp * x**4 - 12.0_dp * x**2 * y**2) * logx
    end function psi_5

    !> Basis function ψ₆ = x⁶ - 12x⁴y² + 8x²y⁴
    pure function psi_6(x, y) result(psi)
        real(dp), intent(in) :: x, y
        real(dp) :: psi
        psi = x**6 - 12.0_dp * x**4 * y**2 + 8.0_dp * x**2 * y**4
    end function psi_6

    !> Basis function ψ₇ = 8y⁶ - 140x²y⁴ + 75x⁴y² + [180x⁴y² - 120x²y⁴ - 15x⁶]ln(x)
    pure function psi_7(x, y) result(psi)
        real(dp), intent(in) :: x, y
        real(dp) :: psi
        real(dp) :: logx
        logx = log(x)
        psi = 8.0_dp * y**6 - 140.0_dp * x**2 * y**4 + 75.0_dp * x**4 * y**2 &
            + (180.0_dp * x**4 * y**2 - 120.0_dp * x**2 * y**4 - 15.0_dp * x**6) * logx
    end function psi_7

    !> Particular solution ψ_p = x⁴/8 + A(x²ln(x)/2 - x⁴/8)
    pure function psi_p(x, y, A_param) result(psi)
        real(dp), intent(in) :: x, y, A_param
        real(dp) :: psi
        psi = x**4 / 8.0_dp + A_param * (x**2 * log(x) / 2.0_dp - x**4 / 8.0_dp)
    end function psi_p

    ! First derivatives of basis functions
    pure subroutine dpsi_1_derivatives(x, y, dpsi_dx, dpsi_dy)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: dpsi_dx, dpsi_dy
        dpsi_dx = 0.0_dp
        dpsi_dy = 0.0_dp
    end subroutine dpsi_1_derivatives

    pure subroutine dpsi_2_derivatives(x, y, dpsi_dx, dpsi_dy)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: dpsi_dx, dpsi_dy
        dpsi_dx = 2.0_dp * x
        dpsi_dy = 0.0_dp
    end subroutine dpsi_2_derivatives

    pure subroutine dpsi_3_derivatives(x, y, dpsi_dx, dpsi_dy)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: dpsi_dx, dpsi_dy
        dpsi_dx = -2.0_dp * x * log(x) - x
        dpsi_dy = 2.0_dp * y
    end subroutine dpsi_3_derivatives

    pure subroutine dpsi_4_derivatives(x, y, dpsi_dx, dpsi_dy)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: dpsi_dx, dpsi_dy
        dpsi_dx = 4.0_dp * x**3 - 8.0_dp * x * y**2
        dpsi_dy = -8.0_dp * x**2 * y
    end subroutine dpsi_4_derivatives

    pure subroutine dpsi_5_derivatives(x, y, dpsi_dx, dpsi_dy)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: dpsi_dx, dpsi_dy
        real(dp) :: logx
        logx = log(x)
        dpsi_dx = -18.0_dp * x * y**2 &
                + 12.0_dp * x**3 * logx + 3.0_dp * x**3 &
                - 24.0_dp * x * y**2 * logx - 12.0_dp * x * y**2
        dpsi_dy = 8.0_dp * y**3 - 18.0_dp * x**2 * y - 24.0_dp * x**2 * y * logx
    end subroutine dpsi_5_derivatives

    pure subroutine dpsi_6_derivatives(x, y, dpsi_dx, dpsi_dy)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: dpsi_dx, dpsi_dy
        dpsi_dx = 6.0_dp * x**5 - 48.0_dp * x**3 * y**2 + 16.0_dp * x * y**4
        dpsi_dy = -24.0_dp * x**4 * y + 32.0_dp * x**2 * y**3
    end subroutine dpsi_6_derivatives

    pure subroutine dpsi_7_derivatives(x, y, dpsi_dx, dpsi_dy)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: dpsi_dx, dpsi_dy
        real(dp) :: logx
        logx = log(x)
        dpsi_dx = -280.0_dp * x * y**4 + 150.0_dp * x**3 * y**2 &
                + 720.0_dp * x**3 * y**2 * logx + 330.0_dp * x**3 * y**2 &
                - 240.0_dp * x * y**4 * logx - 120.0_dp * x * y**4 &
                - 90.0_dp * x**5 * logx - 15.0_dp * x**5
        dpsi_dy = 48.0_dp * y**5 - 560.0_dp * x**2 * y**3 + 150.0_dp * x**4 * y &
                + 360.0_dp * x**4 * y * logx - 480.0_dp * x**2 * y**3 * logx
    end subroutine dpsi_7_derivatives

    ! Particular solution derivatives
    pure function dpsi_p_dx(x, y, A_param) result(dpsi_dx)
        real(dp), intent(in) :: x, y, A_param
        real(dp) :: dpsi_dx
        dpsi_dx = x**3 / 2.0_dp + A_param * (x * log(x) + x / 2.0_dp - x**3 / 2.0_dp)
    end function dpsi_p_dx

    pure function dpsi_p_dy(x, y, A_param) result(dpsi_dy)
        real(dp), intent(in) :: x, y, A_param
        real(dp) :: dpsi_dy
        dpsi_dy = 0.0_dp
    end function dpsi_p_dy

    ! Second derivatives of basis functions
    pure function d2psi_1_dxx(x, y) result(d2psi)
        real(dp), intent(in) :: x, y
        real(dp) :: d2psi
        d2psi = 0.0_dp
    end function d2psi_1_dxx

    pure function d2psi_1_dyy(x, y) result(d2psi)
        real(dp), intent(in) :: x, y
        real(dp) :: d2psi
        d2psi = 0.0_dp
    end function d2psi_1_dyy

    pure function d2psi_2_dxx(x, y) result(d2psi)
        real(dp), intent(in) :: x, y
        real(dp) :: d2psi
        d2psi = 2.0_dp
    end function d2psi_2_dxx

    pure function d2psi_2_dyy(x, y) result(d2psi)
        real(dp), intent(in) :: x, y
        real(dp) :: d2psi
        d2psi = 0.0_dp
    end function d2psi_2_dyy

    pure function d2psi_3_dxx(x, y) result(d2psi)
        real(dp), intent(in) :: x, y
        real(dp) :: d2psi
        d2psi = -2.0_dp * log(x) - 3.0_dp
    end function d2psi_3_dxx

    pure function d2psi_3_dyy(x, y) result(d2psi)
        real(dp), intent(in) :: x, y
        real(dp) :: d2psi
        d2psi = 2.0_dp
    end function d2psi_3_dyy

    pure function d2psi_4_dxx(x, y) result(d2psi)
        real(dp), intent(in) :: x, y
        real(dp) :: d2psi
        d2psi = 12.0_dp * x**2 - 8.0_dp * y**2
    end function d2psi_4_dxx

    pure function d2psi_4_dyy(x, y) result(d2psi)
        real(dp), intent(in) :: x, y
        real(dp) :: d2psi
        d2psi = -8.0_dp * x**2
    end function d2psi_4_dyy

    pure function d2psi_5_dxx(x, y) result(d2psi)
        real(dp), intent(in) :: x, y
        real(dp) :: d2psi
        real(dp) :: logx
        logx = log(x)
        d2psi = -18.0_dp * y**2 + 36.0_dp * x**2 * logx + 21.0_dp * x**2 - 24.0_dp * y**2 * logx - 36.0_dp * y**2
    end function d2psi_5_dxx

    pure function d2psi_5_dyy(x, y) result(d2psi)
        real(dp), intent(in) :: x, y
        real(dp) :: d2psi
        real(dp) :: logx
        logx = log(x)
        d2psi = 24.0_dp * y**2 - 18.0_dp * x**2 - 24.0_dp * x**2 * logx
    end function d2psi_5_dyy

    pure function d2psi_6_dxx(x, y) result(d2psi)
        real(dp), intent(in) :: x, y
        real(dp) :: d2psi
        d2psi = 30.0_dp * x**4 - 144.0_dp * x**2 * y**2 + 16.0_dp * y**4
    end function d2psi_6_dxx

    pure function d2psi_6_dyy(x, y) result(d2psi)
        real(dp), intent(in) :: x, y
        real(dp) :: d2psi
        d2psi = -24.0_dp * x**4 + 96.0_dp * x**2 * y**2
    end function d2psi_6_dyy

    pure function d2psi_7_dxx(x, y) result(d2psi)
        real(dp), intent(in) :: x, y
        real(dp) :: d2psi
        real(dp) :: logx
        logx = log(x)
        d2psi = -280.0_dp * y**4 + 450.0_dp * x**2 * y**2 &
              + 2160.0_dp * x**2 * y**2 * logx + 1710.0_dp * x**2 * y**2 &
              - 240.0_dp * y**4 * logx - 360.0_dp * y**4 &
              - 450.0_dp * x**4 * logx - 165.0_dp * x**4
    end function d2psi_7_dxx

    pure function d2psi_7_dyy(x, y) result(d2psi)
        real(dp), intent(in) :: x, y
        real(dp) :: d2psi
        real(dp) :: logx
        logx = log(x)
        d2psi = 240.0_dp * y**4 - 1680.0_dp * x**2 * y**2 + 150.0_dp * x**4 &
              + 360.0_dp * x**4 * logx - 1440.0_dp * x**2 * y**2 * logx
    end function d2psi_7_dyy

    ! Second derivatives of particular solution
    pure function d2psi_p_dxx(x, y, A_param) result(d2psi)
        real(dp), intent(in) :: x, y, A_param
        real(dp) :: d2psi
        d2psi = 1.5_dp * x**2 + A_param * (log(x) + 1.5_dp - 1.5_dp * x**2)
    end function d2psi_p_dxx

    pure function d2psi_p_dyy(x, y, A_param) result(d2psi)
        real(dp), intent(in) :: x, y, A_param
        real(dp) :: d2psi
        d2psi = 0.0_dp
    end function d2psi_p_dyy

    !> Solve for coefficients c₁-c₇ given plasma shape parameters
    !>
    !> Parameters:
    !>   epsilon - Inverse aspect ratio (a/R₀)
    !>   kappa   - Elongation
    !>   delta   - Triangularity
    !>   A_param - Shafranov parameter
    !>   coeffs  - Output: 7 coefficients
    subroutine solve_coefficients(epsilon, kappa, delta, A_param, coeffs)
        real(dp), intent(in) :: epsilon, kappa, delta, A_param
        real(dp), intent(out) :: coeffs(7)

        real(dp) :: mat(7,7), rhs(7)
        real(dp) :: x_out, x_in, x_high, y_high
        real(dp) :: alpha, N1, N2, N3
        real(dp) :: dpsi_dx, dpsi_dy
        integer :: ipiv(7), info

        ! Boundary points from parametric representation
        x_out = 1.0_dp + epsilon        ! τ=0 (outer equatorial)
        x_in = 1.0_dp - epsilon         ! τ=π (inner equatorial)
        x_high = 1.0_dp - delta*epsilon ! τ=π/2 (high point)
        y_high = kappa*epsilon

        ! Curvature coefficients
        alpha = asin(delta)
        N1 = -(1.0_dp + alpha)**2 / (epsilon * kappa**2)
        N2 = (1.0_dp - alpha)**2 / (epsilon * kappa**2)
        N3 = -kappa / (epsilon * cos(alpha)**2)

        ! Row 1: ψ(x_out, 0) = 0
        mat(1,1) = psi_1(x_out, 0.0_dp)
        mat(1,2) = psi_2(x_out, 0.0_dp)
        mat(1,3) = psi_3(x_out, 0.0_dp)
        mat(1,4) = psi_4(x_out, 0.0_dp)
        mat(1,5) = psi_5(x_out, 0.0_dp)
        mat(1,6) = psi_6(x_out, 0.0_dp)
        mat(1,7) = psi_7(x_out, 0.0_dp)
        rhs(1) = -psi_p(x_out, 0.0_dp, A_param)

        ! Row 2: ψ(x_in, 0) = 0
        mat(2,1) = psi_1(x_in, 0.0_dp)
        mat(2,2) = psi_2(x_in, 0.0_dp)
        mat(2,3) = psi_3(x_in, 0.0_dp)
        mat(2,4) = psi_4(x_in, 0.0_dp)
        mat(2,5) = psi_5(x_in, 0.0_dp)
        mat(2,6) = psi_6(x_in, 0.0_dp)
        mat(2,7) = psi_7(x_in, 0.0_dp)
        rhs(2) = -psi_p(x_in, 0.0_dp, A_param)

        ! Row 3: ψ(x_high, y_high) = 0
        mat(3,1) = psi_1(x_high, y_high)
        mat(3,2) = psi_2(x_high, y_high)
        mat(3,3) = psi_3(x_high, y_high)
        mat(3,4) = psi_4(x_high, y_high)
        mat(3,5) = psi_5(x_high, y_high)
        mat(3,6) = psi_6(x_high, y_high)
        mat(3,7) = psi_7(x_high, y_high)
        rhs(3) = -psi_p(x_high, y_high, A_param)

        ! Row 4: ∂ψ/∂x(x_high, y_high) = 0
        call dpsi_1_derivatives(x_high, y_high, mat(4,1), dpsi_dy)
        call dpsi_2_derivatives(x_high, y_high, mat(4,2), dpsi_dy)
        call dpsi_3_derivatives(x_high, y_high, mat(4,3), dpsi_dy)
        call dpsi_4_derivatives(x_high, y_high, mat(4,4), dpsi_dy)
        call dpsi_5_derivatives(x_high, y_high, mat(4,5), dpsi_dy)
        call dpsi_6_derivatives(x_high, y_high, mat(4,6), dpsi_dy)
        call dpsi_7_derivatives(x_high, y_high, mat(4,7), dpsi_dy)
        rhs(4) = -dpsi_p_dx(x_high, y_high, A_param)

        ! Row 5: ∂²ψ/∂y²(x_out, 0) + N₁·∂ψ/∂x(x_out, 0) = 0
        call dpsi_1_derivatives(x_out, 0.0_dp, dpsi_dx, dpsi_dy)
        mat(5,1) = d2psi_1_dyy(x_out, 0.0_dp) + N1 * dpsi_dx
        call dpsi_2_derivatives(x_out, 0.0_dp, dpsi_dx, dpsi_dy)
        mat(5,2) = d2psi_2_dyy(x_out, 0.0_dp) + N1 * dpsi_dx
        call dpsi_3_derivatives(x_out, 0.0_dp, dpsi_dx, dpsi_dy)
        mat(5,3) = d2psi_3_dyy(x_out, 0.0_dp) + N1 * dpsi_dx
        call dpsi_4_derivatives(x_out, 0.0_dp, dpsi_dx, dpsi_dy)
        mat(5,4) = d2psi_4_dyy(x_out, 0.0_dp) + N1 * dpsi_dx
        call dpsi_5_derivatives(x_out, 0.0_dp, dpsi_dx, dpsi_dy)
        mat(5,5) = d2psi_5_dyy(x_out, 0.0_dp) + N1 * dpsi_dx
        call dpsi_6_derivatives(x_out, 0.0_dp, dpsi_dx, dpsi_dy)
        mat(5,6) = d2psi_6_dyy(x_out, 0.0_dp) + N1 * dpsi_dx
        call dpsi_7_derivatives(x_out, 0.0_dp, dpsi_dx, dpsi_dy)
        mat(5,7) = d2psi_7_dyy(x_out, 0.0_dp) + N1 * dpsi_dx
        rhs(5) = -(d2psi_p_dyy(x_out, 0.0_dp, A_param) + N1 * dpsi_p_dx(x_out, 0.0_dp, A_param))

        ! Row 6: ∂²ψ/∂y²(x_in, 0) + N₂·∂ψ/∂x(x_in, 0) = 0
        call dpsi_1_derivatives(x_in, 0.0_dp, dpsi_dx, dpsi_dy)
        mat(6,1) = d2psi_1_dyy(x_in, 0.0_dp) + N2 * dpsi_dx
        call dpsi_2_derivatives(x_in, 0.0_dp, dpsi_dx, dpsi_dy)
        mat(6,2) = d2psi_2_dyy(x_in, 0.0_dp) + N2 * dpsi_dx
        call dpsi_3_derivatives(x_in, 0.0_dp, dpsi_dx, dpsi_dy)
        mat(6,3) = d2psi_3_dyy(x_in, 0.0_dp) + N2 * dpsi_dx
        call dpsi_4_derivatives(x_in, 0.0_dp, dpsi_dx, dpsi_dy)
        mat(6,4) = d2psi_4_dyy(x_in, 0.0_dp) + N2 * dpsi_dx
        call dpsi_5_derivatives(x_in, 0.0_dp, dpsi_dx, dpsi_dy)
        mat(6,5) = d2psi_5_dyy(x_in, 0.0_dp) + N2 * dpsi_dx
        call dpsi_6_derivatives(x_in, 0.0_dp, dpsi_dx, dpsi_dy)
        mat(6,6) = d2psi_6_dyy(x_in, 0.0_dp) + N2 * dpsi_dx
        call dpsi_7_derivatives(x_in, 0.0_dp, dpsi_dx, dpsi_dy)
        mat(6,7) = d2psi_7_dyy(x_in, 0.0_dp) + N2 * dpsi_dx
        rhs(6) = -(d2psi_p_dyy(x_in, 0.0_dp, A_param) + N2 * dpsi_p_dx(x_in, 0.0_dp, A_param))

        ! Row 7: ∂²ψ/∂x²(x_high, y_high) + N₃·∂ψ/∂y(x_high, y_high) = 0
        call dpsi_1_derivatives(x_high, y_high, dpsi_dx, dpsi_dy)
        mat(7,1) = d2psi_1_dxx(x_high, y_high) + N3 * dpsi_dy
        call dpsi_2_derivatives(x_high, y_high, dpsi_dx, dpsi_dy)
        mat(7,2) = d2psi_2_dxx(x_high, y_high) + N3 * dpsi_dy
        call dpsi_3_derivatives(x_high, y_high, dpsi_dx, dpsi_dy)
        mat(7,3) = d2psi_3_dxx(x_high, y_high) + N3 * dpsi_dy
        call dpsi_4_derivatives(x_high, y_high, dpsi_dx, dpsi_dy)
        mat(7,4) = d2psi_4_dxx(x_high, y_high) + N3 * dpsi_dy
        call dpsi_5_derivatives(x_high, y_high, dpsi_dx, dpsi_dy)
        mat(7,5) = d2psi_5_dxx(x_high, y_high) + N3 * dpsi_dy
        call dpsi_6_derivatives(x_high, y_high, dpsi_dx, dpsi_dy)
        mat(7,6) = d2psi_6_dxx(x_high, y_high) + N3 * dpsi_dy
        call dpsi_7_derivatives(x_high, y_high, dpsi_dx, dpsi_dy)
        mat(7,7) = d2psi_7_dxx(x_high, y_high) + N3 * dpsi_dy
        rhs(7) = -(d2psi_p_dxx(x_high, y_high, A_param) + N3 * dpsi_p_dy(x_high, y_high, A_param))

        ! Solve 7×7 linear system using LAPACK
        call dgesv(7, 1, mat, 7, ipiv, rhs, 7, info)

        if (info /= 0) then
            error stop "analytical_gs_solver: Boundary condition matrix is singular"
        end if

        coeffs = rhs
    end subroutine solve_coefficients

end module analytical_gs_solver
