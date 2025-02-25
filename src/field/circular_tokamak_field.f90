module neo_circular_tokamak_field
    use libneo_kinds, only : dp
    use neo_field_base, only: field_t
    implicit none

    type, extends(field_t) :: circular_tokamak_field_t
        real(dp) :: R_axis, Z_axis
        real(dp) :: B_pol_ampl, B_tor_ampl
        contains
        procedure :: circular_tokamak_field_init
        procedure :: compute_afield
        procedure :: compute_bfield
        procedure :: compute_abfield
    end type circular_tokamak_field_t

    contains

    subroutine circular_tokamak_field_init(self, R_axis, Z_axis, B_pol_ampl, B_tor_ampl)
        class(circular_tokamak_field_t), intent(out) :: self
        real(dp), intent(in), optional :: R_axis, Z_axis
        real(dp), intent(in), optional :: B_pol_ampl, B_tor_ampl

        if (present(R_axis)) then
            self%R_axis = R_axis
        else
            self%R_axis = 1.0_dp
        end if
        if (present(Z_axis)) then
            self%Z_axis = Z_axis
        else
            self%Z_axis = 0.0_dp
        end if
        if (present(B_pol_ampl)) then
            self%B_pol_ampl = B_pol_ampl
        else
            self%B_pol_ampl = 0.314159_dp
        end if
        if (present(B_tor_ampl)) then
            self%B_tor_ampl = B_tor_ampl
        else
            self%B_tor_ampl = 1.0_dp
        end if
    end subroutine circular_tokamak_field_init

    subroutine compute_afield(self, x, A)
        class(circular_tokamak_field_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: A(3)

        A(1) = 0.5d0 * self%B_tor_ampl * x(3)
        A(2) = 0.5d0 * self%B_pol_ampl * &
                       ((x(1) - self%R_axis)**2 - (x(3) - self%Z_axis)**2)
        A(3) = 0.5d0 * self%B_tor_ampl * x(1)
    end subroutine compute_afield

    subroutine compute_bfield(self, x, B)
        class(circular_tokamak_field_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: B(3)

        if (x(1) .gt. 0.0d0) then
            B(1) = - self%B_pol_ampl / x(1) * (x(3) - self%Z_axis)
            B(2) = self%B_tor_ampl
            B(3) = self%B_pol_ampl / x(1) * (x(1) - self%R_axis)
        else
            error stop 'x(1)=R must be positive'
        end if
    end subroutine compute_bfield

    subroutine compute_abfield(self, x, A, B)
        class(circular_tokamak_field_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: A(3), B(3)

        call self%compute_afield(x, A)
        call self%compute_bfield(x, B)
    end subroutine compute_abfield

    end module neo_circular_tokamak_field
