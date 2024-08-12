module neo_example_field
use, intrinsic :: iso_fortran_env, only: dp => real64
use neo_field_base, only: field_t
implicit none

type, extends(field_t) :: example_field_t
    real(dp) :: ampl, ampl2
    contains
    procedure :: example_field_init
    procedure :: compute_afield
    procedure :: compute_bfield
    procedure :: compute_abfield
end type example_field_t

contains

subroutine example_field_init(self, ampl, ampl2)
    class(example_field_t), intent(out) :: self
    real(dp), intent(in), optional :: ampl, ampl2

    if (present(ampl)) then
        self%ampl = ampl
    else
        self%ampl = 1.0d-7
    end if
    if (present(ampl2)) then 
        self%ampl2 = ampl2
    else
        self%ampl2 = 2.0d-6
    end if
end subroutine example_field_init

subroutine compute_afield(self, x, A)
    class(example_field_t), intent(in) :: self
    real(dp), intent(in) :: x(3)
    real(dp), intent(out) :: A(3)

    real(dp) :: R, phi, Z
    real(dp) :: AR, Aphi, AZ

    R = sqrt(x(1)**2 + x(2)**2)
    phi = atan2(x(2), x(1))
    Z = x(3)

    AR = self%ampl * R * cos(phi)
    Aphi = -0.5d0 * self%ampl2 * Z * R
    AZ = -log(R)

    A(1) = AR * cos(phi) - Aphi * sin(phi)
    A(2) = AR * sin(phi) + Aphi * cos(phi)
    A(3) = AZ
end subroutine compute_afield

subroutine compute_bfield(self, x, B)
    class(example_field_t), intent(in) :: self
    real(dp), intent(in) :: x(3)
    real(dp), intent(out) :: B(3)

    real(dp) :: R, phi, Z
    real(dp) :: BR, Bphi, BZ

    R = sqrt(x(1)**2 + x(2)**2)
    phi = atan2(x(2), x(1))
    Z = x(3)

    BR = 0.5d0 * self%ampl2 * R
    Bphi = 1.0d0 / R
    BZ = -self%ampl2 * Z + self%ampl * sin(phi)

    B(1) = BR * cos(phi) - Bphi * sin(phi)
    B(2) = BR * sin(phi) + Bphi * cos(phi)
    B(3) = BZ
end subroutine compute_bfield

subroutine compute_abfield(self, x, A, B)
    class(example_field_t), intent(in) :: self
    real(dp), intent(in) :: x(3)
    real(dp), intent(out) :: A(3), B(3)

    call self%compute_afield(x, A)

    call self%compute_bfield(x, B)
end subroutine compute_abfield

end module neo_example_field
