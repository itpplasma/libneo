module libneo_coordinate_conventions
    !> Coordinate convention constants for magnetic confinement geometry.
    !>
    !> This module defines standard conventions for radial and toroidal angle
    !> coordinates used across different equilibrium representations (VMEC,
    !> Boozer, chartmap, etc.).
    !>
    !> UNKNOWN = -1 is used for any unknown/unspecified convention,
    !> matching SIMPLE convention where TEST=-1 is the special case.
    implicit none

    ! Universal unknown value (used for both zeta and rho conventions)
    integer, parameter :: UNKNOWN = -1    !< Unknown/unspecified convention

    ! Zeta (toroidal angle) convention constants
    integer, parameter :: CYL = 0         !< Cylindrical phi
    integer, parameter :: VMEC = 1        !< VMEC toroidal angle
    integer, parameter :: BOOZER = 2      !< Boozer toroidal angle

    ! Rho (radial) convention constants
    !> RHO_TOR = sqrt(psi_tor/psi_tor_edge) = sqrt(s)
    !> RHO_POL = sqrt(psi_pol/psi_pol_edge)
    !> PSI_TOR_NORM = psi_tor/psi_tor_edge (= s in VMEC)
    !> PSI_POL_NORM = psi_pol/psi_pol_edge
    integer, parameter :: RHO_TOR = 1         !< sqrt(normalized toroidal flux)
    integer, parameter :: RHO_POL = 2         !< sqrt(normalized poloidal flux)
    integer, parameter :: PSI_TOR_NORM = 3    !< Normalized toroidal flux (= s)
    integer, parameter :: PSI_POL_NORM = 4    !< Normalized poloidal flux

end module libneo_coordinate_conventions
