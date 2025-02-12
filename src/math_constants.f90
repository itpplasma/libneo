!> Despite the name this module contains mathematical and physical
!> constants, as it is sometimes hard to distinguish, e.g. for
!> conversion factors.
module math_constants
  use libneo_kinds, only : complex_kind, real_kind

  real(kind=real_kind), parameter :: PI = 3.141592653589793238462643383279502884197_real_kind
  !> \f[ \pi/2 \f]
  real(kind=real_kind), parameter :: PIO2 = 1.57079632679489661923132169163975144209858_real_kind
  !> \f[ 2\pi \f]
  real(kind=real_kind), parameter :: TWOPI = 6.283185307179586476925286766559005768394_real_kind
  !> \f[ \sqrt{\pi} \f]
  real(kind=real_kind), parameter :: SQRT2 = 1.41421356237309504880168872420969807856967_real_kind
  !> Defined as \f$ \gamma = \int\limits_{1}^{\infty} 1/floor(x) - 1/x dx \f$
  !> or alternatively \f$ \sum\limits_{k=1}^{\infty} ( 1/k - \ln ((k+1)/k)\f$
  !> or \f$ \Gamma_prime(1) = - \gamma \f$.
  real(kind=real_kind), parameter :: EULER = 0.5772156649015328606065120900824024310422_real_kind

  !> Imaginary unit:
  complex(kind=complex_kind), parameter :: IMUN = (0.d0,1.d0)

  ! Define physical constants (SI-units)
  real(kind=real_kind), parameter, public :: MU_0 = 1.25663706212e-6_real_kind !< vacuum permeability [N A^-2]

  ! Define physical constants (cgs-units)
  real(kind=real_kind), parameter, public :: C = 2.99792458e10_real_kind   !< speed of light, [cm/s]
  real(kind=real_kind), parameter, public :: E = 4.8032e-10_real_kind      !< elementary charge [StatC]
  real(kind=real_kind), parameter, public :: U = 1.660539040e-24_real_kind !< atomic mass unit [g]

  ! Masses of some particles/atoms. Atoms get designated with m\_ and
  ! symbol (e.g. D for deuterium).
  real(kind=real_kind), parameter, public :: m_e = 9.1093837139e-28 !< Electron mass [g]
  real(kind=real_kind), parameter, public :: m_p = 1.672621637e-24 !< Proton mass [g]
  real(kind=real_kind), parameter, public :: m_a = 6.644656200e-24 !< alpha particle mass [g]
  real(kind=real_kind), parameter, public :: m_D = 3.343583719e-24 !< Deuterium mass [g]
  real(kind=real_kind), parameter, public :: m_C = 19.94406876e-24 !< Carbon mass [g]
  real(kind=real_kind), parameter, public :: m_W = 305.2734971e-24 !< Tungsten mass [g]

  ! Mass ratios.
  real(kind=real_kind), parameter, public :: mc_o_e  = 5.6856793e-8 !< cgs
  real(kind=real_kind), parameter, public :: mc_o_2e = 2.8428397e-8 !< cgs
  real(kind=real_kind), parameter, public :: e_o_mc  = 1.7588048e7  !< cgs

  ! Conversion of units.
  real(kind=real_kind), parameter, public :: magneticfield_si_to_cgs = 1.0e4_real_kind  !< SI to cgs for magnetic fields, i.e. convert from T to Gauss.
  real(kind=real_kind), parameter, public :: magneticflux_si_to_cgs = 1.0e8_real_kind  !< SI to cgs for magnetic flux, i.e. convert from Weber to Maxwell.
  real(kind=real_kind), parameter, public :: current_si_to_cgs = 1.0e-1_real_kind * C  !< SI to cgs, i.e. A to statA
  real(kind=real_kind), parameter, public :: sqg11_convfac = 1.0e6_real_kind  !< SI to cgs
  real(kind=real_kind), parameter, public :: length_si_to_cgs = 1.0e2_real_kind  !< SI to cgs for length, i.e. convert from m to cm.
  real(kind=real_kind), parameter, public :: te_to_vte     = 4.19e7_real_kind !< Thermal energy to thermal velocity, i.e. eV to cm/s.
  real(kind=real_kind), parameter, public :: ev_to_cgs     = 1.6022e-12_real_kind !< eV to gcm^2/s^2, the cgs equivalent of joule.
  real(kind=real_kind), parameter, public :: pressure_si_to_cgs = 1.0e1_real_kind  !< SI to cgs for pressure, i.e. convert from Pa to Barye.

end module math_constants
