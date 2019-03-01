!> \brief Definitions of kind parameters for use in the library and codes that use the library.
!>
!> The module defines kind parameters, that should be used in the libneo
!> library and can be used at code interfacing with libneo.
!>
!> In most cases the *_kind parameters are the ones that should be used.
module libneo_kinds
  ! real kinds
  integer, parameter :: rp = kind(1.0)
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: real_kind=dp

  ! complex kinds
  integer, parameter :: crp = kind( (1.0, 1.0) )
  integer, parameter :: cdp = kind( (1.0d0, 1.0d0) )
  integer, parameter :: complex_kind=cdp
end module libneo_kinds
