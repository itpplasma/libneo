!> \brief Definitions of kind parameters for use in the library and codes that use the library.
!>
!> The module defines kind parameters, that should be used in the libneo
!> library and can be used at code interfacing with libneo.
!>
!> In most cases the *_kind parameters are the ones that should be used.
!>
!> \todo Name of the interger kinds? I think the b comes from byte but
!>   does the byte count need to match the indicated value?
module libneo_kinds
  ! integer kinds
  integer, parameter :: i4b = selected_int_kind(9)
  integer, parameter :: i2b = selected_int_kind(4)
  integer, parameter :: i1b = selected_int_kind(2)
  integer, parameter :: integer_kind=i4b

  ! real kinds
  integer, parameter :: rp = kind(1.0)
  integer, parameter :: dp = kind(1.0d0)

  ! complex kinds
  integer, parameter :: crp = kind( (1.0, 1.0) )
  integer, parameter :: cdp = kind( (1.0d0, 1.0d0) )
  integer, parameter :: complex_kind=cdp
end module libneo_kinds
