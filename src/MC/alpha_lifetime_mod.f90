
module parmot_mod
  use libneo_kinds, only : real_kind

  implicit none

  real(kind=real_kind) :: rmu,ro0
end module parmot_mod

module collis_alp
  use libneo_kinds, only : real_kind

  implicit none

  integer, parameter :: nsorts=3, ns=10000 !original: ns=10000
  integer :: iswmod
  logical :: swcoll=.false.
  real(kind=real_kind), dimension(nsorts)    :: efcolf,velrat,enrat
  real(kind=real_kind), dimension(nsorts,ns) :: efcolf_arr,velrat_arr,enrat_arr
end module collis_alp

module elefie_mod
  use libneo_kinds, only : real_kind

  implicit none

  real(kind=real_kind) :: rbig
  real(kind=real_kind), dimension(:,:), allocatable :: Mtprofile
  real(kind=real_kind), dimension(:,:), allocatable :: plasma
  real(kind=real_kind) :: amb,am1,am2,Zb,Z1,Z2,densi1,densi2,tempi1,tempi2,tempe,v0
  real(kind=real_kind) :: escale=1.0,bscale=1.0
end module elefie_mod

module constants
  use libneo_kinds, only : real_kind

  implicit none

  real(kind=real_kind), parameter :: pi=3.14159265358979d0
  real(kind=real_kind),parameter  :: c=2.9979d10
  real(kind=real_kind),parameter  :: e_charge=4.8032d-10
  real(kind=real_kind),parameter  :: e_mass=9.1094d-28
  real(kind=real_kind),parameter  :: p_mass=1.6726d-24
  real(kind=real_kind),parameter  :: ev=1.6022d-12
end module constants
