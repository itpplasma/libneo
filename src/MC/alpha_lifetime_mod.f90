
module parmot_mod
  use libneo_kinds, only : dp

  implicit none

  real(dp) :: rmu,ro0
end module parmot_mod

module collis_alp
  use libneo_kinds, only : dp

  implicit none

  integer, parameter :: nsorts=3, ns=10000 !original: ns=10000
  integer :: iswmod
  logical :: swcoll=.false.
  real(dp), dimension(nsorts)    :: efcolf,velrat,enrat
  real(dp), dimension(nsorts,ns) :: efcolf_arr,velrat_arr,enrat_arr
end module collis_alp

module elefie_mod
  use libneo_kinds, only : dp

  implicit none

  real(dp) :: rbig
  real(dp), dimension(:,:), allocatable :: Mtprofile
  real(dp), dimension(:,:), allocatable :: plasma
  real(dp) :: amb,am1,am2,Zb,Z1,Z2,densi1,densi2,tempi1,tempi2,tempe,v0
  real(dp) :: escale=1.0,bscale=1.0
end module elefie_mod

module constants
  use libneo_kinds, only : dp

  implicit none

  real(dp), parameter :: pi=3.14159265358979d0
  real(dp),parameter  :: c=2.9979d10
  real(dp),parameter  :: e_charge=4.8032d-10
  real(dp),parameter  :: e_mass=9.1094d-28
  real(dp),parameter  :: p_mass=1.6726d-24
  real(dp),parameter  :: ev=1.6022d-12
end module constants
