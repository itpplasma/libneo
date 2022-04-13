!
  module parmot_mod
    double precision :: rmu,ro0
  end module parmot_mod
!
  module collis_alp
    integer, parameter :: nsorts=3, ns=10000 !original: ns=10000
    integer :: iswmod
    logical :: swcoll=.false.
    double precision, dimension(nsorts)    :: efcolf,velrat,enrat
    double precision, dimension(nsorts,ns) :: efcolf_arr,velrat_arr,enrat_arr
  end module collis_alp
!
  module elefie_mod
    double precision :: rbig
    double precision, dimension(:,:), allocatable :: Mtprofile
    double precision, dimension(:,:), allocatable :: plasma
    double precision :: amb,am1,am2,Zb,Z1,Z2,densi1,densi2,tempi1,tempi2,tempe,v0
    double precision :: escale=1.0,bscale=1.0
  end module elefie_mod

  module constants
    double precision, parameter :: pi=3.14159265358979d0
    double precision,parameter  :: c=2.9979d10
    double precision,parameter  :: e_charge=4.8032d-10
    double precision,parameter  :: e_mass=9.1094d-28
    double precision,parameter  :: p_mass=1.6726d-24
    double precision,parameter  :: ev=1.6022d-12
  end module constants
