module inthecore_mod
  use libneo_kinds, only : dp

  implicit none

  logical :: prop=.true.
  integer :: npoi,ijumpb,ibeg,iend
  real(dp), parameter :: epssep=1.d-6
  real(dp) :: rc,zc,twopi,sig,psi_sep,psi_cut,sigpsi,cutoff
  real(dp), dimension(:), allocatable :: rho2i,theti
  integer          :: incore
  real(dp) :: vacf,dvacdr,dvacdz,d2vacdr2,d2vacdrdz,d2vacdz2
  real(dp) :: plaf,dpladr,dpladz,d2pladr2,d2pladrdz,d2pladz2
end module inthecore_mod
