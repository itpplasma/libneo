module inthecore_mod
  use libneo_kinds, only : real_kind

  implicit none

  logical :: prop=.true.
  integer :: npoi,ijumpb,ibeg,iend
  real(kind=real_kind), parameter :: epssep=1.d-6
  real(kind=real_kind) :: rc,zc,twopi,sig,psi_sep,psi_cut,sigpsi,cutoff
  real(kind=real_kind), dimension(:), allocatable :: rho2i,theti
  integer          :: incore
  real(kind=real_kind) :: vacf,dvacdr,dvacdz,d2vacdr2,d2vacdrdz,d2vacdz2
  real(kind=real_kind) :: plaf,dpladr,dpladz,d2pladr2,d2pladrdz,d2pladz2
end module inthecore_mod
