module inthecore_mod
  logical :: prop=.true.
  integer :: npoi,ijumpb,ibeg,iend
  double precision, parameter :: epssep=1.d-6
  double precision :: rc,zc,twopi,sig,psi_sep,psi_cut,sigpsi,cutoff
  double precision, dimension(:), allocatable :: rho2i,theti
  integer          :: incore
  double precision :: vacf,dvacdr,dvacdz,d2vacdr2,d2vacdrdz,d2vacdz2
  double precision :: plaf,dpladr,dpladz,d2pladr2,d2pladrdz,d2pladz2
end module inthecore_mod
