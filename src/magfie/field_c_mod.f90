module field_c_mod
  use libneo_kinds, only : real_kind

  implicit none

  integer :: icall_c=0
  integer :: ntor=16
  integer :: nr,np,nz
  integer :: icftype
  real(kind=real_kind) :: rmin,pmin,zmin,rmax,pmax,zmax
end module field_c_mod
