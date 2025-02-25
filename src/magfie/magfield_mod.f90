module magfield_mod
  use libneo_kinds, only : dp

  implicit none

  integer :: ierrfield
  integer :: input_format,npmid,nr,np,nz,npoint
  integer,          dimension(:,:,:), allocatable :: ipoint
  real(dp), dimension(:),     allocatable :: rad,phi,zet,Brs,Bzs,Bps
  real(dp), dimension(:,:,:), allocatable :: Bx,By,Bz,Br,Bp
end module
