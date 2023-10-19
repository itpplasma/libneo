module magfield_mod
  use libneo_kinds, only : real_kind

  implicit none

  integer :: ierrfield
  integer :: input_format,npmid,nr,np,nz,npoint
  integer,          dimension(:,:,:), allocatable :: ipoint
  real(kind=real_kind), dimension(:),     allocatable :: rad,phi,zet,Brs,Bzs,Bps
  real(kind=real_kind), dimension(:,:,:), allocatable :: Bx,By,Bz,Br,Bp
end module
