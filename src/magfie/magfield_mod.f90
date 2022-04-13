module magfield_mod
  integer :: ierrfield
  integer :: input_format,npmid,nr,np,nz,npoint
  integer,          dimension(:,:,:), allocatable :: ipoint
  double precision, dimension(:),     allocatable :: rad,phi,zet,Brs,Bzs,Bps
  double precision, dimension(:,:,:), allocatable :: Bx,By,Bz,Br,Bp
end module
