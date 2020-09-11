module bdivfree_mod
  integer :: nr,nz,ntor,icp
  integer, dimension(:,:), allocatable :: ipoint
  double precision :: rmin,zmin,hr,hz,pmin,pfac
  double precision, dimension(:),       allocatable :: rpoi,zpoi
  double precision, dimension(:,:,:),   allocatable :: apav,rbpav_coef
  double precision, dimension(:,:,:,:), allocatable :: aznre,aznim,arnre,arnim
end module bdivfree_mod
