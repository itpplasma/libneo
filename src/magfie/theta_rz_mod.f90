module theta_rz_mod
  integer :: icall=0
  integer :: nsqp,nlab,nthe,icp_pt
  integer, dimension(:,:), allocatable :: ipoint_pt
  real(kind=8) :: hsqpsi,hlabel,htheqt,psiaxis,sigma_qt,raxis,zaxis
  real(kind=8), dimension(:,:),   allocatable :: spllabel
  real(kind=8), dimension(:,:,:), allocatable :: splthet
  real(kind=8), dimension(:),     allocatable :: sqpsi,flab,theqt
end module theta_rz_mod
