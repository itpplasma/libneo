module theta_rz_mod
  use libneo_kinds, only : dp

  implicit none

  integer :: icall=0
  integer :: nsqp,nlab,nthe,icp_pt
  integer, dimension(:,:), allocatable :: ipoint_pt
  real(dp) :: hsqpsi,hlabel,htheqt,psiaxis,sigma_qt,raxis,zaxis
  real(dp), dimension(:,:),   allocatable :: spllabel
  real(dp), dimension(:,:,:), allocatable :: splthet
  real(dp), dimension(:),     allocatable :: sqpsi,flab,theqt
end module theta_rz_mod
