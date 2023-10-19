module theta_rz_mod
  use libneo_kinds, only : real_kind

  implicit none

  integer :: icall=0
  integer :: nsqp,nlab,nthe,icp_pt
  integer, dimension(:,:), allocatable :: ipoint_pt
  real(kind=real_kind) :: hsqpsi,hlabel,htheqt,psiaxis,sigma_qt,raxis,zaxis
  real(kind=real_kind), dimension(:,:),   allocatable :: spllabel
  real(kind=real_kind), dimension(:,:,:), allocatable :: splthet
  real(kind=real_kind), dimension(:),     allocatable :: sqpsi,flab,theqt
end module theta_rz_mod
