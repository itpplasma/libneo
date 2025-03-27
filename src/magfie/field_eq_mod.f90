module field_eq_mod
  use libneo_kinds, only : dp

  implicit none

  logical :: use_fpol = .true.                                      !<=18.12.18
  logical :: skip_read = .false.
  integer :: icall_eq=0
  integer :: nrad,nzet,icp,nwindow_r,nwindow_z
  integer :: ierrfield

  real(dp) :: psib,btf,rtf,hrad,hzet
  real(dp) :: psi_axis,psi_sep,hfpol                            !<=18.12.18
  real(dp), dimension(:,:), allocatable    :: psi, psi0
  real(dp), dimension(:,:), allocatable    :: splfpol           !<=18.12.18
  real(dp), dimension(:,:,:), allocatable  :: splpsi
  real(dp), dimension(:), allocatable      :: rad, zet, xi,f
  integer, dimension(:), allocatable           :: imi,ima,jmi,jma
  integer, dimension(:,:), allocatable         :: ipoint
  real(dp) :: psif,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2

  ! Make temporary variables threadprivate
  !$omp threadprivate(psif,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2)
end module field_eq_mod
