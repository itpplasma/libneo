module field_eq_mod
  use libneo_kinds, only : real_kind

  implicit none

  logical :: use_fpol = .true.                                      !<=18.12.18
  logical :: skip_read = .false.
  integer :: icall_eq=0
  integer :: nrad,nzet,icp,nwindow_r,nwindow_z
  integer :: ierrfield

  real(kind=real_kind) :: psib,btf,rtf,hrad,hzet
  real(kind=real_kind) :: psi_axis,psi_sep,hfpol                            !<=18.12.18
  real(kind=real_kind), dimension(:,:), allocatable    :: psi, psi0
  real(kind=real_kind), dimension(:,:), allocatable    :: splfpol           !<=18.12.18
  real(kind=real_kind), dimension(:,:,:), allocatable  :: splpsi
  real(kind=real_kind), dimension(:), allocatable      :: rad, zet, xi,f
  integer, dimension(:), allocatable           :: imi,ima,jmi,jma
  integer, dimension(:,:), allocatable         :: ipoint
  real(kind=real_kind) :: psif,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2

  ! Make temporary variables threadprivate
  !$omp threadprivate(psif)
end module field_eq_mod
