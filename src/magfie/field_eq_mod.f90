module field_eq_mod
  implicit none

  integer, parameter :: dp = kind(1.0d0)

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
end module field_eq_mod
