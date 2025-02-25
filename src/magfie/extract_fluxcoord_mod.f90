module extract_fluxcoord_mod
  use libneo_kinds, only : dp

  implicit none

  integer :: load_extract_fluxcoord=1
  integer :: nphinorm
  real(dp) :: psif_extract,theta_extract,psifmin,hpsif
  real(dp) :: psifmax,phifmax,sigcos
  real(dp), dimension(:), allocatable :: phinorm_arr
end module extract_fluxcoord_mod
