module extract_fluxcoord_mod
  use libneo_kinds, only : real_kind

  implicit none

  integer :: load_extract_fluxcoord=1
  integer :: nphinorm
  real(kind=real_kind) :: psif_extract,theta_extract,psifmin,hpsif
  real(kind=real_kind) :: psifmax,phifmax,sigcos
  real(kind=real_kind), dimension(:), allocatable :: phinorm_arr
end module extract_fluxcoord_mod
