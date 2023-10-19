module neo_magfie_mod
  use libneo_kinds, only : real_kind

  implicit none

  integer :: magfie_spline
  real(kind=real_kind), dimension(:), allocatable :: magfie_sarray
end module neo_magfie_mod
