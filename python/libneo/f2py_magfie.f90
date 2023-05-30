
subroutine init_vmec(filename, amultharm)
  use new_vmec_stuff_mod, only : netcdffile, multharm
  implicit none

  character(*), intent(in) :: filename
  integer, intent(in) :: amultharm

  netcdffile = filename
  multharm = amultharm

  call spline_vmec_data
end subroutine init_vmec
