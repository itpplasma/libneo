program test_system_utility
  use system_utility, only : get_rusage
  use rusage_type, only : fortran_rusage, write_fortran_rusage

  implicit none

  type(fortran_rusage) :: usage

  usage = get_rusage()

  call write_fortran_rusage(usage)

end program test_system_utility
