module system_utility
  use, intrinsic :: iso_c_binding

  implicit none

  type, bind(c) :: fortran_timeval
    integer(c_long) :: seconds, microseconds
  end type fortran_timeval

  interface
    integer(c_int) function getrusage(who, usage) bind(c)
      use iso_c_binding
      use rusage_type, only : fortran_rusage
      integer(c_int) :: who
      type(fortran_rusage) :: usage(*)
    end function getrusage
  end interface

  private

  public system_mem_usage, get_rusage

contains

  subroutine system_mem_usage(valueRSS)
    implicit none

!~     use ifport, only : getpid !if on intel compiler

    integer, intent(out) :: valueRSS

    character(len=200):: filename=' '
    character(len=80) :: line
    character(len=8)  :: pid_char=' '
    integer :: pid
    logical :: ifxst

    valueRSS=-1    ! return negative number if not found

    !--- get process ID

!~     pid=getpid()
    write(pid_char,'(I8)') pid
    filename='/proc/'//trim(adjustl(pid_char))//'/status'

    !--- read system file

    inquire (file=filename,exist=ifxst)
    if (.not.ifxst) then
      write (*,*) 'system file does not exist'
      return
    end if

    open(unit=100, file=filename, action='read')
    do
      read (100,'(a)',end=120) line
      if (line(1:6).eq.'VmRSS:') then
         read (line(7:),*) valueRSS
         exit
      end if
    end do
    120 continue
    close(100)

  end subroutine system_mem_usage

  function get_rusage()
    use rusage_type, only : fortran_rusage

    type(fortran_rusage) :: get_rusage

    type(fortran_rusage), allocatable :: rusage(:)
    integer :: res

    allocate(rusage(1))

    res = 0

    res = getrusage(0, rusage) ! \todo find better solution: RUSAGE_SELF should be passed.

    if (res /= 0) then
      write(*,*) 'WARNING: resource usage could not be determined.'
    end if

    get_rusage = rusage(1)
  end function get_rusage

end module system_utility
