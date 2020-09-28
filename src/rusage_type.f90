!> Not in system utitlity as it is needed in the module itself and in
!> the interface for the c-function getrusage.
module rusage_type
  use, intrinsic :: iso_c_binding

  implicit none

  type, bind(c) :: fortran_rusage
    integer(c_long) :: ru_utimes !< user CPU time used, seconds
    integer(c_long) :: ru_utimems !< user CPU time used, microseconds
    integer(c_long) :: ru_stimes !< system CPU time used, seconds
    integer(c_long) :: ru_stimems !< system CPU time used, microseconds
    integer(c_long) :: ru_maxrss !< maximum resident set size
    integer(c_long) :: ru_ixrss !< integral shared memory size
    integer(c_long) :: ru_idrss !< integral unshared data size
    integer(c_long) :: ru_isrss !< integral unshared stack size
    integer(c_long) :: ru_minflt !< page reclaims (soft page faults)
    integer(c_long) :: ru_majflt !< page faults (hard page faults)
    integer(c_long) :: ru_nswap !< swaps
    integer(c_long) :: ru_inblock !< block input operations
    integer(c_long) :: ru_oublock !< block output operations
    integer(c_long) :: ru_msgsnd !< IPC messages sent
    integer(c_long) :: ru_msgrcv !< IPC messages received
    integer(c_long) :: ru_nsignals !< signals received
    integer(c_long) :: ru_nvcsw !< voluntary context switches
    integer(c_long) :: ru_nivcsw !< involuntary context switches
  end type fortran_rusage

contains

  subroutine write_fortran_rusage(usage)
    type(fortran_rusage) :: usage

    write(*,*) 'user CPU time used: ', usage%ru_utimes, 's ', usage%ru_utimems, ' us'
    write(*,*) 'system CPU time used: ', usage%ru_stimes, 's ', usage%ru_stimems, ' us'
    write(*,*) 'maximum resident set size: ', usage%ru_maxrss
    write(*,*) 'integral shared memory size: ', usage%ru_ixrss
    write(*,*) 'integral unshared data size: ', usage%ru_idrss
    write(*,*) 'integral unshared stack size: ', usage%ru_isrss
    write(*,*) 'page reclaims (soft page faults): ', usage%ru_minflt
    write(*,*) 'page faults (hard page faults): ', usage%ru_majflt
    write(*,*) 'swaps: ', usage%ru_nswap
    write(*,*) 'block input operations: ', usage%ru_inblock
    write(*,*) 'block output operations: ', usage%ru_oublock
    write(*,*) 'IPC messages sent: ', usage%ru_msgsnd
    write(*,*) 'IPC messages received: ', usage%ru_msgrcv
    write(*,*) 'signals received: ', usage%ru_nsignals
    write(*,*) 'voluntary context switches: ', usage%ru_nvcsw
    write(*,*) 'involuntary context switches: ', usage%ru_nivcsw
  end subroutine write_fortran_rusage

  subroutine write_fortran_rusage_short(usage)
    type(fortran_rusage) :: usage

    write(*,*) 'user CPU time used: ', usage%ru_utimes, 's'
    write(*,*) 'system CPU time used: ', usage%ru_stimes, 's'
    write(*,*) 'maximum resident set size: ', usage%ru_maxrss
  end subroutine write_fortran_rusage_short

end module rusage_type
