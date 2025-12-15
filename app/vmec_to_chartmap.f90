program vmec_to_chartmap
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_chartmap_vmec_generator, only: write_chartmap_from_vmec
    implicit none

    character(len=1024) :: wout_file
    character(len=1024) :: outfile
    character(len=64) :: arg
    integer :: nrho, ntheta, nzeta
    integer :: ierr
    character(len=2048) :: message

    call get_command_argument(1, wout_file)
    call get_command_argument(2, outfile)
    if (len_trim(wout_file) == 0 .or. len_trim(outfile) == 0) then
        print *, "Usage: vmec_to_chartmap.x <wout.nc> <out.chartmap.nc> "// &
            "[nrho ntheta nzeta]"
        error stop 1
    end if

    nrho = 33
    ntheta = 32
    nzeta = 33

    call get_command_argument(3, arg)
    if (len_trim(arg) > 0) read (arg, *) nrho
    call get_command_argument(4, arg)
    if (len_trim(arg) > 0) read (arg, *) ntheta
    call get_command_argument(5, arg)
    if (len_trim(arg) > 0) read (arg, *) nzeta

    call write_chartmap_from_vmec(trim(wout_file), trim(outfile), nrho, ntheta, &
                                  nzeta, ierr, &
                                  message)
    if (ierr /= 0) then
        print *, "vmec_to_chartmap: failed: ", trim(message)
        error stop 2
    end if
end program vmec_to_chartmap
