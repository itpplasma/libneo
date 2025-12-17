program test_vmec_chartmap_generator
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use chartmap_test_utils, only: chartmap_roundtrip_check
    use libneo_chartmap_vmec_generator, only: write_chartmap_from_vmec
    use libneo_coordinates, only: validate_chartmap_file
    implicit none

    integer :: ierr
    integer :: nerrors
    character(len=2048) :: message
    character(len=*), parameter :: wout_file = "wout.nc"
    character(len=*), parameter :: chartmap_file = "wout_vmec.chartmap.nc"

    nerrors = 0

    call write_chartmap_from_vmec(wout_file, chartmap_file, 32, 33, 34, ierr, message)
    if (ierr /= 0) then
        print *, "  FAIL: write_chartmap_from_vmec: ", trim(message)
        nerrors = nerrors + 1
    end if

    call validate_chartmap_file(chartmap_file, ierr, message)
    if (ierr /= 0) then
        print *, "  FAIL: validate_chartmap_file: ", trim(message)
        nerrors = nerrors + 1
    else
        print *, "  PASS: generated chartmap validated"
    end if

    call chartmap_roundtrip_check(chartmap_file, 1.0e-2_dp, nerrors)

    if (nerrors > 0) then
        print *, "FAILED: ", nerrors, " error(s) in VMEC->chartmap generator test"
        error stop 1
    end if

end program test_vmec_chartmap_generator
