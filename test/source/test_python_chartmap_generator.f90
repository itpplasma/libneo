program test_python_chartmap_generator
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use chartmap_test_utils, only: chartmap_roundtrip_check
    use libneo_coordinates, only: validate_chartmap_file
    implicit none

    integer :: ierr
    integer :: nerrors
    character(len=2048) :: message
    character(len=*), parameter :: chartmap_file = "wout_vmec_python.chartmap.nc"

    nerrors = 0

    call validate_chartmap_file(chartmap_file, ierr, message)
    if (ierr /= 0) then
        print *, "  FAIL: validate_chartmap_file: ", trim(message)
        nerrors = nerrors + 1
    else
        print *, "  PASS: python-generated chartmap validated"
    end if

    call chartmap_roundtrip_check(chartmap_file, 5.0e-2_dp, nerrors)

    if (nerrors > 0) then
        print *, "FAILED: ", nerrors, " error(s) in python chartmap generator test"
        error stop 1
    end if

end program test_python_chartmap_generator
