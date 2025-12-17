program test_coils_chartmap_generator
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use chartmap_test_utils, only: chartmap_roundtrip_check
    use libneo_coordinates, only: validate_chartmap_file
    implicit none

    integer :: ierr
    integer :: nerrors
    character(len=2048) :: message
    character(len=*), parameter :: chartmap_file = "coils_offset.chartmap.nc"

    nerrors = 0

    call validate_chartmap_file(chartmap_file, ierr, message)
    if (ierr /= 0) then
        print *, "  FAIL: validate_chartmap_file: ", trim(message)
        nerrors = nerrors + 1
    else
        print *, "  PASS: coils-generated chartmap validated"
    end if

    call chartmap_roundtrip_check(chartmap_file, 1.0e-1_dp, nerrors)

    if (nerrors > 0) then
        print *, "FAILED: ", nerrors, " error(s) in coils->chartmap generator test"
        error stop 1
    end if
end program test_coils_chartmap_generator

