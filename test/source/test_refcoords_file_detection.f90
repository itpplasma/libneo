program test_refcoords_file_detection
    use libneo_coordinates, only: detect_refcoords_file_type, &
                                  refcoords_file_chartmap, refcoords_file_vmec_wout
    implicit none

    integer :: ierr
    integer :: file_type
    character(len=2048) :: message
    integer :: nerrors

    nerrors = 0

    call detect_refcoords_file_type("chartmap.nc", file_type, ierr, message)
    if (ierr /= 0) then
        print *, "  FAIL: detect_refcoords_file_type(chartmap.nc): ", trim(message)
        nerrors = nerrors + 1
    else if (file_type /= refcoords_file_chartmap) then
        print *, "  FAIL: expected CHARTMAP file_type for chartmap.nc, got ", file_type
        nerrors = nerrors + 1
    else
        print *, "  PASS: chartmap.nc detected as CHARTMAP"
    end if

    call detect_refcoords_file_type("wout.nc", file_type, ierr, message)
    if (ierr /= 0) then
        print *, "  FAIL: detect_refcoords_file_type(wout.nc): ", trim(message)
        nerrors = nerrors + 1
    else if (file_type /= refcoords_file_vmec_wout) then
        print *, "  FAIL: expected VMEC_WOUT file_type for wout.nc, got ", file_type
        nerrors = nerrors + 1
    else
        print *, "  PASS: wout.nc detected as VMEC_WOUT"
    end if

    if (nerrors > 0) then
        print *, "FAILED: ", nerrors, " error(s) detected in refcoords file detection"
        error stop 1
    end if

end program test_refcoords_file_detection

