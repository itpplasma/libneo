program cleanup_test_jorek_field
use, intrinsic :: iso_fortran_env, only: dp => real64
use util_for_test, only: print_test, print_ok, print_fail

implicit none


call remove_mockup_jorek_output()


contains


subroutine remove_mockup_jorek_output()
    character(len=*), parameter :: filename = 'jorek_output.h5'
    integer :: stat, unit_id
    logical :: exists

    call print_test("remove_mockup_jorek_output")
    
    open(newunit=unit_id, iostat=stat, file=filename, status='old')
    if (stat == 0) close(unit_id, status='delete')
    inquire(file=trim(filename), exist=exists)
    if (exists) then
        call print_fail
        error stop
    else
        call print_ok
    end if
end subroutine remove_mockup_jorek_output

end program cleanup_test_jorek_field