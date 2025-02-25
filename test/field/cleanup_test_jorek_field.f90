program cleanup_test_jorek_field
use util_for_test, only: print_test, print_ok, print_fail
use util_for_test_jorek_field, only: get_filename, remove_saved_filename

implicit none


call remove_mockup_jorek_output()


contains


subroutine remove_mockup_jorek_output()
    character(len=256) :: filename
    integer :: stat, unit_id
    logical :: exists

    call print_test("remove_mockup_jorek_output")

    call get_filename(filename)
    open(newunit=unit_id, iostat=stat, file=filename, status='old')
    if (stat == 0) close(unit_id, status='delete')
    inquire(file=trim(filename), exist=exists)
    if (exists) then
        call print_fail
        error stop
    else
        call print_ok
    end if
    call remove_saved_filename()
end subroutine remove_mockup_jorek_output

end program cleanup_test_jorek_field
