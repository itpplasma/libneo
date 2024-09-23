module util_for_test_jorek_field
use, intrinsic :: iso_fortran_env, only: dp => real64

implicit none

character(len=*), parameter :: save_location = 'jorek_filename.txt'

contains

subroutine save_filename(filename)
    character(len=*), intent(in) :: filename

    integer :: unit_id

    open(newunit=unit_id, file=save_location, status='new')
    write(unit_id, '(A)') filename
    close(unit_id)
end subroutine save_filename


subroutine get_filename(filename)
    character(len=*), intent(out) :: filename
    
    integer :: unit_id

    open(newunit=unit_id, file=save_location, status='old')
    read(unit_id, '(A)') filename
    close(unit_id)
end subroutine get_filename


subroutine remove_saved_filename()
    logical :: exists
    integer :: unit_id, stat

    inquire(file=save_location, exist=exists)
    if (exists) open(newunit=unit_id, iostat=stat, file=save_location, status='old')
    if (stat == 0) close(unit_id, status='delete')
end subroutine remove_saved_filename    
    
end module util_for_test_jorek_field
    