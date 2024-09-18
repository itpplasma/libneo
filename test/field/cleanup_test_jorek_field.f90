program cleanup_test_jorek_field
use, intrinsic :: iso_fortran_env, only: dp => real64
use util_for_test, only: print_test, print_ok, print_fail
!use hdf5_tools, only: h5_init, h5_deinit

implicit none


call remove_straight_field_jorek_file()


contains


subroutine remove_straight_field_jorek_file()
    call print_test("remove_straight_field_jorek_file")
    
    !call h5_init()
    !call h5_deinit()
    
    call print_ok
end subroutine remove_straight_field_jorek_file

end program cleanup_test_jorek_field