program test_jorek_field
use, intrinsic :: iso_fortran_env, only: dp => real64
use neo_jorek_field, only: jorek_field_t
use util_for_test, only: print_test, print_ok, print_fail
use hdf5_tools, only: hid_t, h5_init, h5_open, h5_get, h5_close, h5_deinit
use util_for_test_jorek_field, only: get_filename

implicit none


call test_jorek_field_init


contains


subroutine test_jorek_field_init
    type(jorek_field_t) :: field
    character(len=100) :: filename
    integer(hid_t) :: file_id
    integer :: dimensions(3), n_var
    real(dp), dimension(:,:,:,:), allocatable :: values

    call print_test("test_jorek_field_init")

    call get_filename(filename)

    call field%jorek_field_init(filename)

    call print_ok
end subroutine test_jorek_field_init


end program test_jorek_field