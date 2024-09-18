program test_jorek_field
use, intrinsic :: iso_fortran_env, only: dp => real64
use util_for_test, only: print_test, print_ok, print_fail
use hdf5_tools, only: hid_t, h5_init, h5_open, h5_get, h5_close, h5_deinit

implicit none


call test_jorek_field_init


contains


subroutine test_jorek_field_init
    character(len=*), parameter :: filename = 'jorek_output.h5'
    integer(hid_t) :: file_id
    integer :: dimensions(3), n_var
    real(dp), dimension(:,:,:,:), allocatable :: values

    call print_test("test_jorek_field_init")

    call h5_init()
    call h5_open(filename, file_id)
    call h5_get(file_id, 'dim', dimensions)
    call h5_get(file_id, 'n_var', n_var)
    allocate(values(dimensions(1), dimensions(2), dimensions(3), n_var))
    call h5_get(file_id, 'values', values)
    call h5_close(file_id)
    call h5_deinit()

    call print_ok
end subroutine test_jorek_field_init


end program test_jorek_field