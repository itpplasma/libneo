program test_jorek_field
use, intrinsic :: iso_fortran_env, only: dp => real64
use util_for_test, only: print_test, print_ok, print_fail

implicit none


call test_jorek_field_init


contains


subroutine test_jorek_field_init

    call print_test("test_jorek_field_init")

    call print_ok
end subroutine test_jorek_field_init


end program test_jorek_field