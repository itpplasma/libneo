program test_libneo_field
use, intrinsic :: iso_fortran_env, only: dp => real64
use test_util, only: print_test, print_ok, print_fail

implicit none


call test_create_field
call test_create_example_field
call test_create_biotsavart_field


contains


subroutine test_create_field
    use libneo_field, only: field_t, example_field_t, create_field

    class(field_t), allocatable :: field

    call print_test("test_create_field")

    field = create_field("example")
    if (.not.allocated(field)) then
        call print_fail
        return
    end if

    call print_ok
end subroutine test_create_field


subroutine test_create_example_field
    use libneo_field, only: field_t, example_field_t, create_example_field

    class(field_t), allocatable :: field

    call print_test("test_create_example_field")

    field = create_example_field()
    if (.not.allocated(field)) then
        call print_fail
        return
    end if

    call print_ok
end subroutine test_create_example_field


subroutine test_create_biotsavart_field
    use libneo_field, only: field_t, create_field, create_biotsavart_field

    class(field_t), allocatable :: field

    call print_test("test_create_biotsavart_field")

    field = create_biotsavart_field()
    if (.not.allocated(field)) then
        call print_fail
        return
    end if

    call print_ok
end subroutine test_create_biotsavart_field
    

end program test_libneo_field

