program test_field
use, intrinsic :: iso_fortran_env, only: dp => real64
use util_for_test, only: print_test, print_ok, print_fail

implicit none


call test_create_field
call test_create_example_field
call test_create_biotsavart_field
call test_create_polylag_field
call test_create_spline_field
call test_create_spline_field_from_mesh
call test_create_jorek_field


contains


subroutine test_create_field
    use neo_field, only: field_t, example_field_t, create_field

    class(field_t), allocatable :: field

    call print_test("test_create_field")

    call create_field(field, "example")
    if (.not.allocated(field)) then
        call print_fail
        error stop
    end if

    call print_ok
end subroutine test_create_field


subroutine test_create_example_field
    use neo_field, only: field_t, example_field_t, create_example_field

    class(field_t), allocatable :: field

    call print_test("test_create_example_field")

    allocate(field, source=create_example_field())
    if (.not.allocated(field)) then
        call print_fail
        error stop
    end if

    call print_ok
end subroutine test_create_example_field


subroutine test_create_biotsavart_field
    use neo_field, only: field_t, create_field, create_biotsavart_field

    class(field_t), allocatable :: field

    call print_test("test_create_biotsavart_field")

    allocate(field, source=create_biotsavart_field())
    if (.not.allocated(field)) then
        call print_fail
        error stop
    end if

    call print_ok
end subroutine test_create_biotsavart_field


subroutine test_create_polylag_field
    use neo_field, only: field_t, create_polylag_field

    class(field_t), allocatable :: field
    real(dp), dimension(3,2) :: limits

    call print_test("test_create_polylag_field")

    limits(1,:) = [1.0_dp, 10.0_dp]
    limits(2,:) = [1.0_dp, 10.0_dp]
    limits(3,:) = [1.0_dp, 10.0_dp]

    allocate(field, source=create_polylag_field(limits))
    if (.not.allocated(field)) then
        call print_fail
        error stop
    end if

    call print_ok
end subroutine test_create_polylag_field


subroutine test_create_spline_field
    use neo_field, only: field_t, create_spline_field

    class(field_t), allocatable :: field
    real(dp), dimension(3,2) :: limits

    call print_test("test_create_spline_field")

    limits(1,:) = [1.0_dp, 10.0_dp]
    limits(2,:) = [1.0_dp, 10.0_dp]
    limits(3,:) = [1.0_dp, 10.0_dp]

    allocate(field, source=create_spline_field(limits))
    if (.not.allocated(field)) then
        call print_fail
        error stop
    end if

    call print_ok
end subroutine test_create_spline_field


subroutine test_create_spline_field_from_mesh
    use neo_field, only: field_t, create_spline_field_from_mesh
    use neo_field, only: field_mesh_t

    class(field_t), allocatable :: field
    type(field_mesh_t) :: field_mesh
    real(dp), dimension(3,2) :: limits

    call print_test("test_create_spline_field_from_mesh")

    limits(1,:) = [1.0_dp, 10.0_dp]
    limits(2,:) = [1.0_dp, 10.0_dp]
    limits(3,:) = [1.0_dp, 10.0_dp]

    call field_mesh%field_mesh_init_with_field(limits)
    allocate(field, source=create_spline_field_from_mesh(field_mesh))
    if (.not.allocated(field)) then
        call print_fail
        error stop
    end if

    call print_ok
end subroutine test_create_spline_field_from_mesh


subroutine test_create_jorek_field
    use neo_field, only: field_t, create_field
    use util_for_test_jorek_field, only: get_filename

    class(field_t), allocatable :: field
    character(len=256) :: filename

    call print_test("test_create_jorek_field")

    call get_filename(filename)
    call create_field(field, "jorek", filename=filename)
    if (.not.allocated(field)) then
        call print_fail
        error stop
    end if

    call print_ok
end subroutine test_create_jorek_field
    

end program test_field

