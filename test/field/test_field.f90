program test_field
use libneo_kinds, only : dp
use util_for_test, only: print_test, print_ok, print_fail
use neo_field, only: field_t, create_field

implicit none


call test_create_example_field
call test_create_biotsavart_field
call test_create_polylag_field
call test_create_spline_field
call test_create_spline_field_from_mesh
call test_create_jorek_field


contains


subroutine test_create_example_field
    class(field_t), allocatable :: field

    call print_test("test_create_example_field")

    call create_field(field, "example")
    if (.not.allocated(field)) then
        call print_fail
        error stop
    end if

    call print_ok
end subroutine test_create_example_field


subroutine test_create_biotsavart_field
    class(field_t), allocatable :: field

    call print_test("test_create_biotsavart_field")

    call create_field(field, "biotsavart")
    if (.not.allocated(field)) then
        call print_fail
        error stop
    end if

    call print_ok
end subroutine test_create_biotsavart_field


subroutine test_create_polylag_field
    class(field_t), allocatable :: field
    real(dp), dimension(3,2) :: limits

    call print_test("test_create_polylag_field")

    limits(1,:) = [1.0_dp, 10.0_dp]
    limits(2,:) = [1.0_dp, 10.0_dp]
    limits(3,:) = [1.0_dp, 10.0_dp]

    call create_field(field, "polylag", limits=limits)
    if (.not.allocated(field)) then
        call print_fail
        error stop
    end if

    call print_ok
end subroutine test_create_polylag_field


subroutine test_create_spline_field
    class(field_t), allocatable :: field
    real(dp), dimension(3,2) :: limits

    call print_test("test_create_spline_field")

    limits(1,:) = [1.0_dp, 10.0_dp]
    limits(2,:) = [1.0_dp, 10.0_dp]
    limits(3,:) = [1.0_dp, 10.0_dp]

    call create_field(field, "spline", limits=limits)
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
    call create_field(field, "spline", field_mesh=field_mesh)
    if (.not.allocated(field)) then
        call print_fail
        error stop
    end if

    call print_ok
end subroutine test_create_spline_field_from_mesh


subroutine test_create_jorek_field
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
