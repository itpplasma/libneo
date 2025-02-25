program test_biotsavart_field
use libneo_kinds, only : dp
use util_for_test, only: print_test, print_ok, print_fail

implicit none


call test_biotsavart_field_init
call test_curla_equal_b


contains


subroutine test_biotsavart_field_init
    use neo_biotsavart_field, only: biotsavart_field_t

    real(dp), parameter :: tol = 1.0e-9_dp

    type(biotsavart_field_t) :: biotsavart_field
    real(dp) :: x(3), A(3), B(3)

    call print_test("test_biotsavart_field_init")

    call biotsavart_field%biotsavart_field_init()

    x = [1.0_dp, 1.0_dp, 1.0_dp]
    call biotsavart_field%compute_abfield(x, A, B)
    if (maxval(abs(A)) > tol) then
        print *, "A != 0 of default biotsavart field"
        call print_fail
        error stop
    end if
    if (maxval(abs(B)) > tol) then
        print *, "B != 0 of default biotsavart field"
        call print_fail
        error stop
    end if

    call print_ok
end subroutine test_biotsavart_field_init


subroutine test_curla_equal_b
    use neo_biotsavart_field, only: biotsavart_field_t
    use util_for_test_field, only: compute_cartesian_curla

    real(dp), parameter :: tol = 1.0e-9_dp

    character(*), parameter :: coils_file = "test_coils.txt"
    type(biotsavart_field_t) :: biotsavart_field
    real(dp) :: x(3), B(3), B_from_A(3)

    call print_test("test_curla_equal_b")

    call create_straight_wire_coils_file(coils_file)
    call biotsavart_field%biotsavart_field_init(coils_file)
    call remove_test_coils_file(coils_file)

    x = [1.0_dp, 1.0_dp, 1.0_dp]

    call biotsavart_field%compute_bfield(x, B)
    B_from_A = compute_cartesian_curla(biotsavart_field, x, tol)

    if (maxval(abs(B - B_from_A)) > tol) then
        print *, "curl A != B"
        call print_fail
        error stop
    end if

    call print_ok
end subroutine test_curla_equal_b


subroutine create_straight_wire_coils_file(filename)
    use neo_biotsavart, only: coils_t, save_coils_to_file

    character(len=*), intent(in) :: filename

    type(coils_t) :: coils

    call init_straight_wire_coils(coils)
    call save_coils_to_file(filename, coils)
end subroutine create_straight_wire_coils_file


subroutine init_straight_wire_coils(coils)
    use neo_biotsavart, only: coils_t, coils_init

    type(coils_t), intent(out) :: coils

    real(dp), dimension(4) :: x, y, z, current
    real(dp), parameter    :: L = 1.0e3_dp

    x = [0.0d0, 0.0d0, 0.0d0, 0.0d0]
    y = [0.0d0, 0.0d0, 0.0d0, 0.0d0]
    z = [-L/2.0d0, -L/2.5d0, 0.0d0, L/2.0d0]
    current = [1.0d0, 1.0d0, 1.0d0, 0.0d0]

    call coils_init(x, y, z, current, coils)
end subroutine init_straight_wire_coils


subroutine remove_test_coils_file(filename)
    character(len=*), intent(in) :: filename

    call system("rm -f " // filename)
end subroutine remove_test_coils_file


end program test_biotsavart_field
