program test_spline_field
use, intrinsic :: iso_fortran_env, only: dp => real64
use util_for_test, only: print_test, print_ok, print_fail
use neo_spline_field, only: spline_field_t
implicit none

call test_spline_field_from_file_init
call test_spline_field_init
call test_against_original_field
call test_curla_equal_b
call test_compute_afield_derivatives

contains

subroutine test_spline_field_from_file_init
    use neo_field, only: create_field
    use neo_field_base, only: field_t

    class(field_t), allocatable :: field
    character(*), parameter :: filename = "field_from_gorilla"
    real(dp), dimension(2) :: deviation
    real(dp), parameter :: tol = 1.0e-7_dp

    call print_test("test_spline_field_from_file_init")

    call create_file(filename)
    call create_field(field, "spline", filename=filename)

    deviation = compare_spline_with_file(field,filename)
    call unlink(filename)

    if (deviation(1).gt.tol .or. deviation(2).gt.tol) then
        print *, "error: tol = ", tol, " , deviation in A-field is = ", deviation(1), &
                                       " , deviation in B-field is = ", deviation(2)
        call print_fail
        error stop
    end if
    
    call print_ok
end subroutine test_spline_field_from_file_init

subroutine create_file(filename)
    use neo_example_field, only: example_field_t
    use neo_field_mesh, only: linspace
    use math_constants, only: pi

    character(*), intent(in) :: filename
    type(example_field_t) :: example_field
    real(dp) :: x(3), A(3), B(3)
    real(dp), dimension(:), allocatable :: x1, x2, x3
    real(dp), dimension(3,2) :: limits !min/max for the three coordinates
    integer, dimension(3) :: n_points !number of vetices per direction
    logical :: is_periodic(3) !describes periodicity in three dimensions
    integer :: unit, i, j, k

    limits(1,:) = (/100.0_dp,200.0_dp/) !R
    limits(2,:) = (/0.0_dp,2*pi/) !phi
    limits(3,:) = (/-50.0_dp,50.0_dp/) !Z
    n_points = 30
    is_periodic = (/.false.,.true.,.false./)

    allocate(x1(n_points(1)), x2(n_points(2)), x3(n_points(3)))
    x1 = linspace(limits(1,1), limits(1,2), n_points(1))
    x2 = linspace(limits(2,1), limits(2,2), n_points(2))
    x3 = linspace(limits(3,1), limits(3,2), n_points(3))

    call example_field%example_field_init(1.0_dp, 2.0_dp)

    call unlink(filename)
    open(newunit = unit, file = filename)
    write(unit,*) n_points
    write(unit,*) limits
    write(unit,*) is_periodic 

    do i = 1, size(x1)
        do j = 1, size(x2)
            do k = 1, size(x3)
                x = [x1(i), x2(j), x3(k)]
                call example_field%compute_abfield(x, A, B)
                write(unit,*) A, B
            end do
        end do
    end do

    close(unit)

end subroutine create_file

function compare_spline_with_file(spline_field, filename) result(deviation)           
    use neo_field_base, only: field_t
    use neo_spline_field, only: spline_field_t

    character(*), intent(in) :: filename
    class(field_t), intent(in) :: spline_field
    real(dp), dimension(2) :: deviation
    integer :: unit
    real(dp), dimension(3) :: A, A_spline, B, B_spline, x
    real(dp), dimension(3,2) :: limits

    open(newunit = unit, file = filename)
    read(unit,*)
    read(unit,*) limits
    read(unit,*)
    read(unit,*) A, B

    x = (/limits(1,1), limits(2,1), limits(3,1)/)

    call spline_field%compute_afield(x, A_spline)
    call spline_field%compute_bfield(x, B_spline)

    deviation(1) = abs(sum(A_spline-A))
    deviation(2) = abs(sum(B_spline-B))

end function compare_spline_with_file


subroutine test_spline_field_init
    type(spline_field_t) :: spline_field

    call print_test("test_spline_field_init")

    call get_example_spline_field(spline_field)

    call print_ok
end subroutine test_spline_field_init


subroutine test_against_original_field
    use neo_field, only: create_field
    use neo_field_base, only: field_t
    use neo_field_mesh, only: field_mesh_t

    real(dp), parameter :: tol = 1.0e-7_dp

    class(field_t), allocatable :: field
    type(field_mesh_t) :: field_mesh
    type(spline_field_t) :: spline_field
    real(dp), dimension(3,2) :: limits
    integer, dimension(3) :: n_points
    real(dp) :: x(3), B(3), B_original(3)

    call print_test("test_against_original_field")

    call create_field(field, "example")

    x = [1.5_dp, 1.5_dp, 1.5_dp]

    limits(1,:) = [1.0_dp, 2.0_dp]
    limits(2,:) = [1.0_dp, 2.0_dp]
    limits(3,:) = [1.0_dp, 2.0_dp]
    n_points = int((limits(:,2) - limits(:,1))/tol**(1.0_dp/5.0_dp)) + 1

    call field_mesh%field_mesh_init_with_field(limits, field, n_points)
    call spline_field%spline_field_init(field_mesh)
    call spline_field%compute_bfield(x, B)
    call field%compute_bfield(x, B_original)

    if (maxval(abs(B - B_original)) > tol) then
        print *, "B != B_original"
        print *, "B = ", B
        print *, "B_original = ", B_original
        call print_fail
        error stop
    end if

    call print_ok
end subroutine test_against_original_field


subroutine test_curla_equal_b
    use util_for_test_field, only: compute_cartesian_curla

    real(dp), parameter :: tol = 1.0e-7_dp

    type(spline_field_t) :: spline_field
    real(dp) :: x(3), B(3), curl_A(3)

    call print_test("test_curla_equal_b")

    x = [1.5_dp, 1.5_dp, 1.5_dp]

    call get_example_spline_field(spline_field, tol)
    call spline_field%compute_bfield(x, B)
    curl_A = compute_cartesian_curla(spline_field, x, tol)

    if (maxval(abs(B - curl_A)) > tol) then
        print *, "curl A != B"
        print *, "B = ", B
        print *, "curl_A = ", curl_A
        call print_fail
        error stop
    end if

    call print_ok
end subroutine test_curla_equal_b


subroutine test_compute_afield_derivatives
    use util_for_test_field, only: compute_cartesian_curla

    real(dp), parameter :: tol = 1.0e-7_dp

    type(spline_field_t) :: spline_field
    real(dp) :: x(3), curla(3), curla_spline_der(3)

    call print_test("test_compute_afield_derivatives")

    call get_example_spline_field(spline_field, tol)

    x = [1.5_dp, 1.5_dp, 1.5_dp]

    curla = compute_cartesian_curla(spline_field, x, tol)
    curla_spline_der = compute_cartesian_curla_from_spline_derivatives(spline_field, x)

    if (maxval(abs(curla - curla_spline_der)) > tol) then
        print *, "curla != curla_spline_der"
        print *, "curla = ", curla
        print *, "curla_spline_der = ", curla_spline_der
        call print_fail
        error stop
    end if

    call print_ok
end subroutine test_compute_afield_derivatives

function compute_cartesian_curla_from_spline_derivatives(spline_field, x) result(curla)
    class(spline_field_t), intent(in) :: spline_field
    real(dp), dimension(3), intent(in) :: x
    real(dp) :: curla(3)

    real(dp) :: dA_dx(3,3)

    call spline_field%compute_afield_derivatives(x, dA_dx)
    curla(1) = dA_dx(3,2) - dA_dx(2,3)
    curla(2) = dA_dx(1,3) - dA_dx(3,1)
    curla(3) = dA_dx(2,1) - dA_dx(1,2)
end function compute_cartesian_curla_from_spline_derivatives


subroutine get_example_spline_field(spline_field, tol)
    use neo_spline_field, only: spline_field_t
    use neo_field, only: create_field
    use neo_field_base, only: field_t
    use neo_field_mesh, only: field_mesh_t

    type(spline_field_t), intent(out) :: spline_field
    real(dp), intent(in), optional :: tol

    class(field_t), allocatable :: field
    type(field_mesh_t) :: field_mesh
    real(dp), dimension(3,2) :: limits
    integer, dimension(3) :: n_points


    limits(1,:) = [1.0_dp, 2.0_dp]
    limits(2,:) = [1.0_dp, 2.0_dp]
    limits(3,:) = [1.0_dp, 2.0_dp]

    call create_field(field, "example")
    if (present(tol)) then
        n_points = int((limits(:,2) - limits(:,1))/tol**(1.0_dp/5.0_dp)) + 1
        call field_mesh%field_mesh_init_with_field(limits, field, n_points)
    else
        call field_mesh%field_mesh_init_with_field(limits, field)
    end if
    call spline_field%spline_field_init(field_mesh)
end subroutine get_example_spline_field

end program test_spline_field