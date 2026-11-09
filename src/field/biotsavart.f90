module neo_biotsavart
implicit none

integer, parameter :: dp = kind(1.0d0)
real(dp), parameter :: clight = 2.99792458d10


!> Filamentary coils as a polygonal chain. Current zero marks end of one coil.
type coils_t
    real(dp), dimension(:), allocatable :: x, y, z, current
end type coils_t


contains


recursive subroutine coils_init(x, y, z, current, coils)
    real(dp), intent(in) :: x(:), y(:), z(:), current(:)
    type(coils_t), intent(out) :: coils

    integer :: n_points

    n_points = size(x)
    call allocate_coils_data(coils, n_points)

    coils%x = x
    coils%y = y
    coils%z = z
    coils%current = current
end subroutine coils_init


recursive subroutine coils_deinit(coils)
    type(coils_t), intent(inout) :: coils

    call deallocate_coils_data(coils)
end subroutine coils_deinit


recursive subroutine load_coils_from_file(filename, coils)
    character(*), intent(in) :: filename
    type(coils_t), intent(out) :: coils

    integer :: unit
    integer :: i, n_points

    open(newunit=unit, file=filename, status="old", action="read")
    read(unit, *) n_points
    call allocate_coils_data(coils, n_points)
    do i = 1, n_points
        read(unit, *) coils%x(i), coils%y(i), coils%z(i), coils%current(i)
    end do
    close(unit)
end subroutine load_coils_from_file


recursive subroutine save_coils_to_file(filename, coils)
    character(*), intent(in) :: filename
    type(coils_t), intent(in) :: coils

    integer :: unit
    integer :: i, n_points

    n_points = size(coils%x)
    open(newunit=unit, file=filename, status="replace", action="write")
    write(unit, *) n_points
    do i = 1, n_points
        write(unit, *) coils%x(i), coils%y(i), coils%z(i), coils%current(i)
    end do
    close(unit)
end subroutine save_coils_to_file


recursive subroutine allocate_coils_data(coils, n_points)
    type(coils_t), intent(out) :: coils
    integer, intent(in) :: n_points

    allocate(coils%x(n_points), coils%y(n_points), coils%z(n_points))
    allocate(coils%current(n_points))
end subroutine allocate_coils_data


recursive subroutine deallocate_coils_data(coils)
    type(coils_t), intent(inout) :: coils

    deallocate(coils%x, coils%y, coils%z, coils%current)
end subroutine deallocate_coils_data


!> Formula of Hanson and Hirshman (2002)
!> A per unit I/c of the straight segment from p_i to p_f at distances
!> d_i = |x - p_i|, d_f = |x - p_f| (Hanson and Hirshman 2002). Consecutive
!> segments share an endpoint distance, so callers pass it instead of
!> recomputing it.
pure function segment_potential(dl, d_i, d_f) result(A)
    real(dp), intent(in) :: dl(3), d_i, d_f
    real(dp) :: A(3), len

    len = norm2(dl)
    A = dl/len*log((d_i + d_f + len)/(d_i + d_f - len))
end function segment_potential

!> B per unit I/c of the straight segment from p_i to p_f, with r_i = x - p_i,
!> r_f = x - p_f and their lengths d_i, d_f.
pure function segment_field(r_i, r_f, d_i, d_f) result(B)
    real(dp), intent(in) :: r_i(3), r_f(3), d_i, d_f
    real(dp) :: B(3), c

    c = (d_i + d_f)/(d_i*d_f*(d_i*d_f + dot_product(r_i, r_f)))
    B(1) = c*(r_i(2)*r_f(3) - r_i(3)*r_f(2))
    B(2) = c*(r_i(3)*r_f(1) - r_i(1)*r_f(3))
    B(3) = c*(r_i(1)*r_f(2) - r_i(2)*r_f(1))
end function segment_field

pure function compute_vector_potential(coils, x) result(A)
    type(coils_t), intent(in) :: coils
    real(dp), intent(in) :: x(3)

    real(dp) :: A(3), r_i(3), r_f(3), d_i, d_f
    integer :: i

    A = 0.0d0
    r_f = x - [coils%x(1), coils%y(1), coils%z(1)]
    d_f = norm2(r_f)
    do i = 1, size(coils%x) - 1
        r_i = r_f
        d_i = d_f
        r_f = x - [coils%x(i + 1), coils%y(i + 1), coils%z(i + 1)]
        d_f = norm2(r_f)
        if (coils%current(i) == 0.0d0) cycle
        A = A + coils%current(i)/clight*segment_potential(r_i - r_f, d_i, d_f)
    end do
end function compute_vector_potential

pure function compute_magnetic_field(coils, x) result(B)
    type(coils_t), intent(in) :: coils
    real(dp), intent(in) :: x(3)

    real(dp) :: B(3), r_i(3), r_f(3), d_i, d_f
    integer :: i

    B = 0.0d0
    r_f = x - [coils%x(1), coils%y(1), coils%z(1)]
    d_f = norm2(r_f)
    do i = 1, size(coils%x) - 1
        r_i = r_f
        d_i = d_f
        r_f = x - [coils%x(i + 1), coils%y(i + 1), coils%z(i + 1)]
        d_f = norm2(r_f)
        if (coils%current(i) == 0.0d0) cycle
        B = B + coils%current(i)/clight*segment_field(r_i, r_f, d_i, d_f)
    end do
end function compute_magnetic_field

!> A, B and their Jacobians dA(i, j) = dA_i/dx_j, dB(i, j) = dB_i/dx_j, summed
!> over the segments with the generated exact segment leaf.
subroutine compute_field_jacobian(coils, x, A, dA, B, dB)
    type(coils_t), intent(in) :: coils
    real(dp), intent(in) :: x(3)
    real(dp), intent(out) :: A(3), dA(3, 3), B(3), dB(3, 3)

    real(dp) :: seg_a(3), seg_b(3), jac_a(3, 3), jac_b(3, 3), c
    integer :: i

    A = 0.0d0
    B = 0.0d0
    dA = 0.0d0
    dB = 0.0d0
    do i = 1, size(coils%x) - 1
        if (coils%current(i) == 0.0d0) cycle
        call biotsavart_segment_jacobian(x(1), x(2), x(3), &
                                         coils%x(i), coils%y(i), coils%z(i), &
                                     coils%x(i + 1), coils%y(i + 1), coils%z(i + 1), &
                         seg_a(1), seg_a(2), seg_a(3), seg_b(1), seg_b(2), seg_b(3), &
       jac_a(1, 1), jac_a(1, 2), jac_a(1, 3), jac_a(2, 1), jac_a(2, 2), jac_a(2, 3), &
                                         jac_a(3, 1), jac_a(3, 2), jac_a(3, 3), &
       jac_b(1, 1), jac_b(1, 2), jac_b(1, 3), jac_b(2, 1), jac_b(2, 2), jac_b(2, 3), &
                                         jac_b(3, 1), jac_b(3, 2), jac_b(3, 3))
        c = coils%current(i)/clight
        A = A + c*seg_a
        B = B + c*seg_b
        dA = dA + c*jac_a
        dB = dB + c*jac_b
    end do
end subroutine compute_field_jacobian

include '../generated/biotsavart_segment_leaves.inc'

end module neo_biotsavart
