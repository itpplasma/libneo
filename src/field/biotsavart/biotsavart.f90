module biotsavart
use, intrinsic :: iso_fortran_env, only: dp => real64
implicit none

real(dp), parameter :: clight = 2.99792458d10


!> Filamentary coils as a polygonal chain. Current zero marks end of one coil.
type coils_t
    real(dp), dimension(:), allocatable :: x, y, z, current
end type coils_t


contains


subroutine coils_init(x, y, z, current, coils)
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


subroutine coils_deinit(coils)
    type(coils_t), intent(inout) :: coils

    call deallocate_coils_data(coils)
end subroutine coils_deinit


subroutine load_coils_from_file(filename, coils)
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


subroutine save_coils_to_file(filename, coils)
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


subroutine allocate_coils_data(coils, n_points)
    type(coils_t), intent(out) :: coils
    integer, intent(in) :: n_points

    allocate(coils%x(n_points), coils%y(n_points), coils%z(n_points))
    allocate(coils%current(n_points))
end subroutine allocate_coils_data


subroutine deallocate_coils_data(coils)
    type(coils_t), intent(inout) :: coils

    deallocate(coils%x, coils%y, coils%z, coils%current)
end subroutine deallocate_coils_data


!> Formula of Hanson and Hirshman (2002)
function compute_vector_potential(coils, x) result(A)
    type(coils_t), intent(in) :: coils
    real(dp), intent(in) :: x(3)

    real(dp) :: A(3), dx_i(3), dx_f(3), dl(3), R_i, R_f, L, eps, log_term
    integer :: i

    A = 0.0d0
    do i = 1, size(coils%x) - 1
        dl = get_segment_vector(coils, i)
        dx_i = get_vector_from_segment_start_to_x(coils, i, x)
        dx_f = get_vector_from_segment_end_to_x(coils, i, x)
        R_i = calc_norm(dx_i)
        R_f = calc_norm(dx_f)
        L = calc_norm(dl)
        eps = L / (R_i + R_f)
        log_term = log((1.0d0 + eps) / (1.0d0 - eps))
        A = A + (coils%current(i) / clight) * (dl / L) * log_term
    end do
end function compute_vector_potential


function compute_magnetic_field(coils, x) result(B)
    type(coils_t), intent(in) :: coils
    real(dp), intent(in) :: x(3)

    real(dp) :: B(3), dx_i(3), dx_f(3), dl(3), R_i, R_f, L, eps, dx_i_hat(3), dl_hat(3)
    integer :: i

    B = 0.0d0
    do i = 1, size(coils%x) - 1
        dl = get_segment_vector(coils, i)
        dx_i = get_vector_from_segment_start_to_x(coils, i, x)
        dx_f = get_vector_from_segment_end_to_x(coils, i, x)
        R_i = calc_norm(dx_i)
        R_f = calc_norm(dx_f)
        L = calc_norm(dl)
        eps = L / (R_i + R_f)
        dx_i_hat = dx_i / R_i
        dl_hat = dl / L
        B = B + (coils%current(i) / clight) * cross_product(dl_hat, dx_i_hat) * &
                (1 / R_f) * (2*eps / (1.0d0 - eps**2))
    end do
end function compute_magnetic_field


function get_segment_vector(coils, i) result(dl)
    type(coils_t), intent(in) :: coils
    integer, intent(in) :: i

    real(dp) :: dl(3)

    dl(1) = coils%x(i+1) - coils%x(i)
    dl(2) = coils%y(i+1) - coils%y(i)
    dl(3) = coils%z(i+1) - coils%z(i)
end function get_segment_vector


function get_vector_from_segment_start_to_x(coils, i, x) result(dx_i)
    type(coils_t), intent(in) :: coils
    integer, intent(in) :: i
    real(dp), intent(in) :: x(3)

    real(dp) :: dx_i(3)

    dx_i(1) = x(1) - coils%x(i)
    dx_i(2) = x(2) - coils%y(i)
    dx_i(3) = x(3) - coils%z(i)
end function get_vector_from_segment_start_to_x


function get_vector_from_segment_end_to_x(coils, i, x) result(dx_f)
    type(coils_t), intent(in) :: coils
    integer, intent(in) :: i
    real(dp), intent(in) :: x(3)

    real(dp) :: dx_f(3)

    dx_f(1) = x(1) - coils%x(i+1)
    dx_f(2) = x(2) - coils%y(i+1)
    dx_f(3) = x(3) - coils%z(i+1)
end function get_vector_from_segment_end_to_x


function calc_norm(v) result(norm)
    real(dp), intent(in) :: v(3)
    real(dp) :: norm

    norm = sqrt(v(1)**2 + v(2)**2 + v(3)**2)
end function calc_norm


function cross_product(a, b) result(c)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c(3)

    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
end function cross_product


end module biotsavart
