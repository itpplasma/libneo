module neo_biotsavart
implicit none

integer, parameter :: dp = kind(1.0d0)
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


subroutine load_coils_from_gpec_file(filename, coils)
    character(*), intent(in) :: filename
    type(coils_t), intent(out) :: coils

    real(dp), parameter :: length_si_to_cgs = 1.0d2
    integer :: unit
    integer :: ncoil, nseg, idum, kc, ks, total_segments, seg_idx, coil_start_idx
    real(dp) :: ddum, nwind

    open(newunit=unit, file=filename, status="old", action="read")
    read(unit, *) ncoil, idum, nseg, ddum
    nseg = nseg - 1  ! Convert from number of points to number of segments
    nwind = ddum

    ! Total: ncoil coils * (nseg points + 1 closing point)
    total_segments = ncoil * (nseg + 1)
    call allocate_coils_data(coils, total_segments)

    seg_idx = 0
    do kc = 1, ncoil
        ! Read nseg points with nonzero current
        do ks = 1, nseg
            seg_idx = seg_idx + 1
            read(unit, *) coils%x(seg_idx), coils%y(seg_idx), coils%z(seg_idx)
            coils%x(seg_idx) = coils%x(seg_idx) * length_si_to_cgs
            coils%y(seg_idx) = coils%y(seg_idx) * length_si_to_cgs
            coils%z(seg_idx) = coils%z(seg_idx) * length_si_to_cgs
            coils%current(seg_idx) = nwind
        end do

        ! Read the closing point (equals first point)
        ! This point has ZERO current so the segment from here to next coil contributes nothing
        seg_idx = seg_idx + 1
        read(unit, *) coils%x(seg_idx), coils%y(seg_idx), coils%z(seg_idx)
        coils%x(seg_idx) = coils%x(seg_idx) * length_si_to_cgs
        coils%y(seg_idx) = coils%y(seg_idx) * length_si_to_cgs
        coils%z(seg_idx) = coils%z(seg_idx) * length_si_to_cgs
        coils%current(seg_idx) = 0.0d0  ! Zero current: marks end of coil
    end do
    close(unit)
end subroutine load_coils_from_gpec_file


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


subroutine compute_magnetic_field_array(coils, points, B)
    type(coils_t), intent(in) :: coils
    real(dp), intent(in) :: points(:, :) !< Cartesian evaluation points
    real(dp), intent(out) :: B(:, :)      !< Cartesian magnetic field

    integer :: n_segments, n_points, iseg, ip
    real(dp), allocatable :: seg_start(:, :), seg_end(:, :)
    real(dp), allocatable :: dl_hat(:, :), seg_length(:), current_fac(:)
    real(dp) :: dx_i(3), dx_f(3), R_i, R_f, eps, denom, factor
    real(dp) :: cross(3)

    if (size(points, 1) /= 3) then
        B(:, :) = 0.0d0
        return
    end if

    if (size(B, 1) /= 3 .or. size(B, 2) /= size(points, 2)) then
        B(:, :) = 0.0d0
        return
    end if

    n_segments = size(coils%x) - 1
    n_points = size(points, 2)

    if (n_segments <= 0 .or. n_points <= 0) then
        if (n_points > 0) B(:, :) = 0.0d0
        return
    end if

    allocate(seg_start(3, n_segments), seg_end(3, n_segments))
    allocate(dl_hat(3, n_segments), seg_length(n_segments))
    allocate(current_fac(n_segments))

    do iseg = 1, n_segments
        seg_start(:, iseg) = [coils%x(iseg), coils%y(iseg), coils%z(iseg)]
        seg_end(:, iseg) = [coils%x(iseg + 1), coils%y(iseg + 1), coils%z(iseg + 1)]
        dl_hat(:, iseg) = seg_end(:, iseg) - seg_start(:, iseg)
        seg_length(iseg) = calc_norm(dl_hat(:, iseg))
        if (seg_length(iseg) > 0.0d0) then
            dl_hat(:, iseg) = dl_hat(:, iseg) / seg_length(iseg)
        else
            dl_hat(:, iseg) = 0.0d0
        end if
        current_fac(iseg) = coils%current(iseg) / clight
    end do

    B(:, :) = 0.0d0

    do iseg = 1, n_segments
        if (current_fac(iseg) == 0.0d0 .or. seg_length(iseg) == 0.0d0) cycle
        do ip = 1, n_points
            dx_i = points(:, ip) - seg_start(:, iseg)
            dx_f = points(:, ip) - seg_end(:, iseg)
            R_i = calc_norm(dx_i)
            R_f = calc_norm(dx_f)
            if (R_i == 0.0d0 .or. R_f == 0.0d0) cycle
            denom = R_i + R_f
            if (denom == 0.0d0) cycle
            eps = seg_length(iseg) / denom
            if (abs(eps) >= 1.0d0) cycle
            cross = cross_product(dl_hat(:, iseg), dx_i / R_i)
            factor = current_fac(iseg) * (2.0d0 * eps / (1.0d0 - eps * eps)) / R_f
            B(:, ip) = B(:, ip) + factor * cross
        end do
    end do

end subroutine compute_magnetic_field_array


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


end module neo_biotsavart
