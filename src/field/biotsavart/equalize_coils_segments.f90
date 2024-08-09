module equalize_coils_segments
use, intrinsic :: iso_fortran_env, only: dp => real64
use biotsavart, only: coils_t
implicit none


contains


subroutine equalize_coils_segments_lengths(coils)
    type(coils_t), intent(inout) :: coils

    real(dp), dimension((size(coils%x) - 1)) :: lengths
    real(dp) :: min_length
    integer, dimension((size(coils%x) - 1)) :: cuts_per_knot

    lengths = compute_coils_segments_lengths(coils)
    min_length = minval(lengths)
    cuts_per_knot = max(floor(lengths / min_length - 1), 0)
    call cut_coils_segments(coils, cuts_per_knot)

end subroutine equalize_coils_segments_lengths


function compute_coils_segments_lengths(coils) result(lenghts)
    type(coils_t), intent(in) :: coils

    real(dp), dimension(size(coils%x) - 1) :: lenghts
    integer :: segment

    do segment = 1, size(coils%x) - 1

        lenghts(segment) = sqrt((coils%x(segment+1) - coils%x(segment))**2 + &
                                (coils%y(segment+1) - coils%y(segment))**2 + &
                                (coils%z(segment+1) - coils%z(segment))**2)
    end do
end function compute_coils_segments_lengths


subroutine cut_coils_segments(coils, cuts_per_knot)
    use biotsavart, only: coils_init, coils_deinit
    
    type(coils_t), intent(inout) :: coils
    integer, dimension(:), intent(in) :: cuts_per_knot
    
    integer :: total_knots
    real(dp), dimension(:), allocatable :: x, y, z, current
    integer :: knot, subknot, global_knot_index
    real(dp), dimension(3) :: subknot_xyz

    total_knots = sum(cuts_per_knot) + size(coils%x)
    allocate(x(total_knots), y(total_knots), &
             z(total_knots), current(total_knots))
    global_knot_index = 0
    do knot = 1, size(cuts_per_knot)
        do subknot = 0, cuts_per_knot(knot)
            subknot_xyz = calc_subknot_xyz(coils, knot, subknot, cuts_per_knot(knot))
            global_knot_index = global_knot_index + 1
            x(global_knot_index) = subknot_xyz(1)
            y(global_knot_index) = subknot_xyz(2)
            z(global_knot_index) = subknot_xyz(3)
            current(global_knot_index) = coils%current(knot)
        end do
    end do
    
    x(total_knots) = coils%x(size(coils%x))
    y(total_knots) = coils%y(size(coils%y))
    z(total_knots) = coils%z(size(coils%z))
    current(total_knots) = 0.0d0
    call coils_deinit(coils)
    call coils_init(x, y, z, current, coils)
    deallocate(x, y, z, current)
end subroutine cut_coils_segments

function calc_subknot_xyz(coils, knot, subknot, cuts_per_knot) result(xyz)
    type(coils_t), intent(in) :: coils
    integer, intent(in) :: knot, subknot
    integer, intent(in) :: cuts_per_knot

    real(dp) :: ratio
    real(dp) :: xyz(3)

    ratio = 1.0d0 / (cuts_per_knot + 1) * subknot
    xyz(1) = coils%x(knot) + (coils%x(knot+1) - coils%x(knot)) * ratio
    xyz(2) = coils%y(knot) + (coils%y(knot+1) - coils%y(knot)) * ratio
    xyz(3) = coils%z(knot) + (coils%z(knot+1) - coils%z(knot)) * ratio
end function calc_subknot_xyz

end module equalize_coils_segments