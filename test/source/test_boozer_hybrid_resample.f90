!> hybrid_grid resamples the source shifts lambda + iota nu and nu by a DFT.
!> On the source nodes the synthesis must return the nodal values exactly.
!> The library is built with -fno-realloc-lhs, so the spectral work array
!> has to be allocated before assignment.
program test_boozer_hybrid_resample
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use boozer_hybrid, only: boozer_hybrid_t, hybrid_surface_t, &
                             boozer_hybrid_init, hybrid_surface, hybrid_grid
    implicit none

    real(dp), parameter :: TOL = 1.0e-10_dp
    integer :: failures

    failures = 0
    call check_nodes('wout.nc', 0.5_dp)
    call check_nodes('wout_li383.nc', 0.3_dp)
    if (failures /= 0) error stop 'test_boozer_hybrid_resample FAILED'
    print '(a)', 'test_boozer_hybrid_resample: all checks passed'

contains

    subroutine check_nodes(filename, s)
        character(*), intent(in) :: filename
        real(dp), intent(in) :: s
        integer, parameter :: REFINE = 3
        type(boozer_hybrid_t) :: h
        type(hybrid_surface_t) :: sv
        real(dp), allocatable :: theta(:), phi(:), values(:, :, :), shifts(:, :, :)
        real(dp), allocatable :: expected(:, :)
        integer :: i

        call boozer_hybrid_init(h, filename, 8, 6)
        call hybrid_surface(h, s, sv, with_spectrum=.false.)
        allocate (theta(REFINE*h%nt), phi(REFINE*h%np))
        theta = [(h%theta(1) + real(i - 1, dp)*(h%theta(2) - h%theta(1))/REFINE, &
                  i = 1, size(theta))]
        phi = [(h%phi(1) + real(i - 1, dp)*(h%phi(2) - h%phi(1))/REFINE, &
                i = 1, size(phi))]
        allocate (values(5, size(theta), size(phi)), shifts(2, size(theta), size(phi)))
        allocate (expected(h%nt, h%np))
        call hybrid_grid(h, sv, theta, phi, values, shifts, .true.)
        expected = reshape(sv%lam + sv%iota*sv%nu, [h%nt, h%np])
        call check(filename//': theta shift on nodes', &
                   shifts(1, 1::REFINE, 1::REFINE), expected)
        expected = reshape(sv%nu, [h%nt, h%np])
        call check(filename//': phi shift on nodes', &
                   shifts(2, 1::REFINE, 1::REFINE), expected)
        if (.not. all(ieee_is_finite(values))) then
            print '(a)', 'FAIL '//filename//': grid values are not finite'
            failures = failures + 1
        end if
    end subroutine check_nodes

    subroutine check(label, actual, expected)
        character(*), intent(in) :: label
        real(dp), intent(in) :: actual(:, :), expected(:, :)
        real(dp) :: error

        error = maxval(abs(actual - expected))/max(1.0_dp, maxval(abs(expected)))
        if (.not. (error <= TOL)) then
            print '(a, 2es14.5)', 'FAIL '//label//': ', error, TOL
            failures = failures + 1
        end if
    end subroutine check

end program test_boozer_hybrid_resample
