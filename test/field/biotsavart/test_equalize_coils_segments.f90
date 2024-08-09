program test_equalize_coils_segments
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use test_util, only: print_test, print_ok, print_fail

    implicit none

    call test_compute_coils_segment_lengths
    call test_calc_subknot_xyz
    call test_cut_coils_segments
    call test_equalize_coils_lenghts


    contains


    subroutine test_compute_coils_segment_lengths
        use biotsavart, only: coils_t, coils_init, coils_deinit
        use equalize_coils_segments, only: compute_coils_segments_lengths

        real(dp), parameter :: tol = 1.0e-9

        type(coils_t) :: coils
        real(dp), dimension(4) :: lengths
        real(dp), dimension(4) :: expected_lengths
        integer :: i

        call print_test("compute_coils_segment_lengths")

        call init_diamond_wire_coils(coils)

        lengths = compute_coils_segments_lengths(coils)
        expected_lengths = [sqrt(3.0d0), sqrt(3.0d0), sqrt(3.0d0), sqrt(3.0d0)]

        do i = 1, 3
            if (abs(lengths(i) - expected_lengths(i)) > tol) then
                print *, "lengths(i) = ", lengths(i)
                print *, "expected_lengths(i) = ", expected_lengths(i)
                call print_fail
                error stop
            end if
        end do

        call coils_deinit(coils)

        call print_ok
    end subroutine test_compute_coils_segment_lengths


    subroutine test_calc_subknot_xyz
        use biotsavart, only: coils_t, coils_init, coils_deinit
        use equalize_coils_segments, only: calc_subknot_xyz

        real(dp), parameter :: tol = 1.0e-9

        type(coils_t) :: coils
        integer :: knot
        real(dp), dimension(3) :: xyz, expected_xyz

        call print_test("calc_subknot_xyz")

        call init_diamond_wire_coils(coils)

        do knot = 1, size(coils%x) - 1
            xyz = calc_subknot_xyz(coils, knot, subknot=2, cuts_per_knot=3)
            expected_xyz = [(coils%x(knot) + coils%x(knot+1))/2.0d0, &
                            (coils%y(knot) + coils%y(knot+1))/2.0d0, &
                            (coils%z(knot) + coils%z(knot+1))/2.0d0]
            if (any(abs(xyz - expected_xyz) > tol)) then
                print *, "xyz = ", xyz
                print *, "expected_xyz = ", expected_xyz
                call print_fail
                error stop
            end if
        end do

        call coils_deinit(coils)
        call print_ok
    end subroutine test_calc_subknot_xyz


    subroutine test_cut_coils_segments
        use biotsavart, only: coils_t, coils_init, coils_deinit
        use equalize_coils_segments, only: cut_coils_segments, &
                                           compute_coils_segments_lengths

        real(dp), parameter :: tol = 1.0e-9

        type(coils_t) :: coils
        integer, dimension(4) :: cuts_per_knot
        real(dp) :: expected_x = 0.5d0, expected_y = 1.0d0, expected_current = 4.0d0

        call print_test("cut_coils_segments")

        call init_diamond_wire_coils(coils)

        cuts_per_knot = [0, 1, 0, 1]

        call cut_coils_segments(coils, cuts_per_knot)
        if (size(coils%x) /= 7) then
            print *, "Coil length mismatch"
            print *, "len(coils%x) = ", size(coils%x)
            call print_fail
            error stop
        end if
        if ((coils%x(3) - expected_x) > tol) then
            print *, "Coil x mismatch"
            print *, "coils%x(3) = ", coils%x(3)
            print *, "expected = ", expected_x
            call print_fail
            error stop
        end if
        if ((coils%y(7)-expected_y) > tol) then
            print *, "Coil y mismatch"
            print *, "coils%y(7) = ", coils%y(3)
            print *, "expected = ", expected_y
            call print_fail
            error stop
        end if
        if ((coils%current(6)-expected_current) > tol) then
            print *, "Coil current mismatch"
            print *, "coils%current(6) = ", coils%current(6)
            print *, "expected = ", expected_current
            call print_fail
            error stop
        end if

        call coils_deinit(coils)
        call print_ok
    end subroutine test_cut_coils_segments


    subroutine test_equalize_coils_lenghts
        use biotsavart, only: coils_t, coils_init, coils_deinit
        use equalize_coils_segments, only: equalize_coils_segments_lengths, &
                                           compute_coils_segments_lengths

        type(coils_t) :: coils, old_coils
        real(dp), dimension(:), allocatable :: lengths
        real(dp) :: min_length

        call print_test("equalize_coils_lenghts")

        call init_diamond_wire_coils(coils)
        call equalize_coils_segments_lengths(coils)
        call init_diamond_wire_coils(old_coils)
        if (.not.(are_coils_equal(coils, old_coils))) then
            print *, "Equalized coil is not unchanged"
            call print_fail
            error stop
        end if
        call coils_deinit(coils)
        call coils_deinit(old_coils)

        call init_diamond_wire_coils(coils)
        coils%x(2) = coils%x(2) + 0.5d0
        coils%y(2) = coils%y(2) - 0.5d0
        coils%z(2) = coils%z(2) - 0.5d0
        min_length = minval(compute_coils_segments_lengths(coils))
        call equalize_coils_segments_lengths(coils)
        allocate(lengths(size(coils%x) - 1))
        lengths = compute_coils_segments_lengths(coils)
        if (any(abs(lengths - min_length) > min_length)) then
            print *, "Coil segments lengths mismatch"
            print *, "lengths = ", lengths
            call print_fail
            error stop
        end if
        deallocate(lengths)

        call coils_deinit(coils)

        call print_ok
    end subroutine test_equalize_coils_lenghts


    function are_coils_equal(coil1, coil2)
        use biotsavart, only: coils_t

        real(dp), parameter :: tol = 1.0e-9

        type(coils_t), intent(in) :: coil1, coil2
        logical :: are_coils_equal
        integer :: i

        are_coils_equal = .true.
        if (size(coil1%x) /= size(coil2%x)) then
            are_coils_equal = .false.
            return
        end if
        do i = 1, size(coil1%x)
            if (abs(coil1%x(i) - coil2%x(i)) > tol) then
                are_coils_equal = .false.
                return
            end if
            if (abs(coil1%y(i) - coil2%y(i)) > tol) then
                are_coils_equal = .false.
                return
            end if
            if (abs(coil1%z(i) - coil2%z(i)) > tol) then
                are_coils_equal = .false.
                return
            end if
            if (abs(coil1%current(i) - coil2%current(i)) > tol) then
                are_coils_equal = .false.
                return
            end if
        end do
    end function are_coils_equal


    subroutine init_diamond_wire_coils(coils)
        use biotsavart, only: coils_t, coils_init

        type(coils_t), intent(out) :: coils

        real(dp), dimension(5) :: x, y, z, current

        x = [-1.0d0, 0.0d0, 1.0d0, 0.0d0, -1.0d0]
        y = [1.0d0, 0.0d0, -1.0d0, 0.0d0, 1.0d0]
        z = [0.0d0, 1.0d0, 0.0d0, -1.0d0, 0.0d0]
        current = [1.0d0, 2.0d0, 3.0d0, 4.0d0, 0.0d0]

        call coils_init(x, y, z, current, coils)
    end subroutine init_diamond_wire_coils


end program test_equalize_coils_segments