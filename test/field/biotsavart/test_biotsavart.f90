program test_biotsavart
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use test_util, only: print_test, print_ok, print_fail

    implicit none

    character(*), parameter :: test_coils_file = "coils.test"
    !character(*), parameter :: coils_file = "/proj/plasma/DATA/W7X/COILS/coils.w7x"
    real(8), parameter :: pi = 3.14159265358979323846d0

    real(dp), parameter :: large_distance = 1.0d3

    call test_load_coils_file
    call test_compute_vector_potential
    call test_compute_vector_potential_circular_loop
    call test_compute_magnetic_field
    call test_compute_magnetic_field_circular_loop

    contains

    subroutine test_load_coils_file
        use biotsavart, only: coils_t, load_coils_from_file, coils_deinit

        type(coils_t) :: coils

        call print_test("load_coils_file")

        call create_straight_wire_coils_file
        call load_coils_from_file(test_coils_file, coils)
        call remove_test_coils_file

        if (size(coils%x) /= 4) then
            print *, "Coil length mismatch"
            print *, "len(coils%x) = ", size(coils%x)
            call print_fail
            error stop
        end if

        call coils_deinit(coils)

        call print_ok
    end subroutine test_load_coils_file


    subroutine test_compute_vector_potential
        use biotsavart, only: coils_t, compute_vector_potential, &
                              coils_deinit, clight

        real(dp), parameter :: tol = 1.0e-9
        integer, parameter :: N_TEST = 3

        type(coils_t) :: coils
        real(dp) :: x_test(3, N_TEST)
        real(dp), dimension(3) :: x, A, A_analytic
        integer :: i

        call print_test("compute_vector_potential")

        x_test(:, 1) = [0.4, 0.3, 0.8]
        x_test(:, 2) = [0.0, 0.2, -0.3]
        x_test(:, 3) = [1.0d2, -1.0d2, 1.0d2]

        call init_straight_wire_coils(coils)

        do i = 1, N_TEST
            x = x_test(:, i)
            A_analytic = vector_potential_straight_wire(x, large_distance, 1.0d0)
            A = compute_vector_potential(coils, x)
            if (any(abs(A - A_analytic)*clight > tol)) then
                print *, "A = ", A(3)*clight
                print *, "A_analytic = ", A_analytic(3)*clight
                print *, "Ratio = ", A(3) / A_analytic(3)
                call print_fail
                error stop
            end if
        end do

        call coils_deinit(coils)

        call print_ok
    end subroutine test_compute_vector_potential


    function vector_potential_straight_wire(x, L, current) result(A)
        use biotsavart, only: clight

        real(dp), dimension(3), intent(in) :: x
        real(dp), intent(in) :: L, current

        real(dp), dimension(3) :: A
        real(dp) :: R, z, L_half, A_z

        R = Rcyl(x)
        z = x(3)
        L_half = L/2.0d0
        A_z = current / clight * log((L_half - z + sqrt((L_half - z)**2 + R**2)) / &
                                     (-L_half - z + sqrt((-L_half - z)**2 + R**2)))
        A = [0.0d0, 0.0d0, A_z]
    end function vector_potential_straight_wire


    subroutine test_compute_vector_potential_circular_loop
        use biotsavart, only: coils_t, compute_vector_potential, &
                              coils_deinit, clight, calc_norm

        real(dp), parameter :: tol = 1.0e-5
        integer, parameter :: N_TEST = 3

        type(coils_t) :: coils
        real(dp) :: x_test(3, N_TEST)
        real(dp), dimension(3) :: x, A, A_analytic
        integer :: number_of_segments, i

        call print_test("compute_vector_potential_circular_loop")

        x_test(:, 1) = [0.0d0, 0.0d0, 0.0d0]
        x_test(:, 2) = [0.0d0, 0.0d0, -2.0d0]
        x_test(:, 3) = [0.0d0, 0.0d0, +1.0d2]

        number_of_segments = int(2*pi/tol) + 2
        call init_circular_loop_coils(coils, number_of_segments)

        do i = 1, N_TEST
            x = x_test(:, i)
            A_analytic = vector_potential_circular_loop_on_axis(x(3), 1.0d0, 1.0d0)
            A = compute_vector_potential(coils, x)
            if (any(abs(A - A_analytic)*clight > tol)) then
                print *, "A = ", calc_norm(A)*clight
                print *, "A_analytic = ", calc_norm(A_analytic)*clight
                call print_fail
                error stop
            end if
        end do

        call coils_deinit(coils)

        call print_ok
    end subroutine test_compute_vector_potential_circular_loop


    function vector_potential_circular_loop_on_axis(z, R0, current) result(A)
        use biotsavart, only: clight

        real(dp), intent(in) :: z, R0, current

        real(dp), dimension(3) :: A

        A = 0.0d0 * current / clight * 2 * pi * R0 / sqrt(R0**2 + z**2)
    end function vector_potential_circular_loop_on_axis


    subroutine test_compute_magnetic_field
        use biotsavart, only: coils_t, compute_magnetic_field, &
                              coils_deinit, clight, calc_norm

        real(dp), parameter :: tol = 1.0e-9
        integer, parameter :: N_TEST = 3

        type(coils_t) :: coils
        real(dp) :: x_test(3, N_TEST), x(3), B(3), B_analytic(3)
        integer :: i

        call print_test("compute_magnetic_field")

        call init_straight_wire_coils(coils)

        x_test(:, 1) = [1.0d0, 0.0d0, 0.0d0]
        x_test(:, 2) = [0.0d0, 0.2d0, -0.3d0]
        x_test(:, 3) = [1.0d2, -1.0d2, 1.0d2]
        
        do i = 1, N_TEST
            x = x_test(:, i)
            B = compute_magnetic_field(coils, x)
            B_analytic = magnetic_field_straight_wire(x, large_distance, 1.0d0)
            if (any(abs(B - B_analytic)*clight > tol)) then
                print *, "B = ", calc_norm(B)*clight
                print *, "B_analytic = ", calc_norm(B_analytic)*clight
                print *, "Ratio = ", calc_norm(B) / calc_norm(B_analytic)
                call print_fail
                error stop
            end if
        end do

        call coils_deinit(coils)

        call print_ok
    end subroutine test_compute_magnetic_field


    function magnetic_field_straight_wire(x, L, current) result(B)
        use biotsavart, only: clight

        real(dp), dimension(3), intent(in) :: x
        real(dp), intent(in) :: L, current

        real(dp), dimension(3) :: B
        real(dp) :: R, z, L_half, B_phi, B_x, B_y, sin_theta_1, sin_theta_2

        R = Rcyl(x)
        z = x(3)
        L_half = L/2.0d0
        sin_theta_1 = (L_half + z) / sqrt((L_half + z)**2 + R**2)
        sin_theta_2 = (L_half - z) / sqrt((L_half - z)**2 + R**2)
        B_phi = current / clight * 1.0d0 / R * (sin_theta_1 + sin_theta_2)
        B_x = -B_phi * x(2) / R
        B_y = B_phi * x(1) / R
        B = [B_x, B_y, 0.0d0]
    end function magnetic_field_straight_wire


    subroutine test_compute_magnetic_field_circular_loop
        use biotsavart, only: coils_t, compute_magnetic_field, &
                              coils_deinit, clight, calc_norm

        real(dp), parameter :: tol = 1.0e-5
        integer, parameter :: N_TEST = 3

        type(coils_t) :: coils
        real(dp) :: x_test(3, N_TEST), x(3), B(3), B_analytic(3)
        integer :: number_of_segments, i

        call print_test("compute_magnetic_field_circular_loop")

        number_of_segments = int(2*pi/tol)*2 + 2
        call init_circular_loop_coils(coils, number_of_segments)

        x_test(:, 1) = [0.0d0, 0.0d0, 0.0d0]
        x_test(:, 2) = [0.0d0, 0.0d0, -2.0d0]
        x_test(:, 3) = [0.0d0, 0.0d0, +1.0d2]

        do i = 1, N_TEST
            x = x_test(:, i)
            B = compute_magnetic_field(coils, x)
            B_analytic = magnetic_field_circular_loop_on_axis(x(3), 1.0d0, 1.0d0)
            if (any(abs(B - B_analytic)*clight > tol)) then
                print *, "B = ", B * clight
                print *, "B_analytic = ", B_analytic * clight
                call print_fail
                error stop
            end if
        end do

        call coils_deinit(coils)

        call print_ok
    end subroutine test_compute_magnetic_field_circular_loop

    
    function magnetic_field_circular_loop_on_axis(z, R0, current) result(B)
        use biotsavart, only: clight

        real(dp), intent(in) :: z, R0, current

        real(dp), dimension(3) :: B
        real(dp) :: B_z

        B_z = current / clight * 2.0d0 * pi * R0**2 / sqrt(R0**2 + z**2)**3
        B = [0.0d0, 0.0d0, B_z]
    end function magnetic_field_circular_loop_on_axis


    function Rcyl(x)
        real(dp), dimension(3), intent(in) :: x
        real(dp) :: Rcyl

        Rcyl = sqrt(x(1)**2 + x(2)**2)
    end function Rcyl


    subroutine create_straight_wire_coils_file
        use biotsavart, only: coils_t, save_coils_to_file

        type(coils_t) :: coils

        call init_straight_wire_coils(coils)
        call save_coils_to_file(test_coils_file, coils)
    end subroutine create_straight_wire_coils_file


    subroutine init_straight_wire_coils(coils)
        use biotsavart, only: coils_t, coils_init

        type(coils_t), intent(out) :: coils

        real(dp), dimension(4) :: x, y, z, current
        real(dp)               :: L

        x = [0.0d0, 0.0d0, 0.0d0, 0.0d0]
        y = [0.0d0, 0.0d0, 0.0d0, 0.0d0]
        L = large_distance
        z = [-L/2.0d0, -L/2.5d0, 0.0d0, L/2.0d0]
        current = [1.0d0, 1.0d0, 1.0d0, 0.0d0]

        call coils_init(x, y, z, current, coils)
    end subroutine init_straight_wire_coils


    subroutine init_circular_loop_coils(coils, number_of_segments)
        use biotsavart, only: coils_t, coils_init

        type(coils_t), intent(out) :: coils
        integer, intent(in) :: number_of_segments

        real(dp), dimension(number_of_segments) :: x, y, z, current
        real(dp)               :: R, theta, dtheta
        integer                :: i

        R = 1.0d0
        dtheta = 2.0d0 * pi / number_of_segments
        do i = 1, number_of_segments
            theta = dtheta * (i - 1)
            x(i) = R * cos(theta)
            y(i) = R * sin(theta)
            z(i) = 0.0d0
            current(i) = 1.0d0
        end do
        current(number_of_segments) = 0.0d0

        call coils_init(x, y, z, current, coils)
    end subroutine init_circular_loop_coils


    subroutine remove_test_coils_file
        call system("rm -f " // test_coils_file)
    end subroutine remove_test_coils_file

    
end program test_biotsavart
