program test_jorek_field_values
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use jorek_field_values, only: evaluate_jorek_variable
    use jorek_restart, only: jorek_restart_t
    use util_for_test, only: print_test, print_ok, print_fail

    implicit none

    real(dp), parameter :: tol = 3.0e-14_dp
    integer :: nfail

    nfail = 0
    call test_toroidal_harmonics
    call test_poloidal_derivatives
    call test_rejected_inputs
    if (nfail > 0) error stop

contains

    subroutine test_toroidal_harmonics
        type(jorek_restart_t) :: data
        real(dp) :: value, derivative(3), phi, expected, expected_phi
        integer :: ierr, nfail_before

        call print_test('variable interpolation follows JOREK Fourier slots')
        nfail_before = nfail
        call make_element(data, 5)
        data%values(:, 1, 1, 1) = 1.0_dp
        data%values(:, 2, 1, 1) = 2.0_dp
        data%values(:, 3, 1, 1) = 3.0_dp
        data%values(:, 4, 1, 1) = 4.0_dp
        data%values(:, 5, 1, 1) = 5.0_dp
        phi = 0.31_dp
        expected = 1.0_dp + 2.0_dp*cos(2.0_dp*phi) &
            + 3.0_dp*sin(2.0_dp*phi) + 4.0_dp*cos(4.0_dp*phi) &
            + 5.0_dp*sin(4.0_dp*phi)
        expected_phi = -4.0_dp*sin(2.0_dp*phi) &
            + 6.0_dp*cos(2.0_dp*phi) - 16.0_dp*sin(4.0_dp*phi) &
            + 20.0_dp*cos(4.0_dp*phi)
        call evaluate_jorek_variable(data, 1, 1, 0.23_dp, 0.71_dp, phi, &
            value, derivative, ierr)
        call check_int('ierr', ierr, 0)
        call check_real('value', value, expected)
        call check_real('s derivative', derivative(1), 0.0_dp)
        call check_real('t derivative', derivative(2), 0.0_dp)
        call check_real('phi derivative', derivative(3), expected_phi)
        call report(nfail_before)
    end subroutine test_toroidal_harmonics

    subroutine test_poloidal_derivatives
        type(jorek_restart_t) :: data
        real(dp) :: value, derivative(3), s, t, expected
        integer :: ierr, nfail_before, node
        real(dp), parameter :: node_s(4) = [0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp]
        real(dp), parameter :: node_t(4) = [0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp]

        call print_test('variable interpolation reproduces bicubic derivatives')
        nfail_before = nfail
        call make_element(data, 1)
        do node = 1, 4
            data%values(node, 1, 1, 1) = polynomial(node_s(node), node_t(node))
            data%values(node, 1, 2, 1) = &
                (2.0_dp + 4.0_dp*node_t(node))/3.0_dp
            data%values(node, 1, 3, 1) = &
                (3.0_dp + 4.0_dp*node_s(node))/3.0_dp
            data%values(node, 1, 4, 1) = 4.0_dp/9.0_dp
        end do
        s = 0.37_dp
        t = 0.61_dp
        expected = polynomial(s, t)
        call evaluate_jorek_variable(data, 1, 1, s, t, 0.4_dp, &
            value, derivative, ierr)
        call check_int('ierr', ierr, 0)
        call check_real('value', value, expected)
        call check_real('s derivative', derivative(1), 2.0_dp + 4.0_dp*t)
        call check_real('t derivative', derivative(2), 3.0_dp + 4.0_dp*s)
        call check_real('phi derivative', derivative(3), 0.0_dp)
        call report(nfail_before)
    end subroutine test_poloidal_derivatives

    subroutine test_rejected_inputs
        type(jorek_restart_t) :: data
        real(dp) :: value, derivative(3)
        integer :: ierr, nfail_before

        call print_test('variable interpolation rejects invalid metadata')
        nfail_before = nfail
        call make_element(data, 1)
        call evaluate_jorek_variable(data, 1, 0, 0.5_dp, 0.5_dp, 0.0_dp, &
            value, derivative, ierr)
        call check_int('variable ierr', ierr, 5)
        data%n_tor = 2
        call evaluate_jorek_variable(data, 1, 1, 0.5_dp, 0.5_dp, 0.0_dp, &
            value, derivative, ierr)
        call check_int('Fourier ierr', ierr, 6)
        data%n_tor = 1
        data%n_period = 0
        call evaluate_jorek_variable(data, 1, 1, 0.5_dp, 0.5_dp, 0.0_dp, &
            value, derivative, ierr)
        call check_int('period ierr', ierr, 6)
        call report(nfail_before)
    end subroutine test_rejected_inputs

    subroutine make_element(data, n_tor)
        type(jorek_restart_t), intent(out) :: data
        integer, intent(in) :: n_tor

        data%n_order = 3
        data%n_degrees = 4
        data%n_tor = n_tor
        data%n_period = 2
        data%n_var = 1
        data%n_nodes = 4
        data%n_elements = 1
        data%n_vertex_max = 4
        allocate (data%values(4, n_tor, 4, 1))
        allocate (data%vertex(1, 4), data%size(1, 4, 4))
        data%values = 0.0_dp
        data%vertex(1, :) = [1, 2, 3, 4]
        data%size = 1.0_dp
    end subroutine make_element

    pure real(dp) function polynomial(s, t)
        real(dp), intent(in) :: s, t

        polynomial = 1.0_dp + 2.0_dp*s + 3.0_dp*t + 4.0_dp*s*t
    end function polynomial

    subroutine check_real(name, actual, expected)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: actual, expected

        if (abs(actual - expected) > tol*max(1.0_dp, abs(expected))) then
            call fail(name//' mismatch')
        end if
    end subroutine check_real

    subroutine check_int(name, actual, expected)
        character(len=*), intent(in) :: name
        integer, intent(in) :: actual, expected

        if (actual /= expected) call fail(name//' mismatch')
    end subroutine check_int

    subroutine fail(message)
        character(len=*), intent(in) :: message

        print *, '    ', message
        nfail = nfail + 1
    end subroutine fail

    subroutine report(nfail_before)
        integer, intent(in) :: nfail_before

        if (nfail > nfail_before) then
            call print_fail
        else
            call print_ok
        end if
    end subroutine report

end program test_jorek_field_values
