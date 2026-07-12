program test_jorek_restart_fixture
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use jorek_bezier, only: evaluate_jorek_geometry
    use jorek_field_values, only: evaluate_jorek_variable
    use jorek_model303_field, only: evaluate_jorek_model303_b
    use jorek_restart, only: jorek_restart_t, load_jorek_restart, &
        free_jorek_restart
    use util_for_test, only: print_test, print_ok, print_fail

    implicit none

    character(len=*), parameter :: fixture_env = 'LIBNEO_JOREK_RESTART_FIXTURE'
    real(dp), parameter :: tol = 2.0e-13_dp
    character(len=1024) :: fixture
    integer :: env_status, nfail

    nfail = 0
    call get_environment_variable(fixture_env, fixture, status=env_status)
    if (env_status /= 0 .or. len_trim(fixture) == 0) then
        print '(A)', 'SKIP: '//fixture_env//' not set'
        stop
    end if

    call test_fixture_facts

    if (nfail > 0) error stop

contains

    subroutine test_fixture_facts
        type(jorek_restart_t) :: data
        integer :: ierr

        call print_test('JOREK restart fixture dimensions and connectivity')

        call load_jorek_restart(trim(fixture), data, ierr)
        call check_int('ierr', ierr, 0)
        if (ierr /= 0) then
            call print_fail
            return
        end if

        call check_int('rst_hdf5_version', data%rst_hdf5_version, 2)
        call check_int('n_nodes', data%n_nodes, 1561)
        call check_int('n_elements', data%n_elements, 1485)
        call check_int('n_degrees', data%n_degrees, 4)
        call check_int('n_dim', data%n_dim, 2)
        call check_int('n_vertex_max', data%n_vertex_max, 4)
        call check_int('n_order', data%n_order, 3)

        call check_array_shapes(data)
        call check_nodal_r_range(data)
        call check_vertex_bounds(data)
        call check_neighbour_bounds(data)
        call check_geometry_oracle(data)
        call check_value_oracle(data)
        call check_model303_field_oracle(data)

        call free_jorek_restart(data)

        if (nfail > 0) then
            call print_fail
        else
            call print_ok
        end if
    end subroutine test_fixture_facts

    subroutine check_array_shapes(data)
        type(jorek_restart_t), intent(in) :: data

        call check_shape('values', shape(data%values), &
            [data%n_nodes, data%n_tor, data%n_degrees, data%n_var])
        call check_shape('x', shape(data%x), &
            [data%n_nodes, data%n_coord_tor, data%n_degrees, &
            data%n_dim])
        call check_shape('vertex', shape(data%vertex), &
            [data%n_elements, data%n_vertex_max])
        call check_shape('size', shape(data%size), &
            [data%n_elements, data%n_vertex_max, data%n_degrees])
        call check_shape('neighbours', shape(data%neighbours), &
            [data%n_elements, data%n_vertex_max])
    end subroutine check_array_shapes

    subroutine check_nodal_r_range(data)
        type(jorek_restart_t), intent(in) :: data

        real(dp) :: r_min, r_max

        r_min = minval(data%x(:, 1, 1, 1))
        r_max = maxval(data%x(:, 1, 1, 1))
        if (r_min <= 0.0_dp .or. r_max <= r_min) then
            call fail('nodal R value DOFs do not span positive radii')
            print *, '    R range:', r_min, r_max
        end if
    end subroutine check_nodal_r_range

    subroutine check_vertex_bounds(data)
        type(jorek_restart_t), intent(in) :: data

        if (minval(data%vertex) < 1 .or. maxval(data%vertex) > data%n_nodes) then
            call fail('vertex indices outside [1, n_nodes]')
        end if
    end subroutine check_vertex_bounds

    subroutine check_neighbour_bounds(data)
        type(jorek_restart_t), intent(in) :: data

        if (minval(data%neighbours) < -1 &
            .or. maxval(data%neighbours) > data%n_elements) then
            call fail('neighbours outside [-1, n_elements]')
        end if
    end subroutine check_neighbour_bounds

    subroutine check_geometry_oracle(data)
        type(jorek_restart_t), intent(in) :: data

        real(dp) :: rz(2), rz_st(2, 2)
        integer :: ierr

        call evaluate_jorek_geometry(data, 1, 0.2_dp, 0.7_dp, rz, rz_st, ierr)
        call check_int('geometry ierr', ierr, 0)
        call check_dp('R', rz(1), 1.7020710470004476_dp)
        call check_dp('Z', rz(2), 0.03437336518385514_dp)
        call check_dp('R_s', rz_st(1, 1), -0.02304197881763608_dp)
        call check_dp('Z_s', rz_st(2, 1), -0.10438992315046813_dp)
        call check_dp('R_t', rz_st(1, 2), 0.0006138917527668017_dp)
        call check_dp('Z_t', rz_st(2, 2), -0.00017457983954904943_dp)
    end subroutine check_geometry_oracle

    subroutine check_value_oracle(data)
        type(jorek_restart_t), intent(in) :: data

        real(dp) :: value, derivative(3)
        integer :: ierr

        call evaluate_jorek_variable(data, 1, 1, 0.2_dp, 0.7_dp, 0.4_dp, &
            value, derivative, ierr)
        call check_int('value ierr', ierr, 0)
        call check_dp('value', value, -0.29533650641295134_dp)
        call check_dp('value_s', derivative(1), 0.006069515481673916_dp)
        call check_dp('value_t', derivative(2), -4.6960480732495265e-8_dp)
        call check_dp('value_phi', derivative(3), 0.0_dp)
    end subroutine check_value_oracle

    subroutine check_model303_field_oracle(data)
        type(jorek_restart_t), intent(in) :: data

        real(dp) :: b_r_z_phi(3)
        integer :: ierr

        call evaluate_jorek_model303_b(data, 1, 0.2_dp, 0.7_dp, 0.4_dp, &
            b_r_z_phi, ierr)
        call check_int('field ierr', ierr, 0)
        call check_dp('B_R', b_r_z_phi(1), -0.03213302883228271_dp)
        call check_dp('B_Z', b_r_z_phi(2), 0.009183002039748268_dp)
        call check_dp('B_phi', b_r_z_phi(3), 2.4859760851386405_dp)
    end subroutine check_model303_field_oracle

    subroutine check_int(name, actual, expected)
        character(len=*), intent(in) :: name
        integer, intent(in) :: actual, expected

        character(len=128) :: message

        if (actual /= expected) then
            write (message, '(A, I0, A, I0)') name//': got ', actual, &
                ', expected ', expected
            call fail(trim(message))
        end if
    end subroutine check_int

    subroutine check_dp(name, actual, expected)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: actual, expected

        if (abs(actual - expected) > tol*max(1.0_dp, abs(expected))) then
            call fail(name//' mismatch')
        end if
    end subroutine check_dp

    subroutine check_shape(name, actual, expected)
        character(len=*), intent(in) :: name
        integer, intent(in) :: actual(:), expected(:)

        if (size(actual) /= size(expected)) then
            call fail(name//': wrong rank')
        else if (any(actual /= expected)) then
            call fail(name//': wrong shape')
        end if
    end subroutine check_shape

    subroutine fail(message)
        character(len=*), intent(in) :: message

        print *, '    ', message
        nfail = nfail + 1
    end subroutine fail

end program test_jorek_restart_fixture
