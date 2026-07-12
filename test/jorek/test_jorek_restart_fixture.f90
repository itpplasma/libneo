program test_jorek_restart_fixture
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use jorek_restart, only: jorek_restart_t, load_jorek_restart, &
        free_jorek_restart
    use util_for_test, only: print_test, print_ok, print_fail

    implicit none

    character(len=*), parameter :: fixture_env = 'LIBNEO_JOREK_RESTART_FIXTURE'
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
        call check_int('n_degrees', data%n_degrees, 4)
        call check_int('n_dim', data%n_dim, 2)
        call check_int('n_vertex_max', data%n_vertex_max, 4)
        call check_int('n_order', data%n_order, 3)

        call check_array_shapes(data)
        call check_nodal_r_range(data)
        call check_vertex_bounds(data)
        call check_neighbour_bounds(data)

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
