program test_jorek_restart
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use hdf5_tools, only: HID_T, h5_init, h5_deinit, h5_create, h5_close, h5_add
    use jorek_restart, only: jorek_restart_t, load_jorek_restart, &
        free_jorek_restart
    use util_for_test, only: print_test, print_ok, print_fail

    implicit none

    integer, parameter :: n_var = 2, n_degrees = 2, n_tor = 1
    integer, parameter :: n_coord_tor = 2, n_dim = 3, n_vertex_max = 4
    integer, parameter :: n_nodes = 3, n_elements = 2
    integer, parameter :: n_order = 1, n_period = 1, jorek_model = 303
    real(dp), parameter :: F0 = 1.5_dp
    real(dp), parameter :: t_norm = 2.5e-7_dp
    real(dp), parameter :: t_now = 0.125_dp

    integer :: nfail

    nfail = 0

    call test_round_trip
    call test_reload_after_free
    call test_reject_unsupported_version
    call test_reject_missing_dataset
    call test_reject_invalid_size
    call test_reject_missing_file

    if (nfail > 0) error stop

contains

    subroutine test_round_trip
        character(len=*), parameter :: filename = 'jorek_restart_synthetic.h5'

        type(jorek_restart_t) :: data
        integer :: ierr, nfail_before

        call print_test('load_jorek_restart round-trips a synthetic file')
        nfail_before = nfail

        call write_synthetic(filename, 2, .false.)
        call load_jorek_restart(filename, data, ierr)
        call check_int('ierr', ierr, 0)
        call check_scalars(data)
        call check_arrays(data)
        call free_jorek_restart(data)
        if (allocated(data%values)) call fail('free left values allocated')
        if (allocated(data%x)) call fail('free left x allocated')
        if (allocated(data%vertex)) call fail('free left vertex allocated')
        if (allocated(data%size)) call fail('free left size allocated')
        if (allocated(data%neighbours)) call fail('free left neighbours allocated')

        call report(nfail_before)
    end subroutine test_round_trip

    subroutine test_reload_after_free
        character(len=*), parameter :: filename = 'jorek_restart_synthetic.h5'

        type(jorek_restart_t) :: data
        integer :: ierr, nfail_before

        call print_test('second load/free cycle succeeds')
        nfail_before = nfail

        call load_jorek_restart(filename, data, ierr)
        call check_int('first ierr', ierr, 0)
        call free_jorek_restart(data)
        call load_jorek_restart(filename, data, ierr)
        call check_int('second ierr', ierr, 0)
        call check_dp('values(3,1,2,2)', data%values(3, 1, 2, 2), &
            values_ref(3, 2, 2))
        call free_jorek_restart(data)

        call report(nfail_before)
    end subroutine test_reload_after_free

    subroutine test_reject_unsupported_version
        character(len=*), parameter :: filename = 'jorek_restart_bad_version.h5'

        type(jorek_restart_t) :: data
        integer :: ierr, nfail_before

        call print_test('load_jorek_restart rejects rst_hdf5_version /= 2')
        nfail_before = nfail

        call write_synthetic(filename, 1, .false.)
        call load_jorek_restart(filename, data, ierr)
        call check_int('ierr', ierr, 2)
        if (allocated(data%values)) call fail('rejected load allocated values')

        call report(nfail_before)
    end subroutine test_reject_unsupported_version

    subroutine test_reject_missing_dataset
        character(len=*), parameter :: filename = 'jorek_restart_no_values.h5'

        type(jorek_restart_t) :: data
        integer :: ierr, nfail_before

        call print_test('load_jorek_restart reports a missing dataset')
        nfail_before = nfail

        call write_synthetic(filename, 2, .true.)
        call load_jorek_restart(filename, data, ierr)
        call check_int('ierr', ierr, 3)

        call report(nfail_before)
    end subroutine test_reject_missing_dataset

    subroutine test_reject_invalid_size
        character(len=*), parameter :: filename = 'jorek_restart_bad_size.h5'

        type(jorek_restart_t) :: data
        integer :: ierr, nfail_before

        call print_test('load_jorek_restart rejects nonpositive dimensions')
        nfail_before = nfail

        call write_synthetic(filename, 2, .false., 0)
        call load_jorek_restart(filename, data, ierr)
        call check_int('ierr', ierr, 4)
        if (allocated(data%values)) call fail('rejected load allocated values')

        call report(nfail_before)
    end subroutine test_reject_invalid_size

    subroutine test_reject_missing_file
        type(jorek_restart_t) :: data
        integer :: ierr, nfail_before

        call print_test('load_jorek_restart reports a missing file')
        nfail_before = nfail

        call load_jorek_restart('jorek_restart_does_not_exist.h5', data, ierr)
        call check_int('ierr', ierr, 1)

        call report(nfail_before)
    end subroutine test_reject_missing_file

    subroutine write_synthetic(filename, version, omit_values, stored_n_nodes)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: version
        logical, intent(in) :: omit_values
        integer, intent(in), optional :: stored_n_nodes

        integer(HID_T) :: h5id
        real(dp) :: values(n_nodes, n_tor, n_degrees, n_var)
        real(dp) :: x(n_nodes, n_coord_tor, n_degrees, n_dim)
        real(dp) :: elem_size(n_elements, n_vertex_max, n_degrees)
        integer :: vertex(n_elements, n_vertex_max)
        integer :: neighbours(n_elements, n_vertex_max)
        integer :: file_n_nodes, i, j, k, l

        file_n_nodes = n_nodes
        if (present(stored_n_nodes)) file_n_nodes = stored_n_nodes

        do l = 1, n_var
            do k = 1, n_degrees
                do i = 1, n_nodes
                    values(i, 1, k, l) = values_ref(i, k, l)
                end do
            end do
        end do
        do l = 1, n_dim
            do k = 1, n_degrees
                do j = 1, n_coord_tor
                    do i = 1, n_nodes
                        x(i, j, k, l) = x_ref(i, j, k, l)
                    end do
                end do
            end do
        end do
        do l = 1, n_degrees
            do k = 1, n_vertex_max
                do i = 1, n_elements
                    elem_size(i, k, l) = size_ref(i, k, l)
                end do
            end do
        end do
        vertex = reshape([1, 2, 2, 3, 3, 1, 1, 2], &
            [n_elements, n_vertex_max])
        neighbours = reshape([2, 1, 0, -1, -1, -1, -1, -1], &
            [n_elements, n_vertex_max])

        call h5_init()
        call h5_create(filename, h5id)
        call h5_add(h5id, 'rst_hdf5_version', version)
        call h5_add(h5id, 'n_var', n_var)
        call h5_add(h5id, 'n_degrees', n_degrees)
        call h5_add(h5id, 'n_tor', n_tor)
        call h5_add(h5id, 'n_coord_tor', n_coord_tor)
        call h5_add(h5id, 'n_dim', n_dim)
        call h5_add(h5id, 'n_nodes', file_n_nodes)
        call h5_add(h5id, 'n_elements', n_elements)
        call h5_add(h5id, 'n_vertex_max', n_vertex_max)
        call h5_add(h5id, 'n_order', n_order)
        call h5_add(h5id, 'n_period', n_period)
        call h5_add(h5id, 'jorek_model', jorek_model)
        call h5_add(h5id, 'F0', F0)
        call h5_add(h5id, 't_norm', t_norm)
        call h5_add(h5id, 't_now', t_now)
        if (.not. omit_values) then
            call h5_add(h5id, 'values', values, lbounds=[1, 1, 1, 1], &
                ubounds=shape(values))
        end if
        call h5_add(h5id, 'x', x, lbounds=[1, 1, 1, 1], ubounds=shape(x))
        call h5_add(h5id, 'vertex', vertex, lbounds=[1, 1], &
            ubounds=shape(vertex))
        call h5_add(h5id, 'size', elem_size, lbounds=[1, 1, 1], &
            ubounds=shape(elem_size))
        call h5_add(h5id, 'neighbours', neighbours, lbounds=[1, 1], &
            ubounds=shape(neighbours))
        call h5_close(h5id)
        call h5_deinit()
    end subroutine write_synthetic

    subroutine check_scalars(data)
        type(jorek_restart_t), intent(in) :: data

        call check_int('rst_hdf5_version', data%rst_hdf5_version, 2)
        call check_int('n_var', data%n_var, n_var)
        call check_int('n_degrees', data%n_degrees, n_degrees)
        call check_int('n_tor', data%n_tor, n_tor)
        call check_int('n_coord_tor', data%n_coord_tor, n_coord_tor)
        call check_int('n_dim', data%n_dim, n_dim)
        call check_int('n_nodes', data%n_nodes, n_nodes)
        call check_int('n_elements', data%n_elements, n_elements)
        call check_int('n_vertex_max', data%n_vertex_max, n_vertex_max)
        call check_int('n_order', data%n_order, n_order)
        call check_int('n_period', data%n_period, n_period)
        call check_int('jorek_model', data%jorek_model, jorek_model)
        call check_dp('F0', data%F0, F0)
        call check_dp('t_norm', data%t_norm, t_norm)
        call check_dp('t_now', data%t_now, t_now)
    end subroutine check_scalars

    subroutine check_arrays(data)
        type(jorek_restart_t), intent(in) :: data

        integer :: i, j, k, l
        character(len=64) :: label

        call check_shape('values', shape(data%values), &
            [n_nodes, n_tor, n_degrees, n_var])
        call check_shape('x', shape(data%x), &
            [n_nodes, n_coord_tor, n_degrees, n_dim])
        call check_shape('vertex', shape(data%vertex), &
            [n_elements, n_vertex_max])
        call check_shape('size', shape(data%size), &
            [n_elements, n_vertex_max, n_degrees])
        call check_shape('neighbours', shape(data%neighbours), &
            [n_elements, n_vertex_max])

        do l = 1, n_var
            do k = 1, n_degrees
                do i = 1, n_nodes
                    write (label, '(A, 3(I0, A))') 'values(', i, ',1,', k, ',', l, ')'
                    call check_dp(trim(label), data%values(i, 1, k, l), &
                        values_ref(i, k, l))
                end do
            end do
        end do
        do l = 1, n_dim
            do k = 1, n_degrees
                do j = 1, n_coord_tor
                    do i = 1, n_nodes
                        write (label, '(A, 4(I0, A))') 'x(', i, ',', j, ',', &
                            k, ',', l, ')'
                        call check_dp(trim(label), data%x(i, j, k, l), &
                            x_ref(i, j, k, l))
                    end do
                end do
            end do
        end do
        do l = 1, n_degrees
            do k = 1, n_vertex_max
                do i = 1, n_elements
                    write (label, '(A, 3(I0, A))') 'size(', i, ',', k, ',', l, ')'
                    call check_dp(trim(label), data%size(i, k, l), &
                        size_ref(i, k, l))
                end do
            end do
        end do
        do k = 1, n_vertex_max
            do i = 1, n_elements
                write (label, '(A, 2(I0, A))') 'vertex(', i, ',', k, ')'
                call check_int(trim(label), data%vertex(i, k), &
                    vertex_ref(i, k))
                write (label, '(A, 2(I0, A))') 'neighbours(', i, ',', k, ')'
                call check_int(trim(label), data%neighbours(i, k), &
                    neighbours_ref(i, k))
            end do
        end do
    end subroutine check_arrays

    pure real(dp) function values_ref(i, k, l)
        integer, intent(in) :: i, k, l

        values_ref = real(i + 10*k + 100*l, dp)/7.0_dp
    end function values_ref

    pure real(dp) function x_ref(i, j, k, l)
        integer, intent(in) :: i, j, k, l

        x_ref = real(i + 10*j + 100*k + 1000*l, dp)/3.0_dp
    end function x_ref

    pure real(dp) function size_ref(i, k, l)
        integer, intent(in) :: i, k, l

        size_ref = real(i + 10*k + 100*l, dp)/11.0_dp
    end function size_ref

    pure integer function vertex_ref(i, k)
        integer, intent(in) :: i, k

        integer, parameter :: ref(n_elements, n_vertex_max) = &
            reshape([1, 2, 2, 3, 3, 1, 1, 2], [n_elements, n_vertex_max])

        vertex_ref = ref(i, k)
    end function vertex_ref

    pure integer function neighbours_ref(i, k)
        integer, intent(in) :: i, k

        integer, parameter :: ref(n_elements, n_vertex_max) = &
            reshape([2, 1, 0, -1, -1, -1, -1, -1], &
            [n_elements, n_vertex_max])

        neighbours_ref = ref(i, k)
    end function neighbours_ref

    subroutine check_shape(name, actual, expected)
        character(len=*), intent(in) :: name
        integer, intent(in) :: actual(:), expected(:)

        if (size(actual) /= size(expected)) then
            call fail(name//': wrong rank')
        else if (any(actual /= expected)) then
            call fail(name//': wrong shape')
        end if
    end subroutine check_shape

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

        character(len=160) :: message

        if (actual /= expected) then
            write (message, '(A, ES23.16, A, ES23.16)') name//': got ', &
                actual, ', expected ', expected
            call fail(trim(message))
        end if
    end subroutine check_dp

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

end program test_jorek_restart
