module jorek_restart
    !> Reader for JOREK restart files (rst_hdf5_version = 2).
    !>
    !> The restart layout is flat: no groups, no attributes, scalars stored as
    !> one-element datasets. Arrays keep the Fortran dimension order of the
    !> JOREK writer, so the C-order file shape (n_var, n_degrees, n_tor,
    !> n_nodes) reads back as (n_nodes, n_tor, n_degrees, n_var). Node indices
    !> in vertex stay 1-based as stored; neighbours keeps nonpositive boundary
    !> sentinels as stored (-1 at a collapsed axis, 0 at a physical boundary).
    !> Values keep JOREK normalized units (t_norm = sqrt(mu0 rho0)).

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use hdf5_tools, only: HID_T, h5_init, h5_deinit, h5_open, h5_close, h5_get, &
        h5_exists

    implicit none
    private

    public :: jorek_restart_t, load_jorek_restart, free_jorek_restart

    integer, parameter :: supported_rst_hdf5_version = 2

    type :: jorek_restart_t
        integer :: rst_hdf5_version = 0
        integer :: n_var = 0
        integer :: n_degrees = 0
        integer :: n_tor = 0
        integer :: n_coord_tor = 0
        integer :: n_dim = 0
        integer :: n_nodes = 0
        integer :: n_elements = 0
        integer :: n_vertex_max = 0
        integer :: n_order = 0
        integer :: n_period = 0
        integer :: jorek_model = 0
        real(dp) :: F0 = 0.0_dp
        real(dp) :: t_norm = 0.0_dp
        real(dp) :: t_now = 0.0_dp
        !> Bezier value DOFs, (n_nodes, n_tor, n_degrees, n_var)
        real(dp), allocatable :: values(:, :, :, :)
        !> Bezier geometry DOFs, (n_nodes, n_coord_tor, n_degrees, n_dim)
        real(dp), allocatable :: x(:, :, :, :)
        !> Element-to-node map, (n_elements, n_vertex_max), 1-based node indices
        integer, allocatable :: vertex(:, :)
        !> Element scalings, (n_elements, n_vertex_max, n_degrees)
        real(dp), allocatable :: size(:, :, :)
        !> Element adjacency; nonpositive values mark boundary sides as stored
        integer, allocatable :: neighbours(:, :)
    end type jorek_restart_t

contains

    subroutine load_jorek_restart(filename, data, ierr)
        character(len=*), intent(in) :: filename
        type(jorek_restart_t), intent(out) :: data
        integer, intent(out) :: ierr

        integer(HID_T) :: h5id
        logical :: file_exists

        ierr = 0
        inquire (file=filename, exist=file_exists)
        if (.not. file_exists) then
            write (*, '(A)') 'load_jorek_restart: file not found: '//trim(filename)
            ierr = 1
            return
        end if

        call h5_init()
        call h5_open(filename, h5id)

        call require_datasets(h5id, ierr)
        if (ierr /= 0) then
            call h5_close(h5id)
            call h5_deinit()
            return
        end if

        call h5_get(h5id, 'rst_hdf5_version', data%rst_hdf5_version)
        if (data%rst_hdf5_version /= supported_rst_hdf5_version) then
            write (*, '(A, I0, A, I0)') &
                'load_jorek_restart: unsupported rst_hdf5_version = ', &
                data%rst_hdf5_version, ', expected ', supported_rst_hdf5_version
            ierr = 2
            call h5_close(h5id)
            call h5_deinit()
            return
        end if

        call read_scalars(h5id, data)
        call validate_sizes(data, ierr)
        if (ierr /= 0) then
            call h5_close(h5id)
            call h5_deinit()
            return
        end if
        call read_arrays(h5id, data)

        call h5_close(h5id)
        call h5_deinit()
    end subroutine load_jorek_restart

    subroutine free_jorek_restart(data)
        type(jorek_restart_t), intent(inout) :: data
        type(jorek_restart_t) :: empty

        data = empty
    end subroutine free_jorek_restart

    subroutine require_datasets(h5id, ierr)
        integer(HID_T), intent(in) :: h5id
        integer, intent(out) :: ierr

        character(len=16), parameter :: required(20) = [character(len=16) :: &
            'rst_hdf5_version', 'n_var', 'n_degrees', 'n_tor', 'n_coord_tor', &
            'n_dim', 'n_nodes', 'n_elements', 'n_vertex_max', 'n_order', &
            'n_period', 'jorek_model', 'F0', 't_norm', 't_now', 'values', 'x', &
            'vertex', 'size', 'neighbours']

        integer :: i

        ierr = 0
        do i = 1, size(required)
            if (.not. h5_exists(h5id, trim(required(i)))) then
                write (*, '(A)') 'load_jorek_restart: missing dataset: ' &
                    //trim(required(i))
                ierr = 3
            end if
        end do
    end subroutine require_datasets

    subroutine read_scalars(h5id, data)
        integer(HID_T), intent(in) :: h5id
        type(jorek_restart_t), intent(inout) :: data

        call h5_get(h5id, 'n_var', data%n_var)
        call h5_get(h5id, 'n_degrees', data%n_degrees)
        call h5_get(h5id, 'n_tor', data%n_tor)
        call h5_get(h5id, 'n_coord_tor', data%n_coord_tor)
        call h5_get(h5id, 'n_dim', data%n_dim)
        call h5_get(h5id, 'n_nodes', data%n_nodes)
        call h5_get(h5id, 'n_elements', data%n_elements)
        call h5_get(h5id, 'n_vertex_max', data%n_vertex_max)
        call h5_get(h5id, 'n_order', data%n_order)
        call h5_get(h5id, 'n_period', data%n_period)
        call h5_get(h5id, 'jorek_model', data%jorek_model)
        call h5_get(h5id, 'F0', data%F0)
        call h5_get(h5id, 't_norm', data%t_norm)
        call h5_get(h5id, 't_now', data%t_now)
    end subroutine read_scalars

    subroutine validate_sizes(data, ierr)
        type(jorek_restart_t), intent(in) :: data
        integer, intent(out) :: ierr

        integer :: sizes(8)

        sizes = [data%n_var, data%n_degrees, data%n_tor, data%n_coord_tor, &
            data%n_dim, data%n_nodes, data%n_elements, data%n_vertex_max]
        ierr = 0
        if (any(sizes <= 0)) then
            write (*, '(A, 8(1X, I0))') &
                'load_jorek_restart: nonpositive array dimension:', sizes
            ierr = 4
        end if
    end subroutine validate_sizes

    subroutine read_arrays(h5id, data)
        integer(HID_T), intent(in) :: h5id
        type(jorek_restart_t), intent(inout) :: data

        allocate (data%values(data%n_nodes, data%n_tor, data%n_degrees, &
            data%n_var))
        allocate (data%x(data%n_nodes, data%n_coord_tor, data%n_degrees, &
            data%n_dim))
        allocate (data%vertex(data%n_elements, data%n_vertex_max))
        allocate (data%size(data%n_elements, data%n_vertex_max, &
            data%n_degrees))
        allocate (data%neighbours(data%n_elements, data%n_vertex_max))

        call h5_get(h5id, 'values', data%values)
        call h5_get(h5id, 'x', data%x)
        call h5_get(h5id, 'vertex', data%vertex)
        call h5_get(h5id, 'size', data%size)
        call h5_get(h5id, 'neighbours', data%neighbours)
    end subroutine read_arrays

end module jorek_restart
