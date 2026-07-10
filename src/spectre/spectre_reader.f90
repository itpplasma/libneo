module spectre_reader
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use hdf5_tools, only: HID_T, h5_init, h5_deinit, h5_open, h5_close, h5_get

    implicit none
    private

    public :: spectre_volume_t, spectre_data_t, load_spectre, free_spectre

    integer, parameter :: igeometry_toroidal = 3

    type :: spectre_volume_t
        real(dp), allocatable :: Ate(:, :), Aze(:, :), Ato(:, :), Azo(:, :)
    end type spectre_volume_t

    type :: spectre_data_t
        integer :: Igeometry = 0
        integer :: Nfp = 0
        integer :: Nvol = 0
        integer :: Mvol = 0
        integer :: Mpol = 0
        integer :: Ntor = 0
        integer :: mn = 0
        integer, allocatable :: im(:), in(:), Lrad(:)
        real(dp), allocatable :: mu(:), pressure(:), tflux(:), pflux(:)
        real(dp), allocatable :: Rbc(:, :), Zbs(:, :), Rbs(:, :), Zbc(:, :)
        type(spectre_volume_t), allocatable :: vol(:)
    end type spectre_data_t

contains

    subroutine load_spectre(filename, data, ierr)
        character(len=*), intent(in) :: filename
        type(spectre_data_t), intent(out) :: data
        integer, intent(out) :: ierr

        integer(HID_T) :: h5id
        logical :: file_exists

        ierr = 0
        inquire (file=filename, exist=file_exists)
        if (.not. file_exists) then
            write (*, '(A)') 'load_spectre: file not found: '//trim(filename)
            ierr = 1
            return
        end if

        call h5_init()
        call h5_open(filename, h5id)
        call h5_get(h5id, 'input/physics/Igeometry', data%Igeometry)
        if (data%Igeometry /= igeometry_toroidal) then
            write (*, '(A, I0)') 'load_spectre: unsupported Igeometry = ', &
                data%Igeometry
            ierr = 2
            call h5_close(h5id)
            call h5_deinit()
            return
        end if

        call read_sizes(h5id, data)
        call read_profiles(h5id, data)
        call read_interface_geometry(h5id, data)
        call read_vector_potential(h5id, data)

        call h5_close(h5id)
        call h5_deinit()
    end subroutine load_spectre

    subroutine free_spectre(data)
        type(spectre_data_t), intent(inout) :: data

        data = spectre_data_t()
    end subroutine free_spectre

    subroutine read_sizes(h5id, data)
        integer(HID_T), intent(in) :: h5id
        type(spectre_data_t), intent(inout) :: data

        call h5_get(h5id, 'input/physics/Nfp', data%Nfp)
        call h5_get(h5id, 'input/physics/Nvol', data%Nvol)
        call h5_get(h5id, 'input/physics/Mpol', data%Mpol)
        call h5_get(h5id, 'input/physics/Ntor', data%Ntor)
        call h5_get(h5id, 'output/Mvol', data%Mvol)
        call h5_get(h5id, 'output/mn', data%mn)
    end subroutine read_sizes

    subroutine read_profiles(h5id, data)
        integer(HID_T), intent(in) :: h5id
        type(spectre_data_t), intent(inout) :: data

        allocate (data%Lrad(data%Mvol))
        allocate (data%im(data%mn), data%in(data%mn))
        allocate (data%mu(data%Mvol))
        allocate (data%pressure(data%Nvol))
        allocate (data%tflux(data%Mvol), data%pflux(data%Mvol))

        call h5_get(h5id, 'input/physics/Lrad', data%Lrad)
        call h5_get(h5id, 'output/im', data%im)
        call h5_get(h5id, 'output/in', data%in)
        call h5_get(h5id, 'output/mu', data%mu)
        call h5_get(h5id, 'input/physics/pressure', data%pressure)
        call h5_get(h5id, 'input/physics/tflux', data%tflux)
        call h5_get(h5id, 'input/physics/pflux', data%pflux)
    end subroutine read_profiles

    subroutine read_interface_geometry(h5id, data)
        integer(HID_T), intent(in) :: h5id
        type(spectre_data_t), intent(inout) :: data

        allocate (data%Rbc(data%mn, 0:data%Mvol), data%Zbs(data%mn, 0:data%Mvol))
        allocate (data%Rbs(data%mn, 0:data%Mvol), data%Zbc(data%mn, 0:data%Mvol))

        call h5_get(h5id, 'output/Rbc', data%Rbc)
        call h5_get(h5id, 'output/Zbs', data%Zbs)
        call h5_get(h5id, 'output/Rbs', data%Rbs)
        call h5_get(h5id, 'output/Zbc', data%Zbc)
    end subroutine read_interface_geometry

    subroutine read_vector_potential(h5id, data)
        integer(HID_T), intent(in) :: h5id
        type(spectre_data_t), intent(inout) :: data

        integer :: lvol
        real(dp), allocatable :: buf(:, :)

        allocate (data%vol(data%Mvol))
        allocate (buf(sum(data%Lrad + 1), data%mn))

        ! The writer concatenates per-volume radial blocks of Lrad(lvol)+1 rows
        ! along the first dimension, in volume order.
        call h5_get(h5id, 'vector_potential/Ate', buf)
        do lvol = 1, data%Mvol
            call unpack_volume_block(buf, data%Lrad, lvol, data%vol(lvol)%Ate)
        end do
        call h5_get(h5id, 'vector_potential/Aze', buf)
        do lvol = 1, data%Mvol
            call unpack_volume_block(buf, data%Lrad, lvol, data%vol(lvol)%Aze)
        end do
        call h5_get(h5id, 'vector_potential/Ato', buf)
        do lvol = 1, data%Mvol
            call unpack_volume_block(buf, data%Lrad, lvol, data%vol(lvol)%Ato)
        end do
        call h5_get(h5id, 'vector_potential/Azo', buf)
        do lvol = 1, data%Mvol
            call unpack_volume_block(buf, data%Lrad, lvol, data%vol(lvol)%Azo)
        end do
    end subroutine read_vector_potential

    subroutine unpack_volume_block(buf, Lrad, lvol, coeff)
        real(dp), intent(in) :: buf(:, :)
        integer, intent(in) :: Lrad(:)
        integer, intent(in) :: lvol
        real(dp), allocatable, intent(out) :: coeff(:, :)

        integer :: row

        row = sum(Lrad(1:lvol - 1) + 1)
        allocate (coeff(0:Lrad(lvol), size(buf, 2)))
        coeff = buf(row + 1:row + Lrad(lvol) + 1, :)
    end subroutine unpack_volume_block

end module spectre_reader
