module neo_jorek_field
use, intrinsic :: iso_fortran_env, only: dp => real64
use neo_field_mesh, only: field_mesh_t
use neo_spline_field, only: spline_field_t

implicit none


type, extends(spline_field_t) :: jorek_field_t
    contains
        procedure :: jorek_field_init
end type jorek_field_t


contains


subroutine jorek_field_init(self, jorek_filename)

    class(jorek_field_t), intent(out) :: self
    character(*), intent(in), optional :: jorek_filename

    type(field_mesh_t) :: field_mesh

    call load_field_mesh_from_jorek(jorek_filename, field_mesh)
    call self%spline_field_init(field_mesh)

end subroutine jorek_field_init

subroutine load_field_mesh_from_jorek(jorek_filename, field_mesh)
    use neo_field_mesh, only: linspace

    character(*), intent(in) :: jorek_filename
    type(field_mesh_t), intent(out) :: field_mesh

    integer :: n_R, n_Z, n_phi, idx_R
    integer :: dims(3)
    real(dp) :: Rmin, Rmax, Zmin, Zmax, phimin, phimax
    real(dp), dimension(:), allocatable :: R, Z, phi
    real(dp), dimension(:,:,:,:), allocatable :: values
    real(dp), dimension(:,:,:), allocatable :: A_R, A_Z, A_phi, A_3, B_R, B_Z, B_phi
    logical :: is_periodic(3)

    
    call read_dims_and_values_from_jorek(jorek_filename, dims, values)
    n_R = dims(1)
    n_Z = dims(2)
    n_phi = dims(3)
    allocate(R(n_R), Z(n_Z), phi(n_phi))

    call get_ranges_from_filename(Rmin, Rmax, Zmin, Zmax, phimin, phimax, &
                                  jorek_filename)
    R = linspace(Rmin, Rmax, n_R)
    Z = linspace(Zmin, Zmax, n_Z)
    phi = linspace(phimin, phimax, n_phi)

    allocate(A_R(n_R, n_Z, n_phi), A_Z(n_R, n_Z, n_phi), &
             A_phi(n_R, n_Z, n_phi), A_3(n_R, n_Z, n_phi), &
             B_R(n_R, n_Z, n_phi), B_Z(n_R, n_Z, n_phi), B_phi(n_R, n_Z, n_phi))


    A_R = values(1,:,:,:)
    A_Z = values(2,:,:,:)
    A_3 = values(3,:,:,:)
    do idx_R = 1, n_R
        A_phi(idx_R,:,:) = A_3(idx_R,:,:) * R(idx_R)
    end do
    B_R = values(12,:,:,:)
    B_Z = values(13,:,:,:)
    B_phi = values(14,:,:,:)

    is_periodic = [.false., .false., .true.]
    call field_mesh%A1%mesh_init(R, Z, phi, A_R, is_periodic)
    call field_mesh%A2%mesh_init(R, Z, phi, A_Z, is_periodic)
    call field_mesh%A3%mesh_init(R, Z, phi, A_phi, is_periodic)
    call field_mesh%B1%mesh_init(R, Z, phi, B_R, is_periodic)
    call field_mesh%B2%mesh_init(R, Z, phi, B_Z, is_periodic)
    call field_mesh%B3%mesh_init(R, Z, phi, B_phi, is_periodic)

    deallocate(A_R, A_Z, A_phi, A_3, B_R, B_Z, B_phi, values)

end subroutine load_field_mesh_from_jorek

subroutine read_dims_and_values_from_jorek(jorek_filename, dims, values)
    use hdf5_tools, only: hid_t, h5_init, h5_open, h5_get, h5_close, h5_deinit

    character(*), intent(in) :: jorek_filename
    integer, intent(out) :: dims(3)
    real(dp), dimension(:,:,:,:), allocatable, intent(out) :: values

    integer, parameter :: n_variables = 17
    integer(hid_t) :: file_id
    integer :: flipped_dims(3)
    real(dp), dimension(:,:,:,:), allocatable :: flipped_values

    call h5_init()
    call h5_open(jorek_filename, file_id)
    call h5_get(file_id, 'dim', flipped_dims)
    allocate(flipped_values(flipped_dims(1), flipped_dims(2), &
                            flipped_dims(3), n_variables))
    call h5_get(file_id, 'values', flipped_values)
    call h5_close(file_id)
    call h5_deinit()
    call unflip_array_dims(flipped_values, values)
    dims(1) = flipped_dims(3)
    dims(2) = flipped_dims(2)
    dims(3) = flipped_dims(1)
end subroutine read_dims_and_values_from_jorek

subroutine unflip_array_dims(flipped_array, array)
    real(dp), dimension(:,:,:,:), intent(in) :: flipped_array
    real(dp), dimension(:,:,:,:), allocatable :: array

    integer :: dims(4)
    integer :: i, j, k, l

    dims = shape(flipped_array)
    allocate(array(dims(4), dims(3), dims(2), dims(1)))

    do i = 1, dims(4)
        do j = 1, dims(3)
            do k = 1, dims(2)
                do l = 1, dims(1)
                    array(i,j,k,l) = flipped_array(l,k,j,i)
                end do
            end do
        end do
    end do
end subroutine unflip_array_dims

subroutine get_ranges_from_filename(Rmin, Rmax, Zmin, Zmax, phimin, phimax, filename)
    character(*), intent(in) :: filename
    real(dp), intent(out) :: Rmin, Rmax, Zmin, Zmax, phimin, phimax

    integer, parameter :: ndigits = 4
    integer :: pos

    pos = index(filename, 'Rmin')
    read(filename(pos+len('Rmin'):pos+len('Rmin')+ndigits), *) Rmin

    pos = index(filename, 'Rmax')
    read(filename(pos+len('Rmax'):pos+len('Rmax')+ndigits), *) Rmax

    pos = index(filename, 'Zmin')
    read(filename(pos+len('Zmin'):pos+len('Zmin')+ndigits), *) Zmin

    pos = index(filename, 'Zmax')
    read(filename(pos+len('Zmax'):pos+len('Zmax')+ndigits), *) Zmax

    pos = index(filename, 'phimin')
    read(filename(pos+len('phimin'):pos+len('phimin')+ndigits), *) phimin

    pos = index(filename, 'phimax')
    read(filename(pos+len('phimax'):pos+len('phimax')+ndigits), *) phimax
end subroutine get_ranges_from_filename


end module neo_jorek_field