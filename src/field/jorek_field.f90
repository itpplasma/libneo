module neo_jorek_field
use, intrinsic :: iso_fortran_env, only: dp => real64
use neo_field_mesh, only: field_mesh_t
use neo_spline_field, only: spline_field_t
use interpolate, only: SplineData3D
use neo_mesh, only: mesh_t

implicit none


type, extends(spline_field_t) :: jorek_field_t
    type(SplineData3D) :: fluxfunction_spline
    contains
        procedure :: jorek_field_init
end type jorek_field_t


contains


subroutine jorek_field_init(self, jorek_filename)

    class(jorek_field_t), intent(out) :: self
    character(*), intent(in), optional :: jorek_filename

    type(field_mesh_t) :: field_mesh
    type(mesh_t) :: fluxfunction_mesh

    call load_field_mesh_from_jorek(jorek_filename, field_mesh)
    call load_fluxfunction_mesh_from_jorek(jorek_filename, fluxfunction_mesh)
    call self%spline_field_init(field_mesh)
end subroutine jorek_field_init

subroutine load_field_mesh_from_jorek(jorek_filename, field_mesh)
    character(*), intent(in) :: jorek_filename
    type(field_mesh_t), intent(out) :: field_mesh

    integer :: n_R, n_Z, n_phi, idx_R
    real(dp), dimension(:), allocatable :: R, Z, phi
    real(dp), dimension(:,:,:,:), allocatable :: values
    real(dp), dimension(:,:,:), allocatable :: A_R, A_Z, A_phi, A_3, B_R, B_Z, B_phi
    logical :: is_periodic(3)

    call get_grid_and_values_from_jorek(jorek_filename, R, phi, Z, values, n_R, n_phi, n_Z)
    allocate(A_R(n_R, n_phi, n_Z), A_Z(n_R, n_phi, n_Z), &
             A_phi(n_R, n_phi, n_Z), A_3(n_R, n_phi, n_Z), &
             B_R(n_R, n_phi, n_Z), B_Z(n_R, n_phi, n_Z), B_phi(n_R, n_phi, n_Z))

    A_R = values(1,:,:,:)
    A_Z = values(2,:,:,:)
    A_3 = values(3,:,:,:)
    do idx_R = 1, n_R
        A_phi(idx_R,:,:) = A_3(idx_R,:,:) / R(idx_R)
    end do
    B_R = values(12,:,:,:)
    B_Z = values(13,:,:,:)
    B_phi = values(14,:,:,:)

    is_periodic = [.false., .true., .false.]
    call field_mesh%A1%mesh_init(R, phi, Z, A_R, is_periodic)
    call field_mesh%A2%mesh_init(R, phi, Z, -A_phi, is_periodic)
    call field_mesh%A3%mesh_init(R, phi, Z, A_Z, is_periodic)
    call field_mesh%B1%mesh_init(R, phi, Z, B_R, is_periodic)
    call field_mesh%B2%mesh_init(R, phi, Z, -B_phi, is_periodic)
    call field_mesh%B3%mesh_init(R, phi, Z, B_Z, is_periodic)

    deallocate(A_R, A_Z, A_phi, A_3, B_R, B_Z, B_phi, values, R, Z, phi)
end subroutine load_field_mesh_from_jorek

subroutine load_fluxfunction_mesh_from_jorek(jorek_filename, fluxfunction_mesh)
    character(*), intent(in) :: jorek_filename
    type(mesh_t), intent(out) :: fluxfunction_mesh

    integer :: n_R, n_Z, n_phi
    real(dp), dimension(:), allocatable :: R, Z, phi
    real(dp), dimension(:,:,:,:), allocatable :: values
    real(dp), dimension(:,:,:), allocatable :: fluxfunction
    logical :: is_periodic(3)

    call get_grid_and_values_from_jorek(jorek_filename, R, phi, Z, values, n_R, n_phi, n_Z)
    allocate(fluxfunction(n_R, n_phi, n_Z))
             
    fluxfunction = values(11,:,:,:)

    is_periodic = [.false., .true., .false.]
    call fluxfunction_mesh%mesh_init(R, phi, Z, fluxfunction, is_periodic)
    deallocate(fluxfunction, values, R, Z, phi)
end subroutine load_fluxfunction_mesh_from_jorek


subroutine get_grid_and_values_from_jorek(jorek_filename, R, phi, Z, values, & 
                                                          n_R, n_phi, n_Z)
    use neo_field_mesh, only: linspace

    character(*), intent(in) :: jorek_filename
    real(dp), dimension(:), allocatable, intent(out) :: R, phi, Z
    real(dp), dimension(:,:,:,:), allocatable, intent(out) :: values
    integer, intent(out) :: n_R, n_Z, n_phi

    real(dp) :: Rmin, Rmax, Zmin, Zmax, phimin, phimax
    
    call read_dims_and_values_from_jorek(jorek_filename, n_R, n_phi, n_Z, values)
    allocate(R(n_R), phi(n_phi), Z(n_Z))
    call get_ranges_from_filename(Rmin, Rmax, Zmin, Zmax, phimin, phimax, &
                                jorek_filename)
    R = linspace(Rmin, Rmax, n_R)
    Z = linspace(Zmin, Zmax, n_Z)
    phi = linspace(phimin, phimax, n_phi)
end subroutine get_grid_and_values_from_jorek

subroutine read_dims_and_values_from_jorek(jorek_filename, n_R, n_phi, n_Z, values)
    use hdf5_tools, only: hid_t, h5_init, h5_open, h5_get, h5_close, h5_deinit

    character(*), intent(in) :: jorek_filename
    integer, intent(out) :: n_R, n_phi, n_Z
    real(dp), dimension(:,:,:,:), allocatable, intent(out) :: values

    integer, parameter :: n_var = 17
    integer(hid_t) :: file_id
    integer :: flipped_dims(3)
    real(dp), dimension(:,:,:,:), allocatable :: flipped_values

    call h5_init()
    call h5_open(jorek_filename, file_id)
    call h5_get(file_id, 'dim', flipped_dims)
    n_R = flipped_dims(3)
    n_Z = flipped_dims(2)
    n_phi = flipped_dims(1)
    allocate(flipped_values(n_phi,n_Z, n_R, n_var))
    call h5_get(file_id, 'values', flipped_values)
    call h5_close(file_id)
    call h5_deinit()
    call phi_Z_R_var2var_R_phi_Z(flipped_values, values)
end subroutine read_dims_and_values_from_jorek

subroutine phi_Z_R_var2var_R_phi_Z(flipped_array, array)
    real(dp), dimension(:,:,:,:), intent(in) :: flipped_array
    real(dp), dimension(:,:,:,:), allocatable :: array

    integer :: dims(4)
    integer :: n_var, n_R, n_phi, n_Z
    integer :: var_idx, R_idx, phi_idx, Z_idx

    dims = shape(flipped_array)
    n_phi = dims(1)
    n_Z = dims(2)
    n_R = dims(3)
    n_var = dims(4)

    allocate(array(n_var, n_R, n_phi, n_Z))

    do var_idx = 1, n_var
        do R_idx = 1, n_R
            do phi_idx = 1, n_phi
                do Z_idx = 1, n_Z
                    array(var_idx,R_idx,phi_idx,Z_idx) = &
                                            flipped_array(phi_idx,Z_idx,R_idx,var_idx)
                end do
            end do
        end do
    end do
end subroutine phi_Z_R_var2var_R_phi_Z

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