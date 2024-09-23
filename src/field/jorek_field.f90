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

    type(field_mesh_t) :: mesh

    call load_field_mesh_from_jorek(jorek_filename, mesh)
    call self%spline_field_init(mesh)

end subroutine jorek_field_init

subroutine load_field_mesh_from_jorek(jorek_filename, mesh)
    use hdf5_tools, only: hid_t, h5_init, h5_open, h5_get, h5_deinit

    character(*), intent(in) :: jorek_filename
    type(field_mesh_t), intent(out) :: mesh

    integer(hid_t) :: file_id
    integer :: n_R, n_Z, n_phi
    integer :: dims(3)
    real(dp), dimension(:,:,:,:), allocatable :: values
    real(dp), dimension(:,:,:), allocatable :: Ar, Az, Aphi, A3, Br, Bz, Bphi

    call h5_init()
    call h5_open(jorek_filename, file_id)
    call h5_get(file_id, 'dim', dims)
    call h5_get(file_id, 'values', values)

    n_phi = dims(1)
    n_Z = dims(2)
    n_R = dims(3)

    allocate(Ar(n_phi, n_Z, n_R), Az(n_phi, n_Z, n_R), Aphi(n_phi, n_Z, n_R), &
             A3(n_phi, n_Z, n_R), Br(n_phi, n_Z, n_R), Bz(n_phi, n_Z, n_R), &
             Bphi(n_phi, n_Z, n_R))

    Ar = values(:,:,:,1)
    Az = values(:,:,:,2)
    A3 = values(:,:,:,3)
    Br = values(:,:,:,12)
    Bz = values(:,:,:,13)
    Bphi = values(:,:,:,14)
end subroutine load_field_mesh_from_jorek

subroutine get_ranges_from_filename(Rmin, Rmax, Zmin, Zmax, phimin, phimax, filename)
    character(*), intent(in) :: filename
    real(dp), intent(out) :: Rmin, Rmax, Zmin, Zmax, phimin, phimax

    integer, parameter :: ndigits = 4
    character(len=20) :: temp
    integer :: pos

    pos = index(filename, 'Rmin')
    read(filename(pos+len('Rmin')+1:pos+len('Rmin')+1+ndigits), *) Rmin

    pos = index(filename, 'Rmax')
    read(filename(pos+len('Rmax')+1:pos+len('Rmax')+1+ndigits), *) Rmax

    pos = index(filename, 'Zmin')
    read(filename(pos+len('Zmin')+1:pos+len('Zmin')+1+ndigits), *) Zmin

    pos = index(filename, 'Zmax')
    read(filename(pos+len('Zmax')+1:pos+len('Zmax')+1+ndigits), *) Zmax

    pos = index(filename, 'phimin')
    read(filename(pos+len('phimin')+1:pos+len('phimin')+1+ndigits), *) phimin

    pos = index(filename, 'phimax')
    read(filename(pos+len('phimax')+1:pos+len('phimax')+1+ndigits), *) phimax
end subroutine get_ranges_from_filename

end module neo_jorek_field