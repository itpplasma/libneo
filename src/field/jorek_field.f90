module neo_jorek_field
use, intrinsic :: iso_fortran_env, only: dp => real64
use neo_field_mesh, only: field_mesh_t
use neo_field_base, only: field_t
use neo_spline_field, only: spline_field_t, make_spline_from_mesh
use interpolate, only: SplineData3D, SplineData2D
use interpolate, only: evaluate_splines_3d, evaluate_splines_3d_der, evaluate_splines_2d
use neo_mesh, only: mesh_t
use neo_mesh_2d, only: mesh_2d_t

implicit none


type, extends(field_t) :: jorek_field_t
    type(SplineData2D) :: fluxfunction_spline
    type(SplineData3D) :: A1_spline
    type(SplineData3D) :: A2_spline
    type(SplineData3D) :: A3_spline
    contains
        procedure :: jorek_field_init
        procedure :: compute_abfield
        procedure :: compute_afield
        procedure :: compute_bfield
        procedure :: compute_curla
        procedure :: compute_afield_derivatives
        procedure :: compute_fluxfunction
end type jorek_field_t


contains

subroutine jorek_field_init(self, jorek_filename)

    class(jorek_field_t), intent(out) :: self
    character(*), intent(in), optional :: jorek_filename

    type(field_mesh_t) :: field_mesh
    type(mesh_2d_t) :: fluxfunction_mesh

    call load_field_mesh_from_jorek(jorek_filename, field_mesh)
    call make_spline_from_mesh(field_mesh%A1, self%A1_spline)
    call make_spline_from_mesh(field_mesh%A2, self%A2_spline)
    call make_spline_from_mesh(field_mesh%A3, self%A3_spline)

    call load_fluxfunction_mesh_from_jorek(jorek_filename, fluxfunction_mesh)
    call make_spline_from_mesh_2d(fluxfunction_mesh, self%fluxfunction_spline)
end subroutine jorek_field_init

subroutine compute_abfield(self, x, A, B)
    class(jorek_field_t), intent(in) :: self
    real(dp), intent(in) :: x(3)
    real(dp), intent(out) :: A(3), B(3)

    call self%compute_afield(x, A)
    call self%compute_bfield(x, B)
end subroutine compute_abfield

subroutine compute_afield(self, x, A)
    class(jorek_field_t), intent(in) :: self
    real(dp), intent(in) :: x(3)
    real(dp), intent(out) :: A(3)

    call evaluate_splines_3d(self%A1_spline, x, A(1))
    call evaluate_splines_3d(self%A2_spline, x, A(2))
    call evaluate_splines_3d(self%A3_spline, x, A(3))
end subroutine compute_afield

subroutine compute_bfield(self, x, B)
    class(jorek_field_t), intent(in) :: self
    real(dp), intent(in) :: x(3)
    real(dp), intent(out) :: B(3)

    real(dp) :: curla(3), fluxfunction, R

    call self%compute_curla(x, curla)
    call self%compute_fluxfunction(x, fluxfunction)
    R = x(1)

    B(1) = curla(1)
    B(2) = curla(2) + fluxfunction / R
    B(3) = curla(3)
end subroutine compute_bfield

subroutine compute_curla(self, x, curla)
    class(jorek_field_t), intent(in) :: self
    real(dp), intent(in) :: x(3)
    real(dp), intent(out) :: curla(3)

    real(dp) :: dA_dx(3,3)
    real(dp) :: dAphi_dZ, dAZ_dphi, dAR_dphi, dAphi_dR, dAZ_dR, dAR_dZ
    real(dp) :: A(3), Aphi, R

    call self%compute_afield_derivatives(x, dA_dx)
    dAR_dphi = dA_dx(1,2)
    dAR_dZ = dA_dx(1,3)
    dAphi_dR = dA_dx(2,1)
    dAphi_dZ = dA_dx(2,3)
    dAZ_dR = dA_dx(3,1)
    dAZ_dphi = dA_dx(3,2)
    call self%compute_afield(x, A)
    Aphi = A(2)
    R = x(1)

    curla(1) = dAZ_dphi / R - dAphi_dZ
    curla(2) = dAR_dZ - dAZ_dR
    curla(3) = dAphi_dR + Aphi / R - dAR_dphi / R
end subroutine compute_curla

subroutine compute_afield_derivatives(self, x, dA_dx)
    class(jorek_field_t), intent(in) :: self
    real(dp), intent(in) :: x(3)
    real(dp), intent(out) :: dA_dx(3,3)

    real(dp), dimension(3) :: dA1_dx, dA2_dx, dA3_dx
    real(dp) :: dummy

    call evaluate_splines_3d_der(self%A1_spline, x, dummy, dA1_dx)
    call evaluate_splines_3d_der(self%A2_spline, x, dummy, dA2_dx)
    call evaluate_splines_3d_der(self%A3_spline, x, dummy, dA3_dx)

    dA_dx(1,:) = dA1_dx
    dA_dx(2,:) = dA2_dx
    dA_dx(3,:) = dA3_dx
end subroutine compute_afield_derivatives

subroutine compute_fluxfunction(self, x, fluxfunction)
    class(jorek_field_t), intent(in) :: self
    real(dp), intent(in) :: x(3)
    real(dp), intent(out) :: fluxfunction

    real(dp) :: x_2d(2)

    x_2d(1) = x(1)
    x_2d(2) = x(3)
    call evaluate_splines_2d(self%fluxfunction_spline, x_2d, fluxfunction)
end subroutine compute_fluxfunction


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

    A_R = values(2,:,:,:)
    A_Z = values(3,:,:,:)
    A_3 = values(4,:,:,:)
    do idx_R = 1, n_R
        A_phi(idx_R,:,:) = A_3(idx_R,:,:) / R(idx_R)
    end do
    B_R = values(12,:,:,:)
    B_Z = values(13,:,:,:)
    B_phi = values(14,:,:,:)

    call switch_phi_orientation(A_R)
    call switch_phi_orientation(A_Z)
    call switch_phi_orientation(A_phi)
    A_phi = -A_phi
    call switch_phi_orientation(B_R)
    call switch_phi_orientation(B_Z)
    call switch_phi_orientation(B_phi)
    B_phi = -B_phi

    is_periodic = [.false., .true., .false.]
    call field_mesh%A1%mesh_init(R, phi, Z, A_R, is_periodic)
    call field_mesh%A2%mesh_init(R, phi, Z, A_phi, is_periodic)
    call field_mesh%A3%mesh_init(R, phi, Z, A_Z, is_periodic)
    call field_mesh%B1%mesh_init(R, phi, Z, B_R, is_periodic)
    call field_mesh%B2%mesh_init(R, phi, Z, B_phi, is_periodic)
    call field_mesh%B3%mesh_init(R, phi, Z, B_Z, is_periodic)

    deallocate(A_R, A_Z, A_phi, A_3, B_R, B_Z, B_phi, values, R, Z, phi)
end subroutine load_field_mesh_from_jorek

subroutine load_fluxfunction_mesh_from_jorek(jorek_filename, fluxfunction_mesh)
    character(*), intent(in) :: jorek_filename
    type(mesh_2d_t), intent(out) :: fluxfunction_mesh

    integer :: n_R, n_Z, n_phi
    real(dp), dimension(:), allocatable :: R, Z, phi
    real(dp), dimension(:,:,:,:), allocatable :: values
    real(dp), dimension(:,:), allocatable :: fluxfunction
    logical :: is_periodic(2)

    call get_grid_and_values_from_jorek(jorek_filename, R, phi, Z, values, n_R, n_phi, n_Z)
    allocate(fluxfunction(n_R, n_Z))
             
    fluxfunction = values(11,:,1,:)
    fluxfunction = -fluxfunction

    is_periodic = [.false., .false.]
    call fluxfunction_mesh%mesh_init(R, Z, fluxfunction, is_periodic)
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
    n_phi = flipped_dims(1)
    n_R = flipped_dims(2)
    n_Z = flipped_dims(3)
    allocate(flipped_values(n_phi, n_R, n_Z, n_var))
    call h5_get(file_id, 'values', flipped_values)
    call h5_close(file_id)
    call h5_deinit()
    call phi_R_Z_var2var_R_phi_Z(flipped_values, values)
end subroutine read_dims_and_values_from_jorek

subroutine phi_R_Z_var2var_R_phi_Z(flipped_array, array)
    real(dp), dimension(:,:,:,:), intent(in) :: flipped_array
    real(dp), dimension(:,:,:,:), allocatable :: array

    integer :: dims(4)
    integer :: n_var, n_R, n_phi, n_Z
    integer :: var_idx, R_idx, phi_idx, Z_idx

    dims = shape(flipped_array)
    n_phi = dims(1)
    n_R = dims(2)
    n_Z = dims(3)
    n_var = dims(4)

    allocate(array(n_var, n_R, n_phi, n_Z))

    do var_idx = 1, n_var
        do R_idx = 1, n_R
            do phi_idx = 1, n_phi
                do Z_idx = 1, n_Z
                    array(var_idx,R_idx,phi_idx,Z_idx) = &
                                            flipped_array(phi_idx,R_idx,Z_idx,var_idx)
                end do
            end do
        end do
    end do
end subroutine phi_R_Z_var2var_R_phi_Z

subroutine get_ranges_from_filename(Rmin, Rmax, Zmin, Zmax, phimin, phimax, filename)
    character(*), intent(in) :: filename
    real(dp), intent(out) :: Rmin, Rmax, Zmin, Zmax, phimin, phimax

    call get_value_from_filename(filename, 'Rmin', Rmin)
    call get_value_from_filename(filename, 'Rmax', Rmax)
    call get_value_from_filename(filename, 'Zmin', Zmin)
    call get_value_from_filename(filename, 'Zmax', Zmax)
    call get_value_from_filename(filename, 'phimin', phimin)
    call get_value_from_filename(filename, 'phimax', phimax)
end subroutine get_ranges_from_filename

subroutine get_value_from_filename(filename, keyword, value)
    character(*), intent(in) :: filename
    character (*), intent(in) :: keyword
    real(dp), intent(out) :: value

    integer, parameter :: ndigits_max = 10
    character(ndigits_max) :: value_str
    integer :: keyword_start, value_start, value_end

    keyword_start = index(filename, keyword)
    value_start = len(keyword) + keyword_start
    value_str = filename(value_start:)
    value_end = index(value_str, '_') - 1
    if (value_end < 0) then
        value_end = len(value_str)
    end if
    value_str = value_str(:value_end)
    read(value_str, *) value
end subroutine get_value_from_filename


subroutine switch_phi_orientation(array)
    real(dp), dimension(:,:,:), intent(inout) :: array

    integer :: n_phi

    n_phi = size(array, 2)
    array = array(:, n_phi:1:-1, :)
end subroutine switch_phi_orientation


subroutine make_spline_from_mesh_2d(mesh, spline, order_in)
    use interpolate, only: construct_splines_2d

    class(mesh_2d_t), intent(in) :: mesh
    type(SplineData2D), intent(out) :: spline
    integer, optional :: order_in(2)

    real(dp) :: x_min(2), x_max(2)
    integer :: order(2)
        
    if (present(order_in)) then
        order = order_in
    else
        order = 3
    end if
    x_min = [mesh%x1(1), mesh%x2(1)]
    x_max = [mesh%x1(mesh%n1), mesh%x2(mesh%n2)]
    call construct_splines_2d(x_min, x_max, mesh%value, order, mesh%is_periodic, spline)
end subroutine make_spline_from_mesh_2d

end module neo_jorek_field