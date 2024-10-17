module neo_jorek_field
use, intrinsic :: iso_fortran_env, only: dp => real64
use neo_field_mesh, only: field_mesh_t
use neo_spline_field, only: spline_field_t, make_spline_from_mesh
use interpolate, only: SplineData3D
use neo_mesh, only: mesh_t

implicit none


type, extends(spline_field_t) :: jorek_field_t
    type(SplineData3D) :: fluxfunction_spline
    contains
        procedure :: jorek_field_init
        procedure :: compute_fluxfunction
end type jorek_field_t


contains

subroutine jorek_field_init(self, jorek_filename)

    class(jorek_field_t), intent(out) :: self
    character(*), intent(in), optional :: jorek_filename

    type(field_mesh_t) :: field_mesh
    type(mesh_t) :: fluxfunction_mesh

    call load_field_mesh_from_jorek(jorek_filename, field_mesh)
    call load_fluxfunction_mesh_from_jorek(jorek_filename, fluxfunction_mesh)
    call set_b_mesh_to_curla_plus_fluxfunction(field_mesh, fluxfunction_mesh)
    call self%spline_field_init(field_mesh)
    call make_spline_from_mesh(fluxfunction_mesh, self%fluxfunction_spline)
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
    type(mesh_t), intent(out) :: fluxfunction_mesh

    integer :: n_R, n_Z, n_phi
    real(dp), dimension(:), allocatable :: R, Z, phi
    real(dp), dimension(:,:,:,:), allocatable :: values
    real(dp), dimension(:,:,:), allocatable :: fluxfunction
    logical :: is_periodic(3)

    call get_grid_and_values_from_jorek(jorek_filename, R, phi, Z, values, n_R, n_phi, n_Z)
    allocate(fluxfunction(n_R, n_phi, n_Z))
             
    fluxfunction = values(11,:,:,:)
    call switch_phi_orientation(fluxfunction)
    fluxfunction = -fluxfunction

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


subroutine set_b_mesh_to_curla_plus_fluxfunction(field_mesh, fluxfunction_mesh)
    type(field_mesh_t), intent(inout) :: field_mesh
    type(mesh_t), intent(in) :: fluxfunction_mesh

    integer :: idx_R, idx_phi, idx_Z, knote(3)
    real(dp) :: R
    real(dp) :: dAphi_dZ, dAZ_dphi, dAR_dphi, dAphi_dR, dAZ_dR, dAR_dZ, dummy
    real(dp) :: Aphi, fluxfunction
    real(dp) :: BR, Bphi, BZ

    do idx_R = 1, field_mesh%A1%n1
        R = field_mesh%A1%x1(idx_R)
        do idx_phi = 1, field_mesh%A1%n2
            do idx_Z = 1, field_mesh%A1%n3
                knote = [idx_R, idx_phi, idx_Z]
                call compute_mesh_derivs(field_mesh%A1, knote, dummy, dAR_dphi, dAR_dZ)
                call compute_mesh_derivs(field_mesh%A2, knote, dAphi_dR, dummy, dAphi_dZ)
                call compute_mesh_derivs(field_mesh%A3, knote, dAZ_dR, dAZ_dphi, dummy)
                Aphi = field_mesh%A2%value(idx_R, idx_phi, idx_Z)
                fluxfunction = fluxfunction_mesh%value(idx_R, idx_phi, idx_Z)

                BR = dAZ_dphi / R - dAphi_dZ
                Bphi = (dAR_dZ - dAZ_dR) + fluxfunction / R
                BZ = dAphi_dR + Aphi / R - dAR_dphi / R

                field_mesh%B1%value(idx_R, idx_phi, idx_Z) = BR
                field_mesh%B2%value(idx_R, idx_phi, idx_Z) = Bphi
                field_mesh%B3%value(idx_R, idx_phi, idx_Z) = BZ
            end do
        end do
    end do
end subroutine set_b_mesh_to_curla_plus_fluxfunction

subroutine compute_mesh_derivs(mesh, knote, df_dx1, df_dx2, df_dx3)
    type(mesh_t), intent(in) :: mesh
    integer, intent(in) :: knote(3)
    real(dp), intent(out) :: df_dx1, df_dx2, df_dx3

    integer :: idx1, idx2, idx3

    idx1 = knote(1)
    idx2 = knote(2)
    idx3 = knote(3)

    if (idx1 == 1) then
        df_dx1 = (-0.5_dp * mesh%value(idx1 + 2, idx2, idx3) + &
                  2_dp * mesh%value(idx1 + 1, idx2, idx3) - &
                  1.5_dp * mesh%value(idx1, idx2, idx3)) / mesh%dx1
    elseif (idx1 == mesh%n1) then
        df_dx1 = (1.5_dp * mesh%value(idx1, idx2, idx3) - &
                  2_dp * mesh%value(idx1 - 1, idx2, idx3) + &
                  0.5_dp * mesh%value(idx1 - 2, idx2, idx3)) / mesh%dx1
    else
        df_dx1 = (mesh%value(idx1+1, idx2, idx3) - &
                  mesh%value(idx1-1, idx2, idx3)) / &
                 (mesh%x1(idx1+1) - mesh%x1(idx1-1))
    end if

    if (idx2 == 1) then
        df_dx2 = (-0.5_dp * mesh%value(idx1, idx2 + 2, idx3) + &
                  2_dp * mesh%value(idx1, idx2 + 1, idx3) - &
                  1.5_dp * mesh%value(idx1, idx2, idx3)) / mesh%dx2
    elseif (idx2 == mesh%n2) then
        df_dx2 = (1.5_dp * mesh%value(idx1, idx2, idx3) - &
                  2_dp * mesh%value(idx1, idx2 - 1, idx3) + &
                  0.5_dp * mesh%value(idx1, idx2 - 2, idx3)) / mesh%dx2
    else
        df_dx2 = (mesh%value(idx1, idx2+1, idx3) - &
                  mesh%value(idx1, idx2-1, idx3)) / &
                 (mesh%x2(idx2+1) - mesh%x2(idx2-1))
    endif


    if (idx3 == 1) then
        df_dx3 = (-0.5_dp * mesh%value(idx1, idx2, idx3 + 2) + &
                  2_dp * mesh%value(idx1, idx2, idx3 + 1) - &
                  1.5_dp * mesh%value(idx1, idx2, idx3)) / mesh%dx3
    elseif (idx3 == mesh%n3) then
        df_dx3 = (1.5_dp * mesh%value(idx1, idx2, idx3) - &
                  2_dp * mesh%value(idx1, idx2, idx3 - 1) + &
                  0.5_dp * mesh%value(idx1, idx2, idx3 - 2)) / mesh%dx3
    else
        df_dx3 = (mesh%value(idx1, idx2, idx3+1) - &
                  mesh%value(idx1, idx2, idx3-1)) / &
                 (mesh%x3(idx3+1) - mesh%x3(idx3-1))
    endif
end subroutine compute_mesh_derivs

subroutine compute_fluxfunction(self, x, fluxfunction)
    use interpolate, only: evaluate_splines_3d

    class(jorek_field_t), intent(in) :: self
    real(dp), intent(in) :: x(3)
    real(dp), intent(out) :: fluxfunction

    call evaluate_splines_3d(self%fluxfunction_spline, x, fluxfunction)
end subroutine compute_fluxfunction

end module neo_jorek_field