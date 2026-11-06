!> BEAMS3D field grid of a libneo field.
!>
!> BEAMS3D writes the magnetic field on its cylindrical grid to
!> beams3d_<ext>.h5 as B_R, B_PHI, B_Z (T) of shape (nr, nphi, nz) on raxis
!> (m), phiaxis (rad, one field period including its end point) and zaxis
!> (m).  Restarted with -continue_grid it reads only these with nr, nphi,
!> nz, S_ARR and U_ARR (LIBSTELL read_beams3d_mag).  write_beams3d_field
!> writes such a file with the grid and the flux labels S_ARR, U_ARR of a
!> BEAMS3D file and B_R, B_PHI, B_Z of the given field at the grid nodes;
!> the BEAMS3D run restarted from it writes the full file for ASCOT5.  The
!> field is evaluated in cm and G.
module neo_beams3d_field_grid
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use hdf5_tools, only: hid_t, h5_init, h5_deinit, h5_open, h5_create, &
                          h5_close, h5_get, h5_add
    use neo_field_base, only: field_t
    use neo_frozen_current_leaves, only: cart_to_cyl
    use neo_equilibrium_source, only: equilibrium_source_t
    use neo_flux_label, only: boundary_section_t, boundary_section, &
                              continued_flux_label
    implicit none
    private

    public :: beams3d_grid_t, read_beams3d_grid, beams3d_grid_field
    public :: write_beams3d_field

    type :: beams3d_grid_t
        real(dp), allocatable :: r(:), phi(:), z(:) !! m, rad, m
    end type beams3d_grid_t

contains

    subroutine read_beams3d_grid(filename, grid)
        character(*), intent(in) :: filename
        type(beams3d_grid_t), intent(out) :: grid
        integer(hid_t) :: h5id
        integer :: nr, nphi, nz

        call h5_init()
        call h5_open(filename, h5id)
        call h5_get(h5id, 'nr', nr)
        call h5_get(h5id, 'nphi', nphi)
        call h5_get(h5id, 'nz', nz)
        allocate (grid%r(nr), grid%phi(nphi), grid%z(nz))
        call h5_get(h5id, 'raxis', grid%r)
        call h5_get(h5id, 'phiaxis', grid%phi)
        call h5_get(h5id, 'zaxis', grid%z)
        call h5_close(h5id)
        call h5_deinit()
    end subroutine read_beams3d_grid

    !> (B_R, B_phi, B_Z) in T at the grid nodes, shape (nr, nphi, nz) each.
    subroutine beams3d_grid_field(field, grid, br, bphi, bz)
        class(field_t), intent(in) :: field
        type(beams3d_grid_t), intent(in) :: grid
        real(dp), allocatable, intent(out) :: br(:, :, :), bphi(:, :, :), bz(:, :, :)

        real(dp), parameter :: m_to_cm = 100.0_dp, gauss_to_tesla = 1.0e-4_dp
        real(dp) :: x(3), b(3), c, s
        integer :: i, j, k, nr, nphi, nz

        nr = size(grid%r)
        nphi = size(grid%phi)
        nz = size(grid%z)
        allocate (br(nr, nphi, nz), bphi(nr, nphi, nz), bz(nr, nphi, nz))
        !$omp parallel do collapse(3) schedule(dynamic, 16) default(shared) &
        !$omp private(i, j, k, x, b, c, s)
        do k = 1, nz
            do j = 1, nphi
                do i = 1, nr
                    c = cos(grid%phi(j))
                    s = sin(grid%phi(j))
                    x = [grid%r(i)*c, grid%r(i)*s, grid%z(k)]*m_to_cm
                    call field%compute_bfield(x, b)
                    call cart_to_cyl(c, s, b(1), b(2), br(i, j, k), bphi(i, j, k))
                    bz(i, j, k) = b(3)
                end do
            end do
        end do
        !$omp end parallel do
        br = br*gauss_to_tesla
        bphi = bphi*gauss_to_tesla
        bz = bz*gauss_to_tesla
    end subroutine beams3d_grid_field

    !> Grid and flux labels of template, field of field, written to output.
    !> With src, S_ARR outside the boundary of src is replaced by the
    !> continued flux label of neo_flux_label (inside it is kept, at most 1).
    subroutine write_beams3d_field(template, output, field, src)
        character(*), intent(in) :: template, output
        class(field_t), intent(in) :: field
        class(equilibrium_source_t), intent(in), optional :: src

        type(beams3d_grid_t) :: grid
        real(dp), allocatable :: br(:, :, :), bphi(:, :, :), bz(:, :, :)
        real(dp), allocatable :: s_arr(:, :, :), u_arr(:, :, :)
        integer(hid_t) :: h5id
        integer :: nr, nphi, nz

        call read_beams3d_grid(template, grid)
        nr = size(grid%r)
        nphi = size(grid%phi)
        nz = size(grid%z)
        allocate (s_arr(nr, nphi, nz), u_arr(nr, nphi, nz))
        call h5_init()
        call h5_open(template, h5id)
        call h5_get(h5id, 'S_ARR', s_arr)
        call h5_get(h5id, 'U_ARR', u_arr)
        call h5_close(h5id)
        if (present(src)) call continue_labels(src, grid, s_arr)
        call beams3d_grid_field(field, grid, br, bphi, bz)
        call h5_create(output, h5id)
        call h5_add(h5id, 'nr', nr)
        call h5_add(h5id, 'nphi', nphi)
        call h5_add(h5id, 'nz', nz)
        call h5_add(h5id, 'raxis', grid%r, lbound(grid%r), ubound(grid%r))
        call h5_add(h5id, 'phiaxis', grid%phi, lbound(grid%phi), ubound(grid%phi))
        call h5_add(h5id, 'zaxis', grid%z, lbound(grid%z), ubound(grid%z))
        call h5_add(h5id, 'B_R', br, lbound(br), ubound(br), &
                    'Radial magnetic field', 'T')
        call h5_add(h5id, 'B_PHI', bphi, lbound(bphi), ubound(bphi), &
                    'Toroidal magnetic field', 'T')
        call h5_add(h5id, 'B_Z', bz, lbound(bz), ubound(bz), &
                    'Vertical magnetic field', 'T')
        call h5_add(h5id, 'S_ARR', s_arr, lbound(s_arr), ubound(s_arr))
        call h5_add(h5id, 'U_ARR', u_arr, lbound(u_arr), ubound(u_arr))
        call h5_close(h5id)
        call h5_deinit()
    end subroutine write_beams3d_field

    subroutine continue_labels(src, grid, s_arr)
        class(equilibrium_source_t), intent(in) :: src
        type(beams3d_grid_t), intent(in) :: grid
        real(dp), intent(inout) :: s_arr(:, :, :)

        real(dp), parameter :: m_to_cm = 100.0_dp
        integer, parameter :: ntheta = 1024
        type(boundary_section_t) :: sec
        integer :: i, j, k

        do j = 1, size(grid%phi)
            call boundary_section(src, grid%phi(j), ntheta, sec)
            !$omp parallel do collapse(2) default(shared) private(i, k)
            do k = 1, size(grid%z)
                do i = 1, size(grid%r)
                    s_arr(i, j, k) = continued_flux_label(sec, grid%r(i)*m_to_cm, &
                                                          grid%z(k)*m_to_cm, &
                                                          s_arr(i, j, k))
                end do
            end do
            !$omp end parallel do
        end do
    end subroutine continue_labels

end module neo_beams3d_field_grid
