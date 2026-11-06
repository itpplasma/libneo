!> BEAMS3D field grid export: from a file in BEAMS3D's layout (with a
!> rank-4 dataset as BEAMS3D writes) a grid file with B_R, B_PHI, B_Z of a
!> field at its nodes, in T, and the template's flux labels.  The field is that of N rectangular toroidal-field coils, so
!> B_phi = -2 N I/(c R) between their legs (Ampere) and B_R = B_Z = 0 up to
!> the ripple, (R/R_out)^N with R_out the outer legs.
program test_beams3d_field_grid
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use hdf5_tools, only: hid_t, h5_init, h5_deinit, h5_create, h5_open, h5_close, &
                          h5_add, h5_get
    use neo_biotsavart, only: coils_init, clight
    use neo_biotsavart_field, only: biotsavart_field_t
    use neo_beams3d_field_grid, only: write_beams3d_field
    implicit none

    real(dp), parameter :: pi = 3.141592653589793238462643383279502884_dp
    integer, parameter :: ncoil = 48, nr = 5, nphi = 4, nz = 3
    real(dp), parameter :: current = 1.0e5_dp !! A
    character(*), parameter :: template = 'beams3d_template_test.h5'
    character(*), parameter :: output = 'beams3d_frozen_test.h5'

    type(biotsavart_field_t) :: field
    real(dp) :: raxis(nr), phiaxis(nphi), zaxis(nz), s_arr(nr, nphi, nz)
   real(dp) :: br(nr, nphi, nz), bphi(nr, nphi, nz), bz(nr, nphi, nz), s_out(nr, nphi, nz)
    real(dp) :: ref, dist(2, 3, 2, 2)
    integer(hid_t) :: h5id
    integer :: i

    raxis = [(2.5_dp + 0.25_dp*(i - 1), i=1, nr)]
    phiaxis = [(0.4_dp*pi*(i - 1)/(nphi - 1), i=1, nphi)]
    zaxis = [-0.5_dp, 0.0_dp, 0.5_dp]
    s_arr = 1.5_dp
    s_arr(3, :, 2) = 0.25_dp
    call h5_init()
    call h5_create(template, h5id)
    call h5_add(h5id, 'nr', nr)
    call h5_add(h5id, 'nphi', nphi)
    call h5_add(h5id, 'nz', nz)
    call h5_add(h5id, 'raxis', raxis, lbound(raxis), ubound(raxis))
    call h5_add(h5id, 'phiaxis', phiaxis, lbound(phiaxis), ubound(phiaxis))
    call h5_add(h5id, 'zaxis', zaxis, lbound(zaxis), ubound(zaxis))
    br = 0.0_dp
    call h5_add(h5id, 'B_R', br, lbound(br), ubound(br))
    call h5_add(h5id, 'B_PHI', br, lbound(br), ubound(br))
    call h5_add(h5id, 'B_Z', br, lbound(br), ubound(br))
    call h5_add(h5id, 'S_ARR', s_arr, lbound(s_arr), ubound(s_arr))
    call h5_add(h5id, 'U_ARR', s_arr, lbound(s_arr), ubound(s_arr))
    dist = 1.0_dp
    call h5_add(h5id, 'dist_prof', dist, lbound(dist), ubound(dist))
    call h5_close(h5id)
    call h5_deinit()

    call tf_coils(field)
    call write_beams3d_field(template, output, field)

    call h5_init()
    call h5_open(output, h5id)
    call h5_get(h5id, 'B_R', br)
    call h5_get(h5id, 'B_PHI', bphi)
    call h5_get(h5id, 'B_Z', bz)
    call h5_get(h5id, 'S_ARR', s_out)
    call h5_close(h5id)
    call h5_deinit()

    do i = 1, nr
        ref = -2.0e-7_dp*ncoil*current/raxis(i)
        call check(all(abs(bphi(i, :, :)/ref - 1.0_dp) < 1.0e-6_dp), "B_PHI(R)")
    end do
    call check(maxval(abs(br)) + maxval(abs(bz)) < 1.0e-6_dp*maxval(abs(bphi)), &
               "B_R, B_Z")
    call check(all(s_out == s_arr), "flux labels kept")
    print *, "PASS"

contains

    !> N coils, legs at R = 1 m and 6 m, |Z| <= 3 m, in cm and statA.
    subroutine tf_coils(field)
        type(biotsavart_field_t), intent(out) :: field
        real(dp) :: x(5*ncoil), y(5*ncoil), z(5*ncoil), c(5*ncoil), rz(2, 5), phi
        integer :: k, j, m

        rz = reshape([100.0_dp, -300.0_dp, 600.0_dp, -300.0_dp, 600.0_dp, 300.0_dp, &
                      100.0_dp, 300.0_dp, 100.0_dp, -300.0_dp], [2, 5])
        m = 0
        do k = 1, ncoil
            phi = 2.0_dp*pi*(k - 0.5_dp)/ncoil
            do j = 1, 5
                m = m + 1
                x(m) = rz(1, j)*cos(phi)
                y(m) = rz(1, j)*sin(phi)
                z(m) = rz(2, j)
                c(m) = merge(0.0_dp, current*clight/10.0_dp, j == 5)
            end do
        end do
        call coils_init(x, y, z, c, field%coils)
    end subroutine tf_coils

    subroutine check(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what

        if (.not. ok) then
            print *, "FAIL: ", what
            error stop 1
        end if
    end subroutine check

end program test_beams3d_field_grid
