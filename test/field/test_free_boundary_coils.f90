!> Coils of a free-boundary run: N planar rectangular toroidal-field coils in
!> two groups give B_phi = -2 I_enc/(c R) between their legs (Ampere), with
!> I_enc = (N/2)(e_1 I_1 + e_2 I_2) in raw mgrid mode.  The groups alternate
!> and carry different currents, so the ripple at R = 3 m from the legs at
!> 5 m is of order (3/5)^(N/2) = 8e-8.  The filament field agrees with an mgrid
!> holding that field, disagrees when EXTCUR is wrong, and has no normal
!> component on a torus between the legs.
program test_free_boundary_coils
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use netcdf
    use neo_biotsavart, only: coils_t
    use neo_free_boundary_coils, only: free_boundary_coils, coil_field_cyl, &
                                       mgrid_consistency_t, check_mgrid_consistency, &
                                       normal_coil_field_t, normal_coil_field, &
                                       subtract_coil_field
    use neo_equilibrium_source, only: source_points_t
    implicit none

    real(dp), parameter :: pi = 3.141592653589793238462643383279502884_dp
    real(dp), parameter :: mu0_over_2pi = 2.0e-7_dp
    integer, parameter :: ncoil = 64
    real(dp), parameter :: r_in = 1.0_dp, r_out = 5.0_dp, h = 4.0_dp
    real(dp), parameter :: file_current(2) = [108.0_dp, 36.0_dp]
    real(dp), parameter :: extcur(2) = [1.2e4_dp, -3.0e3_dp]
    character(*), parameter :: coils_file = "coils.test_free_boundary"
    character(*), parameter :: mgrid_file = "mgrid_test_free_boundary.nc"

    type(coils_t) :: coils
    type(mgrid_consistency_t) :: mg
    type(normal_coil_field_t) :: nf
    type(source_points_t) :: pts
    real(dp) :: b(3), i_enc, b_ref

    call write_tf_coils()
    i_enc = 0.5_dp*ncoil*sum(extcur*file_current)

    call free_boundary_coils(coils_file, extcur, 'R', coils)
    b = coil_field_cyl(coils, [300.0_dp, 0.37_dp, 50.0_dp])
    b_ref = -mu0_over_2pi*i_enc/3.0_dp*1.0e4_dp
    print *, "raw mode B_phi relative error", b(2)/b_ref - 1.0_dp
    call check(abs(b(2)/b_ref - 1.0_dp) < 1.0e-6_dp, "raw-mode B_phi")
    call check(abs(b(1)) + abs(b(3)) < 1.0e-6_dp*abs(b_ref), "raw-mode B_R, B_Z")

    call free_boundary_coils(coils_file, extcur, 'S', coils)
    b = coil_field_cyl(coils, [300.0_dp, 0.37_dp, 50.0_dp])
    b_ref = -mu0_over_2pi*0.5_dp*ncoil*sum(extcur)/3.0_dp*1.0e4_dp
    call check(abs(b(2)/b_ref - 1.0_dp) < 1.0e-6_dp, "scaled-mode B_phi")

    call write_mgrid()
    call free_boundary_coils(coils_file, extcur, 'R', coils)
    call check_mgrid_consistency(coils, mgrid_file, extcur, mg)
    print *, "mgrid nodes", mg%nodes, " rms relative difference", mg%rms_rel
    call check(mg%rms_rel < 1.0e-4_dp, "filament field agrees with the mgrid")
    call check_mgrid_consistency(coils, mgrid_file, 1.01_dp*extcur, mg)
    call check(mg%rms_rel > 5.0e-3_dp, "a wrong EXTCUR is detected")

    allocate (pts%x(3, 2), pts%b(3, 2))
    pts%x = reshape([300.0_dp, 0.1_dp, 0.0_dp, 280.0_dp, 1.3_dp, 20.0_dp], [3, 2])
    pts%b(:, 1) = coil_field_cyl(coils, pts%x(:, 1)) + [1.0_dp, 2.0_dp, 3.0_dp]
    pts%b(:, 2) = coil_field_cyl(coils, pts%x(:, 2))
    call subtract_coil_field(coils, pts)
    call check(norm2(pts%b(:, 1) - [1.0_dp, 2.0_dp, 3.0_dp]) < 1.0e-9_dp*abs(b_ref) &
               .and. norm2(pts%b(:, 2)) < 1.0e-9_dp*abs(b_ref), "B_p = B - B_coils")

    call torus_normal_field(nf)
    print *, "normal coil field on the torus", nf%max_rel
    call check(nf%max_rel < 1.0e-4_dp, "no normal field of a toroidal field")
    print *, "PASS"

contains

    subroutine write_tf_coils()
        integer :: unit, c, g, k
        real(dp) :: phi, cp, sp, rz(2, 5), cur

        rz = reshape([r_in, -h, r_out, -h, r_out, h, r_in, h, r_in, -h], [2, 5])
        open (newunit=unit, file=coils_file, status='replace', action='write')
        write (unit, '(a)') 'periods 1'
        write (unit, '(a)') 'begin filament'
        write (unit, '(a)') 'mirror NIL'
        do c = 1, ncoil
            phi = 2.0_dp*pi*(c - 0.5_dp)/ncoil
            cp = cos(phi)
            sp = sin(phi)
            g = 1 + mod(c, 2)
            do k = 1, 5
                cur = file_current(g)
                if (k == 5) cur = 0.0_dp
                if (k < 5) then
                    write (unit, '(4es24.16)') rz(1, k)*cp, rz(1, k)*sp, rz(2, k), cur
                else
                    write (unit, '(4es24.16,i4,a)') rz(1, k)*cp, rz(1, k)*sp, &
                        rz(2, k), cur, g, ' tf'
                end if
            end do
        end do
        write (unit, '(a)') 'end'
        close (unit)
    end subroutine write_tf_coils

    !> mgrid in raw mode with the group fields B_phi = -mu0 (N/2) I_g/(2 pi R).
    subroutine write_mgrid()
        integer, parameter :: ir = 5, jz = 4, kp = 3
        real(dp), parameter :: rmin = 2.5_dp, rmax = 3.5_dp, zmin = -1.0_dp, zmax = 1.0_dp
        real(dp) :: bp(ir, jz, kp), zero(ir, jz, kp)
        character(len=3), parameter :: comp(3) = ['br_', 'bp_', 'bz_']
        integer :: ncid, d(3), vid, g, i
        character(len=8) :: name

        zero = 0.0_dp
        call ok(nf90_create(mgrid_file, NF90_CLOBBER, ncid))
        call ok(nf90_def_dim(ncid, 'rad', ir, d(1)))
        call ok(nf90_def_dim(ncid, 'zee', jz, d(2)))
        call ok(nf90_def_dim(ncid, 'phi', kp, d(3)))
        call def_scalar(ncid, 'nfp', NF90_INT)
        call def_scalar(ncid, 'nextcur', NF90_INT)
        call def_scalar(ncid, 'rmin', NF90_DOUBLE)
        call def_scalar(ncid, 'rmax', NF90_DOUBLE)
        call def_scalar(ncid, 'zmin', NF90_DOUBLE)
        call def_scalar(ncid, 'zmax', NF90_DOUBLE)
        do g = 1, 2
            do i = 1, 3
                write (name, '(a,i3.3)') comp(i), g
                call ok(nf90_def_var(ncid, trim(name), NF90_DOUBLE, d, vid))
            end do
        end do
        call ok(nf90_enddef(ncid))
        call put_int(ncid, 'nfp', 1)
        call put_int(ncid, 'nextcur', 2)
        call put_dbl(ncid, 'rmin', rmin)
        call put_dbl(ncid, 'rmax', rmax)
        call put_dbl(ncid, 'zmin', zmin)
        call put_dbl(ncid, 'zmax', zmax)
        do g = 1, 2
            do i = 1, ir
                bp(i, :, :) = -mu0_over_2pi*0.5_dp*ncoil*file_current(g)/ &
                              (rmin + (i - 1)*(rmax - rmin)/(ir - 1))
            end do
            write (name, '(a,i3.3)') 'br_', g
            call put_field(ncid, name, zero)
            write (name, '(a,i3.3)') 'bp_', g
            call put_field(ncid, name, bp)
            write (name, '(a,i3.3)') 'bz_', g
            call put_field(ncid, name, zero)
        end do
        call ok(nf90_close(ncid))
    end subroutine write_mgrid

    subroutine torus_normal_field(res)
        type(normal_coil_field_t), intent(out) :: res
        integer, parameter :: nt = 16, np = 8
        real(dp) :: x(3, nt*np), area(3, nt*np), b(3, nt*np), th, ph
        integer :: i, j, m

        do j = 1, np
            do i = 1, nt
                m = i + nt*(j - 1)
                th = 2.0_dp*pi*(i - 1)/nt
                ph = 2.0_dp*pi*(j - 1)/np
                x(:, m) = [300.0_dp + 40.0_dp*cos(th), ph, 40.0_dp*sin(th)]
                area(:, m) = [cos(th), 0.0_dp, sin(th)]
                b(:, m) = coil_field_cyl(coils, x(:, m))
            end do
        end do
        call normal_coil_field(coils, x, area, b, res)
    end subroutine torus_normal_field

    subroutine def_scalar(ncid, name, xtype)
        integer, intent(in) :: ncid, xtype
        character(*), intent(in) :: name
        integer :: vid

        call ok(nf90_def_var(ncid, name, xtype, vid))
    end subroutine def_scalar

    subroutine put_int(ncid, name, v)
        integer, intent(in) :: ncid, v
        character(*), intent(in) :: name
        integer :: vid

        call ok(nf90_inq_varid(ncid, name, vid))
        call ok(nf90_put_var(ncid, vid, v))
    end subroutine put_int

    subroutine put_dbl(ncid, name, v)
        integer, intent(in) :: ncid
        character(*), intent(in) :: name
        real(dp), intent(in) :: v
        integer :: vid

        call ok(nf90_inq_varid(ncid, name, vid))
        call ok(nf90_put_var(ncid, vid, v))
    end subroutine put_dbl

    subroutine put_field(ncid, name, v)
        integer, intent(in) :: ncid
        character(*), intent(in) :: name
        real(dp), intent(in) :: v(:, :, :)
        integer :: vid

        call ok(nf90_inq_varid(ncid, trim(name), vid))
        call ok(nf90_put_var(ncid, vid, v))
    end subroutine put_field

    subroutine ok(status)
        integer, intent(in) :: status

        if (status /= NF90_NOERR) then
            print *, trim(nf90_strerror(status))
            error stop 1
        end if
    end subroutine ok

    subroutine check(cond, what)
        logical, intent(in) :: cond
        character(*), intent(in) :: what

        if (.not. cond) then
            print *, "FAIL: ", what
            error stop 1
        end if
    end subroutine check

end program test_free_boundary_coils
