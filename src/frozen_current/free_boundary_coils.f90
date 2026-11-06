!> Coils of a free-boundary equilibrium and their filament field.
!>
!> VMEC's vacuum field is sum_g EXTCUR(g) b_g with b_g the coil-group field
!> stored in the mgrid.  MAKEGRID stores b_g for the currents of the coils
!> file in raw mode ('R') and divides it by the group's file current in
!> scaled mode ('S').  The physical current of a filament point is therefore
!> EXTCUR(g) I_file in raw mode and EXTCUR(g) I_file / I_g in scaled mode,
!> with I_g the file current of the group's first point.  The field is the
!> exact line integral over the filaments (neo_biotsavart, Hanson and
!> Hirshman), in cm and G; the coils file is in m and A.
module neo_free_boundary_coils
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neo_biotsavart, only: coils_t, compute_magnetic_field, &
                              compute_vector_potential, clight
    use neo_makegrid_coils, only: makegrid_coils_t, read_makegrid_coils, &
                                  makegrid_to_biotsavart
    use neo_frozen_current_leaves, only: cart_to_cyl
    implicit none
    private

    public :: free_boundary_coils, read_wout_coil_currents
    public :: coil_field_cyl, coil_field_at_points, subtract_coil_field
    public :: mgrid_consistency_t
    public :: check_mgrid_consistency, normal_coil_field_t, normal_coil_field

    real(dp), parameter :: m_to_cm = 100.0_dp
    real(dp), parameter :: a_to_stata = clight/10.0_dp

    type :: mgrid_consistency_t
        integer :: nodes = 0
        real(dp) :: max_abs = 0.0_dp !! max |B_fil - B_mgrid| (G)
        real(dp) :: rms_rel = 0.0_dp !! rms |B_fil - B_mgrid| / rms |B_mgrid|
        real(dp) :: max_rel = 0.0_dp !! max |B_fil - B_mgrid| / max |B_mgrid|
    end type mgrid_consistency_t

    type :: normal_coil_field_t
        real(dp) :: max_rel = 0.0_dp !! max |B_coils . n| / |B|
        real(dp) :: rms_rel = 0.0_dp !! area-weighted rms of |B_coils . n| / |B|
    end type normal_coil_field_t

contains

    !> Filaments of the coils file with the physical currents of the
    !> free-boundary run, in cm and statA.
    subroutine free_boundary_coils(coils_file, extcur, mgrid_mode, coils)
        character(*), intent(in) :: coils_file
        real(dp), intent(in) :: extcur(:)
        character(1), intent(in) :: mgrid_mode
        type(coils_t), intent(out) :: coils

        type(makegrid_coils_t) :: mc
        real(dp) :: scale(size(extcur))
        integer :: ic, g

        call read_makegrid_coils(coils_file, mc)
        scale = extcur
        if (mgrid_mode == 'S') then
            do g = 1, size(extcur)
                if (any(mc%group == g)) scale(g) = extcur(g)/group_file_current(mc, g)
            end do
        end if
        do ic = 1, mc%n_coils
            g = mc%group(ic)
            if (g < 1 .or. g > size(extcur)) error stop &
                'free_boundary_coils: coil group without EXTCUR entry'
            mc%current(mc%first(ic):mc%last(ic)) = &
                mc%current(mc%first(ic):mc%last(ic))*scale(g)
        end do
        call makegrid_to_biotsavart(mc, coils)
        coils%x = coils%x*m_to_cm
        coils%y = coils%y*m_to_cm
        coils%z = coils%z*m_to_cm
        coils%current = coils%current*a_to_stata
    end subroutine free_boundary_coils

    real(dp) function group_file_current(mc, g) result(i_g)
        type(makegrid_coils_t), intent(in) :: mc
        integer, intent(in) :: g
        integer :: ic

        do ic = 1, mc%n_coils
            if (mc%group(ic) == g) then
                i_g = mc%current(mc%first(ic))
                return
            end if
        end do
        error stop 'free_boundary_coils: empty coil group'
    end function group_file_current

    !> EXTCUR and mgrid mode of a free-boundary wout file.
    subroutine read_wout_coil_currents(wout_file, extcur, mgrid_mode)
        use netcdf
        character(*), intent(in) :: wout_file
        real(dp), allocatable, intent(out) :: extcur(:)
        character(1), intent(out) :: mgrid_mode

        integer :: ncid, vid, nextcur, lfreeb, st

        call nc_ok(nf90_open(wout_file, NF90_NOWRITE, ncid), 'open '//wout_file)
        call nc_ok(nf90_inq_varid(ncid, 'lfreeb__logical__', vid), 'lfreeb')
        call nc_ok(nf90_get_var(ncid, vid, lfreeb), 'lfreeb')
        if (lfreeb /= 1) error stop 'read_wout_coil_currents: fixed-boundary wout'
        call nc_ok(nf90_inq_varid(ncid, 'nextcur', vid), 'nextcur')
        call nc_ok(nf90_get_var(ncid, vid, nextcur), 'nextcur')
        allocate (extcur(nextcur))
        call nc_ok(nf90_inq_varid(ncid, 'extcur', vid), 'extcur')
        call nc_ok(nf90_get_var(ncid, vid, extcur), 'extcur')
        mgrid_mode = 'R'
        st = nf90_inq_varid(ncid, 'mgrid_mode', vid)
        if (st == NF90_NOERR) call nc_ok(nf90_get_var(ncid, vid, mgrid_mode), &
                                         'mgrid_mode')
        call nc_ok(nf90_close(ncid), 'close')
    end subroutine read_wout_coil_currents

    !> Coil field (B_R, B_phi, B_Z) at the cylindrical point x = (R, phi, Z).
    function coil_field_cyl(coils, x) result(b)
        type(coils_t), intent(in) :: coils
        real(dp), intent(in) :: x(3)
        real(dp) :: b(3), bc(3), c, s

        c = cos(x(2))
        s = sin(x(2))
        bc = compute_magnetic_field(coils, [x(1)*c, x(1)*s, x(3)])
        call cart_to_cyl(c, s, bc(1), bc(2), b(1), b(2))
        b(3) = bc(3)
    end function coil_field_cyl

    !> Coil field at cylindrical points x(3, n).
    subroutine coil_field_at_points(coils, x, b)
        type(coils_t), intent(in) :: coils
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: b(:, :)
        integer :: k

        !$omp parallel do schedule(dynamic, 64) default(shared) private(k)
        do k = 1, size(x, 2)
            b(:, k) = coil_field_cyl(coils, x(:, k))
        end do
        !$omp end parallel do
    end subroutine coil_field_at_points

    !> Plasma field B_p = B - B_coils at equilibrium source points.
    subroutine subtract_coil_field(coils, pts)
        use neo_equilibrium_source, only: source_points_t
        type(coils_t), intent(in) :: coils
        type(source_points_t), intent(inout) :: pts
        real(dp), allocatable :: bc(:, :)

        allocate (bc(3, size(pts%x, 2)))
        call coil_field_at_points(coils, pts%x, bc)
        pts%b = pts%b - bc
    end subroutine subtract_coil_field

    !> Compare the filament field with sum_g EXTCUR(g) b_g of the mgrid at
    !> its nodes.  mgrid fields are in T and positions in m.
    subroutine check_mgrid_consistency(coils, mgrid_file, extcur, res)
        use netcdf
        type(coils_t), intent(in) :: coils
        character(*), intent(in) :: mgrid_file
        real(dp), intent(in) :: extcur(:)
        type(mgrid_consistency_t), intent(out) :: res

        real(dp), parameter :: tesla_to_gauss = 1.0e4_dp
        real(dp), parameter :: twopi = 6.283185307179586476925286766559_dp
        real(dp), allocatable :: br(:, :, :), bp(:, :, :), bz(:, :, :), g3(:, :, :)
        real(dp), allocatable :: x(:, :), bf(:, :), bm(:, :)
        real(dp) :: rmin, rmax, zmin, zmax, d2, m2, dmax, mmax, dk, mk
        integer :: ncid, vid, ir, jz, kp, nfp, nextcur, g, i, j, k, m
        character(len=8) :: name

        call nc_ok(nf90_open(mgrid_file, NF90_NOWRITE, ncid), 'open '//mgrid_file)
        call dim(ncid, 'rad', ir)
        call dim(ncid, 'zee', jz)
        call dim(ncid, 'phi', kp)
        call get_int(ncid, 'nfp', nfp)
        call get_int(ncid, 'nextcur', nextcur)
        call get_dbl(ncid, 'rmin', rmin)
        call get_dbl(ncid, 'rmax', rmax)
        call get_dbl(ncid, 'zmin', zmin)
        call get_dbl(ncid, 'zmax', zmax)
        if (nextcur /= size(extcur)) error stop &
            'check_mgrid_consistency: EXTCUR and mgrid group counts differ'
        allocate (br(ir, jz, kp), bp(ir, jz, kp), bz(ir, jz, kp), g3(ir, jz, kp))
        br = 0.0_dp
        bp = 0.0_dp
        bz = 0.0_dp
        do g = 1, nextcur
            if (extcur(g) == 0.0_dp) cycle
            write (name, '(a,i3.3)') 'br_', g
            call get_group(ncid, name, g3)
            br = br + extcur(g)*g3
            write (name, '(a,i3.3)') 'bp_', g
            call get_group(ncid, name, g3)
            bp = bp + extcur(g)*g3
            write (name, '(a,i3.3)') 'bz_', g
            call get_group(ncid, name, g3)
            bz = bz + extcur(g)*g3
        end do
        call nc_ok(nf90_close(ncid), 'close')

        allocate (x(3, ir*jz*kp), bf(3, ir*jz*kp), bm(3, ir*jz*kp))
        do k = 1, kp
            do j = 1, jz
                do i = 1, ir
                    m = i + ir*((j - 1) + jz*(k - 1))
                    x(:, m) = [rmin + (i - 1)*(rmax - rmin)/(ir - 1), &
                               (k - 1)*twopi/(nfp*kp), &
                         zmin + (j - 1)*(zmax - zmin)/(jz - 1)]*[m_to_cm, 1.0_dp, m_to_cm]
                    bm(:, m) = [br(i, j, k), bp(i, j, k), bz(i, j, k)]*tesla_to_gauss
                end do
            end do
        end do
        call coil_field_at_points(coils, x, bf)
        d2 = 0.0_dp
        m2 = 0.0_dp
        dmax = 0.0_dp
        mmax = 0.0_dp
        do m = 1, size(x, 2)
            dk = norm2(bf(:, m) - bm(:, m))
            mk = norm2(bm(:, m))
            d2 = d2 + dk**2
            m2 = m2 + mk**2
            dmax = max(dmax, dk)
            mmax = max(mmax, mk)
        end do
        res%nodes = size(x, 2)
        res%max_abs = dmax
        res%rms_rel = sqrt(d2/m2)
        res%max_rel = dmax/mmax
    end subroutine check_mgrid_consistency

    !> Normal coil field on the boundary points of an equilibrium source,
    !> relative to the equilibrium |B| there.
    subroutine normal_coil_field(coils, x, area, b, res)
        type(coils_t), intent(in) :: coils
        real(dp), intent(in) :: x(:, :), area(:, :), b(:, :)
        type(normal_coil_field_t), intent(out) :: res

        real(dp), allocatable :: bc(:, :)
        real(dp) :: da, bn, a2, s2
        integer :: k

        allocate (bc(3, size(x, 2)))
        call coil_field_at_points(coils, x, bc)
        a2 = 0.0_dp
        s2 = 0.0_dp
        do k = 1, size(x, 2)
            da = norm2(area(:, k))
            bn = dot_product(bc(:, k), area(:, k))/(da*norm2(b(:, k)))
            res%max_rel = max(res%max_rel, abs(bn))
            s2 = s2 + bn**2*da
            a2 = a2 + da
        end do
        res%rms_rel = sqrt(s2/a2)
    end subroutine normal_coil_field

    subroutine get_group(ncid, name, v)
        use netcdf
        integer, intent(in) :: ncid
        character(*), intent(in) :: name
        real(dp), intent(out) :: v(:, :, :)
        integer :: vid

        call nc_ok(nf90_inq_varid(ncid, trim(name), vid), trim(name))
        call nc_ok(nf90_get_var(ncid, vid, v), trim(name))
    end subroutine get_group

    subroutine dim(ncid, name, n)
        use netcdf
        integer, intent(in) :: ncid
        character(*), intent(in) :: name
        integer, intent(out) :: n
        integer :: did

        call nc_ok(nf90_inq_dimid(ncid, name, did), name)
        call nc_ok(nf90_inquire_dimension(ncid, did, len=n), name)
    end subroutine dim

    subroutine get_int(ncid, name, v)
        use netcdf
        integer, intent(in) :: ncid
        character(*), intent(in) :: name
        integer, intent(out) :: v
        integer :: vid

        call nc_ok(nf90_inq_varid(ncid, name, vid), name)
        call nc_ok(nf90_get_var(ncid, vid, v), name)
    end subroutine get_int

    subroutine get_dbl(ncid, name, v)
        use netcdf
        integer, intent(in) :: ncid
        character(*), intent(in) :: name
        real(dp), intent(out) :: v
        integer :: vid

        call nc_ok(nf90_inq_varid(ncid, name, vid), name)
        call nc_ok(nf90_get_var(ncid, vid, v), name)
    end subroutine get_dbl

    subroutine nc_ok(status, what)
        use netcdf
        integer, intent(in) :: status
        character(*), intent(in) :: what

        if (status /= NF90_NOERR) then
            print *, 'neo_free_boundary_coils: ', trim(what), ': ', &
                trim(nf90_strerror(status))
            error stop 1
        end if
    end subroutine nc_ok

end module neo_free_boundary_coils
