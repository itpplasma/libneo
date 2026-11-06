!> Frozen-current field of a free-boundary VMEC equilibrium.
!>
!> Usage: frozen_current.x <namelist file>
!>
!> Reads the namelist &frozen_current (lengths in m), builds the coil
!> filaments with the wout file's EXTCUR, optionally checks them against the
!> mgrid, samples B - B_coils on the plasma volume and boundary, solves the
!> per-mode problem and reports the coil normal field on the boundary, the
!> per-mode energy and gauge fraction, and the interior residual
!> |curl A_p - B_p| relative to |B_p|.  With a BEAMS3D file as template it
!> writes a copy whose B_R, B_PHI, B_Z are the frozen-current field
!> (export_part 'total'), its plasma part ('plasma') or its coil part
!> ('coils').  extcur and mgrid_mode, when given, replace the wout values,
!> for example to evaluate another code's coil set.
program frozen_current_app
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use omp_lib, only: omp_get_wtime
    use neo_equilibrium_source, only: source_points_t, source_resolution_t
    use neo_vmec_equilibrium_source, only: vmec_equilibrium_source_t, &
                                           vmec_equilibrium_source_init
    use neo_free_boundary_coils, only: free_boundary_coils, read_wout_coil_currents, &
                                       subtract_coil_field, mgrid_consistency_t, &
                                       check_mgrid_consistency, normal_coil_field_t, &
                                       normal_coil_field
    use neo_frozen_current_field, only: frozen_current_patch_t, &
                                        frozen_current_field_t, solve_frozen_current, &
                                        frozen_current_plasma_cyl
    use neo_beams3d_field_grid, only: write_beams3d_field
    use neo_biotsavart_field, only: biotsavart_field_t
    implicit none

    character(len=1024) :: wout_file = 'wout.nc', coils_file = 'coils', mgrid_file = ''
    character(len=1024) :: beams3d_template = '', beams3d_output = ''
    character(len=1024) :: input
    character(len=16) :: export_part = 'total'
    character(len=1) :: mgrid_mode = ' '
    real(dp) :: extcur_in(64) = huge(1.0_dp)
    integer :: ns = 64, ntheta = 256, nzeta = 64, nr = 96, nz = 96, degree = 3
    integer :: kmax = 12, panels_per_span = 1
    real(dp) :: rmin = 0.0_dp, rmax = 0.0_dp, zmin = 0.0_dp, zmax = 0.0_dp
    namelist /frozen_current/ wout_file, coils_file, mgrid_file, beams3d_template, &
        beams3d_output, ns, ntheta, nzeta, nr, nz, degree, kmax, panels_per_span, &
        rmin, rmax, zmin, zmax, export_part, mgrid_mode, extcur_in

    real(dp), parameter :: m_to_cm = 100.0_dp
    type(vmec_equilibrium_source_t) :: src
    type(source_resolution_t) :: res
    type(source_points_t) :: vol, surf
    type(frozen_current_patch_t) :: patch
    type(frozen_current_field_t) :: field, plasma_only
    type(biotsavart_field_t) :: coils_only
    type(mgrid_consistency_t) :: mg
    type(normal_coil_field_t) :: nf
    real(dp), allocatable :: extcur(:), btot(:, :)
    character(len=1) :: mode
    real(dp) :: t0, t1, bp2, b2, d2, a(3), b(3)
    integer :: unit, k, status

    if (command_argument_count() /= 1) error stop 'usage: frozen_current.x <namelist>'
    call get_command_argument(1, input)
    open (newunit=unit, file=trim(input), status='old', action='read')
    read (unit, nml=frozen_current)
    close (unit)

    t0 = omp_get_wtime()
    call read_wout_coil_currents(trim(wout_file), extcur, mode)
    if (mgrid_mode /= ' ') mode = mgrid_mode
    if (extcur_in(1) /= huge(1.0_dp)) then
        deallocate (extcur)
        allocate (extcur(count(extcur_in /= huge(1.0_dp))))
        extcur = pack(extcur_in, extcur_in /= huge(1.0_dp))
    end if
    print '(a,a,a,*(es12.4))', ' mgrid mode ', mode, ', EXTCUR', extcur
    call free_boundary_coils(trim(coils_file), extcur, mode, field%coils)
    if (len_trim(mgrid_file) > 0) then
        call check_mgrid_consistency(field%coils, trim(mgrid_file), extcur, mg)
        print '(a,i0,a,es10.3,a,es10.3,a,es10.3,a)', ' mgrid check: ', mg%nodes, &
            ' nodes, rms ', mg%rms_rel, ', max ', mg%max_rel, ' (rel.), max ', &
            mg%max_abs*1.0e-4_dp, ' T'
    end if

    if (trim(export_part) == 'coils') then
        coils_only%coils = field%coils
        call write_beams3d_field(trim(beams3d_template), trim(beams3d_output), &
                                 coils_only)
        print '(a,a)', ' wrote the coil field to ', trim(beams3d_output)
        stop
    end if

    call vmec_equilibrium_source_init(src, trim(wout_file))
    res%ns = ns
    res%ntheta = ntheta
    res%nzeta = nzeta
    call src%sample_volume(res, vol)
    call src%sample_surface(res, surf)
    allocate (btot(3, size(vol%x, 2)))
    btot = vol%b
    call normal_coil_field(field%coils, surf%x, surf%area, surf%b, nf)
    print '(a,es10.3,a,es10.3)', ' boundary B_coils.n/|B|: max ', nf%max_rel, &
        ', rms ', nf%rms_rel
    call subtract_coil_field(field%coils, vol)
    call subtract_coil_field(field%coils, surf)
    bp2 = sum(vol%w*sum(vol%b**2, 1))
    b2 = sum(vol%w*sum(btot**2, 1))
    print '(a,es10.3)', ' volume rms |B_p|/|B|: ', sqrt(bp2/b2)
    t1 = omp_get_wtime()
    print '(a,i0,a,i0,a,f8.1,a)', ' source: ', size(vol%x, 2), ' volume and ', &
        size(surf%x, 2), ' boundary points, ', t1 - t0, ' s'

    patch%r0 = rmin*m_to_cm
    patch%r1 = rmax*m_to_cm
    patch%z0 = zmin*m_to_cm
    patch%z1 = zmax*m_to_cm
    patch%nr = nr
    patch%nz = nz
    patch%degree = degree
    patch%kmax = kmax
    patch%panels_per_span = panels_per_span
    call solve_frozen_current(patch, vol, surf, field%plasma, status)
    if (status /= 0) error stop 'frozen_current: solve failed'
    t0 = omp_get_wtime()
    print '(a,i0,a,f8.1,a)', ' solve: ', kmax + 1, ' modes, ', t0 - t1, ' s'
    do k = 0, kmax
        print '(a,i4,a,es10.3,a,es10.3)', ' n = ', k*src%nfp, '  energy ', &
            field%plasma%energy(k), '  gauge ', field%plasma%gauge(k)
    end do

    d2 = 0.0_dp
    !$omp parallel do default(shared) private(k, a, b, status) reduction(+:d2)
    do k = 1, size(vol%x, 2)
        call frozen_current_plasma_cyl(field%plasma, vol%x(:, k), a, b, status)
        d2 = d2 + vol%w(k)*sum((b - vol%b(:, k))**2)
    end do
    !$omp end parallel do
    t1 = omp_get_wtime()
    print '(a,es10.3,a,es10.3,a,f8.1,a)', ' interior |curl A_p - B_p|: ', &
        sqrt(d2/bp2), ' of |B_p|, ', sqrt(d2/b2), ' of |B| (', t1 - t0, ' s)'

    if (len_trim(beams3d_template) > 0) then
        select case (trim(export_part))
        case ('plasma')
            plasma_only%plasma = field%plasma
            call write_beams3d_field(trim(beams3d_template), trim(beams3d_output), &
                                     plasma_only)
        case default
            call write_beams3d_field(trim(beams3d_template), trim(beams3d_output), &
                                     field)
        end select
        print '(a,a,a,f8.1,a)', ' wrote ', trim(beams3d_output), ' (', &
            omp_get_wtime() - t1, ' s)'
    end if
end program frozen_current_app
