program near_axis_field_diag
    !> Near-axis VMEC field diagnostic. Splines the field with a chosen axis
    !> healing variant and dumps a radial cut at fixed (theta, varphi) so the
    !> continuity and rho^|m| behaviour can be inspected. One row per rho point:
    !>   rho s theta A_theta A_phi dA_theta_ds dA_phi_ds Bmod Bcov_theta Bcov_phi
    !> Radial derivatives are taken in post from the fine rho grid.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use new_vmec_stuff_mod, only: netcdffile, ns_s, ns_tp, multharm, &
                                  old_axis_healing, old_axis_healing_boundary, &
                                  axis_healing_power_law, rho_axis_heal, &
                                  axis_healing_polyfit, axis_healing_polyfit_degree
    use spline_vmec_sub, only: spline_vmec_data, vmec_field
    implicit none

    character(len=1024) :: wout_file, variant, outfile, arg
    integer :: nrho, irho, itheta
    real(dp) :: rho, s, theta, varphi, rho_min, rho_max
    real(dp) :: A_theta, A_phi, dA_theta_ds, dA_phi_ds, aiota
    real(dp) :: sqg, alam, dl_ds, dl_dt, dl_dp
    real(dp) :: Bctr_th, Bctr_ph, Bcov_s, Bcov_th, Bcov_ph, Bmod2, Bmod
    real(dp), parameter :: pi = 3.14159265358979d0
    real(dp) :: thetas(3)

    call get_command_argument(1, wout_file)
    call get_command_argument(2, variant)
    call get_command_argument(3, outfile)
    if (len_trim(wout_file) == 0 .or. len_trim(variant) == 0 .or. len_trim(outfile) == 0) then
        print *, "Usage: near_axis_field_diag.x <wout.nc> <power|oldheal|oldbnd> <out.dat> [rho_axis_heal]"
        error stop 1
    end if

    ! Geometry/resolution matching the benchmark runs (mh5 ns3).
    ns_s = 3
    ns_tp = 3
    multharm = 5
    rho_axis_heal = 0.1d0
    call get_command_argument(4, arg)
    if (len_trim(arg) > 0) read (arg, *) rho_axis_heal

    select case (trim(variant))
    case ('polyfit')
        axis_healing_polyfit = .True.
        axis_healing_polyfit_degree = 3
        axis_healing_power_law = .False.
    case ('power')
        axis_healing_power_law = .True.
        old_axis_healing = .True.
        old_axis_healing_boundary = .True.
    case ('oldheal')
        axis_healing_power_law = .False.
        old_axis_healing = .True.
        old_axis_healing_boundary = .False.
    case ('oldbnd')
        axis_healing_power_law = .False.
        old_axis_healing = .True.
        old_axis_healing_boundary = .True.
    case default
        print *, "unknown variant: ", trim(variant)
        error stop 2
    end select

    netcdffile = trim(wout_file)
    call spline_vmec_data

    thetas = [0.0d0, 0.5d0*pi, pi]
    varphi = 0.0d0
    nrho = 4000
    rho_min = 1.0d-4
    rho_max = 0.40d0

    open (unit=21, file=trim(outfile), status='replace', action='write')
    write (21, '(A)') '# rho s theta A_theta A_phi dA_theta_ds dA_phi_ds Bmod Bcov_theta Bcov_phi'
    do itheta = 1, 3
        theta = thetas(itheta)
        do irho = 0, nrho
            rho = rho_min + (rho_max - rho_min)*dble(irho)/dble(nrho)
            s = rho*rho
            call vmec_field(s, theta, varphi, A_theta, A_phi, dA_theta_ds, dA_phi_ds, &
                            aiota, sqg, alam, dl_ds, dl_dt, dl_dp, Bctr_th, Bctr_ph, &
                            Bcov_s, Bcov_th, Bcov_ph)
            Bmod2 = Bctr_th*Bcov_th + Bctr_ph*Bcov_ph
            Bmod = sqrt(max(Bmod2, 0.0d0))
            write (21, '(10ES20.11)') rho, s, theta, A_theta, A_phi, dA_theta_ds, &
                dA_phi_ds, Bmod, Bcov_th, Bcov_ph
        end do
        write (21, '(A)') ''
    end do
    close (21)
    print *, 'wrote ', trim(outfile), ' variant=', trim(variant), ' rho_axis_heal=', rho_axis_heal
end program near_axis_field_diag
