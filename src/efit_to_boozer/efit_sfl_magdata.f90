!> Contour construction (module efit_sfl) as the source of the efit_to_boozer_mod
!> grids that spline_magdata_in_symfluxcoord and magdata_in_symfluxcoord_ext
!> interpolate, in place of the field-line integration flint_for_Boozer. Shared
!> by efit_to_boozer.x, hamada_to_boozer.x and the f2py module _efit_to_boozer.
!>
!> The grids keep the flint_for_Boozer layout: nlabel surfaces, here at
!> s = (i/nlabel)^2 up to the edge s = 1, and ntheta points theta_j = 2 pi j/ntheta
!> of the symmetry-flux (PEST) angle. G is nu of the configured route (Boozer
!> by default), evaluated from its spectrum over theta without re-interpolation.
module efit_sfl_magdata
    use efit_sfl, only: efit_sfl_t, efit_sfl_config_t, efit_sfl_status_t, &
                        efit_sfl_surface_t, contour_t, efit_sfl_init, efit_sfl_symflux, &
                        efit_sfl_x_of_s, EFIT_SFL_OK
    implicit none
    private

    public :: efit_to_boozer_read_options, efit_sfl_load, efit_sfl_magdata_init, &
              efit_sfl_symflux_eval

    integer, parameter :: dp = kind(1.0d0)

contains

    !> Optional lines of efit_to_boozer.inp after psimax (defaults if absent):
    !> construction (1 contour, 2 legacy field lines), route,
    !> allow_poor_equilibrium, override_reason, n_theta of the contours, write_gs_map.
    subroutine efit_to_boozer_read_options(unit, construction, cfg, write_gs_map)
        integer, intent(in) :: unit
        integer, intent(out) :: construction, write_gs_map
        type(efit_sfl_config_t), intent(inout) :: cfg
        integer :: ios

        construction = 1
        write_gs_map = 0
        read (unit, *, iostat=ios) construction
        if (ios /= 0) construction = 1
        if (ios == 0) read (unit, *, iostat=ios) cfg%route
        if (ios == 0) read (unit, *, iostat=ios) cfg%allow_poor_equilibrium
        if (ios == 0) read (unit, *, iostat=ios) cfg%override_reason
        if (ios == 0) read (unit, *, iostat=ios) cfg%n_theta
        if (ios == 0) read (unit, *, iostat=ios) write_gs_map
    end subroutine efit_to_boozer_read_options

    !> Contour construction of the G-EQDSK named in field_divB0.inp with the
    !> edge at psimax of efit_to_boozer.inp (if set). A failed check returns in status.
    subroutine efit_sfl_load(sfl, cfg, status)
        use field_sub, only: read_field_input
        use input_files, only: gfile
        use efit_to_boozer_mod, only: psimax
        type(efit_sfl_t), intent(inout) :: sfl
        type(efit_sfl_config_t), intent(in) :: cfg
        type(efit_sfl_status_t), intent(out) :: status
        type(efit_sfl_config_t) :: c

        call read_field_input
        c = cfg
        if (psimax < 1.d299) c%psi_edge_rel = psimax
        call efit_sfl_init(sfl, trim(gfile), c, status)
    end subroutine efit_sfl_load

    !> efit_sfl_load and the filled grids (allocated by the caller as for
    !> flint_for_Boozer). A failed check or contour returns in status.
    subroutine efit_sfl_magdata_init(sfl, cfg, sigma, status)
        type(efit_sfl_t), intent(inout) :: sfl
        type(efit_sfl_config_t), intent(in) :: cfg
        real(dp), intent(out) :: sigma
        type(efit_sfl_status_t), intent(out) :: status

        call efit_sfl_load(sfl, cfg, status)
        sigma = sfl%sig_theta
        if (status%code /= EFIT_SFL_OK .and. .not. status%overridden) return
        call fill_grids(sfl, status)
    end subroutine efit_sfl_magdata_init

    subroutine fill_grids(sfl, status)
        use efit_to_boozer_mod, only: nlabel, ntheta, twopi, rmn, rmx, zmn, zmx, &
            raxis, zaxis, rbeg, rsmall, qsaf, psi_pol, psi_tor_vac, psi_tor, C_const, &
            R_spl, Z_spl, bmod_spl, sqgnorm_spl, Gfunc_spl
        type(efit_sfl_t), intent(in) :: sfl
        type(efit_sfl_status_t), intent(inout) :: status
        type(efit_sfl_surface_t) :: sv
        type(contour_t) :: ct
        real(dp) :: s, psi_unit, r, z, nu, bmod, dnu, dr_dth(sfl%nt)
        integer :: i, j, mpol

        mpol = sfl%nt/4
        ! Poloidal flux in the units of the file (efit_sfl works per radian).
        psi_unit = merge(twopi, 1.0_dp, sfl%diag%exp_bpol == 1)
        do i = 1, nlabel
            s = (real(i, dp)/real(nlabel, dp))**2
            call efit_sfl_symflux(sfl, efit_sfl_x_of_s(sfl, s), mpol, sv, ct, status)
            if (status%code /= EFIT_SFL_OK) then
                status%s = s
                status%isurf = i
                return
            end if
            rbeg(i) = ct%r(1)
            rsmall(i) = sqrt(sum(ct%rho**2)/real(sfl%nt, dp))
            qsaf(i) = 1.0_dp/sv%iota
            psi_pol(i) = (ct%psit - sfl%psi_ax)*psi_unit
            psi_tor(i) = s*sfl%psitor_max
            C_const(i) = sv%curr_g + sv%iota*sv%curr_i
            ! -(1/2 pi) oint R Z B_R dphi along the field line = -B_phi R <Z dR/dtheta_g / R>
            dr_dth = -ct%rho*sfl%sig_theta*((ct%pzraw*ct%ca - ct%praw*ct%sa)/ct%prho*ct%ca + ct%sa)
            psi_tor_vac(i) = -sv%curr_g*sum(ct%z*dr_dth/ct%r)/real(sfl%nt, dp)
            do j = 1, ntheta
                call efit_sfl_symflux_eval(sv, twopi*real(j, dp)/real(ntheta, dp), &
                                           r, z, nu, bmod, dnu)
                R_spl(0, j, i) = r
                Z_spl(0, j, i) = z
                bmod_spl(0, j, i) = bmod
                sqgnorm_spl(0, j, i) = r**2/ct%f
                Gfunc_spl(0, j, i) = nu
            end do
        end do
        rmn = sfl%g%R_eqd(1)
        rmx = sfl%g%R_eqd(sfl%nw)
        zmn = sfl%g%Z_eqd(1)
        zmx = sfl%g%Z_eqd(sfl%nh)
        raxis = sfl%r0
        zaxis = sfl%z0
        open (1, file='btor_rbig.dat')
        write (1, *) sfl%g%bcentr, sfl%g%rcentr
        close (1)
    end subroutine fill_grids

    !> R, Z, nu, |B| and d nu/d theta at the symmetry-flux angle theta from the
    !> spectra of efit_sfl_symflux (v = -nu/(2 pi)).
    pure subroutine efit_sfl_symflux_eval(sv, theta, r, z, nu, bmod, dnu_dtheta)
        type(efit_sfl_surface_t), intent(in) :: sv
        real(dp), intent(in) :: theta
        real(dp), intent(out) :: r, z, nu, bmod, dnu_dtheta
        real(dp), parameter :: twopi = 8.0_dp*atan(1.0_dp)
        real(dp) :: c, s, v, dv
        complex(dp) :: e, em
        integer :: m

        r = 0; z = 0; v = 0; bmod = 0; dv = 0
        e = exp(cmplx(0.0_dp, theta, dp))
        em = (1.0_dp, 0.0_dp)
        do m = 0, ubound(sv%rmnc, 1)
            if (m > 0) em = em*e
            c = real(em)
            s = aimag(em)
            r = r + sv%rmnc(m)*c + sv%rmns(m)*s
            z = z + sv%zmnc(m)*c + sv%zmns(m)*s
            v = v + sv%vmnc(m)*c + sv%vmns(m)*s
            bmod = bmod + sv%bmnc(m)*c + sv%bmns(m)*s
            dv = dv + m*(sv%vmns(m)*c - sv%vmnc(m)*s)
        end do
        nu = -twopi*v
        dnu_dtheta = -twopi*dv
    end subroutine efit_sfl_symflux_eval

end module efit_sfl_magdata
