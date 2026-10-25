!> Boozer spectrum of a VMEC wout file with tangent or adjoint derivatives.
!>
!> Usage: boozer_spectrum_ad.x primal|jvp|vjp <wout.nc> <in.nc> <out.nc>
!>
!> in.nc holds s(nsurf), the integers mboz, nboz and grid_refinement, optionally
!> the integer construction (1 = hybrid, the default; 2 = legacy), the integers
!> n_theta and n_phi (hybrid grid override, default 0 = automatic), the
!> integer route (hybrid only: 1 sfl_covariant = Boozer, the default,
!> 2 sfl_jacobian_b2, 3 sfl_jacobian_const = Hamada I, 4 sfl_straight_current
!> = Hamada II, 5 sfl_pest), the integers mnu and nnu (nu truncation of the
!> routes 2 to 4, default 0 = 4 (mpol - 1), 4 ntor), the integer diagnostics
!> (primal only, 1 writes the route diagnostics, see below), the integer
!> lambda_cubic (legacy only: 1 selects lambda_half_to_full = 'cubic', default
!> 0 = 'linear'), and for
!> jvp the wout-shaped tangents rmnc, zmns, lmns (ns x mnmax), iotaf (ns, legacy)
!> and iotas (ns, hybrid)
!> and the scalars phiedge (edge toroidal flux, scaling the phi profile) and
!> rmajor (Rmajor_p), for vjp the cotangents bmnc_b, bmns_b, numns_b, numnc_b,
!> gradsmnc_b, gradsmns_b (nsurf x mnboz) and iota, G, I, psi_edge, R (nsurf).
!> Absent arrays are zero. out.nc receives xm_b, xn_b, the primal spectrum
!> (with |grad s| = gradsmnc_b, psi_edge and R as rabe uses them) and either
!> its tangent (prefix d_) or the wout cotangents rmnc, zmns, lmns, iotaf,
!> iotas, phiedge and rmajor. Units: m, T, T m, 1/m and Wb. With diagnostics
!> = 1, out.nc also holds per surface diag_d_jac, diag_d_jacb2, diag_d_cov,
!> diag_d_cur (relative variations of J*, J* B^2, of the covariant components
!> and of J* J^phi*), diag_d_cur_n (ptp of J* J^phi* over |I'| + |G'|),
!> diag_defect (the property the route does not enforce),
!> diag_eps_fb (averaged radial force balance), diag_pprime, diag_Iprime,
!> diag_Gprime, diag_n_zeroed (exactly resonant modes set to zero) and the
!> near-resonant modes of the Jacobian routes res_surface, res_m, res_n,
!> res_dist (iota m - n N), res_f (|f_mn|) and res_amp (|f_mn/(iota m - n N)|).
program boozer_spectrum_ad
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use netcdf
    use new_vmec_stuff_mod, only: lambda_half_to_full
    use boozer_vmec_ad, only: boozer_ad_context_t, boozer_wout_tangent_t, &
                              boozer_spectrum_t, boozer_ad_init, boozer_ad_primal, &
                              boozer_ad_jvp, boozer_ad_vjp, wout_tangent_alloc, &
                              spectrum_alloc, BOOZER_HYBRID, SFL_COVARIANT, &
                              sfl_route_name, hybrid_diag_t, boozer_ad_diagnostics
    implicit none

    character(len=1024) :: mode, wout, fin, fout
    type(boozer_ad_context_t) :: ctx
    type(boozer_wout_tangent_t) :: t
    type(boozer_spectrum_t) :: y, yd, yb
    real(dp), allocatable :: s(:)
    integer :: ncid, mboz, nboz, refine, nsurf, construction, nth, nph, route, mnu, nnu
    integer :: want_diag
    type(hybrid_diag_t), allocatable :: diag(:)
    integer(int64) :: c0, c1, rate

    if (command_argument_count() /= 4) then
        print *, 'usage: boozer_spectrum_ad.x primal|jvp|vjp wout.nc in.nc out.nc'
        error stop 2
    end if
    call get_command_argument(1, mode)
    call get_command_argument(2, wout)
    call get_command_argument(3, fin)
    call get_command_argument(4, fout)

    call check(nf90_open(trim(fin), NF90_NOWRITE, ncid))
    nsurf = dim_len(ncid, 's')
    allocate (s(nsurf))
    call get_var(ncid, 's', s)
    mboz = get_int(ncid, 'mboz')
    nboz = get_int(ncid, 'nboz')
    refine = get_int(ncid, 'grid_refinement')
    lambda_half_to_full = 'linear'
    if (get_int_opt(ncid, 'lambda_cubic', 0) == 1) lambda_half_to_full = 'cubic'

    construction = get_int_opt(ncid, 'construction', BOOZER_HYBRID)
    nth = get_int_opt(ncid, 'n_theta', 0)
    nph = get_int_opt(ncid, 'n_phi', 0)
    route = get_int_opt(ncid, 'route', SFL_COVARIANT)
    mnu = get_int_opt(ncid, 'mnu', 0)
    nnu = get_int_opt(ncid, 'nnu', 0)
    want_diag = get_int_opt(ncid, 'diagnostics', 0)
    call boozer_ad_init(ctx, trim(wout), mboz, nboz, grid_refinement=refine, &
                        construction=construction, n_theta=nth, n_phi=nph, route=route, &
                        mnu=mnu, nnu=nnu)
    call system_clock(c0, rate)
    select case (trim(mode))
    case ('primal')
        call boozer_ad_primal(ctx, s, y)
        if (want_diag == 1) call boozer_ad_diagnostics(ctx, s, diag)
    case ('jvp')
        call wout_tangent_alloc(ctx, t)
        call get_opt2(ncid, 'rmnc', t%rmnc)
        call get_opt2(ncid, 'zmns', t%zmns)
        call get_opt2(ncid, 'lmns', t%lmns)
        call get_opt1(ncid, 'iotaf', t%iotaf)
        call get_opt1(ncid, 'iotas', t%iotas)
        call get_opt0(ncid, 'phiedge', t%phiedge)
        call get_opt0(ncid, 'rmajor', t%rmajor)
        call boozer_ad_jvp(ctx, s, t, y, yd)
    case ('vjp')
        call spectrum_alloc(ctx, nsurf, yb)
        call get_opt2(ncid, 'bmnc_b', yb%bmnc)
        call get_opt2(ncid, 'bmns_b', yb%bmns)
        call get_opt2(ncid, 'numns_b', yb%numns)
        call get_opt2(ncid, 'numnc_b', yb%numnc)
        call get_opt1(ncid, 'iota', yb%iota)
        call get_opt1(ncid, 'G', yb%g)
        call get_opt1(ncid, 'I', yb%i)
        call get_opt2(ncid, 'gradsmnc_b', yb%gradsmnc)
        call get_opt2(ncid, 'gradsmns_b', yb%gradsmns)
        call get_opt1(ncid, 'psi_edge', yb%psi_edge)
        call get_opt1(ncid, 'R', yb%rmajor)
        call boozer_ad_vjp(ctx, s, yb, t, y)
    case default
        error stop 'mode must be primal, jvp or vjp'
    end select
    call system_clock(c1)
    call check(nf90_close(ncid))
    call write_output()

contains

    subroutine write_output()
        integer :: dsurf, dmn, dns, dmnmax, v(40), k

        call check(nf90_create(trim(fout), NF90_CLOBBER, ncid))
        call check(nf90_def_dim(ncid, 'surface', nsurf, dsurf))
        call check(nf90_def_dim(ncid, 'mnboz', ctx%mnboz, dmn))
        call check(nf90_put_att(ncid, NF90_GLOBAL, 'mode', trim(mode)))
        call check(nf90_put_att(ncid, NF90_GLOBAL, 'seconds', &
                                real(c1 - c0, dp)/real(rate, dp)))
        call check(nf90_put_att(ncid, NF90_GLOBAL, 'n_theta', ctx%grid%n_theta))
        call check(nf90_put_att(ncid, NF90_GLOBAL, 'n_phi', ctx%grid%n_phi))
        call check(nf90_put_att(ncid, NF90_GLOBAL, 'construction', construction))
        call check(nf90_put_att(ncid, NF90_GLOBAL, 'route', route))
        call check(nf90_put_att(ncid, NF90_GLOBAL, 'route_name', trim(sfl_route_name(route))))
        if (construction == BOOZER_HYBRID) then
            call check(nf90_put_att(ncid, NF90_GLOBAL, 'mnu', ctx%hyb%mnu))
            call check(nf90_put_att(ncid, NF90_GLOBAL, 'nnu', ctx%hyb%nnu))
        end if
        call check(nf90_put_att(ncid, NF90_GLOBAL, 'lambda_half_to_full', &
                                trim(lambda_half_to_full)))
        k = 0
        call def(v, k, 'xm_b', NF90_INT, [dmn])
        call def(v, k, 'xn_b', NF90_INT, [dmn])
        call def(v, k, 's', NF90_DOUBLE, [dsurf])
        call def_spectrum(v, k, '', dsurf, dmn)
        if (trim(mode) == 'jvp') call def_spectrum(v, k, 'd_', dsurf, dmn)
        if (trim(mode) == 'vjp') then
            call check(nf90_def_dim(ncid, 'radius', ctx%ns, dns))
            call check(nf90_def_dim(ncid, 'mn_mode', ctx%nstrm, dmnmax))
            call def(v, k, 'rmnc', NF90_DOUBLE, [dmnmax, dns])
            call def(v, k, 'zmns', NF90_DOUBLE, [dmnmax, dns])
            call def(v, k, 'lmns', NF90_DOUBLE, [dmnmax, dns])
            call def(v, k, 'iotaf', NF90_DOUBLE, [dns])
            call def0(v, k, 'phiedge')
            call def0(v, k, 'rmajor')
            call def(v, k, 'iotas', NF90_DOUBLE, [dns])
        end if
        call check(nf90_enddef(ncid))
        call check(nf90_put_var(ncid, v(1), ctx%xm_b))
        call check(nf90_put_var(ncid, v(2), ctx%xn_b))
        call check(nf90_put_var(ncid, v(3), s))
        call put_spectrum(v(4:14), y)
        if (trim(mode) == 'jvp') call put_spectrum(v(15:25), yd)
        if (trim(mode) == 'vjp') then
            call check(nf90_put_var(ncid, v(15), t%rmnc))
            call check(nf90_put_var(ncid, v(16), t%zmns))
            call check(nf90_put_var(ncid, v(17), t%lmns))
            call check(nf90_put_var(ncid, v(18), t%iotaf))
            call check(nf90_put_var(ncid, v(19), t%phiedge))
            call check(nf90_put_var(ncid, v(20), t%rmajor))
            call check(nf90_put_var(ncid, v(21), t%iotas))
        end if
        if (allocated(diag)) call write_diagnostics(dsurf)
        call check(nf90_close(ncid))
    end subroutine write_output

    !> Route diagnostics, appended to the output file.
    subroutine write_diagnostics(dsurf)
        integer, intent(in) :: dsurf
        integer :: dres, v(17), k, nres, j, i0, n
        integer, allocatable :: isurf(:), rm(:), rn(:)
        real(dp), allocatable :: rd(:), rf(:), ra(:)

        nres = 0
        do j = 1, nsurf
            nres = nres + size(diag(j)%res_m)
        end do
        allocate (isurf(nres), rm(nres), rn(nres), rd(nres), rf(nres), ra(nres))
        i0 = 0
        do j = 1, nsurf
            n = size(diag(j)%res_m)
            isurf(i0 + 1:i0 + n) = j
            rm(i0 + 1:i0 + n) = diag(j)%res_m
            rn(i0 + 1:i0 + n) = diag(j)%res_n
            rd(i0 + 1:i0 + n) = diag(j)%res_dist
            rf(i0 + 1:i0 + n) = diag(j)%res_f
            ra(i0 + 1:i0 + n) = diag(j)%res_amp
            i0 = i0 + n
        end do
        call check(nf90_redef(ncid))
        call check(nf90_def_dim(ncid, 'resonance', max(nres, 1), dres))
        k = 0
        call def(v, k, 'diag_d_jac', NF90_DOUBLE, [dsurf])
        call def(v, k, 'diag_d_jacb2', NF90_DOUBLE, [dsurf])
        call def(v, k, 'diag_d_cov', NF90_DOUBLE, [dsurf])
        call def(v, k, 'diag_d_cur', NF90_DOUBLE, [dsurf])
        call def(v, k, 'diag_defect', NF90_DOUBLE, [dsurf])
        call def(v, k, 'diag_eps_fb', NF90_DOUBLE, [dsurf])
        call def(v, k, 'diag_pprime', NF90_DOUBLE, [dsurf])
        call def(v, k, 'diag_Iprime', NF90_DOUBLE, [dsurf])
        call def(v, k, 'diag_Gprime', NF90_DOUBLE, [dsurf])
        call def(v, k, 'diag_n_zeroed', NF90_INT, [dsurf])
        call def(v, k, 'res_surface', NF90_INT, [dres])
        call def(v, k, 'res_m', NF90_INT, [dres])
        call def(v, k, 'res_n', NF90_INT, [dres])
        call def(v, k, 'res_dist', NF90_DOUBLE, [dres])
        call def(v, k, 'res_f', NF90_DOUBLE, [dres])
        call def(v, k, 'res_amp', NF90_DOUBLE, [dres])
        call def(v, k, 'diag_d_cur_n', NF90_DOUBLE, [dsurf])
        call check(nf90_put_att(ncid, NF90_GLOBAL, 'n_resonance', nres))
        call check(nf90_enddef(ncid))
        call check(nf90_put_var(ncid, v(1), diag%d_jac))
        call check(nf90_put_var(ncid, v(2), diag%d_jacb2))
        call check(nf90_put_var(ncid, v(3), diag%d_cov))
        call check(nf90_put_var(ncid, v(4), diag%d_cur))
        call check(nf90_put_var(ncid, v(5), diag%defect))
        call check(nf90_put_var(ncid, v(6), diag%eps_fb))
        call check(nf90_put_var(ncid, v(7), diag%pprime))
        call check(nf90_put_var(ncid, v(8), diag%cip))
        call check(nf90_put_var(ncid, v(9), diag%cgp))
        call check(nf90_put_var(ncid, v(10), diag%n_zeroed))
        if (nres > 0) then
            call check(nf90_put_var(ncid, v(11), isurf))
            call check(nf90_put_var(ncid, v(12), rm))
            call check(nf90_put_var(ncid, v(13), rn))
            call check(nf90_put_var(ncid, v(14), rd))
            call check(nf90_put_var(ncid, v(15), rf))
            call check(nf90_put_var(ncid, v(16), ra))
        end if
        call check(nf90_put_var(ncid, v(17), diag%d_cur_n))
    end subroutine write_diagnostics

    subroutine def(v, k, name, xtype, dims)
        integer, intent(inout) :: v(:), k
        character(*), intent(in) :: name
        integer, intent(in) :: xtype, dims(:)

        k = k + 1
        call check(nf90_def_var(ncid, name, xtype, dims, v(k)))
    end subroutine def

    subroutine def0(v, k, name)
        integer, intent(inout) :: v(:), k
        character(*), intent(in) :: name

        k = k + 1
        call check(nf90_def_var(ncid, name, NF90_DOUBLE, v(k)))
    end subroutine def0

    subroutine def_spectrum(v, k, prefix, dsurf, dmn)
        integer, intent(inout) :: v(:), k
        character(*), intent(in) :: prefix
        integer, intent(in) :: dsurf, dmn

        call def(v, k, prefix//'bmnc_b', NF90_DOUBLE, [dmn, dsurf])
        call def(v, k, prefix//'bmns_b', NF90_DOUBLE, [dmn, dsurf])
        call def(v, k, prefix//'numns_b', NF90_DOUBLE, [dmn, dsurf])
        call def(v, k, prefix//'numnc_b', NF90_DOUBLE, [dmn, dsurf])
        call def(v, k, prefix//'iota', NF90_DOUBLE, [dsurf])
        call def(v, k, prefix//'G', NF90_DOUBLE, [dsurf])
        call def(v, k, prefix//'I', NF90_DOUBLE, [dsurf])
        call def(v, k, prefix//'gradsmnc_b', NF90_DOUBLE, [dmn, dsurf])
        call def(v, k, prefix//'gradsmns_b', NF90_DOUBLE, [dmn, dsurf])
        call def(v, k, prefix//'psi_edge', NF90_DOUBLE, [dsurf])
        call def(v, k, prefix//'R', NF90_DOUBLE, [dsurf])
    end subroutine def_spectrum

    subroutine put_spectrum(v, sp)
        integer, intent(in) :: v(11)
        type(boozer_spectrum_t), intent(in) :: sp

        call check(nf90_put_var(ncid, v(1), sp%bmnc))
        call check(nf90_put_var(ncid, v(2), sp%bmns))
        call check(nf90_put_var(ncid, v(3), sp%numns))
        call check(nf90_put_var(ncid, v(4), sp%numnc))
        call check(nf90_put_var(ncid, v(5), sp%iota))
        call check(nf90_put_var(ncid, v(6), sp%g))
        call check(nf90_put_var(ncid, v(7), sp%i))
        call check(nf90_put_var(ncid, v(8), sp%gradsmnc))
        call check(nf90_put_var(ncid, v(9), sp%gradsmns))
        call check(nf90_put_var(ncid, v(10), sp%psi_edge))
        call check(nf90_put_var(ncid, v(11), sp%rmajor))
    end subroutine put_spectrum

    integer function dim_len(id, name)
        integer, intent(in) :: id
        character(*), intent(in) :: name
        integer :: varid, dimids(1)

        call check(nf90_inq_varid(id, name, varid))
        call check(nf90_inquire_variable(id, varid, dimids=dimids))
        call check(nf90_inquire_dimension(id, dimids(1), len=dim_len))
    end function dim_len

    integer function get_int(id, name)
        integer, intent(in) :: id
        character(*), intent(in) :: name
        integer :: varid

        call check(nf90_inq_varid(id, name, varid))
        call check(nf90_get_var(id, varid, get_int))
    end function get_int

    integer function get_int_opt(id, name, default)
        integer, intent(in) :: id, default
        character(*), intent(in) :: name
        integer :: varid

        get_int_opt = default
        if (nf90_inq_varid(id, name, varid) /= NF90_NOERR) return
        call check(nf90_get_var(id, varid, get_int_opt))
    end function get_int_opt

    subroutine get_var(id, name, a)
        integer, intent(in) :: id
        character(*), intent(in) :: name
        real(dp), intent(out) :: a(:)
        integer :: varid

        call check(nf90_inq_varid(id, name, varid))
        call check(nf90_get_var(id, varid, a))
    end subroutine get_var

    subroutine get_opt0(id, name, a)
        integer, intent(in) :: id
        character(*), intent(in) :: name
        real(dp), intent(inout) :: a
        integer :: varid

        if (nf90_inq_varid(id, name, varid) /= NF90_NOERR) return
        call check(nf90_get_var(id, varid, a))
    end subroutine get_opt0

    subroutine get_opt1(id, name, a)
        integer, intent(in) :: id
        character(*), intent(in) :: name
        real(dp), intent(inout) :: a(:)
        integer :: varid

        if (nf90_inq_varid(id, name, varid) /= NF90_NOERR) return
        call check(nf90_get_var(id, varid, a))
    end subroutine get_opt1

    subroutine get_opt2(id, name, a)
        integer, intent(in) :: id
        character(*), intent(in) :: name
        real(dp), intent(inout) :: a(:, :)
        integer :: varid

        if (nf90_inq_varid(id, name, varid) /= NF90_NOERR) return
        call check(nf90_get_var(id, varid, a))
    end subroutine get_opt2

    subroutine check(status)
        integer, intent(in) :: status

        if (status /= NF90_NOERR) then
            print *, trim(nf90_strerror(status))
            error stop 1
        end if
    end subroutine check

end program boozer_spectrum_ad
