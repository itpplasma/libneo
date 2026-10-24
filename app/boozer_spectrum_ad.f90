!> Boozer spectrum of a VMEC wout file with tangent or adjoint derivatives.
!>
!> Usage: boozer_spectrum_ad.x primal|jvp|vjp <wout.nc> <in.nc> <out.nc>
!>
!> in.nc holds s(nsurf), the integers mboz, nboz and grid_refinement, optionally
!> the integer lambda_cubic (1 selects lambda_half_to_full = 'cubic', default 0 =
!> 'linear'), and for
!> jvp the wout-shaped tangents rmnc, zmns, lmns (ns x mnmax) and iotaf (ns),
!> for vjp the cotangents bmnc_b, bmns_b, numns_b, numnc_b (nsurf x mnboz)
!> and iota, G, I (nsurf). Absent arrays are zero. out.nc receives xm_b, xn_b,
!> the primal spectrum and either its tangent (prefix d_) or the wout
!> cotangents rmnc, zmns, lmns and iotaf. Units: m, T and T m.
program boozer_spectrum_ad
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use netcdf
    use new_vmec_stuff_mod, only: lambda_half_to_full
    use boozer_vmec_ad, only: boozer_ad_context_t, boozer_wout_tangent_t, &
                              boozer_spectrum_t, boozer_ad_init, boozer_ad_primal, &
                              boozer_ad_jvp, boozer_ad_vjp, wout_tangent_alloc, &
                              spectrum_alloc
    implicit none

    character(len=1024) :: mode, wout, fin, fout
    type(boozer_ad_context_t) :: ctx
    type(boozer_wout_tangent_t) :: t
    type(boozer_spectrum_t) :: y, yd, yb
    real(dp), allocatable :: s(:)
    integer :: ncid, mboz, nboz, refine, nsurf
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

    call boozer_ad_init(ctx, trim(wout), mboz, nboz, grid_refinement=refine)
    call system_clock(c0, rate)
    select case (trim(mode))
    case ('primal')
        call boozer_ad_primal(ctx, s, y)
    case ('jvp')
        call wout_tangent_alloc(ctx, t)
        call get_opt2(ncid, 'rmnc', t%rmnc)
        call get_opt2(ncid, 'zmns', t%zmns)
        call get_opt2(ncid, 'lmns', t%lmns)
        call get_opt1(ncid, 'iotaf', t%iotaf)
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
        call boozer_ad_primal(ctx, s, y)
        call boozer_ad_vjp(ctx, s, yb, t)
    case default
        error stop 'mode must be primal, jvp or vjp'
    end select
    call system_clock(c1)
    call check(nf90_close(ncid))
    call write_output()

contains

    subroutine write_output()
        integer :: dsurf, dmn, dns, dmnmax, v(32), k

        call check(nf90_create(trim(fout), NF90_CLOBBER, ncid))
        call check(nf90_def_dim(ncid, 'surface', nsurf, dsurf))
        call check(nf90_def_dim(ncid, 'mnboz', ctx%mnboz, dmn))
        call check(nf90_put_att(ncid, NF90_GLOBAL, 'mode', trim(mode)))
        call check(nf90_put_att(ncid, NF90_GLOBAL, 'seconds', &
                                real(c1 - c0, dp)/real(rate, dp)))
        call check(nf90_put_att(ncid, NF90_GLOBAL, 'n_theta', ctx%grid%n_theta))
        call check(nf90_put_att(ncid, NF90_GLOBAL, 'n_phi', ctx%grid%n_phi))
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
        end if
        call check(nf90_enddef(ncid))
        call check(nf90_put_var(ncid, v(1), ctx%xm_b))
        call check(nf90_put_var(ncid, v(2), ctx%xn_b))
        call check(nf90_put_var(ncid, v(3), s))
        call put_spectrum(v(4:10), y)
        if (trim(mode) == 'jvp') call put_spectrum(v(11:17), yd)
        if (trim(mode) == 'vjp') then
            call check(nf90_put_var(ncid, v(11), t%rmnc))
            call check(nf90_put_var(ncid, v(12), t%zmns))
            call check(nf90_put_var(ncid, v(13), t%lmns))
            call check(nf90_put_var(ncid, v(14), t%iotaf))
        end if
        call check(nf90_close(ncid))
    end subroutine write_output

    subroutine def(v, k, name, xtype, dims)
        integer, intent(inout) :: v(:), k
        character(*), intent(in) :: name
        integer, intent(in) :: xtype, dims(:)

        k = k + 1
        call check(nf90_def_var(ncid, name, xtype, dims, v(k)))
    end subroutine def

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
    end subroutine def_spectrum

    subroutine put_spectrum(v, sp)
        integer, intent(in) :: v(7)
        type(boozer_spectrum_t), intent(in) :: sp

        call check(nf90_put_var(ncid, v(1), sp%bmnc))
        call check(nf90_put_var(ncid, v(2), sp%bmns))
        call check(nf90_put_var(ncid, v(3), sp%numns))
        call check(nf90_put_var(ncid, v(4), sp%numnc))
        call check(nf90_put_var(ncid, v(5), sp%iota))
        call check(nf90_put_var(ncid, v(6), sp%g))
        call check(nf90_put_var(ncid, v(7), sp%i))
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
