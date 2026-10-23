!> Equilibrium-quality checks of boozer_hybrid_quality on the coarse public
!> test wouts (circular tokamak ns 17, li383 ns 49):
!>  1. metrics against an independent reference (boozer-bench
!>     reports/scripts/libneo_quality/reference.py: lambda from the stationary
!>     magnetic energy by jax.hessian instead of the assembled weak form, numpy
!>     FFT derivatives, pyrabe spline rows) on the check surfaces 1, 8, 16;
!>     circular surface 1 uses an independent s**(m/2+j), j=0..2, endpoint
!>     solve matching the pyrabe value and first two derivatives below the
!>     first positive source knot. The other reference surfaces are unchanged;
!>  2. the default tolerances reject both files with the expected first code
!>     and an actionable message (surface, metric, cause, resonance, action,
!>     override), returned in status instead of stopping;
!>  3. allow_poor_equilibrium keeps the code with overridden = .true.;
!>  4. negative tolerances make every metric report-only (code SFL_OK);
!>  5. boozer_ad_init records the check in ctx%quality.
program test_sfl_quality
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use boozer_hybrid_quality, only: sfl_quality_config_t, sfl_quality_t, sfl_status_t, &
                                     sfl_quality_check, SFL_OK, SFL_ERR_RESOLUTION, &
                                     SFL_ERR_FORCE_BALANCE
    use boozer_vmec_ad, only: boozer_ad_context_t, boozer_ad_init
    implicit none

    logical :: failed
    integer, parameter :: KS(3) = [1, 8, 16]

    failed = .false.
    call check_reference('wout_circular_tokamak.nc', 7.368343302587206e-06_dp, &
                         [7.878128988601874e-05_dp, 2.285497830151195e-05_dp, &
                          2.975357900882880e-04_dp], &
                         [8.410217911893529e-04_dp, 3.074045040945442e-05_dp, &
                          3.491247588011836e-04_dp], &
                         [1.663730685320132e-03_dp, 1.151512650679340e-03_dp, &
                          1.312227811239188e-03_dp])
    call check_reference('wout_li383.nc', 3.319349848838740e-05_dp, &
                         [2.418495592286882e-03_dp, 2.161437712747278e-04_dp, &
                          8.351875815271113e-05_dp], &
                         [1.598540744261418e-02_dp, 1.136380705172979e-02_dp, &
                          4.684173852325216e-02_dp], &
                         [3.619139612569459e-02_dp, 1.243284843547301e-04_dp, &
                          1.206248881889407e-03_dp])
    call check_status('wout_circular_tokamak.nc', SFL_ERR_FORCE_BALANCE, 'eps_fb')
    call check_status('wout_li383.nc', SFL_ERR_RESOLUTION, 'eps_rad')
    call check_override('wout_li383.nc')
    call check_report_only('wout_li383.nc')
    call check_context('wout_circular_tokamak.nc')
    if (failed) error stop 'test_sfl_quality failed'
    print *, 'test_sfl_quality: all checks passed'

contains

    subroutine report(name, err, tol)
        character(*), intent(in) :: name
        real(dp), intent(in) :: err, tol

        if (err <= tol) then
            write (*, '(a, a, es10.3, a, es10.3)') '  ok   ', name//': ', err, ' <= ', tol
        else
            write (*, '(a, a, es10.3, a, es10.3)') '  FAIL ', name//': ', err, ' > ', tol
            failed = .true.
        end if
    end subroutine report

    subroutine expect(name, cond)
        character(*), intent(in) :: name
        logical, intent(in) :: cond

        if (cond) then
            print '(a, a)', '  ok   ', name
        else
            print '(a, a)', '  FAIL ', name
            failed = .true.
        end if
    end subroutine expect

    !> eps_fb uses central differences in s on both sides, so its ratio agrees
    !> to the cancellation of its three terms only (1e-4 relative).
    subroutine check_reference(wout, rad, lam, js, fb)
        character(*), intent(in) :: wout
        real(dp), intent(in) :: rad, lam(3), js(3), fb(3)
        type(sfl_quality_t) :: q
        type(sfl_status_t) :: st

        print '(a)', wout//': metrics against the independent reference'
        call sfl_quality_check(wout, q, status=st)
        call report('eps_rad', abs(q%eps_rad/rad - 1.0_dp), 1.0e-10_dp)
        call report('eps_lam', maxval(abs(q%eps_lam(KS)/lam - 1.0_dp)), 1.0e-6_dp)
        call report('eps_js', maxval(abs(q%eps_js(KS)/js - 1.0_dp)), 1.0e-7_dp)
        call report('eps_fb', maxval(abs(q%eps_fb(KS)/fb - 1.0_dp)), 1.0e-3_dp)
    end subroutine check_reference

    subroutine check_status(wout, code, metric)
        character(*), intent(in) :: wout, metric
        integer, intent(in) :: code
        type(sfl_quality_t) :: q
        type(sfl_status_t) :: st

        print '(a)', wout//': default tolerances'
        call sfl_quality_check(wout, q, status=st)
        call expect('first failing code', st%code == code)
        call expect('not overridden', .not. st%overridden)
        call expect('message names the metric', index(st%message, metric//' ') > 0)
        call expect('message names the surface', index(st%message, 'at s = ') > 0)
        call expect('message names the nearest resonance', &
                    index(st%message, 'nearest resonance (m,n)') > 0)
        call expect('message names the action', index(st%message, 'Raise ') > 0)
        call expect('message names the override', &
                    index(st%message, 'allow_poor_equilibrium') > 0)
        call expect('status carries surface and value', st%s > 0.0_dp .and. &
                    st%value > st%tol)
    end subroutine check_status

    subroutine check_override(wout)
        character(*), intent(in) :: wout
        type(sfl_quality_config_t) :: cfg
        type(sfl_quality_t) :: q

        print '(a)', wout//': override'
        cfg%allow_poor_equilibrium = .true.
        cfg%override_reason = 'test_sfl_quality'
        call sfl_quality_check(wout, q, cfg)
        call expect('returns without status, code kept', q%status%code /= SFL_OK)
        call expect('overridden', q%overridden .and. q%status%overridden)
        call expect('reason kept', trim(q%override_reason) == 'test_sfl_quality')
    end subroutine check_override

    subroutine check_report_only(wout)
        character(*), intent(in) :: wout
        type(sfl_quality_config_t) :: cfg
        type(sfl_quality_t) :: q

        print '(a)', wout//': report-only tolerances'
        cfg%tol = -1.0_dp
        call sfl_quality_check(wout, q, cfg)
        call expect('code SFL_OK', q%status%code == SFL_OK)
        call expect('metrics still computed', q%eps_rad > 0.0_dp .and. all(q%eps_lam > 0.0_dp))
    end subroutine check_report_only

    subroutine check_context(wout)
        character(*), intent(in) :: wout
        type(boozer_ad_context_t) :: ctx
        type(sfl_quality_config_t) :: cfg
        type(sfl_status_t) :: st

        print '(a)', wout//': boozer_ad_init'
        call boozer_ad_init(ctx, wout, 4, 0, status=st)
        call expect('status returned by boozer_ad_init', st%code == SFL_ERR_FORCE_BALANCE)
        call expect('record in ctx%quality', ctx%quality%status%code == st%code .and. &
                    size(ctx%quality%eps_lam) == 16)
        cfg%allow_poor_equilibrium = .true.
        cfg%override_reason = 'test_sfl_quality'
        call boozer_ad_init(ctx, wout, 4, 0, quality=cfg)
        call expect('override recorded in ctx%quality', ctx%quality%overridden)
    end subroutine check_context

end program test_sfl_quality
