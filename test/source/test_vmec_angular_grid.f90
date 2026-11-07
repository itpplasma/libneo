!> The VMEC spline grid spans one field period in phi, so it must resolve the
!> per-period toroidal mode number n_max/nper as finely as the poloidal grid
!> resolves m_max. Sampling n_max (which includes nper) per period multiplies
!> the spline arrays by nper.
program test_vmec_angular_grid
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use new_vmec_stuff_mod, only: netcdffile, multharm, nper, n_theta, n_phi, sR
    use spline_vmec_sub, only: vmec_angular_grid, spline_vmec_data
    implicit none

    real(dp), parameter :: TWOPI = 2.0_dp*acos(-1.0_dp)
    integer :: failures

    failures = 0
    call check_grid('five periods', m_max=12, n_max=60, nper_in=5, mult=3, &
                    n_phi_expected=37)
    call check_grid('one period', m_max=8, n_max=6, nper_in=1, mult=5, &
                    n_phi_expected=31)
    call check_grid('axisymmetric', m_max=8, n_max=0, nper_in=1, mult=5, &
                    n_phi_expected=11)
    call check_spline_arrays()

    if (failures /= 0) error stop 'test_vmec_angular_grid FAILED'
    print '(a)', 'test_vmec_angular_grid: passed'

contains

    subroutine check_grid(label, m_max, n_max, nper_in, mult, n_phi_expected)
        character(len=*), intent(in) :: label
        integer, intent(in) :: m_max, n_max, nper_in, mult, n_phi_expected
        integer :: nth, nph
        real(dp) :: h_t, h_p

        call vmec_angular_grid(m_max, n_max, mult, nper_in, 5, nth, nph, h_t, h_p)
        call expect(label//': n_theta', nth == m_max*mult + 1)
        call expect(label//': n_phi', nph == n_phi_expected)
        call expect(label//': grid spans one period', &
                    abs(h_p*(nph - 1)*nper_in - TWOPI) < 1.0e-12_dp)
    end subroutine check_grid

    !> Public LandremanPaul2021 QA equilibrium: nfp = 2, max xm = 7, max xn = 16.
    subroutine check_spline_arrays()
        netcdffile = 'wout.nc'
        multharm = 3
        call spline_vmec_data
        call expect('QA nper', nper == 2)
        call expect('QA n_theta', n_theta == 7*3 + 1)
        call expect('QA n_phi per field period', n_phi == (16/2)*3 + 1)
        call expect('QA spline array size', size(sR, 6) == n_phi)
    end subroutine check_spline_arrays

    subroutine expect(label, ok)
        character(len=*), intent(in) :: label
        logical, intent(in) :: ok

        if (.not. ok) then
            failures = failures + 1
            print '(2a)', 'FAIL: ', label
        end if
    end subroutine expect

end program test_vmec_angular_grid
