!> The optional cubic half-to-full mesh map of the VMEC lambda harmonics is
!> exact for cubic polynomials in s, while the default two-point average
!> carries an O(hs^2) error. Synthetic data, no wout needed.
program test_lambda_half_to_full
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use vmecin_sub, only: lambda_half_to_full_cubic
    implicit none

    integer, parameter :: kpar = 20, nmodes = 2
    real(dp) :: lhalf(nmodes, 0:kpar), lfull(nmodes, 0:kpar), hs, s, err
    integer :: i

    hs = 1.0_dp/kpar
    lhalf = 0.0_dp
    do i = 1, kpar
        s = (i - 0.5_dp)*hs
        lhalf(1, i) = poly(s)
        lhalf(2, i) = -2.0_dp*poly(s)
    end do
    lfull = 0.0_dp
    call lambda_half_to_full_cubic(lhalf, lfull, kpar)

    err = 0.0_dp
    do i = 1, kpar
        s = i*hs
        err = max(err, abs(lfull(1, i) - poly(s)), abs(lfull(2, i) + 2.0_dp*poly(s)))
    end do
    print '(A,ES10.3)', 'max error of cubic half-to-full map: ', err
    if (err > 1.0e-13_dp .or. lfull(1, 0) /= 0.0_dp) then
        print *, 'FAIL'
        error stop 1
    end if
    print *, 'PASS'

contains

    pure real(dp) function poly(s)
        real(dp), intent(in) :: s
        poly = 0.3_dp - 1.1_dp*s + 0.7_dp*s**2 + 2.3_dp*s**3
    end function poly

end program test_lambda_half_to_full
