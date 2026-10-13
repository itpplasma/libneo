program test_magfie_cyl_status
    ! magfie_with_status reports points outside the LCFS / convex boundary via
    ! its optional status argument, prints nothing, and leaves the returned
    ! values identical to the legacy magfie entry point.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_set_flag, &
                                             ieee_set_halting_mode, ieee_invalid, &
                                             ieee_divide_by_zero
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use input_files, only: convexfile
    use field_mod, only: icall, iequil, ipert
    use field_eq_mod, only: psi_axis, psi_sep, ierrfield, allow_sol
    use field_sub, only: psif
    use magfie_cyl_mod, only: magfie_with_status, MAGFIE_OK, MAGFIE_OUTSIDE_LCFS, &
                              MAGFIE_OUTSIDE_BOUNDARY

    implicit none

    integer :: u, i, nerr, status
    real(dp) :: x(3)

    nerr = 0

    ! The zero field below gives 0/0 in magfie; debug builds trap it otherwise.
    call ieee_set_halting_mode(ieee_invalid, .false.)
    call ieee_set_halting_mode(ieee_divide_by_zero, .false.)

    ! Circular convex wall: radius 1 around (R,Z) = (2,0).
    convexfile = 'test_magfie_cyl_convex.dat'
    open (newunit=u, file=convexfile, status='replace', action='write')
    do i = 0, 71
        write (u, *) 2.0_dp + cos(2.0_dp*acos(-1.0_dp)*i/72.0_dp), &
            sin(2.0_dp*acos(-1.0_dp)*i/72.0_dp)
    end do
    close (u)

    ! Zero equilibrium field; only the status logic is exercised.
    icall = 1
    iequil = 0
    ipert = 0
    allow_sol = .false.
    psi_axis = 0.0_dp
    psi_sep = 1.0_dp

    psif = 0.5_dp
    call evaluate([2.0_dp, 0.0_dp, 0.0_dp], MAGFIE_OK, 0)
    psif = 2.0_dp
    call evaluate([2.3_dp, 0.0_dp, 0.1_dp], MAGFIE_OUTSIDE_LCFS, 1)
    allow_sol = .true.
    call evaluate([2.3_dp, 0.0_dp, 0.1_dp], MAGFIE_OK, 0)
    allow_sol = .false.
    call evaluate([4.0_dp, 0.0_dp, 0.0_dp], MAGFIE_OUTSIDE_BOUNDARY, 1)

    ! Do not leave the expected flags raised at program exit.
    call ieee_set_flag([ieee_invalid, ieee_divide_by_zero], .false.)

    if (nerr /= 0) error stop 'test_magfie_cyl_status failed'

contains

    subroutine evaluate(xin, expected_status, expected_ierr)
        real(dp), intent(in) :: xin(3)
        integer, intent(in) :: expected_status, expected_ierr
        real(dp) :: b1, s1, d1(3), c1(3), t1(3), r1(3)
        real(dp) :: b2, s2, d2(3), c2(3), t2(3), r2(3)

        x = xin
        status = -99
        call magfie_with_status(x, b1, s1, d1, c1, t1, r1, status)
        call check(status == expected_status, 'wrong status')
        call check(ierrfield == expected_ierr, 'wrong ierrfield')
        ierrfield = -99
        call magfie_with_status(x, b2, s2, d2, c2, t2, r2)
        call check(ierrfield == expected_ierr, 'ierrfield without status')
        call check(same(b1, b2) .and. same(s1, s2), 'bmod/sqrtg differ')
        call check(all([(same(d1(i), d2(i)), i=1, 3)]), 'bder differs')
        call check(all([(same(c1(i), c2(i)), i=1, 3)]), 'hcovar differs')
        call check(all([(same(t1(i), t2(i)), i=1, 3)]), 'hctrvr differs')
        call check(all([(same(r1(i), r2(i)), i=1, 3)]), 'hcurl differs')
        call check(same(s1, xin(1)), 'sqrtg must equal R')
    end subroutine evaluate

    pure logical function same(a, b)
        real(dp), intent(in) :: a, b
        same = (a == b) .or. (ieee_is_nan(a) .and. ieee_is_nan(b))
    end function same

    subroutine check(ok, msg)
        logical, intent(in) :: ok
        character(*), intent(in) :: msg

        if (.not. ok) then
            nerr = nerr + 1
            write (*, '(2a)') 'FAIL: ', msg
        end if
    end subroutine check

end program test_magfie_cyl_status
