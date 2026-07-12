program test_spectre_basis
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use spectre_reader, only: spectre_data_t, load_spectre, free_spectre
    use spectre_basis, only: spectre_vecpot_t, eval_spectre_vector_potential
    use util_for_test, only: print_test, print_ok, print_fail

    implicit none

    character(len=1024) :: fixture, reffile
    type(spectre_data_t) :: data
    integer :: nfail, ierr

    nfail = 0
    call get_command_argument(1, fixture)
    call get_command_argument(2, reffile)
    if (len_trim(fixture) == 0) error stop 'usage: test_spectre_basis.x <h5> <ref.dat>'
    if (len_trim(reffile) == 0) error stop 'usage: test_spectre_basis.x <h5> <ref.dat>'

    call load_spectre(trim(fixture), data, ierr)
    if (ierr /= 0) error stop 'load_spectre failed'

    call test_reference_values
    call test_fd_convergence
    call test_polynomial_extension
    call test_axis_regularity

    call free_spectre(data)
    if (nfail > 0) error stop

contains

    subroutine test_reference_values
        type(spectre_vecpot_t) :: av
        integer :: u, io, lvol, nrow, nfail_before
        real(dp) :: s, th, ze, ath_ref, azt_ref, rest(5)
        character(len=256) :: header

        call print_test('vector potential matches SPECTRE reference to 1e-12')
        nfail_before = nfail
        nrow = 0

        open (newunit=u, file=trim(reffile), status='old', action='read')
        read (u, '(A)') header
        do
            read (u, *, iostat=io) lvol, s, th, ze, ath_ref, azt_ref, rest
            if (io /= 0) exit
            call eval_spectre_vector_potential(data, lvol, s, th, ze, av)
            call check_close('Ath', av%Ath, ath_ref, 1.0e-12_dp)
            call check_close('Azt', av%Azt, azt_ref, 1.0e-12_dp)
            nrow = nrow + 1
        end do
        close (u)

        if (nrow < 1) call fail('no reference rows read')
        call report(nfail_before)
    end subroutine test_reference_values

    subroutine test_fd_convergence
        integer :: nfail_before

        call print_test('analytic derivatives match central finite differences')
        nfail_before = nfail

        call check_point(1, -1.0_dp + 1.0e-3_dp, 0.3_dp, 0.15_dp)
        call check_point(1, -0.4_dp, 2.1_dp, 0.9_dp)
        call check_point(1, 0.3_dp, 1.0_dp, 0.5_dp)
        call check_point(2, -0.5_dp, 0.7_dp, 0.4_dp)
        call check_point(2, 0.1_dp, 2.0_dp, 1.0_dp)
        call check_point(2, 0.8_dp, 4.0_dp, 0.2_dp)
        call check_point(3, -0.6_dp, 0.5_dp, 0.3_dp)
        call check_point(3, 0.2_dp, 3.0_dp, 0.8_dp)
        call check_point(3, 0.7_dp, 1.5_dp, 0.6_dp)

        call report(nfail_before)
    end subroutine test_fd_convergence

    subroutine check_point(lvol, s, th, ze)
        integer, intent(in) :: lvol
        real(dp), intent(in) :: s, th, ze

        type(spectre_vecpot_t) :: av
        integer :: dirs(2, 3), k

        call eval_spectre_vector_potential(data, lvol, s, th, ze, av)

        ! first derivatives: slot dir maps to coordinate dir
        do k = 1, 3
            call check_slot(1, [av%dAth(k), av%dAzt(k)], &
                            fd1(lvol, s, th, ze, k, 1.0e-4_dp), &
                            fd1(lvol, s, th, ze, k, 1.0e-5_dp), &
                            fd1(lvol, s, th, ze, k, 1.0e-6_dp))
        end do

        ! pure second derivatives ss, tt, zz -> slots 1, 4, 6
        call check_second_pure(lvol, s, th, ze, 1, [av%d2Ath(1), av%d2Azt(1)])
        call check_second_pure(lvol, s, th, ze, 2, [av%d2Ath(4), av%d2Azt(4)])
        call check_second_pure(lvol, s, th, ze, 3, [av%d2Ath(6), av%d2Azt(6)])

        ! mixed second derivatives st, sz, tz -> slots 2, 3, 5
        dirs(:, 1) = [1, 2]
        dirs(:, 2) = [1, 3]
        dirs(:, 3) = [2, 3]
        call check_mixed(lvol, s, th, ze, dirs(:, 1), [av%d2Ath(2), av%d2Azt(2)])
        call check_mixed(lvol, s, th, ze, dirs(:, 2), [av%d2Ath(3), av%d2Azt(3)])
        call check_mixed(lvol, s, th, ze, dirs(:, 3), [av%d2Ath(5), av%d2Azt(5)])
    end subroutine check_point

    subroutine check_second_pure(lvol, s, th, ze, dir, analytic)
        integer, intent(in) :: lvol, dir
        real(dp), intent(in) :: s, th, ze, analytic(2)

        call check_slot(2, analytic, fd2(lvol, s, th, ze, dir, 1.0e-4_dp), &
                        fd2(lvol, s, th, ze, dir, 1.0e-5_dp), &
                        fd2(lvol, s, th, ze, dir, 1.0e-6_dp))
    end subroutine check_second_pure

    subroutine check_mixed(lvol, s, th, ze, dir, analytic)
        integer, intent(in) :: lvol, dir(2)
        real(dp), intent(in) :: s, th, ze, analytic(2)

        call check_slot(2, analytic, fdm(lvol, s, th, ze, dir, 1.0e-4_dp), &
                        fdm(lvol, s, th, ze, dir, 1.0e-5_dp), &
                        fdm(lvol, s, th, ze, dir, 1.0e-6_dp))
    end subroutine check_mixed

    ! order 1: first derivative (central diff is O(h^2), truncation-limited at
    ! these steps, so error decreases from h=1e-4 to h=1e-5).
    ! order 2: second derivative (central diff is roundoff-limited below the
    ! near-optimal step h~1e-4, so it is checked tightly at h=1e-4 and only
    ! bounded at h=1e-5 rather than required to keep decreasing).
    subroutine check_slot(order, analytic, num1, num2, num3)
        integer, intent(in) :: order
        real(dp), intent(in) :: analytic(2), num1(2), num2(2), num3(2)

        real(dp) :: e1, e2, e3, scale
        integer :: c

        do c = 1, 2
            e1 = abs(analytic(c) - num1(c))
            e2 = abs(analytic(c) - num2(c))
            e3 = abs(analytic(c) - num3(c))
            scale = max(1.0_dp, abs(analytic(c)))
            if (order == 1) then
                if (e2 >= 1.0e-7_dp*scale) call fail('order-1 error too large')
                if (e3 >= 1.0e-7_dp*scale) call fail('order-1 h=1e-6 too large')
                if (e2 >= e1) then
                    if (e2 >= 1.0e-10_dp) call fail('order-1 not decreasing')
                end if
            else
                if (e1 >= 1.0e-6_dp*scale) call fail('order-2 error too large')
                if (e2 >= 1.0e-5_dp*scale) call fail('order-2 h=1e-5 too large')
            end if
        end do
    end subroutine check_slot

    function aval(lvol, s, th, ze) result(a)
        integer, intent(in) :: lvol
        real(dp), intent(in) :: s, th, ze
        real(dp) :: a(2)
        type(spectre_vecpot_t) :: av

        call eval_spectre_vector_potential(data, lvol, s, th, ze, av)
        a = [av%Ath, av%Azt]
    end function aval

    function shifted(lvol, s, th, ze, dir, delta) result(a)
        integer, intent(in) :: lvol, dir
        real(dp), intent(in) :: s, th, ze, delta
        real(dp) :: a(2), cs, ct, cz

        cs = s
        ct = th
        cz = ze
        select case (dir)
        case (1)
            cs = cs + delta
        case (2)
            ct = ct + delta
        case (3)
            cz = cz + delta
        end select
        a = aval(lvol, cs, ct, cz)
    end function shifted

    function shifted2(lvol, s, th, ze, dir, d1v, d2v) result(a)
        integer, intent(in) :: lvol, dir(2)
        real(dp), intent(in) :: s, th, ze, d1v, d2v
        real(dp) :: a(2), cs, ct, cz

        cs = s
        ct = th
        cz = ze
        call bump(cs, ct, cz, dir(1), d1v)
        call bump(cs, ct, cz, dir(2), d2v)
        a = aval(lvol, cs, ct, cz)
    end function shifted2

    subroutine bump(cs, ct, cz, dir, delta)
        real(dp), intent(inout) :: cs, ct, cz
        integer, intent(in) :: dir
        real(dp), intent(in) :: delta

        select case (dir)
        case (1)
            cs = cs + delta
        case (2)
            ct = ct + delta
        case (3)
            cz = cz + delta
        end select
    end subroutine bump

    function fd1(lvol, s, th, ze, dir, h) result(d)
        integer, intent(in) :: lvol, dir
        real(dp), intent(in) :: s, th, ze, h
        real(dp) :: d(2)

        d = (shifted(lvol, s, th, ze, dir, h) &
             - shifted(lvol, s, th, ze, dir, -h))/(2.0_dp*h)
    end function fd1

    function fd2(lvol, s, th, ze, dir, h) result(d)
        integer, intent(in) :: lvol, dir
        real(dp), intent(in) :: s, th, ze, h
        real(dp) :: d(2)

        d = (shifted(lvol, s, th, ze, dir, h) - 2.0_dp*aval(lvol, s, th, ze) &
             + shifted(lvol, s, th, ze, dir, -h))/h**2
    end function fd2

    function fdm(lvol, s, th, ze, dir, h) result(d)
        integer, intent(in) :: lvol, dir(2)
        real(dp), intent(in) :: s, th, ze, h
        real(dp) :: d(2)

        d = (shifted2(lvol, s, th, ze, dir, h, h) &
             - shifted2(lvol, s, th, ze, dir, h, -h) &
             - shifted2(lvol, s, th, ze, dir, -h, h) &
             + shifted2(lvol, s, th, ze, dir, -h, -h))/(4.0_dp*h**2)
    end function fdm

    subroutine test_polynomial_extension
        type(spectre_vecpot_t) :: av
        real(dp) :: outside(2), inside(2), num(2)
        integer :: nfail_before

        call print_test('interior polynomial extension at s>1 stays finite')
        nfail_before = nfail

        call eval_spectre_vector_potential(data, 2, 1.05_dp, 0.6_dp, 0.3_dp, av)
        outside = [av%Ath, av%Azt]
        if (.not. all(ieee_is_finite(outside))) call fail('non-finite at s=1.05')
        if (.not. ieee_is_finite(av%dAth(1))) call fail('non-finite dAth at s=1.05')
        if (.not. ieee_is_finite(av%d2Ath(1))) call fail('non-finite d2Ath at s=1.05')

        ! central FD around s=1.05 uses points on both sides of the interface;
        ! matching the analytic radial derivative there proves no clamping.
        num = fd1(2, 1.05_dp, 0.6_dp, 0.3_dp, 1, 1.0e-5_dp)
        call check_close('ext dAth/ds', av%dAth(1), num(1), 1.0e-6_dp)
        call check_close('ext dAzt/ds', av%dAzt(1), num(2), 1.0e-6_dp)

        inside = aval(2, 0.999_dp, 0.6_dp, 0.3_dp)
        if (any(abs(outside - inside) > 1.0_dp)) call fail('extension jumps across s=1')

        call report(nfail_before)
    end subroutine test_polynomial_extension

    subroutine test_axis_regularity
        type(spectre_vecpot_t) :: av
        real(dp) :: vals(20)
        integer :: nfail_before

        call print_test('axis volume finite at s = -1 + 1e-12')
        nfail_before = nfail

        call eval_spectre_vector_potential(data, 1, -1.0_dp + 1.0e-12_dp, &
                                           0.4_dp, 0.7_dp, av)
        vals = [av%Ath, av%Azt, av%dAth, av%dAzt, av%d2Ath, av%d2Azt]
        if (.not. all(ieee_is_finite(vals))) call fail('non-finite output on axis')

        call report(nfail_before)
    end subroutine test_axis_regularity

    subroutine check_close(name, got, want, tol)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: got, want, tol

        if (abs(got - want) > tol) then
            print *, '    ', name, ' got ', got, ' want ', want, &
                ' diff ', abs(got - want)
            nfail = nfail + 1
        end if
    end subroutine check_close

    subroutine fail(message)
        character(len=*), intent(in) :: message

        print *, '    ', message
        nfail = nfail + 1
    end subroutine fail

    subroutine report(nfail_before)
        integer, intent(in) :: nfail_before

        if (nfail == nfail_before) then
            call print_ok
        else
            call print_fail
        end if
    end subroutine report

end program test_spectre_basis
