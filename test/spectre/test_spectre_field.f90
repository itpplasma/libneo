program test_spectre_field
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_spectre, only: spectre_field_t, create_spectre_field
    use util_for_test, only: print_test, print_ok, print_fail

    implicit none

    real(dp), parameter :: twopi = 8.0_dp*atan(1.0_dp)

    character(len=1024) :: fixture, reffile
    type(spectre_field_t) :: field
    integer :: ierr, Mvol, nfail

    nfail = 0
    call get_command_argument(1, fixture)
    call get_command_argument(2, reffile)
    if (len_trim(fixture) == 0 .or. len_trim(reffile) == 0) then
        print *, 'usage: test_spectre_field.x <spectre.h5> <reference_eval.dat>'
        error stop
    end if

    call create_spectre_field(field, trim(fixture), ierr)
    if (ierr /= 0) then
        print *, 'create_spectre_field failed, ierr =', ierr
        error stop
    end if
    Mvol = field%data%Mvol

    call test_reference_field
    call test_div_b_vanishes
    call test_interfaces_are_flux_surfaces
    call test_axis_toroidal_flux_map
    call test_beltrami_residual
    call test_reload_consistency

    if (nfail > 0) error stop

contains

    subroutine eval_bctr(x, Bctr, Bcov, Bmod, sqrtg, g)
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: Bctr(3), Bcov(3), Bmod, sqrtg, g(3, 3)
        real(dp) :: Acov(3), hcov(3), sqgBctr(3), ginv(3, 3)

        call field%evaluate(x, Acov, hcov, Bmod, sqgBctr)
        call field%coords%metric_tensor(x, g, ginv, sqrtg)
        Bctr = sqgBctr/sqrtg
        Bcov = hcov*Bmod
    end subroutine eval_bctr

    subroutine test_reference_field
        integer :: iunit, io, lvol, nrows
        real(dp) :: col(11), s, theta, zeta, rho_g, x(3)
        real(dp) :: Bctr(3), Bcov(3), Bmod, sqrtg, g(3, 3)
        real(dp) :: worst_axis, worst_int, rel(4), w
        character(len=256) :: line

        call print_test('field matches SPECTRE reference (Bmod, B^theta/zeta/rho)')
        open (newunit=iunit, file=trim(reffile), status='old', action='read')
        read (iunit, '(A)') line
        nrows = 0
        worst_axis = 0.0_dp
        worst_int = 0.0_dp
        do
            read (iunit, *, iostat=io) col
            if (io /= 0) exit
            nrows = nrows + 1
            lvol = nint(col(1))
            s = col(2); theta = col(3); zeta = col(4)
            rho_g = real(lvol - 1, dp) + 0.5_dp*(s + 1.0_dp)
            x = [rho_g, theta, zeta]

            call eval_bctr(x, Bctr, Bcov, Bmod, sqrtg, g)

            rel(1) = abs(Bmod - col(11))/abs(col(11))
            rel(2) = abs(Bctr(2) - col(8))/abs(col(8))
            rel(3) = abs(Bctr(3) - col(9))/abs(col(9))
            rel(4) = abs(Bctr(1) - 0.5_dp*col(7))/abs(0.5_dp*col(7))
            w = maxval(rel)
            if (lvol == 1) then
                worst_axis = max(worst_axis, w)
            else
                worst_int = max(worst_int, w)
            end if
        end do
        close (iunit)

        print *, '    nrows =', nrows, ' worst_interior =', worst_int, &
            ' worst_axis =', worst_axis
        if (nrows == 90 .and. worst_int < 1.0e-10_dp &
            .and. worst_axis < 1.0e-10_dp) then
            call print_ok
        else
            nfail = nfail + 1
            call print_fail
        end if
    end subroutine test_reference_field

    subroutine test_div_b_vanishes
        real(dp), parameter :: h = 1.0e-4_dp
        real(dp), parameter :: frac(5) = [0.15_dp, 0.35_dp, 0.5_dp, 0.65_dp, 0.85_dp]
        integer :: iv, ip
        real(dp) :: x(3), d(3), scale, res, worst

        call print_test('div B = 0 (FD of sqrtg*B^i vanishes)')
        worst = 0.0_dp
        do iv = 1, Mvol
            do ip = 1, 5
                x(1) = real(iv - 1, dp) + frac(ip)
                x(2) = 0.7_dp + 0.31_dp*real(ip, dp)
                x(3) = 0.2_dp + 0.17_dp*real(ip, dp)
                d(1) = fd_sqgb(x, 1, 1, h)
                d(2) = fd_sqgb(x, 2, 2, h)
                d(3) = fd_sqgb(x, 3, 3, h)
                scale = abs(d(1)) + abs(d(2)) + abs(d(3))
                res = abs(d(1) + d(2) + d(3))/max(scale, tiny(1.0_dp))
                worst = max(worst, res)
            end do
        end do

        print *, '    worst relative div B =', worst
        if (worst < 1.0e-6_dp) then
            call print_ok
        else
            nfail = nfail + 1
            call print_fail
        end if
    end subroutine test_div_b_vanishes

    function fd_sqgb(x, comp, dir, h) result(deriv)
        real(dp), intent(in) :: x(3)
        integer, intent(in) :: comp, dir
        real(dp), intent(in) :: h
        real(dp) :: deriv, xp(3), xm(3), sp(3), sm(3)

        xp = x; xm = x
        xp(dir) = x(dir) + h
        xm(dir) = x(dir) - h
        call eval_sqgb(xp, sp)
        call eval_sqgb(xm, sm)
        deriv = (sp(comp) - sm(comp))/(2.0_dp*h)
    end function fd_sqgb

    subroutine eval_sqgb(x, sqgBctr)
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: sqgBctr(3)
        real(dp) :: Acov(3), hcov(3), Bmod

        call field%evaluate(x, Acov, hcov, Bmod, sqgBctr)
    end subroutine eval_sqgb

    subroutine test_interfaces_are_flux_surfaces
        real(dp), parameter :: eps = 1.0e-6_dp
        integer :: k, it, iz, side
        real(dp) :: x(3), Bctr(3), Bcov(3), Bmod, sqrtg, g(3, 3)
        real(dp) :: ratio, worst

        call print_test('interfaces are flux surfaces (B^rho ~ 0 both sides)')
        worst = 0.0_dp
        do k = 1, Mvol - 1
            do side = 1, 2
                do it = 0, 5
                    do iz = 0, 2
                        x(2) = twopi*real(it, dp)/6.0_dp
                        x(3) = twopi*real(iz, dp)/3.0_dp
                        if (side == 1) then
                            x(1) = real(k, dp) - eps
                        else
                            x(1) = real(k, dp) + eps
                        end if
                        call eval_bctr(x, Bctr, Bcov, Bmod, sqrtg, g)
                        ratio = abs(Bctr(1))/max(abs(Bctr(2)), tiny(1.0_dp))
                        worst = max(worst, ratio)
                    end do
                end do
            end do
        end do

        print *, '    worst |B^rho|/|B^theta| at interfaces =', worst
        if (worst < 1.0e-3_dp) then
            call print_ok
        else
            nfail = nfail + 1
            call print_fail
        end if
    end subroutine test_interfaces_are_flux_surfaces

    subroutine test_axis_toroidal_flux_map
        integer :: i, jerr
        real(dp) :: label, previous, rho, rho_back, target, upper

        call print_test('axis toroidal-flux label is monotone and invertible')
        call field%axis_toroidal_flux_label(0.0_dp, label, jerr)
        if (jerr /= 0 .or. abs(label) > 1.0e-14_dp) then
            nfail = nfail + 1
            call print_fail
            return
        end if
        call field%axis_toroidal_flux_label(1.0_dp, upper, jerr)
        if (jerr /= 0 .or. abs(upper - field%data%tflux(1)/ &
            field%data%tflux(Mvol)) > 1.0e-14_dp) then
            nfail = nfail + 1
            call print_fail
            return
        end if

        previous = label
        do i = 1, 100
            rho = real(i, dp)/100.0_dp
            call field%axis_toroidal_flux_label(rho, label, jerr)
            if (jerr /= 0 .or. label <= previous) then
                nfail = nfail + 1
                call print_fail
                return
            end if
            previous = label
        end do
        do i = 0, 10
            target = upper*real(i, dp)/10.0_dp
            call field%axis_rho_from_toroidal_flux(target, rho, jerr)
            call field%axis_toroidal_flux_label(rho, label, jerr)
            call field%axis_rho_from_toroidal_flux(label, rho_back, jerr)
            if (jerr /= 0 .or. abs(label - target) > 1.0e-13_dp .or. &
                abs(rho_back - rho) > 1.0e-13_dp) then
                print *, '    target, rho, label, rho_back, ierr =', &
                    target, rho, label, rho_back, jerr
                nfail = nfail + 1
                call print_fail
                return
            end if
        end do
        call print_ok
    end subroutine test_axis_toroidal_flux_map

    subroutine test_beltrami_residual
        real(dp), parameter :: h = 1.0e-4_dp
        real(dp), parameter :: frac(3) = [0.3_dp, 0.5_dp, 0.7_dp]
        integer :: iv, ip
        real(dp) :: x(3), Bctr(3), Bcov(3), Bmod, sqrtg, g(3, 3)
        real(dp) :: curlB(3), res(3), mu, rr, worst

        call print_test('Beltrami residual |curl B - mu B|/|B|')
        worst = 0.0_dp
        do iv = 1, Mvol
            mu = field%data%mu(iv)
            do ip = 1, 3
                x(1) = real(iv - 1, dp) + frac(ip)
                x(2) = 0.9_dp + 0.4_dp*real(ip, dp)
                x(3) = 0.3_dp + 0.2_dp*real(ip, dp)
                call eval_bctr(x, Bctr, Bcov, Bmod, sqrtg, g)
                call curl_bcov(x, h, sqrtg, curlB)
                res = curlB - mu*Bctr
                rr = sqrt(dot_product(res, matmul(g, res)))/Bmod
                worst = max(worst, rr)
            end do
        end do

        print *, '    worst Beltrami residual =', worst
        if (worst < 1.0e-2_dp) then
            call print_ok
        else
            nfail = nfail + 1
            call print_fail
        end if
    end subroutine test_beltrami_residual

    subroutine curl_bcov(x, h, sqrtg, curlB)
        real(dp), intent(in) :: x(3), h, sqrtg
        real(dp), intent(out) :: curlB(3)
        real(dp) :: dBcov(3, 3)
        integer :: dir

        do dir = 1, 3
            dBcov(:, dir) = fd_bcov(x, dir, h)
        end do
        curlB(1) = (dBcov(3, 2) - dBcov(2, 3))/sqrtg
        curlB(2) = (dBcov(1, 3) - dBcov(3, 1))/sqrtg
        curlB(3) = (dBcov(2, 1) - dBcov(1, 2))/sqrtg
    end subroutine curl_bcov

    function fd_bcov(x, dir, h) result(deriv)
        real(dp), intent(in) :: x(3)
        integer, intent(in) :: dir
        real(dp), intent(in) :: h
        real(dp) :: deriv(3), xp(3), xm(3)
        real(dp) :: Bctr(3), Bcp(3), Bcm(3), Bmod, sqrtg, g(3, 3)

        xp = x; xm = x
        xp(dir) = x(dir) + h
        xm(dir) = x(dir) - h
        call eval_bctr(xp, Bctr, Bcp, Bmod, sqrtg, g)
        call eval_bctr(xm, Bctr, Bcm, Bmod, sqrtg, g)
        deriv = (Bcp - Bcm)/(2.0_dp*h)
    end function fd_bcov

    subroutine test_reload_consistency
        type(spectre_field_t) :: field2
        integer :: jerr
        real(dp) :: x(3), a1(3), h1(3), b1, s1(3), a2(3), h2(3), b2, s2(3)
        real(dp) :: dmax

        call print_test('evaluate threadsafe and consistent after reload')
        x = [1.4_dp, 0.6_dp, 0.25_dp]

        call field%evaluate(x, a1, h1, b1, s1)
        call field%evaluate(x, a2, h2, b2, s2)
        dmax = max(maxval(abs(a1 - a2)), maxval(abs(h1 - h2)), abs(b1 - b2), &
                   maxval(abs(s1 - s2)))

        call create_spectre_field(field2, trim(fixture), jerr)
        if (jerr /= 0) then
            print *, '    reload failed, ierr =', jerr
            nfail = nfail + 1
            call print_fail
            return
        end if
        call field2%evaluate(x, a2, h2, b2, s2)
        dmax = max(dmax, maxval(abs(a1 - a2)), maxval(abs(h1 - h2)), &
                   abs(b1 - b2), maxval(abs(s1 - s2)))

        print *, '    max deviation across calls/reload =', dmax
        if (dmax == 0.0_dp) then
            call print_ok
        else
            nfail = nfail + 1
            call print_fail
        end if
    end subroutine test_reload_consistency

end program test_spectre_field
