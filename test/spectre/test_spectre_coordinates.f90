program test_spectre_coordinates
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use spectre_reader, only: spectre_data_t, load_spectre, free_spectre
    use libneo_coordinates, only: coordinate_system_t, &
                                  spectre_coordinate_system_t, &
                                  make_spectre_coordinate_system
    use util_for_test, only: print_test, print_ok, print_fail

    implicit none

    real(dp), parameter :: twopi = 8.0_dp*atan(1.0_dp)

    character(len=1024) :: fixture, reffile
    type(spectre_data_t) :: data
    class(coordinate_system_t), allocatable :: cs
    integer :: ierr, Mvol, nfail

    nfail = 0
    call get_command_argument(1, fixture)
    call get_command_argument(2, reffile)
    if (len_trim(fixture) == 0 .or. len_trim(reffile) == 0) then
        print *, 'usage: test_spectre_coordinates.x <spectre.h5> <reference_eval.dat>'
        error stop
    end if

    call load_spectre(trim(fixture), data, ierr)
    if (ierr /= 0) then
        print *, 'load_spectre failed, ierr =', ierr
        error stop
    end if
    Mvol = data%Mvol
    call make_spectre_coordinate_system(cs, data)
    call free_spectre(data)

    call test_interface_continuity
    call test_roundtrip
    call test_warm_roundtrip
    call test_warm_interface_ownership
    call test_jacobian_vs_reference
    call test_metric_vs_fd
    call test_axis_single_signed

    if (nfail > 0) error stop

contains

    subroutine test_interface_continuity
        integer :: k
        real(dp) :: u_lo(3), u_hi(3), x_lo(3), x_hi(3), dmax

        call print_test('geometry is C0 across interfaces (rho_g = 1, 2)')
        dmax = 0.0_dp
        do k = 1, Mvol - 1
            u_lo = [real(k, dp) - 1.0e-13_dp, 0.7_dp, 0.4_dp]
            u_hi = [real(k, dp) + 1.0e-13_dp, 0.7_dp, 0.4_dp]
            call cs%evaluate_cyl(u_lo, x_lo)
            call cs%evaluate_cyl(u_hi, x_hi)
            dmax = max(dmax, abs(x_lo(1) - x_hi(1)), abs(x_lo(3) - x_hi(3)))
        end do
        if (dmax < 1.0e-11_dp) then
            call print_ok
        else
            print *, '    max interface jump =', dmax
            nfail = nfail + 1
            call print_fail
        end if
    end subroutine test_interface_continuity

    subroutine test_roundtrip
        integer, parameter :: npt = 200
        integer :: i, seed_size
        integer, allocatable :: seed(:)
        real(dp) :: r(3), u(3), xcyl(3), u2(3), dtheta
        real(dp) :: err_rho, err_theta, worst_rho, worst_theta
        integer :: jerr, lvol_ref, lvol_got, bad_vol

        call print_test('from_cyl(evaluate_cyl(u)) = u for 200 interior points')
        call random_seed(size=seed_size)
        allocate (seed(seed_size))
        seed = 20260711
        call random_seed(put=seed)

        worst_rho = 0.0_dp
        worst_theta = 0.0_dp
        bad_vol = 0
        do i = 1, npt
            call random_number(r)
            u(1) = 0.05_dp + r(1)*(real(Mvol, dp) - 0.05_dp - 1.0e-3_dp)
            u(2) = twopi*r(2)
            u(3) = twopi*r(3)

            call cs%evaluate_cyl(u, xcyl)
            call cs%from_cyl(xcyl, u2, jerr)
            if (jerr /= 0) then
                bad_vol = bad_vol + 1
                cycle
            end if

            err_rho = abs(u2(1) - u(1))
            dtheta = modulo(u2(2) - u(2) + twopi/2.0_dp, twopi) - twopi/2.0_dp
            err_theta = abs(dtheta)
            worst_rho = max(worst_rho, err_rho)
            worst_theta = max(worst_theta, err_theta)

            lvol_ref = min(int(u(1)) + 1, Mvol)
            lvol_got = min(int(u2(1)) + 1, Mvol)
            if (lvol_ref /= lvol_got) bad_vol = bad_vol + 1
        end do

        if (bad_vol == 0 .and. worst_rho < 1.0e-8_dp &
            .and. worst_theta < 1.0e-8_dp) then
            call print_ok
        else
            print *, '    bad_vol =', bad_vol, ' worst_rho =', worst_rho, &
                ' worst_theta =', worst_theta
            nfail = nfail + 1
            call print_fail
        end if
    end subroutine test_roundtrip

    subroutine test_warm_roundtrip
        integer :: lvol, jerr
        real(dp) :: u_ref(3), u_warm(3), xcyl(3), dtheta
        real(dp) :: worst_rho, worst_theta

        call print_test('warm inverse converges from perturbed seeds in owned volume')
        worst_rho = 0.0_dp
        worst_theta = 0.0_dp
        do lvol = 1, Mvol
            u_ref = [real(lvol, dp) - 0.45_dp, 0.7_dp, 0.4_dp]
            call cs%evaluate_cyl(u_ref, xcyl)
            u_warm = [u_ref(1) + 0.1_dp, u_ref(2) + 0.15_dp, -1.0_dp]
            select type (spectre => cs)
            type is (spectre_coordinate_system_t)
                call spectre%from_cyl_warm(xcyl, u_warm, lvol, jerr)
            class default
                jerr = 1
            end select
            if (jerr /= 0) exit
            worst_rho = max(worst_rho, abs(u_warm(1) - u_ref(1)))
            dtheta = modulo(u_warm(2) - u_ref(2) + 0.5_dp*twopi, twopi) &
                - 0.5_dp*twopi
            worst_theta = max(worst_theta, abs(dtheta))
            if (u_warm(3) /= xcyl(2)) jerr = 1
        end do

        if (jerr == 0 .and. worst_rho < 1.0e-8_dp &
            .and. worst_theta < 1.0e-8_dp) then
            call print_ok
        else
            print *, '    ierr =', jerr, ' worst_rho =', worst_rho, &
                ' worst_theta =', worst_theta
            nfail = nfail + 1
            call print_fail
        end if
    end subroutine test_warm_roundtrip

    subroutine test_warm_interface_ownership
        integer :: k, jerr_lo, jerr_hi
        real(dp) :: u_ref(3), u_lo(3), u_hi(3), xcyl(3), worst

        call print_test('warm inverse accepts exact shared interface from each owner')
        worst = 0.0_dp
        jerr_lo = 0
        jerr_hi = 0
        do k = 1, Mvol - 1
            u_ref = [real(k, dp), 0.7_dp, 0.4_dp]
            call cs%evaluate_cyl(u_ref, xcyl)
            u_lo = [real(k, dp) - 0.05_dp, 0.8_dp, 0.0_dp]
            u_hi = [real(k, dp) + 0.05_dp, 0.6_dp, 0.0_dp]
            select type (spectre => cs)
            type is (spectre_coordinate_system_t)
                call spectre%from_cyl_warm(xcyl, u_lo, k, jerr_lo)
                call spectre%from_cyl_warm(xcyl, u_hi, k + 1, jerr_hi)
            class default
                jerr_lo = 1
                jerr_hi = 1
            end select
            worst = max(worst, abs(u_lo(1) - real(k, dp)))
            worst = max(worst, abs(u_hi(1) - real(k, dp)))
            if (jerr_lo /= 0 .or. jerr_hi /= 0) exit
        end do

        if (jerr_lo == 0 .and. jerr_hi == 0 .and. worst < 1.0e-8_dp) then
            call print_ok
        else
            print *, '    ierr_lo =', jerr_lo, ' ierr_hi =', jerr_hi, &
                ' worst_rho =', worst
            nfail = nfail + 1
            call print_fail
        end if
    end subroutine test_warm_interface_ownership

    subroutine test_jacobian_vs_reference
        integer :: iunit, io, lvol, nrows, sign_bad
        real(dp) :: col(11), s, theta, zeta, jac, rho_g, u(3)
        real(dp) :: g(3, 3), ginv(3, 3), sqrtg, e_cov(3, 3)
        real(dp) :: relerr, worst
        character(len=256) :: line

        call print_test('|sqrtg| = 2*|reference jac| at all reference rows')
        open (newunit=iunit, file=trim(reffile), status='old', action='read')
        read (iunit, '(A)') line
        nrows = 0
        worst = 0.0_dp
        sign_bad = 0
        do
            read (iunit, *, iostat=io) col
            if (io /= 0) exit
            nrows = nrows + 1
            lvol = nint(col(1))
            s = col(2)
            theta = col(3)
            zeta = col(4)
            jac = col(10)
            rho_g = real(lvol - 1, dp) + 0.5_dp*(s + 1.0_dp)
            u = [rho_g, theta, zeta]

            call cs%metric_tensor(u, g, ginv, sqrtg)
            relerr = abs(sqrtg - 2.0_dp*abs(jac))/(2.0_dp*abs(jac))
            worst = max(worst, relerr)

            call cs%covariant_basis(u, e_cov)
            if (det3(e_cov)*jac <= 0.0_dp) sign_bad = sign_bad + 1
        end do
        close (iunit)

        if (worst < 1.0e-10_dp .and. sign_bad == 0 .and. nrows == 90) then
            call print_ok
        else
            print *, '    nrows =', nrows, ' worst relerr =', worst, &
                ' sign_bad =', sign_bad
            nfail = nfail + 1
            call print_fail
        end if
    end subroutine test_jacobian_vs_reference

    subroutine test_metric_vs_fd
        integer :: iv, ip, i, j
        real(dp) :: u(3), g(3, 3), ginv(3, 3), sqrtg, gfd(3, 3)
        real(dp) :: worst
        real(dp), parameter :: theta_p(3) = [0.4_dp, 2.7_dp, 5.1_dp]

        call print_test('analytic metric matches finite diff of evaluate_cart')
        worst = 0.0_dp
        do iv = 1, Mvol
            do ip = 1, 3
                u(1) = real(iv, dp) - 0.5_dp
                u(2) = theta_p(ip)
                u(3) = 0.3_dp + 0.2_dp*real(ip, dp)
                call cs%metric_tensor(u, g, ginv, sqrtg)
                call metric_fd(u, gfd)
                do i = 1, 3
                    do j = 1, 3
                        worst = max(worst, abs(g(i, j) - gfd(i, j)) &
                                    /max(1.0_dp, abs(g(i, j))))
                    end do
                end do
            end do
        end do

        if (worst < 1.0e-6_dp) then
            call print_ok
        else
            print *, '    worst metric-FD relerr =', worst
            nfail = nfail + 1
            call print_fail
        end if
    end subroutine test_metric_vs_fd

    subroutine metric_fd(u, gfd)
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: gfd(3, 3)

        real(dp), parameter :: h = 1.0e-6_dp
        real(dp) :: ecol(3, 3), up(3), um(3), xp(3), xm(3)
        integer :: k, i, j

        do k = 1, 3
            up = u
            um = u
            up(k) = u(k) + h
            um(k) = u(k) - h
            call cs%evaluate_cart(up, xp)
            call cs%evaluate_cart(um, xm)
            ecol(:, k) = (xp - xm)/(2.0_dp*h)
        end do
        do i = 1, 3
            do j = 1, 3
                gfd(i, j) = dot_product(ecol(:, i), ecol(:, j))
            end do
        end do
    end subroutine metric_fd

    subroutine test_axis_single_signed
        integer :: ir, it
        real(dp) :: rho_vals(3), u(3), g(3, 3), ginv(3, 3), sqrtg, e_cov(3, 3)
        real(dp) :: det, ref_sign
        logical :: ok

        call print_test('axis Jacobian finite and single-signed as rho_g -> 0')
        rho_vals = [1.0e-6_dp, 1.0e-3_dp, 0.1_dp]
        ok = .true.
        ref_sign = 0.0_dp
        do ir = 1, 3
            do it = 0, 15
                u(1) = rho_vals(ir)
                u(2) = twopi*real(it, dp)/16.0_dp
                u(3) = 0.3_dp
                call cs%metric_tensor(u, g, ginv, sqrtg)
                call cs%covariant_basis(u, e_cov)
                det = det3(e_cov)
                if (.not. (sqrtg == sqrtg)) ok = .false.
                if (sqrtg <= 0.0_dp) ok = .false.
                if (ref_sign == 0.0_dp) ref_sign = sign(1.0_dp, det)
                if (det*ref_sign <= 0.0_dp) ok = .false.
            end do
        end do

        if (ok) then
            call print_ok
        else
            nfail = nfail + 1
            call print_fail
        end if
    end subroutine test_axis_single_signed

    pure function det3(a) result(d)
        real(dp), intent(in) :: a(3, 3)
        real(dp) :: d

        d = a(1, 1)*(a(2, 2)*a(3, 3) - a(2, 3)*a(3, 2)) &
            - a(1, 2)*(a(2, 1)*a(3, 3) - a(2, 3)*a(3, 1)) &
            + a(1, 3)*(a(2, 1)*a(3, 2) - a(2, 2)*a(3, 1))
    end function det3

end program test_spectre_coordinates
