module spectre_basis
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use spectre_reader, only: spectre_data_t

    implicit none
    private

    public :: spectre_vecpot_t, get_cheby_d2, get_zernike_d2, &
              eval_spectre_vector_potential, eval_spectre_toroidal_flux_fraction

    ! Covariant vector potential and its s/theta/zeta derivatives at one point.
    ! Second-derivative slots are ordered ss, st, sz, tt, tz, zz.
    type :: spectre_vecpot_t
        real(dp) :: Ath = 0.0_dp
        real(dp) :: Azt = 0.0_dp
        real(dp) :: dAth(3) = 0.0_dp
        real(dp) :: dAzt(3) = 0.0_dp
        real(dp) :: d2Ath(6) = 0.0_dp
        real(dp) :: d2Azt(6) = 0.0_dp
    end type spectre_vecpot_t

contains

    ! Chebyshev radial basis on [-1, 1] with SPEC recombination and 1/(l+1)
    ! scaling; columns hold value, first and second derivative in s.
    pure subroutine get_cheby_d2(s, lrad, cheby)
        real(dp), intent(in) :: s
        integer, intent(in) :: lrad
        real(dp), intent(out) :: cheby(0:lrad, 0:2)

        integer :: ll

        cheby = 0.0_dp
        cheby(0, 0:2) = [1.0_dp, 0.0_dp, 0.0_dp]
        if (lrad >= 1) cheby(1, 0:2) = [s, 1.0_dp, 0.0_dp]
        do ll = 2, lrad
            cheby(ll, 0) = 2.0_dp*s*cheby(ll - 1, 0) - cheby(ll - 2, 0)
            cheby(ll, 1) = 2.0_dp*cheby(ll - 1, 0) + 2.0_dp*s*cheby(ll - 1, 1) &
                           - cheby(ll - 2, 1)
            cheby(ll, 2) = 4.0_dp*cheby(ll - 1, 1) + 2.0_dp*s*cheby(ll - 1, 2) &
                           - cheby(ll - 2, 2)
        end do

        do ll = 1, lrad
            cheby(ll, 0) = cheby(ll, 0) - real((-1)**ll, dp)
        end do
        do ll = 0, lrad
            cheby(ll, :) = cheby(ll, :)/real(ll + 1, dp)
        end do
    end subroutine get_cheby_d2

    ! Radial Zernike basis in sbar with SPEC's m<=1 axis-vanishing corrections
    ! and 1/(n+1) scaling; columns hold value, first and second derivative in
    ! sbar. The corrections make the m=0 basis vanish and the m=1 basis have
    ! zero slope at sbar=0, keeping the field regular on the axis.
    pure subroutine get_zernike_d2(r, lrad, mpol, zernike)
        real(dp), intent(in) :: r
        integer, intent(in) :: lrad, mpol
        real(dp), intent(out) :: zernike(0:lrad, 0:mpol, 0:2)

        real(dp) :: rm, rm1, rm2
        real(dp) :: f1, f2, f3, f4
        integer :: m, n

        rm = 1.0_dp
        rm1 = 0.0_dp
        rm2 = 0.0_dp
        zernike = 0.0_dp
        do m = 0, mpol
            if (lrad >= m) then
                zernike(m, m, 0) = rm
                zernike(m, m, 1) = real(m, dp)*rm1
                zernike(m, m, 2) = real(m*(m - 1), dp)*rm2
            end if

            if (lrad >= m + 2) then
                zernike(m + 2, m, 0) = real(m + 2, dp)*rm*r**2 &
                                       - real(m + 1, dp)*rm
                zernike(m + 2, m, 1) = real((m + 2)**2, dp)*rm*r &
                                       - real((m + 1)*m, dp)*rm1
                zernike(m + 2, m, 2) = real((m + 2)**2*(m + 1), dp)*rm &
                                       - real((m + 1)*m*(m - 1), dp)*rm2
            end if

            do n = m + 4, lrad, 2
                f1 = real(n, dp)/real(n**2 - m**2, dp)
                f2 = real(4*(n - 1), dp)
                f3 = real((n - 2 + m)**2, dp)/real(n - 2, dp) &
                     + real((n - m)**2, dp)/real(n, dp)
                f4 = real((n - 2)**2 - m**2, dp)/real(n - 2, dp)

                zernike(n, m, 0) = f1*((f2*r**2 - f3)*zernike(n - 2, m, 0) &
                                       - f4*zernike(n - 4, m, 0))
                zernike(n, m, 1) = f1*(2.0_dp*f2*r*zernike(n - 2, m, 0) &
                                       + (f2*r**2 - f3)*zernike(n - 2, m, 1) &
                                       - f4*zernike(n - 4, m, 1))
                zernike(n, m, 2) = f1*(2.0_dp*f2*(2.0_dp*r*zernike(n - 2, m, 1) &
                                                  + zernike(n - 2, m, 0)) &
                                       + (f2*r**2 - f3)*zernike(n - 2, m, 2) &
                                       - f4*zernike(n - 4, m, 2))
            end do

            rm2 = rm1
            rm1 = rm
            rm = rm*r
        end do

        do n = 2, lrad, 2
            zernike(n, 0, 0) = zernike(n, 0, 0) - real((-1)**(n/2), dp)
        end do
        if (mpol >= 1) then
            do n = 3, lrad, 2
                zernike(n, 1, 0) = zernike(n, 1, 0) &
                                   - real((-1)**((n - 1)/2), dp) &
                                   *real((n + 1)/2, dp)*r
                zernike(n, 1, 1) = zernike(n, 1, 1) &
                                   - real((-1)**((n - 1)/2), dp) &
                                   *real((n + 1)/2, dp)
            end do
        end if

        do m = 0, mpol
            do n = m, lrad, 2
                zernike(n, m, :) = zernike(n, m, :)/real(n + 1, dp)
            end do
        end do
    end subroutine get_zernike_d2

    ! Evaluate the covariant vector potential (A_theta, A_zeta) and all first
    ! and second (s, theta, zeta) derivatives at (lvol, s, theta, zeta).
    ! Interior volumes use the Chebyshev basis in s; the axis volume (lvol==1)
    ! uses the Zernike basis in sbar=(1+s)/2 with chain-rule factors 1/2 and
    ! 1/4 for the first and second radial derivative. The basis is a smooth
    ! polynomial, so |s| slightly outside [-1, 1] returns the analytic
    ! continuation without clamping.
    pure subroutine eval_spectre_vector_potential(data, lvol, s, theta, zeta, av)
        type(spectre_data_t), intent(in) :: data
        integer, intent(in) :: lvol
        real(dp), intent(in) :: s, theta, zeta
        type(spectre_vecpot_t), intent(out) :: av

        integer :: ii, ll, lrad, mi
        logical :: axis
        real(dp) :: sbar, arg, carg, sarg, rmi, rni
        real(dp) :: cheby(0:data%Lrad(lvol), 0:2)
        real(dp) :: zernike(0:data%Lrad(lvol), 0:data%Mpol, 0:2)
        real(dp) :: tt(0:2)

        lrad = data%Lrad(lvol)
        axis = lvol == 1

        if (axis) then
            sbar = 0.5_dp*(s + 1.0_dp)
            call get_zernike_d2(sbar, lrad, data%Mpol, zernike)
        else
            call get_cheby_d2(s, lrad, cheby)
        end if

        do ii = 1, data%mn
            mi = data%im(ii)
            rmi = real(data%im(ii), dp)
            rni = real(data%in(ii), dp)
            arg = rmi*theta - rni*zeta
            carg = cos(arg)
            sarg = sin(arg)
            do ll = 0, lrad
                if (axis) then
                    tt = [zernike(ll, mi, 0), 0.5_dp*zernike(ll, mi, 1), &
                          0.25_dp*zernike(ll, mi, 2)]
                else
                    tt = cheby(ll, 0:2)
                end if
                call accumulate_mode(data%vol(lvol)%Ate(ll, ii), &
                                     data%vol(lvol)%Ato(ll, ii), carg, sarg, &
                                     rmi, rni, tt, av%Ath, av%dAth, av%d2Ath)
                call accumulate_mode(data%vol(lvol)%Aze(ll, ii), &
                                     data%vol(lvol)%Azo(ll, ii), carg, sarg, &
                                     rmi, rni, tt, av%Azt, av%dAzt, av%d2Azt)
            end do
        end do
    end subroutine eval_spectre_vector_potential

    pure subroutine eval_spectre_toroidal_flux_fraction(data, lvol, s, fraction, &
                                                         derivative, ierr)
        type(spectre_data_t), intent(in) :: data
        integer, intent(in) :: lvol
        real(dp), intent(in) :: s
        real(dp), intent(out) :: fraction, derivative
        integer, intent(out) :: ierr

        real(dp) :: inner, outer, value, slope, scale

        call eval_axisymmetric_atheta(data, lvol, -1.0_dp, inner, slope, ierr)
        if (ierr /= 0) return
        call eval_axisymmetric_atheta(data, lvol, 1.0_dp, outer, slope, ierr)
        if (ierr /= 0) return
        call eval_axisymmetric_atheta(data, lvol, s, value, slope, ierr)
        if (ierr /= 0) return

        scale = max(abs(inner), abs(outer), 1.0_dp)
        if (abs(outer - inner) <= 64.0_dp*epsilon(scale)*scale) then
            ierr = 2
            return
        end if
        fraction = (value - inner)/(outer - inner)
        derivative = slope/(outer - inner)
    end subroutine eval_spectre_toroidal_flux_fraction

    pure subroutine eval_axisymmetric_atheta(data, lvol, s, value, derivative, ierr)
        type(spectre_data_t), intent(in) :: data
        integer, intent(in) :: lvol
        real(dp), intent(in) :: s
        real(dp), intent(out) :: value, derivative
        integer, intent(out) :: ierr

        integer :: ii, ll
        real(dp) :: cheby(0:data%Lrad(lvol), 0:2)
        real(dp) :: zernike(0:data%Lrad(lvol), 0:data%Mpol, 0:2)

        ii = 0
        do ll = 1, data%mn
            if (data%im(ll) == 0 .and. data%in(ll) == 0) then
                ii = ll
                exit
            end if
        end do
        if (ii == 0) then
            ierr = 1
            return
        end if

        if (lvol == 1) then
            call get_zernike_d2(0.5_dp*(s + 1.0_dp), data%Lrad(lvol), &
                                data%Mpol, zernike)
            value = dot_product(data%vol(lvol)%Ate(:, ii), zernike(:, 0, 0))
            derivative = 0.5_dp*dot_product(data%vol(lvol)%Ate(:, ii), &
                                            zernike(:, 0, 1))
        else
            call get_cheby_d2(s, data%Lrad(lvol), cheby)
            value = dot_product(data%vol(lvol)%Ate(:, ii), cheby(:, 0))
            derivative = dot_product(data%vol(lvol)%Ate(:, ii), cheby(:, 1))
        end if
        ierr = 0
    end subroutine eval_axisymmetric_atheta

    ! Add one (mode, radial-degree) contribution to a covariant component and
    ! its derivatives. cos_c/sin_c weight even/odd Fourier coefficients; the
    ! angular derivatives follow from d/dtheta = +m, d/dzeta = -n on arg.
    pure subroutine accumulate_mode(coeff_e, coeff_o, carg, sarg, rmi, rni, &
                                    tt, val, d1, d2)
        real(dp), intent(in) :: coeff_e, coeff_o, carg, sarg, rmi, rni, tt(0:2)
        real(dp), intent(inout) :: val, d1(3), d2(6)

        real(dp) :: p, q

        p = coeff_e*carg + coeff_o*sarg
        q = coeff_o*carg - coeff_e*sarg

        val = val + p*tt(0)
        d1(1) = d1(1) + p*tt(1)
        d1(2) = d1(2) + rmi*q*tt(0)
        d1(3) = d1(3) - rni*q*tt(0)
        d2(1) = d2(1) + p*tt(2)
        d2(2) = d2(2) + rmi*q*tt(1)
        d2(3) = d2(3) - rni*q*tt(1)
        d2(4) = d2(4) - rmi*rmi*p*tt(0)
        d2(5) = d2(5) + rmi*rni*p*tt(0)
        d2(6) = d2(6) - rni*rni*p*tt(0)
    end subroutine accumulate_mode

end module spectre_basis
