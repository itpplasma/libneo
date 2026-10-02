program test_geqdsk_derivative_cache
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use geqdsk_tools, only: geqdsk_t, geqdsk_read, geqdsk_write, &
        geqdsk_classify, geqdsk_standardise, geqdsk_deinit, geqdsk_check_consistency
    use field_eq_mod, only: reset_field_eq_state
    implicit none
    real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)
    real(dp), parameter :: r0 = 620.0_dp, a = 0.1_dp, d = 2.0e5_dp
    real(dp), parameter :: f0 = r0*53000.0_dp, psi_edge = 5.0e7_dp
    character(len=*), parameter :: fixture = 'geqdsk_derivative_cache_fixture.eqdsk'
    type(geqdsk_t) :: input_eq, eq
    real(dp) :: flux_scale, err_cache, err_ampere, max_cache, max_ampere
    integer :: sign_psi, sign_f, full_flux, correction, failures, cases, unit

    failures = 0
    cases = 0
    max_cache = 0.0_dp
    max_ampere = 0.0_dp
    call check_public_no_cache()
    do sign_psi = -1, 1, 2
        do sign_f = -1, 1, 2
            do full_flux = 0, 1
                flux_scale = (2.0_dp*pi)**full_flux
                do correction = 0, 3
                    call make_fixture(input_eq, sign_psi, sign_f, flux_scale, correction)
                    call geqdsk_write(input_eq, fixture)
                    call geqdsk_read(eq, fixture)
                    call geqdsk_classify(eq)
                    call geqdsk_standardise(eq)
                    call check_analytic_cache(eq, sign_psi, sign_f, err_cache)
                    call initialize_field(eq)
                    call check_native_ampere(eq, sign_psi, err_ampere)
                    if (.not. ieee_is_finite(err_cache)) error stop 'nonfinite cache error'
                    if (.not. ieee_is_finite(err_ampere)) error stop 'nonfinite curl error'
                    max_cache = max(max_cache, err_cache)
                    max_ampere = max(max_ampere, err_ampere)
                    if (err_cache > 2.0e-8_dp .or. err_ampere > 2.0e-6_dp) then
                        failures = failures + 1
                        write(*, '(a,4i4,2es15.6)') 'FAILED signs, flux, correction, ' // &
                            'analytic-cache and native-Ampere errors:', &
                            sign_psi, sign_f, full_flux, correction, err_cache, err_ampere
                    end if
                    cases = cases + 1
                    call reset_field_eq_state()
                    call geqdsk_deinit(eq)
                    call geqdsk_deinit(input_eq)
                end do
            end do
        end do
    end do
    open(newunit=unit, file=fixture, status='old')
    close(unit, status='delete')
    write(*, '(a,2i5,2es16.7)') 'CACHE_GS_ORACLE cases, failures, maxima:', &
        cases, failures, max_cache, max_ampere
    if (failures /= 0) error stop 'corrected profiles violate analytic and Ampere oracles'

contains

    subroutine check_public_no_cache()
        type(geqdsk_t) :: profile_eq
        integer :: i

        profile_eq%nw = 17
        profile_eq%cocos%sgn_Btor = 1
        profile_eq%cocos%sgn_F = 1
        allocate(profile_eq%psi_eqd(17), profile_eq%pres(17), profile_eq%pprime(17))
        allocate(profile_eq%fpol(17), profile_eq%ffprim(17))
        profile_eq%psi_eqd = [(real(i - 1, dp), i = 1, 17)]
        profile_eq%pres = 1.0_dp + 0.2_dp*profile_eq%psi_eqd
        profile_eq%pprime = 0.2_dp
        profile_eq%fpol = 2.0_dp + 0.1_dp*profile_eq%psi_eqd
        profile_eq%ffprim = 0.1_dp*profile_eq%fpol
        if (allocated(profile_eq%fprime)) error stop 'manual fixture has a cache'
        call geqdsk_check_consistency(profile_eq)
        if (allocated(profile_eq%fprime)) error stop 'public check allocated a cache'
        if (.not. all(ieee_is_finite(profile_eq%ffprim))) &
            error stop 'public check returned a nonfinite profile'
        if (maxval(abs(profile_eq%ffprim - 0.1_dp*profile_eq%fpol)) > 1.0e-12_dp) &
            error stop 'public check changed coherent F derivative'
        if (maxval(abs(profile_eq%pprime - 0.2_dp)) > 1.0e-12_dp) &
            error stop 'public check changed coherent pressure derivative'
        write(*, '(a)') 'PUBLIC_NO_CACHE_ORACLE passed without cache allocation'
    end subroutine check_public_no_cache

    pure real(dp) function flux(r, z)
        real(dp), intent(in) :: r, z
        flux = a/8.0_dp*(r**2 - r0**2)**2 + d/2.0_dp*z**2
    end function flux

    subroutine make_fixture(eq, sign_psi, sign_f, flux_scale, correction)
        type(geqdsk_t), intent(out) :: eq
        integer, intent(in) :: sign_psi, sign_f, correction
        real(dp), intent(in) :: flux_scale
        real(dp) :: x, r, z, theta, delta_r2, height, loop
        real(dp) :: br, bz, dr_dtheta, dz_dtheta
        integer :: i, j

        eq%nw = 65
        eq%nh = 65
        eq%nbbbs = 128
        eq%limitr = 128
        eq%header = 'Exact GS derivative-cache regression'
        eq%rdim = 160.0_dp
        eq%zdim = 160.0_dp
        eq%rcentr = r0
        eq%rleft = r0 - 80.0_dp
        eq%zmid = 0.0_dp
        eq%rmaxis = r0
        eq%zmaxis = 0.0_dp
        eq%simag = 0.0_dp
        eq%sibry = sign_psi*flux_scale*psi_edge
        eq%bcentr = sign_f*f0/r0
        allocate(eq%fpol(eq%nw), eq%pres(eq%nw), eq%ffprim(eq%nw))
        allocate(eq%pprime(eq%nw), eq%qpsi(eq%nw), eq%psirz(eq%nw, eq%nh))
        allocate(eq%rbbbs(eq%nbbbs), eq%zbbbs(eq%nbbbs))
        allocate(eq%rlim(eq%limitr), eq%zlim(eq%limitr))
        do i = 1, eq%nw
            x = psi_edge*real(i - 1, dp)/real(eq%nw - 1, dp)
            eq%fpol(i) = sign_f*sqrt(f0**2 - 2.0_dp*d*x)
            eq%pres(i) = 1.0e6_dp - a*x/(4.0_dp*pi)
            eq%qpsi(i) = sign_psi*eq%fpol(i)/(sqrt(a*d)* &
                sqrt(r0**4 - 8.0_dp*x/a))
        end do
        eq%ffprim = -sign_psi*d/flux_scale
        eq%pprime = -sign_psi*a/(4.0_dp*pi*flux_scale)
        select case (correction)
        case (1)
            eq%fpol = -eq%fpol
        case (2)
            eq%ffprim = -eq%ffprim
        case (3)
            eq%ffprim = eq%ffprim/(2.0_dp*pi)
        end select
        do j = 1, eq%nh
            z = -80.0_dp + 160.0_dp*real(j - 1, dp)/real(eq%nh - 1, dp)
            do i = 1, eq%nw
                r = eq%rleft + eq%rdim*real(i - 1, dp)/real(eq%nw - 1, dp)
                eq%psirz(i, j) = sign_psi*flux_scale*flux(r, z)
            end do
        end do
        delta_r2 = sqrt(8.0_dp*psi_edge/a)
        height = sqrt(2.0_dp*psi_edge/d)
        loop = 0.0_dp
        do i = 1, eq%nbbbs
            theta = 2.0_dp*pi*real(i - 1, dp)/real(eq%nbbbs, dp)
            r = sqrt(r0**2 + delta_r2*cos(theta))
            z = height*sin(theta)
            eq%rbbbs(i) = r
            eq%zbbbs(i) = z
            eq%rlim(i) = r
            eq%zlim(i) = z
            br = -sign_psi*d*z/r
            bz = sign_psi*a/2.0_dp*(r**2 - r0**2)
            dr_dtheta = -delta_r2*sin(theta)/(2.0_dp*r)
            dz_dtheta = height*cos(theta)
            loop = loop + br*dr_dtheta + bz*dz_dtheta
        end do
        eq%current = -2.99792458e10_dp*loop/(2.0_dp*real(eq%nbbbs, dp))
    end subroutine make_fixture

    subroutine check_analytic_cache(eq, sign_psi, sign_f, err)
        type(geqdsk_t), intent(in) :: eq
        integer, intent(in) :: sign_psi, sign_f
        real(dp), intent(out) :: err
        real(dp) :: x, expected
        integer :: i

        err = 0.0_dp
        do i = 1, eq%nw
            x = psi_edge*real(i - 1, dp)/real(eq%nw - 1, dp)
            expected = -sign_psi*d/(sign_f*sqrt(f0**2 - 2.0_dp*d*x))
            err = max(err, abs((eq%fprime(i) - expected)/expected))
        end do
    end subroutine check_analytic_cache

    subroutine initialize_field(eq)
        use field_eq_mod, only: use_fpol, skip_read, icall_eq, nrad, nzet, &
            nwindow_r, nwindow_z, psi_axis, psi_sep, btf, rtf, splfpol, &
            rad, zet, psi, psi0
        use field_sub, only: field_eq
        type(geqdsk_t), intent(in) :: eq
        real(dp) :: b(3), db_dr(3), db_dz(3), db_dphi(3)

        use_fpol = .true.
        skip_read = .true.
        icall_eq = -1
        nwindow_r = 0
        nwindow_z = 0
        nrad = eq%nw
        nzet = eq%nh
        allocate(rad(nrad), zet(nzet), psi0(nrad, nzet), psi(nrad, nzet))
        allocate(splfpol(0:5, nrad))
        psi_axis = eq%simag*1.0e-8_dp
        psi_sep = eq%sibry*1.0e-8_dp
        btf = eq%bcentr*1.0e-4_dp
        rtf = eq%rcentr*1.0e-2_dp
        splfpol(0, :) = eq%fpol*1.0e-6_dp
        psi = eq%psirz*1.0e-8_dp
        rad = eq%R_eqd*1.0e-2_dp
        zet = eq%Z_eqd*1.0e-2_dp
        call field_eq(r0, 0.0_dp, 0.0_dp, b(1), b(2), b(3), &
            db_dr(1), db_dphi(1), db_dz(1), db_dr(2), db_dphi(2), &
            db_dz(2), db_dr(3), db_dphi(3), db_dz(3))
    end subroutine initialize_field

    subroutine check_native_ampere(eq, sign_psi, err)
        use field_sub, only: field_eq
        type(geqdsk_t), intent(in) :: eq
        integer, intent(in) :: sign_psi
        real(dp), intent(out) :: err
        real(dp) :: x, r, z, b(3), db_dr(3), db_dz(3), db_dphi(3)
        real(dp) :: ampere_pol(2), profile_pol(2), norm, expected_phi
        integer :: i, direction

        err = 0.0_dp
        do i = 9, eq%nw - 8, 8
            x = psi_edge*real(i - 1, dp)/real(eq%nw - 1, dp)
            do direction = 1, 2
                if (direction == 1) then
                    r = sqrt(r0**2 + sqrt(8.0_dp*x/a))
                    z = 0.0_dp
                else
                    r = r0
                    z = sqrt(2.0_dp*x/d)
                end if
                call field_eq(r, 0.0_dp, z, b(1), b(2), b(3), &
                    db_dr(1), db_dphi(1), db_dz(1), db_dr(2), db_dphi(2), &
                    db_dz(2), db_dr(3), db_dphi(3), db_dz(3))
                ampere_pol(1) = -db_dz(2)
                ampere_pol(2) = db_dr(2) + b(2)/r
                profile_pol(1) = eq%fprime(i)*b(1)
                profile_pol(2) = eq%fprime(i)*b(3)
                norm = sqrt(sum(ampere_pol**2))
                err = max(err, sqrt(sum((profile_pol - ampere_pol)**2))/norm)
                expected_phi = -sign_psi*(a*r + d/r)
                err = max(err, abs((db_dz(1) - db_dr(3) - &
                    expected_phi)/expected_phi))
            end do
        end do
    end subroutine check_native_ampere
end program test_geqdsk_derivative_cache
