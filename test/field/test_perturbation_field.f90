! Module procedures, not internal ones: passing an internal procedure as an
! actual argument needs an executable-stack trampoline, which macOS arm64 forbids.
module perturbation_field_test_amplitudes
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    complex(dp), parameter :: c1 = (1.0_dp, 0.5_dp), c2 = (-0.7_dp, 0.2_dp)
    complex(dp), parameter :: c3 = (0.3_dp, -1.1_dp), cg = (0.4_dp, 0.9_dp)
    complex(dp), parameter :: iu = (0.0_dp, 1.0_dp)

contains

    ! Smooth complex amplitudes, physical cylindrical components. The hand-derived
    ! curl below is the independent oracle:
    !   dB_R   = (i n/R) a_Z - d_Z a_phi
    !   dB_phi = d_Z a_R - d_R a_Z
    !   dB_Z   = a_phi/R + d_R a_phi - (i n/R) a_R
    subroutine amp_smooth(n, R, Z, dA, dPhi)
        integer, intent(in) :: n
        real(dp), intent(in) :: R, Z
        complex(dp), intent(out) :: dA(3), dPhi

        dA(1) = c1*R*sin(Z)
        dA(2) = c2*Z*cos(R)
        dA(3) = c3*R**2*cos(Z)
        dPhi = real(n + 1, dp)*c2*R*Z
    end subroutine amp_smooth

    subroutine curl_smooth(n, R, Z, dB)
        integer, intent(in) :: n
        real(dp), intent(in) :: R, Z
        complex(dp), intent(out) :: dB(3)

        dB(1) = iu*n*c3*R*cos(Z) - c2*cos(R)
        dB(2) = c1*R*cos(Z) - 2.0_dp*c3*R*cos(Z)
        dB(3) = c2*Z*cos(R)/R - c2*Z*sin(R) - iu*n*c1*sin(Z)
    end subroutine curl_smooth

    ! amp_smooth plus grad chi with chi = cg sin(R) cos(Z) exp(i n phi):
    ! grad chi = (d_R chi, i n chi / R, d_Z chi).
    subroutine amp_gauged(n, R, Z, dA, dPhi)
        integer, intent(in) :: n
        real(dp), intent(in) :: R, Z
        complex(dp), intent(out) :: dA(3), dPhi

        call amp_smooth(n, R, Z, dA, dPhi)
        dA(1) = dA(1) + cg*cos(R)*cos(Z)
        dA(2) = dA(2) + iu*n*cg*sin(R)*cos(Z)/R
        dA(3) = dA(3) - cg*sin(R)*sin(Z)
    end subroutine amp_gauged

end module perturbation_field_test_amplitudes

program test_perturbation_field
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use perturbation_field_test_amplitudes, only: amp_smooth, curl_smooth, amp_gauged
    use neo_perturbation_field, only: perturbation_field_t, UNITS_SI, UNITS_GAUSSIAN
    use neo_perturbation_field_netcdf, only: write_perturbation_field_netcdf, &
        read_perturbation_field_netcdf
    implicit none

    integer, parameter :: dp = kind(1.0d0)
    complex(dp), parameter :: iu = (0.0_dp, 1.0_dp)
    real(dp), parameter :: R_lo = 1.0_dp, R_hi = 2.0_dp, Z_lo = -0.5_dp, Z_hi = 0.5_dp
    integer, parameter :: n_pts = 20
    integer :: n_fail = 0

    call fix_random_seed
    call test_curl_matches_analytic
    call test_real_field_sums_modes
    call test_divergence_free
    call test_gauge_invariance
    call test_fail_closed
    call test_netcdf_round_trip
    call test_netcdf_rejects_bad_header

    if (n_fail > 0) then
        print *, 'test_perturbation_field: FAILED checks:', n_fail
        error stop 1
    end if
    print *, 'test_perturbation_field: all checks passed'

contains

    subroutine fix_random_seed
        integer, allocatable :: seed(:)
        integer :: n, i

        call random_seed(size=n)
        allocate (seed(n))
        seed = [(12345 + 7*i, i=1, n)]
        call random_seed(put=seed)
    end subroutine fix_random_seed

    subroutine check(ok, label)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: label

        if (ok) then
            print '(a,a)', 'PASS ', label
        else
            print '(a,a)', 'FAIL ', label
            n_fail = n_fail + 1
        end if
    end subroutine check

    subroutine make_field(field, ntor, amp, nR, nZ, with_potential)
        type(perturbation_field_t), intent(out) :: field
        integer, intent(in) :: ntor(:)
        interface
            subroutine amp(n, R, Z, dA, dPhi)
                import :: dp
                integer, intent(in) :: n
                real(dp), intent(in) :: R, Z
                complex(dp), intent(out) :: dA(3), dPhi
            end subroutine amp
        end interface
        integer, intent(in) :: nR, nZ
        logical, intent(in) :: with_potential

        integer :: ierr

        call field%init_analytic(R_lo, R_hi, nR, Z_lo, Z_hi, nZ, ntor, amp, &
            UNITS_SI, ierr, with_potential=with_potential)
        if (ierr /= 0) error stop 'init_analytic failed'
    end subroutine make_field

    subroutine random_point(x)
        real(dp), intent(out) :: x(3)

        call random_number(x)
        x(1) = R_lo + 0.05_dp + 0.9_dp*(R_hi - R_lo)*x(1)
        x(2) = 8.0_dp*atan(1.0_dp)*x(2)
        x(3) = Z_lo + 0.05_dp + 0.9_dp*(Z_hi - Z_lo)*x(3)
    end subroutine random_point

    subroutine test_curl_matches_analytic
        type(perturbation_field_t) :: field
        complex(dp) :: dA(3, 2), dB(3, 2), dA_ref(3), dB_ref(3), dPhi(2), dPhi_ref
        real(dp) :: x(3), err_a, err_b, err_p
        integer :: i, k, ierr

        call make_field(field, [0, 3], amp_smooth, 81, 81, .true.)
        err_a = 0; err_b = 0; err_p = 0
        do i = 1, n_pts
            call random_point(x)
            call field%eval_modes(x(1), x(3), dA, dB, ierr, dPhi)
            if (ierr /= 0) error stop 'eval_modes failed inside grid'
            do k = 1, 2
                call amp_smooth(field%ntor(k), x(1), x(3), dA_ref, dPhi_ref)
                call curl_smooth(field%ntor(k), x(1), x(3), dB_ref)
                err_a = max(err_a, maxval(abs(dA(:, k) - dA_ref)))
                err_b = max(err_b, maxval(abs(dB(:, k) - dB_ref)))
                err_p = max(err_p, abs(dPhi(k) - dPhi_ref))
            end do
        end do
        print '(a,3es10.2)', '  max err dA, dB, dPhi:', err_a, err_b, err_p
        call check(err_a < 1.0e-9_dp, 'spline dA matches analytic dA')
        call check(err_b < 1.0e-7_dp, 'curl dA matches hand-derived curl')
        call check(err_p < 1.0e-9_dp, 'spline dPhi matches analytic dPhi')
    end subroutine test_curl_matches_analytic

    subroutine test_real_field_sums_modes
        type(perturbation_field_t) :: field
        complex(dp) :: dB_ref(3)
        real(dp) :: x(3), dA(3), dB(3), B0(3), dbm, expect(3), err, err_bmod
        integer :: i, k, ierr
        integer, parameter :: ntor(2) = [0, 3]

        call make_field(field, ntor, amp_smooth, 81, 81, .false.)
        err = 0; err_bmod = 0
        B0 = [0.1_dp, -2.0_dp, 0.3_dp]
        do i = 1, n_pts
            call random_point(x)
            expect = 0
            do k = 1, 2
                call curl_smooth(ntor(k), x(1), x(3), dB_ref)
                expect = expect + real(dB_ref*exp(iu*ntor(k)*x(2)), dp)
            end do
            call field%eval(x, dA, dB, ierr)
            if (ierr /= 0) error stop 'eval failed inside grid'
            err = max(err, maxval(abs(dB - expect)))
            call field%dbmod(x, B0, dbm, ierr)
            err_bmod = max(err_bmod, abs(dbm - dot_product(B0, expect)/norm2(B0)))
        end do
        call check(err < 1.0e-7_dp, 'real dB = Re(sum_n dB_n exp(i n phi))')
        call check(err_bmod < 1.0e-7_dp, 'Eulerian d|B| = b0 . dB')
    end subroutine test_real_field_sums_modes

    ! Divergence of the real field from 4th-order central differences in R, phi
    ! and Z: (1/R) d_R(R B_R) + (1/R) d_phi B_phi + d_Z B_Z. The phi derivative
    ! checks the i n / R terms independently of the implementation.
    real(dp) function fd_div(field, x, use_a) result(div)
        type(perturbation_field_t), intent(in) :: field
        real(dp), intent(in) :: x(3)
        logical, intent(in) :: use_a

        real(dp), parameter :: h = 1.0e-3_dp, w(4) = [1, -8, 8, -1]/12.0_dp
        real(dp), parameter :: s(4) = [-2, -1, 1, 2]
        real(dp) :: xs(3), v(3)
        integer :: j, d

        div = 0
        do d = 1, 3
            do j = 1, 4
                xs = x
                xs(d) = x(d) + s(j)*h
                v = field_vector(field, xs, use_a)
                if (d == 1) div = div + w(j)*xs(1)*v(1)/(h*x(1))
                if (d == 2) div = div + w(j)*v(2)/(h*x(1))
                if (d == 3) div = div + w(j)*v(3)/h
            end do
        end do
    end function fd_div

    function field_vector(field, x, use_a) result(v)
        type(perturbation_field_t), intent(in) :: field
        real(dp), intent(in) :: x(3)
        logical, intent(in) :: use_a
        real(dp) :: v(3)

        real(dp) :: dA(3), dB(3)
        integer :: ierr

        call field%eval(x, dA, dB, ierr)
        if (ierr /= 0) error stop 'eval failed in fd_div'
        v = dB
        if (use_a) v = dA
    end function field_vector

    subroutine test_divergence_free
        type(perturbation_field_t) :: field
        real(dp) :: x(3), dA(3), dB(3), div_b, div_a, scale
        integer :: i, ierr

        call make_field(field, [0, 1, 3, 7], amp_smooth, 41, 37, .false.)
        div_b = 0; div_a = 0; scale = 0
        do i = 1, n_pts
            call random_point(x)
            call field%eval(x, dA, dB, ierr)
            scale = max(scale, maxval(abs(dB)))
            div_b = max(div_b, abs(fd_div(field, x, .false.)))
            div_a = max(div_a, abs(fd_div(field, x, .true.)))
        end do
        print '(a,3es10.2)', '  max |div dB|, |div dA| (control), |dB|:', &
            div_b, div_a, scale
        call check(div_b < 1.0e-9_dp*scale, 'div dB = 0 at random points')
        call check(div_a > 1.0e-2_dp*scale, 'control: div dA is detected nonzero')
    end subroutine test_divergence_free

    subroutine test_gauge_invariance
        type(perturbation_field_t) :: f1, f2
        complex(dp) :: dA1(3, 3), dB1(3, 3), dA2(3, 3), dB2(3, 3)
        real(dp) :: x(3), err_b, diff_a
        integer :: i, ierr1, ierr2

        call make_field(f1, [0, 2, 5], amp_smooth, 101, 101, .false.)
        call make_field(f2, [0, 2, 5], amp_gauged, 101, 101, .false.)
        err_b = 0; diff_a = 0
        do i = 1, n_pts
            call random_point(x)
            call f1%eval_modes(x(1), x(3), dA1, dB1, ierr1)
            call f2%eval_modes(x(1), x(3), dA2, dB2, ierr2)
            if (ierr1 /= 0 .or. ierr2 /= 0) error stop 'eval_modes failed'
            err_b = max(err_b, maxval(abs(dB1 - dB2)))
            diff_a = max(diff_a, maxval(abs(dA1 - dA2)))
        end do
        print '(a,2es10.2)', '  gauge: max |dB1-dB2|, |dA1-dA2|:', err_b, diff_a
        call check(diff_a > 0.1_dp, 'gauge: potentials differ')
        call check(err_b < 1.0e-7_dp, 'gauge: dA + grad chi gives the same dB')
    end subroutine test_gauge_invariance

    logical function all_nan_c(a) result(ok)
        complex(dp), intent(in) :: a(:, :)

        ok = all(ieee_is_nan(a%re)) .and. all(ieee_is_nan(a%im))
    end function all_nan_c

    subroutine test_fail_closed
        type(perturbation_field_t) :: field, bad
        complex(dp) :: dA(3, 1), dB(3, 1)
        real(dp) :: rA(3), rB(3), dbm, outside(2, 5)
        integer :: i, ierr

        call make_field(field, [2], amp_smooth, 21, 21, .false.)
        outside(:, 1) = [R_lo - 1.0e-9_dp, 0.0_dp]
        outside(:, 2) = [R_hi + 1.0e-9_dp, 0.0_dp]
        outside(:, 3) = [1.5_dp, Z_lo - 1.0e-9_dp]
        outside(:, 4) = [1.5_dp, Z_hi + 1.0e-9_dp]
        outside(:, 5) = [ieee_nan(), 0.0_dp]
        do i = 1, 5
            call field%eval_modes(outside(1, i), outside(2, i), dA, dB, ierr)
            call check(ierr /= 0, 'fail closed: eval_modes outside grid errors')
            call check(all_nan_c(dA) .and. all_nan_c(dB), &
                'fail closed: outputs are NaN, not zero')
            call field%eval([outside(1, i), 0.3_dp, outside(2, i)], rA, rB, ierr)
            call check(ierr /= 0, 'fail closed: eval outside grid errors')
            call check(all(ieee_is_nan(rB)), 'fail closed: real dB is NaN')
            call field%dbmod([outside(1, i), 0.3_dp, outside(2, i)], &
                [0.0_dp, 1.0_dp, 0.0_dp], dbm, ierr)
            call check(ierr /= 0 .and. ieee_is_nan(dbm), 'fail closed: dbmod')
        end do
        call field%eval_modes(R_lo, Z_hi, dA, dB, ierr)
        call check(ierr == 0, 'grid corner is inside')
        call field%eval([1.5_dp, 0.0_dp, 0.0_dp], rA, rB, ierr)
        call field%dbmod([1.5_dp, 0.0_dp, 0.0_dp], [0.0_dp, 0.0_dp, 0.0_dp], dbm, ierr)
        call check(ierr /= 0, 'fail closed: zero B0 in dbmod errors')

        call bad%init_analytic(-0.1_dp, 1.0_dp, 11, -1.0_dp, 1.0_dp, 11, [1], &
            amp_smooth, UNITS_SI, ierr)
        call check(ierr /= 0, 'reject grid with R_min <= 0')
        call bad%init_analytic(1.0_dp, 2.0_dp, 11, -1.0_dp, 1.0_dp, 11, [1, 1], &
            amp_smooth, UNITS_SI, ierr)
        call check(ierr /= 0, 'reject duplicate toroidal mode numbers')
        call bad%init_analytic(1.0_dp, 2.0_dp, 11, -1.0_dp, 1.0_dp, 11, [1], &
            amp_smooth, 'furlongs', ierr)
        call check(ierr /= 0, 'reject unknown unit system')
    end subroutine test_fail_closed

    real(dp) function ieee_nan()
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan

        ieee_nan = ieee_value(1.0_dp, ieee_quiet_nan)
    end function ieee_nan

    subroutine test_netcdf_round_trip
        character(len=*), parameter :: fname = 'test_perturbation_field.nc'
        type(perturbation_field_t) :: f1, f2
        complex(dp) :: dA1(3, 2), dB1(3, 2), dA2(3, 2), dB2(3, 2), p1(2), p2(2)
        real(dp) :: x(3)
        integer :: i, ierr, ierr1, ierr2
        logical :: same

        call f1%init_analytic(R_lo, R_hi, 23, Z_lo, Z_hi, 19, [-1, 4], amp_smooth, &
            UNITS_GAUSSIAN, ierr, with_potential=.true.)
        if (ierr /= 0) error stop 'init failed'
        call write_perturbation_field_netcdf(f1, fname, ierr)
        call check(ierr == 0, 'netcdf write')
        call read_perturbation_field_netcdf(fname, f2, ierr)
        call check(ierr == 0, 'netcdf read')
        if (ierr /= 0) return
        call check(all(f2%ntor == [-1, 4]), 'round trip: mode numbers')
        call check(f2%units == UNITS_GAUSSIAN, 'round trip: unit system')
        call check(f2%has_potential, 'round trip: potential present')
        call check(all(f2%R == f1%R) .and. all(f2%Z == f1%Z), 'round trip: grid')
        same = all(f2%dA == f1%dA) .and. all(f2%dPhi == f1%dPhi)
        call check(same, 'round trip: amplitudes bit-identical')
        same = .true.
        do i = 1, n_pts
            call random_point(x)
            call f1%eval_modes(x(1), x(3), dA1, dB1, ierr1, p1)
            call f2%eval_modes(x(1), x(3), dA2, dB2, ierr2, p2)
            same = same .and. ierr1 == 0 .and. ierr2 == 0
            same = same .and. all(dB1 == dB2) .and. all(p1 == p2)
        end do
        call check(same, 'round trip: evaluation identical')
        call delete_file(fname)
    end subroutine test_netcdf_round_trip

    subroutine test_netcdf_rejects_bad_header
        use netcdf, only: nf90_create, nf90_def_dim, nf90_put_att, nf90_enddef, &
            nf90_close, NF90_CLOBBER, NF90_GLOBAL
        character(len=*), parameter :: fname = 'test_perturbation_field_bad.nc'
        type(perturbation_field_t) :: f
        integer :: ncid, dimid, ierr, st

        st = nf90_create(fname, NF90_CLOBBER, ncid)
        st = nf90_def_dim(ncid, 'R', 4, dimid)
        st = nf90_put_att(ncid, NF90_GLOBAL, 'format', 'something_else')
        st = nf90_enddef(ncid)
        st = nf90_close(ncid)
        call read_perturbation_field_netcdf(fname, f, ierr)
        call check(ierr /= 0, 'netcdf read rejects foreign file')
        call read_perturbation_field_netcdf('does_not_exist.nc', f, ierr)
        call check(ierr /= 0, 'netcdf read rejects missing file')
        call delete_file(fname)
    end subroutine test_netcdf_rejects_bad_header

    subroutine delete_file(fname)
        character(len=*), intent(in) :: fname

        integer :: u

        open (newunit=u, file=fname, status='old')
        close (u, status='delete')
    end subroutine delete_file

end program test_perturbation_field
