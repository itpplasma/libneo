program test_cyl_grid_field
    use util_for_test, only: print_test, print_ok, print_fail
    use neo_cyl_grid_field, only: cyl_grid_field_t
    use neo_mfbe2001, only: read_mfbe2001, write_mfbe2001
    implicit none

    integer, parameter :: dp = kind(1.0d0)
    integer, parameter :: nfp = 5
    real(dp), parameter :: pi = 3.14159265358979323846_dp
    real(dp), parameter :: r0 = 5.5_dp, b0 = 2.5_dp

    call test_interpolation_order()
    call test_div_b()
    call test_periodicity_and_outside()
    call test_mfbe2001_round_trip()
    call test_archived_file()

contains

    !> Divergence-free test field B = B0 R0/R e_phi + curl A with
    !> A_R = a R sin(Z) cos(k phi), A_phi = c sin(R) cos(Z),
    !> A_Z = d cos(R) sin(Z) sin(k phi), k = nfp.
    subroutine analytic_b(x, B)
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: B(3)
        real(dp), parameter :: a = 0.1_dp, c = 0.2_dp, d = 0.1_dp
        real(dp) :: R, p, Z, k

        R = x(1)
        p = x(2)
        Z = x(3)
        k = real(nfp, dp)
        B(1) = d*k*cos(R)*sin(Z)*cos(k*p)/R + c*sin(R)*sin(Z)
        B(2) = b0*r0/R + a*R*cos(Z)*cos(k*p) + d*sin(R)*sin(Z)*sin(k*p)
        B(3) = (c*(sin(R) + R*cos(R))*cos(Z) + a*k*R*sin(Z)*sin(k*p))/R
    end subroutine analytic_b

    subroutine make_field(n, nphi, field)
        integer, intent(in) :: n, nphi
        type(cyl_grid_field_t), intent(out) :: field
        real(dp), allocatable :: r(:), z(:), b(:, :, :, :)
        real(dp) :: x(3)
        integer :: i, j, k

        allocate (r(n), z(n), b(3, n, n, nphi))
        do i = 1, n
            r(i) = 4.0_dp + 3.0_dp*real(i - 1, dp)/real(n - 1, dp)
            z(i) = -1.5_dp + 3.0_dp*real(i - 1, dp)/real(n - 1, dp)
        end do
        do k = 1, nphi
            do j = 1, n
                do i = 1, n
                    x = [r(i), 2.0_dp*pi*real(k - 1, dp)/real(nfp*nphi, dp), z(j)]
                    call analytic_b(x, b(:, i, j, k))
                end do
            end do
        end do
        call field%cyl_grid_field_init(r, z, nfp, b)
    end subroutine make_field

    function sample_point(m) result(x)
        integer, intent(in) :: m
        real(dp) :: x(3)
        x(1) = 4.6_dp + 1.8_dp*fract(0.6180339887_dp*m)
        x(2) = 2.0_dp*pi*fract(0.7548776662_dp*m)
        x(3) = -0.9_dp + 1.8_dp*fract(0.5698402910_dp*m)
    end function sample_point

    real(dp) function fract(y)
        real(dp), intent(in) :: y
        fract = y - floor(y)
    end function fract

    real(dp) function max_error(field)
        type(cyl_grid_field_t), intent(in) :: field
        real(dp) :: x(3), B(3), Bx(3)
        integer :: m

        max_error = 0.0_dp
        do m = 1, 400
            x = sample_point(m)
            call field%compute_bfield(x, B)
            call analytic_b(x, Bx)
            max_error = max(max_error, maxval(abs(B - Bx)))
        end do
    end function max_error

    subroutine test_interpolation_order()
        type(cyl_grid_field_t) :: coarse, fine
        real(dp) :: e1, e2, order

        call print_test("test_interpolation_order")
        call make_field(20, 16, coarse)
        call make_field(39, 32, fine)
        e1 = max_error(coarse)
        e2 = max_error(fine)
        ! grid spacing ratio (20-1)/(39-1) = 1/2 in R, Z; phi planes 16 -> 32
        order = log(e1/e2)/log(2.0_dp)
        print *, "  max error coarse/fine:", e1, e2, " observed order:", order
        if (e2 > 1.0e-5_dp .or. order < 4.5_dp) then
            call print_fail
            error stop "interpolation error order below fifth order"
        end if
        call print_ok
    end subroutine test_interpolation_order

    subroutine test_div_b()
        type(cyl_grid_field_t) :: field
        real(dp) :: x(3), B(3), dB(3, 3), div, scale, worst
        integer :: m

        call print_test("test_div_b")
        call make_field(39, 32, field)
        worst = 0.0_dp
        do m = 1, 400
            x = sample_point(m)
            call field%compute_bfield_jacobian(x, B, dB)
            div = B(1)/x(1) + dB(1, 1) + dB(2, 2)/x(1) + dB(3, 3)
            scale = abs(B(1)/x(1)) + abs(dB(1, 1)) + abs(dB(2, 2)/x(1)) + abs(dB(3, 3))
            worst = max(worst, abs(div)/scale)
        end do
        print *, "  max |div B| / sum |terms|:", worst
        if (worst > 1.0e-4_dp) then
            call print_fail
            error stop "interpolated field is not divergence free"
        end if
        call print_ok
    end subroutine test_div_b

    subroutine test_periodicity_and_outside()
        type(cyl_grid_field_t) :: field
        real(dp) :: x(3), B1(3), B2(3), period

        call print_test("test_periodicity_and_outside")
        call make_field(39, 32, field)
        period = 2.0_dp*pi/real(nfp, dp)
        x = [5.3_dp, 0.37_dp, 0.21_dp]
        call field%compute_bfield(x, B1)
        x(2) = x(2) + 3.0_dp*period
        call field%compute_bfield(x, B2)
        if (maxval(abs(B1 - B2)) > 1.0e-12_dp) then
            call print_fail
            error stop "field not periodic in phi"
        end if
        x = [3.0_dp, 0.1_dp, 0.0_dp]
        call field%compute_bfield(x, B1)
        if (field%contains_point(x) .or. any(B1 /= 0.0_dp)) then
            call print_fail
            error stop "outside point not flagged"
        end if
        if (.not. field%contains_point([5.5_dp, 0.0_dp, 0.0_dp])) then
            call print_fail
            error stop "inside point flagged as outside"
        end if
        call print_ok
    end subroutine test_periodicity_and_outside

    subroutine same_field(a, b, what)
        type(cyl_grid_field_t), intent(in) :: a, b
        character(*), intent(in) :: what
        if (a%nr /= b%nr .or. a%nz /= b%nz .or. a%nphi /= b%nphi .or. &
            a%nfp /= b%nfp) error stop what//": grid shape differs"
        if (abs(a%r_min - b%r_min) > 1.0e-12_dp .or. abs(a%dr - b%dr) > 1.0e-12_dp .or. &
            abs(a%z_min - b%z_min) > 1.0e-12_dp .or. abs(a%dz - b%dz) > 1.0e-12_dp) &
            error stop what//": grid geometry differs"
        if (any(a%b /= b%b) .or. any(a%valid .neqv. b%valid)) &
            error stop what//": field data differ"
    end subroutine same_field

    subroutine test_mfbe2001_round_trip()
        type(cyl_grid_field_t) :: field, back
        character(len=*), parameter :: fname = "test_cyl_grid_field_rt.mfbe"
        integer :: unit

        call print_test("test_mfbe2001_round_trip (synthetic)")
        call make_field(24, 12, field)
        field%b(:, 1:2, :, :) = 0.0_dp     ! outside-domain fill
        field%valid(1:2, :, :) = .false.
        call write_mfbe2001(fname, field)
        call read_mfbe2001(fname, back)
        call same_field(field, back, "synthetic")
        open (newunit=unit, file=fname, status="old")
        close (unit, status="delete")
        call print_ok
    end subroutine test_mfbe2001_round_trip

    !> Optional: file named by HINT2_MFBE2001 (an archived HINT2 field).
    subroutine test_archived_file()
        type(cyl_grid_field_t) :: field, back
        character(len=1024) :: path
        character(len=*), parameter :: fname = "test_cyl_grid_field_arch.mfbe"
        integer :: unit, stat
        real(dp) :: x(3), B(3), bmod

        call get_environment_variable("HINT2_MFBE2001", path, status=stat)
        if (stat /= 0) then
            print *, "test_archived_file: HINT2_MFBE2001 not set, skipped"
            return
        end if
        call print_test("test_mfbe2001_round_trip (archived file)")
        call read_mfbe2001(trim(path), field)
        print *, "  nr nz nphi nfp:", field%nr, field%nz, field%nphi, field%nfp
        print *, "  R_min dR Z_min dZ:", field%r_min, field%dr, field%z_min, field%dz
        call write_mfbe2001(fname, field)
        call read_mfbe2001(fname, back)
        call same_field(field, back, "archived")
        open (newunit=unit, file=fname, status="old")
        close (unit, status="delete")
        ! node values are reproduced exactly by the interpolation
        x = [field%r_min + 70.0_dp*field%dr, 0.0_dp, field%z_min + 75.0_dp*field%dz]
        call field%compute_bfield(x, B)
        bmod = sqrt(sum(B**2))
        print *, "  |B| at (node 71, plane 1, node 76):", bmod
        if (any(abs(B - field%b(:, 71, 76, 1)) > 1.0e-10_dp*max(bmod, 1.0_dp))) &
            error stop "interpolation does not reproduce nodes"
        call print_ok
    end subroutine test_archived_file

end program test_cyl_grid_field
