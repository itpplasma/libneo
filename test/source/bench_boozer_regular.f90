!> Cost of the regular near-axis Boozer evaluation against today's spline
!> evaluation on the LandremanPaul2021 QA field (internal VMEC->Boozer,
!> quintic splines, grid refinement 5).
!>
!> Per evaluation, at random points in s, theta, phi:
!>   spline3d      one 3D |B| spline value and gradient (the bare kernel)
!>   splint m0/m2  splint_boozer_coord with mode_secders 0 and 2
!>   regular m0..2 eval_boozer_regular modes 0, 1, 2
!> Points are either uniformly random (cache-cold field data) or orbit-like
!> (1000 tracks of small steps in s, theta, phi, as during tracing). Each
!> timing is repeated over the thread counts given on the command line
!> (default 1); points are split statically across threads and the wall
!> time per evaluation is reported.
!> Usage: bench_boozer_regular.x [n_points] [random|orbit] [threads ...]
program bench_boozer_regular
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use omp_lib, only: omp_get_wtime, omp_set_num_threads
    use boozer_coordinates_mod, only: use_B_r
    use boozer_sub, only: get_boozer_coordinates, splint_boozer_coord
    use boozer_regular, only: regular_field_t, set_boozer_regular_enabled, &
                              eval_boozer_regular, boozer_regular_field
    use interpolate, only: evaluate_batch_splines_3d_der
    implicit none

    real(dp), parameter :: TWOPI = 6.283185307179586476925286766559_dp
    integer, parameter :: N_REPEAT = 5
    integer :: n_points, n_threads, i_arg, k
    real(dp), allocatable :: s(:), theta(:), phi(:), x(:), y(:)
    character(len=32) :: arg, pattern
    character(len=12), parameter :: names(6) = [character(12) :: 'spline3d', &
        'splint m0', 'splint m2', 'regular m0', 'regular m1', 'regular m2']
    real(dp) :: t_ns(6), sink

    sink = 0.0_dp
    n_points = 1000000
    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg)
        read (arg, *) n_points
    end if

    use_B_r = .false.
    call set_boozer_regular_enabled(.true.)
    call get_boozer_coordinates('wout.nc', radial_spline_order=5, &
                                angular_spline_order=5, grid_refinment=5)

    allocate (s(n_points), theta(n_points), phi(n_points), x(n_points), &
              y(n_points))
    pattern = 'random'
    if (command_argument_count() >= 2) call get_command_argument(2, pattern)
    call random_seed(put=[(12345 + k, k=1, 64)])
    call random_number(s)
    call random_number(theta)
    call random_number(phi)
    if (trim(pattern) == 'orbit') call orbit_points()
    s = 1.0e-6_dp + 0.98_dp*s
    theta = TWOPI*theta
    phi = TWOPI*phi
    x = sqrt(s)*cos(theta)
    y = sqrt(s)*sin(theta)

    print '(a,i0,3a)', 'points: ', n_points, ' (', trim(pattern), &
        '); wall ns per evaluation (best of 5)'
    print '(a8,6a12)', 'threads', names
    i_arg = 3
    do
        n_threads = 1
        if (command_argument_count() >= i_arg) then
            call get_command_argument(i_arg, arg)
            read (arg, *) n_threads
        else if (i_arg > 3) then
            exit
        end if
        call omp_set_num_threads(n_threads)
        do k = 1, 6
            t_ns(k) = best_time(k)*1.0e9_dp/real(n_points, dp)
        end do
        print '(i8,6f12.1,a,f6.3,a,f6.3)', n_threads, t_ns, &
            '   m1/spline3d=', t_ns(5)/t_ns(1), ' m1/splint0=', t_ns(5)/t_ns(2)
        i_arg = i_arg + 1
    end do
    print *, 'checksum', sink

contains

    !> 1000 tracks of consecutive points; s, theta, phi in units of their
    !> ranges, advanced by small steps from random starts.
    subroutine orbit_points()
        integer :: i, track_length

        track_length = max(1, n_points/1000)
        do i = 1, n_points
            if (mod(i - 1, track_length) == 0) cycle
            s(i) = modulo(s(i - 1) + 1.0e-4_dp*(s(i) - 0.5_dp), 1.0_dp)
            theta(i) = theta(i - 1) + 2.0e-3_dp
            phi(i) = phi(i - 1) + 1.0e-3_dp
        end do
    end subroutine orbit_points

    !> Best wall time of N_REPEAT sweeps; reported as wall ns per evaluation.
    real(dp) function best_time(kind) result(best)
        integer, intent(in) :: kind
        real(dp) :: t0, acc
        integer :: rep

        best = huge(1.0_dp)
        do rep = 1, N_REPEAT
            t0 = omp_get_wtime()
            acc = sweep(kind)
            best = min(best, omp_get_wtime() - t0)
            sink = sink + acc
        end do
    end function best_time

    real(dp) function sweep(kind) result(acc)
        integer, intent(in) :: kind
        type(regular_field_t) :: f
        real(dp) :: q(3), b(1), db(3, 1)
        real(dp) :: ath, aph, dath, daph, d2aph, d3aph, bth, dbth, d2bth
        real(dp) :: bph, dbph, d2bph, bm, dbm(3), d2bm(6), br, dbr(3), d2br(6)
        integer :: i

        acc = 0.0_dp
        select case (kind)
        case (1)
            !$omp parallel do schedule(static) private(q, b, db) reduction(+:acc)
            do i = 1, n_points
                q = [sqrt(s(i)), theta(i), phi(i)]
                call evaluate_batch_splines_3d_der(boozer_regular_field%remainder, &
                                                   q, b, db)
                acc = acc + b(1) + db(1, 1)
            end do
        case (2, 3)
            !$omp parallel do schedule(static) reduction(+:acc) &
            !$omp private(ath, aph, dath, daph, d2aph, d3aph, bth, dbth, d2bth) &
            !$omp private(bph, dbph, d2bph, bm, dbm, d2bm, br, dbr, d2br)
            do i = 1, n_points
                call splint_boozer_coord(s(i), theta(i), phi(i), 2*(kind - 2), &
                                         ath, aph, dath, daph, d2aph, d3aph, &
                                         bth, dbth, d2bth, bph, dbph, d2bph, &
                                         bm, dbm, d2bm, br, dbr, d2br)
                acc = acc + bm + dbm(1) + bth + aph + d2bm(1)
            end do
        case default
            !$omp parallel do schedule(static) private(f) reduction(+:acc)
            do i = 1, n_points
                call eval_boozer_regular(x(i), y(i), phi(i), kind - 4, f)
                acc = acc + f%Bmod + f%dBmod(1) + f%hth + f%Aph + f%d2Bmod(1)
            end do
        end select
    end function sweep

end program bench_boozer_regular
