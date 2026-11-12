program test_lagrange_speed
    ! The generated cubic point kernel behind polylag_3::plag3d must agree with
    ! the dense evaluation that forms every node product, and be at least twice
    ! as fast as it in the same binary, so a slower schedule fails here instead
    ! of slowing every field evaluation in the consumers.
    use libneo_kinds, only: dp
    use polylag_3, only: plag3d, mp
    use plag_coeff_sub, only: plag_coeff_uniform
    implicit none
    integer, parameter :: ncall = 400000
    real(dp) :: fp(mp, mp, mp), xp(mp), yp(mp), zp(mp), r(4), d(4), s1, s2
    real(dp) :: t0, t1, t_gen, t_dense, x, worst
    integer :: i, j, rep
    character(16) :: arg

    xp = 0.3_dp + 0.25_dp*[0, 1, 2, 3]
    yp = 2.0_dp + 0.1_dp*[0, 1, 2, 3]
    zp = -1.0_dp + 0.4_dp*[0, 1, 2, 3]
    do i = 1, mp
        do j = 1, mp
            fp(:, i, j) = sin(1.3_dp*[1, 2, 3, 4] + 0.7_dp*i - 0.4_dp*j)
        end do
    end do
    worst = 0
    do i = 1, 50
        x = xp(1) + 0.75_dp*i/50
        call plag3d(x, 2.13_dp, -0.6_dp, fp, 4.0_dp, 10.0_dp, 2.5_dp, xp, yp, zp, &
                    r(1), r(2), r(3), r(4))
        call dense(x, 2.13_dp, -0.6_dp, d)
        worst = max(worst, maxval(abs(r - d))/maxval(abs(d)))
    end do
    print '(a,es10.3)', 'generated vs dense, largest relative difference', worst
    if (worst > 1.0e-13_dp) error stop 'generated kernel disagrees'

    t_gen = huge(1.0_dp)
    t_dense = huge(1.0_dp)
    s1 = 0
    s2 = 0
    do rep = 1, 3
        call cpu_time(t0)
        do i = 1, ncall
            x = xp(1) + 0.75_dp*mod(i, 997)/997
            call plag3d(x, 2.13_dp, -0.6_dp, fp, 4.0_dp, 10.0_dp, 2.5_dp, xp, yp, zp, &
                        r(1), r(2), r(3), r(4))
            s1 = s1 + sum(r)
        end do
        call cpu_time(t1)
        t_gen = min(t_gen, t1 - t0)
        call cpu_time(t0)
        do i = 1, ncall
            x = xp(1) + 0.75_dp*mod(i, 997)/997
            call dense(x, 2.13_dp, -0.6_dp, d)
            s2 = s2 + sum(d)
        end do
        call cpu_time(t1)
        t_dense = min(t_dense, t1 - t0)
    end do
    print '(a,2f8.1,a,es12.4)', 'ns per call generated, dense', t_gen/ncall*1e9_dp, &
        t_dense/ncall*1e9_dp, '  checksum', s1 - s2
    ! Without optimisation the ratio says nothing about the kernel schedule.
    call get_command_argument(1, arg)
    if (arg == '--no-speed') stop
    if (t_dense < 2*t_gen) error stop 'generated kernel not twice as fast as dense'
contains
    subroutine dense(x, y, z, out)
        real(dp), intent(in) :: x, y, z
        real(dp), intent(out) :: out(4)
        real(dp) :: cx(0:1, mp), cy(0:1, mp), cz(0:1, mp), w
        integer :: a, b, c

        call plag_coeff_uniform(mp, 1, x, xp(1), 4.0_dp, cx)
        call plag_coeff_uniform(mp, 1, y, yp(1), 10.0_dp, cy)
        call plag_coeff_uniform(mp, 1, z, zp(1), 2.5_dp, cz)
        out = 0
        do c = 1, mp
            do b = 1, mp
                do a = 1, mp
                    w = fp(a, b, c)
                    out(1) = out(1) + w*cx(0, a)*cy(0, b)*cz(0, c)
                    out(2) = out(2) + w*cx(1, a)*cy(0, b)*cz(0, c)
                    out(3) = out(3) + w*cx(0, a)*cy(1, b)*cz(0, c)
                    out(4) = out(4) + w*cx(0, a)*cy(0, b)*cz(1, c)
                end do
            end do
        end do
    end subroutine
end program test_lagrange_speed
