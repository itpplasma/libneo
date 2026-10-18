!> Pointwise mathematics of the regular near-axis Boozer field.
!>
!> |B| is represented by its second-order expansion about the magnetic axis
!> plus a cubic remainder,
!>     B = F0 + X F1c + Y F1s + s F20 + (X**2 - Y**2) F2c + 2 X Y F2s
!>         + rho**3 R(rho, theta, phi),
!> with axis functions F*(phi), X = rho cos(theta), Y = rho sin(theta) and
!> s = rho**2. Flux functions are I = s i(rho), G = G0 + s g(rho),
!> A_phi = A0 + s a(rho) and A_theta = torflux s. The polar chain rule uses
!> d/dX = cos d/drho - sin/rho d/dtheta and d/dY = sin d/drho +
!> cos/rho d/dtheta; every 1/rho cancels after expansion, so the emitted
!> leaves contain no division by rho and are finite on the axis. The quotient
!> leaves form h_theta = I/B and h_phi = G/B together, sharing 1/B.
program gen_boozer_regular
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, func, cos, sin, operator(+), &
                            operator(-), operator(*), operator(/), operator(**)
    use fortsym_diff, only: diff
    use fortsym_subs, only: subs_many
    use fortsym_kernel, only: kernel_binding_t
    use libneo_gen_support, only: leaf_writer_t, fortsym_revision, binding
    implicit none
    character(8), parameter :: axis_names(6) = [character(8) :: 'f0', 'f1c', &
        'f1s', 'f20', 'f2c', 'f2s']
    type(leaf_writer_t) :: w
    type(arena_t), target :: arena
    type(expr_t) :: rho, th, ph, x, y, c, sn, g0, ap0, torflux, ct, st
    type(expr_t) :: bfield, iflux, gflux, aphi, atheta, denom
    type(expr_t) :: numer(2), quot(2)
    type(expr_t) :: axis(6), polar_old(2), polar_new(2)
    type(expr_t) :: trig_old(3), trig_new(3)
    type(kernel_binding_t), allocatable :: polar_bind(:), quot_bind(:)
    character(len=256) :: path
    integer :: k

    path = 'src/generated/boozer_regular_leaves.inc'
    if (command_argument_count() >= 1) call get_command_argument(1, path)

    call arena%init()
    rho = sym(arena, 'rho')
    th = sym(arena, 'theta')
    ph = sym(arena, 'phi')
    x = sym(arena, 'x')
    y = sym(arena, 'y')
    c = sym(arena, 'c')
    sn = sym(arena, 'sn')
    g0 = sym(arena, 'g0')
    ap0 = sym(arena, 'ap0')
    torflux = sym(arena, 'torflux')
    ct = cos(th)
    st = sin(th)
    polar_old = [ct, st]
    polar_new = [c, sn]
    trig_old = [sn**2, sn**3, sn**4]
    trig_new = [1 - c**2, sn*(1 - c**2), (1 - c**2)**2]

    allocate (polar_bind(0))
    do k = 1, size(axis_names)
        axis(k) = func(upper(axis_names(k)), [ph])
        polar_bind = [polar_bind, &
                      binding(upper(axis_names(k)), [integer ::], &
                              trim(axis_names(k))), &
                      binding(upper(axis_names(k)), [1], &
                              trim(axis_names(k))//'_p'), &
                      binding(upper(axis_names(k)), [1, 1], &
                              trim(axis_names(k))//'_pp')]
    end do
    bfield = axis(1) + rho*ct*axis(2) + rho*st*axis(3) + &
             rho**2*(axis(4) + (ct**2 - st**2)*axis(5) + 2*ct*st*axis(6)) + &
             rho**3*func('R', [rho, th, ph])
    iflux = rho**2*func('FI', [rho])
    gflux = g0 + rho**2*func('FG', [rho])
    aphi = ap0 + rho**2*func('FA', [rho])
    atheta = torflux*rho**2
    polar_bind = [polar_bind, &
                  binding('R', [integer ::], 'r'), binding('R', [1], 'r_r'), &
                  binding('R', [2], 'r_t'), binding('R', [3], 'r_p'), &
                  binding('R', [1, 1], 'r_rr'), binding('R', [1, 2], 'r_rt'), &
                  binding('R', [1, 3], 'r_rp'), binding('R', [2, 2], 'r_tt'), &
                  binding('R', [2, 3], 'r_tp'), binding('R', [3, 3], 'r_pp'), &
                  binding('FI', [integer ::], 'fi'), binding('FI', [1], 'fi_r'), &
                  binding('FI', [1, 1], 'fi_rr'), &
                  binding('FG', [integer ::], 'fg'), binding('FG', [1], 'fg_r'), &
                  binding('FG', [1, 1], 'fg_rr'), &
                  binding('FA', [integer ::], 'fa'), binding('FA', [1], 'fa_r'), &
                  binding('FA', [1, 1], 'fa_rr')]

    numer = [func('NI', [x, y]), func('NG', [x, y])]
    denom = func('B', [x, y, ph])
    quot = [numer(1)/denom, numer(2)/denom]
    quot_bind = [binding('NI', [integer ::], 'i'), binding('NI', [1], 'i_x'), &
                 binding('NI', [2], 'i_y'), binding('NI', [1, 1], 'i_xx'), &
                 binding('NI', [1, 2], 'i_xy'), binding('NI', [2, 2], 'i_yy'), &
                 binding('NG', [integer ::], 'g'), binding('NG', [1], 'g_x'), &
                 binding('NG', [2], 'g_y'), binding('NG', [1, 1], 'g_xx'), &
                 binding('NG', [1, 2], 'g_xy'), binding('NG', [2, 2], 'g_yy'), &
                 binding('B', [integer ::], 'b'), binding('B', [1], 'b_x'), &
                 binding('B', [2], 'b_y'), binding('B', [3], 'b_p'), &
                 binding('B', [1, 1], 'b_xx'), binding('B', [1, 2], 'b_xy'), &
                 binding('B', [1, 3], 'b_xp'), binding('B', [2, 2], 'b_yy'), &
                 binding('B', [2, 3], 'b_yp'), binding('B', [3, 3], 'b_pp')]

    call w%open(arena, 'gen_boozer_regular', trim(path))
    call w%line('! Generated by gen/gen_boozer_regular.f90; do not edit.')
    call w%line('! FortSym revision: '//fortsym_revision)
    call w%line('! Regular near-axis Boozer field leaves, included by')
    call w%line('! src/boozer_regular.f90. Inputs are axis jets in phi and polar')
    call w%line('! spline jets at (rho, theta, phi) with c = cos(theta) and')
    call w%line('! sn = sin(theta); outputs are values and derivatives in')
    call w%line('! (X, Y, phi). Second-derivative order: xx, xy, xp, yy, yp, pp.')
    call emit_polar_leaves()
    call emit_quotient_leaves()
    call w%close()

contains

    pure function upper(name) result(u)
        character(*), intent(in) :: name
        character(len_trim(name)) :: u
        integer :: j

        u = trim(name)
        do j = 1, len(u)
            if (u(j:j) >= 'a' .and. u(j:j) <= 'z') &
                u(j:j) = achar(iachar(u(j:j)) - 32)
        end do
    end function

    function dx(e) result(r)
        type(expr_t), intent(in) :: e
        type(expr_t) :: r
        r = w%expand(ct*diff(e, rho) - st*diff(e, th)/rho)
    end function

    function dy(e) result(r)
        type(expr_t), intent(in) :: e
        type(expr_t) :: r
        r = w%expand(st*diff(e, rho) + ct*diff(e, th)/rho)
    end function

    function dp(e) result(r)
        type(expr_t), intent(in) :: e
        type(expr_t) :: r
        r = w%expand(diff(e, ph))
    end function

    !> Replace cos/sin(theta) by the inputs c, sn and reduce powers of sn with
    !> sn**2 = 1 - c**2; the leaves are at most quartic in (c, sn).
    function cart(e) result(r)
        type(expr_t), intent(in) :: e
        type(expr_t) :: r
        r = w%expand(subs_many(e, polar_old, polar_new))
        r = w%expand(subs_many(r, trig_old, trig_new))
    end function

    subroutine emit_polar_leaves()
        character(8), parameter :: base(4) = [character(8) :: 'rho', 'c', 'sn', &
                                               'torflux']
        character(8), parameter :: jet0(6) = [character(8) :: 'r', 'fi', 'g0', &
                                              'fg', 'ap0', 'fa']
        character(8), parameter :: jet1(6) = [character(8) :: 'r_r', 'r_t', &
                                              'r_p', 'fi_r', 'fg_r', 'fa_r']
        character(8), parameter :: jet2(9) = [character(8) :: 'r_rr', 'r_rt', &
            'r_rp', 'r_tt', 'r_tp', 'r_pp', 'fi_rr', 'fg_rr', 'fa_rr']
        character(8) :: axis1(6), axis2(6)
        type(expr_t) :: b1(3), i1(2), g1(2), a1(2), t1(2)
        integer :: j

        do j = 1, 6
            axis1(j) = trim(axis_names(j))//'_p'
            axis2(j) = trim(axis_names(j))//'_pp'
        end do

        call w%emit('boozer_regular_polar_d0', [base, axis_names, jet0], &
                    [character(8) :: 'b', 'i', 'g', 'a', 'at'], &
                    [cart(bfield), cart(iflux), cart(gflux), cart(aphi), &
                     cart(atheta)], polar_bind)

        b1 = [dx(bfield), dy(bfield), dp(bfield)]
        i1 = [dx(iflux), dy(iflux)]
        g1 = [dx(gflux), dy(gflux)]
        a1 = [dx(aphi), dy(aphi)]
        t1 = [dx(atheta), dy(atheta)]
        call w%emit('boozer_regular_polar_d1', [base, axis_names, jet0, axis1, &
                                                jet1], &
                    [character(8) :: 'b', 'b_x', 'b_y', 'b_p', 'i', 'i_x', &
                     'i_y', 'g', 'g_x', 'g_y', 'a', 'a_x', 'a_y', 'at', 'at_x', &
                     'at_y'], &
                    [cart(bfield), cart(b1(1)), cart(b1(2)), cart(b1(3)), &
                     cart(iflux), cart(i1(1)), cart(i1(2)), &
                     cart(gflux), cart(g1(1)), cart(g1(2)), &
                     cart(aphi), cart(a1(1)), cart(a1(2)), &
                     cart(atheta), cart(t1(1)), cart(t1(2))], polar_bind)

        call w%emit('boozer_regular_polar_d2', [base, axis_names, jet0, axis1, &
                                                jet1, axis2, jet2], &
                    [character(8) :: 'b', 'b_x', 'b_y', 'b_p', 'b_xx', 'b_xy', &
                     'b_xp', 'b_yy', 'b_yp', 'b_pp', 'i', 'i_x', 'i_y', 'i_xx', &
                     'i_xy', 'i_yy', 'g', 'g_x', 'g_y', 'g_xx', 'g_xy', 'g_yy', &
                     'a', 'a_x', 'a_y', 'a_xx', 'a_xy', 'a_yy', 'at', 'at_x', &
                     'at_y', 'at_xx', 'at_xy', 'at_yy'], &
                    [cart(bfield), cart(b1(1)), cart(b1(2)), cart(b1(3)), &
                     cart(dx(b1(1))), cart(dy(b1(1))), cart(dp(b1(1))), &
                     cart(dy(b1(2))), cart(dp(b1(2))), cart(dp(b1(3))), &
                     cart(iflux), cart(i1(1)), cart(i1(2)), cart(dx(i1(1))), &
                     cart(dy(i1(1))), cart(dy(i1(2))), &
                     cart(gflux), cart(g1(1)), cart(g1(2)), cart(dx(g1(1))), &
                     cart(dy(g1(1))), cart(dy(g1(2))), &
                     cart(aphi), cart(a1(1)), cart(a1(2)), cart(dx(a1(1))), &
                     cart(dy(a1(1))), cart(dy(a1(2))), &
                     cart(atheta), cart(t1(1)), cart(t1(2)), cart(dx(t1(1))), &
                     cart(dy(t1(1))), cart(dy(t1(2)))], polar_bind)
    end subroutine

    subroutine emit_quotient_leaves()
        character(4), parameter :: ins1(3) = [character(4) :: 'i', 'i_x', 'i_y']
        character(4), parameter :: ins2(3) = [character(4) :: 'i_xx', 'i_xy', 'i_yy']
        character(4), parameter :: bs1(4) = [character(4) :: 'b', 'b_x', 'b_y', &
                                             'b_p']
        character(4), parameter :: bs2(6) = [character(4) :: 'b_xx', 'b_xy', &
            'b_xp', 'b_yy', 'b_yp', 'b_pp']
        character(4), parameter :: d1(4) = [character(4) :: '', '_x', '_y', '_p']
        character(4), parameter :: d2(6) = [character(4) :: '_xx', '_xy', '_xp', &
            '_yy', '_yp', '_pp']
        type(expr_t) :: q1(3, 2), roots1(8), roots2(20)
        character(8) :: out1(8), out2(20)
        integer :: k, j

        do k = 1, 2
            q1(:, k) = [diff(quot(k), x), diff(quot(k), y), diff(quot(k), ph)]
            roots1(4*k - 3:4*k) = [quot(k), q1(:, k)]
            roots2(10*k - 9:10*k) = [quot(k), q1(:, k), diff(q1(1, k), x), &
                                     diff(q1(1, k), y), diff(q1(1, k), ph), &
                                     diff(q1(2, k), y), diff(q1(2, k), ph), &
                                     diff(q1(3, k), ph)]
            do j = 1, 4
                out1(4*k - 4 + j) = trim(merge('hth', 'hph', k == 1))//trim(d1(j))
                out2(10*k - 10 + j) = out1(4*k - 4 + j)
            end do
            do j = 1, 6
                out2(10*k - 6 + j) = trim(merge('hth', 'hph', k == 1))//trim(d2(j))
            end do
        end do
        call w%emit('boozer_regular_quotient_d1', &
                    [ins1, gname(ins1), bs1], out1, roots1, quot_bind)
        call w%emit('boozer_regular_quotient_d2', &
                    [ins1, ins2, gname(ins1), gname(ins2), bs1, bs2], out2, roots2, &
                    quot_bind)
    end subroutine

    pure function gname(names) result(g)
        character(4), intent(in) :: names(:)
        character(4) :: g(size(names))
        integer :: j

        do j = 1, size(names)
            g(j) = 'g'//names(j)(2:)
        end do
    end function
end program gen_boozer_regular
