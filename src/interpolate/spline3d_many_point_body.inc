if (periodic(1)) then
    t = x(1, ipt) - x_min(1)
    k_wrap = int(t/period(1))
    w = t - real(k_wrap, dp)*period(1)
    if (w < 0.0_dp) w = w + period(1)
    xj1 = w + x_min(1)
else
    xj1 = x(1, ipt)
end if
if (periodic(2)) then
    t = x(2, ipt) - x_min(2)
    k_wrap = int(t/period(2))
    w = t - real(k_wrap, dp)*period(2)
    if (w < 0.0_dp) w = w + period(2)
    xj2 = w + x_min(2)
else
    xj2 = x(2, ipt)
end if
if (periodic(3)) then
    t = x(3, ipt) - x_min(3)
    k_wrap = int(t/period(3))
    w = t - real(k_wrap, dp)*period(3)
    if (w < 0.0_dp) w = w + period(3)
    xj3 = w + x_min(3)
else
    xj3 = x(3, ipt)
end if

x_norm1 = (xj1 - x_min(1))/h_step(1)
x_norm2 = (xj2 - x_min(2))/h_step(2)
x_norm3 = (xj3 - x_min(3))/h_step(3)

i1 = max(0, min(num_points(1) - 2, int(x_norm1)))
i2 = max(0, min(num_points(2) - 2, int(x_norm2)))
i3 = max(0, min(num_points(3) - 2, int(x_norm3)))

x_local1 = (x_norm1 - real(i1, dp))*h_step(1)
x_local2 = (x_norm2 - real(i2, dp))*h_step(2)
x_local3 = (x_norm3 - real(i3, dp))*h_step(3)

do iq = 1, nq
    ! k3 = order3 term (2D Horner in x1,x2)
    v = spl%coeff(iq, order1, order2, order3, i1 + 1, i2 + 1, i3 + 1)
    do k1 = order1 - 1, 0, -1
        v = spl%coeff(iq, k1, order2, order3, i1 + 1, i2 + 1, i3 + 1) + x_local1*v
    end do
    w2 = v
    do k2 = order2 - 1, 0, -1
        v = spl%coeff(iq, order1, k2, order3, i1 + 1, i2 + 1, i3 + 1)
        do k1 = order1 - 1, 0, -1
            v = spl%coeff(iq, k1, k2, order3, i1 + 1, i2 + 1, i3 + 1) + x_local1*v
        end do
        w2 = v + x_local2*w2
    end do
    yq = w2

    do k3 = order3 - 1, 0, -1
        v = spl%coeff(iq, order1, order2, k3, i1 + 1, i2 + 1, i3 + 1)
        do k1 = order1 - 1, 0, -1
            v = spl%coeff(iq, k1, order2, k3, i1 + 1, i2 + 1, i3 + 1) + x_local1*v
        end do
        w2 = v
        do k2 = order2 - 1, 0, -1
            v = spl%coeff(iq, order1, k2, k3, i1 + 1, i2 + 1, i3 + 1)
            do k1 = order1 - 1, 0, -1
                v = spl%coeff(iq, k1, k2, k3, i1 + 1, i2 + 1, i3 + 1) + x_local1*v
            end do
            w2 = v + x_local2*w2
        end do
        yq = w2 + x_local3*yq
    end do

    y_batch(iq, ipt) = yq
end do

