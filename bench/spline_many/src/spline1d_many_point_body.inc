base = (ipt - 1) * num_quantities

if (periodic_int == 1) then
    t = x(ipt) - x_min
    k_wrap = int(t / period)
    w = t - real(k_wrap, dp) * period
    if (w < 0.0_dp) w = w + period
    xj = w + x_min
else
    xj = x(ipt)
end if

x_norm = (xj - x_min) / h_step
idx = max(0, min(num_points - 2, int(x_norm)))
x_local = (x_norm - real(idx, dp)) * h_step

do iq = 1, num_quantities
    y(base + iq) = coeff(iq, order, idx + 1)
end do

do k_power = order - 1, 0, -1
    do iq = 1, num_quantities
        y(base + iq) = coeff(iq, k_power, idx + 1) + x_local * y(base + iq)
    end do
end do

