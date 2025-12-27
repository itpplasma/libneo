    do j = 1, 3
        if (spl%periodic(j)) then
            xj = x(j)
            if (xj < x_min(j) .or. xj >= x_min(j) + period(j)) then
                xj = modulo(xj - x_min(j), period(j)) + x_min(j)
            end if
        else
            xj = x(j)
        end if
        x_norm(j) = (xj - x_min(j))*inv_h_step(j)
        interval_index(j) = max(0, min(spl%num_points(j) - 2, int(x_norm(j))))
        x_local(j) = (x_norm(j) - real(interval_index(j), dp))*h_step(j)
    end do

    x1 = x_local(1)
    x2 = x_local(2)
    x3 = x_local(3)

    i1 = interval_index(1) + 1
    i2 = interval_index(2) + 1
    i3 = interval_index(3) + 1
