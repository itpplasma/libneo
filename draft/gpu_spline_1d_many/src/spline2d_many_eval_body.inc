#if defined(SPLINE2D_BACKEND_OPENACC)
        !$acc parallel loop present(coeff, x, y) &
        !$acc& private(ipt, iq, k1, k2, i1, i2, base, xj1, xj2, x_norm1, x_norm2, &
        !$acc& x_local1, x_local2, t, w, v, yq, k_wrap) gang vector vector_length(256)
#elif defined(SPLINE2D_BACKEND_OMP_TARGET)
!$omp target teams distribute parallel do thread_limit(256) &
!$omp& private(ipt, iq, k1, k2, i1, i2, base, xj1, xj2, x_norm1, x_norm2, x_local1, &
!$omp& x_local2, t, w, v, yq, k_wrap)
#endif

        do ipt = 1, npts
            base = (ipt - 1) * num_quantities

            if (periodic_int(1) == 1) then
                t = x(1, ipt) - x_min(1)
                k_wrap = int(t / period(1))
                w = t - real(k_wrap, dp) * period(1)
                if (w < 0.0_dp) w = w + period(1)
                xj1 = w + x_min(1)
            else
                xj1 = x(1, ipt)
            end if
            if (periodic_int(2) == 1) then
                t = x(2, ipt) - x_min(2)
                k_wrap = int(t / period(2))
                w = t - real(k_wrap, dp) * period(2)
                if (w < 0.0_dp) w = w + period(2)
                xj2 = w + x_min(2)
            else
                xj2 = x(2, ipt)
            end if

            x_norm1 = (xj1 - x_min(1)) / h_step(1)
            x_norm2 = (xj2 - x_min(2)) / h_step(2)

            i1 = max(0, min(num_points(1) - 2, int(x_norm1)))
            i2 = max(0, min(num_points(2) - 2, int(x_norm2)))

            x_local1 = (x_norm1 - real(i1, dp)) * h_step(1)
            x_local2 = (x_norm2 - real(i2, dp)) * h_step(2)

            do iq = 1, num_quantities
                v = coeff(iq, order1, order2, i1 + 1, i2 + 1)
                do k1 = order1 - 1, 0, -1
                    v = coeff(iq, k1, order2, i1 + 1, i2 + 1) + x_local1 * v
                end do
                yq = v
                do k2 = order2 - 1, 0, -1
                    v = coeff(iq, order1, k2, i1 + 1, i2 + 1)
                    do k1 = order1 - 1, 0, -1
                        v = coeff(iq, k1, k2, i1 + 1, i2 + 1) + x_local1 * v
                    end do
                    yq = v + x_local2 * yq
                end do
                y(base + iq) = yq
            end do
        end do

#if defined(SPLINE2D_BACKEND_OPENACC)
        !$acc end parallel loop
#elif defined(SPLINE2D_BACKEND_OMP_TARGET)
!$omp end target teams distribute parallel do
#endif
