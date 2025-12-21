#if defined(SPLINE1D_BACKEND_OPENACC)
        !$acc parallel loop present(coeff, x, y) &
        !$acc& private(ipt, iq, k_power, idx, base, xj, x_norm, x_local, t, w) &
        !$acc& gang vector vector_length(256)
#elif defined(SPLINE1D_BACKEND_OMP_TARGET)
!$omp target teams distribute parallel do thread_limit(256) &
!$omp& private(ipt, iq, k_power, idx, base, xj, x_norm, x_local, t, w)
#endif

        do ipt = 1, size(x)
            base = (ipt - 1) * num_quantities

            if (periodic_int == 1) then
                t = x(ipt) - x_min
                w = t - floor(t / period) * period
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
        end do

#if defined(SPLINE1D_BACKEND_OPENACC)
        !$acc end parallel loop
#elif defined(SPLINE1D_BACKEND_OMP_TARGET)
!$omp end target teams distribute parallel do
#endif

