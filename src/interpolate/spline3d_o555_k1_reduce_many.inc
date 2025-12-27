    ! First reduction (k1): o555 specialization, NQ >= 1 (vectorized over iq)
    ! Computes value, dx1, dx1x1 polynomials over x1.
    !
    ! Required symbols:
    ! - spl, x1, i1,i2,i3, N2,N3, NQ
    ! - coeff_23, coeff_23_dx1, coeff_23_dx1x1, c
    ! - k2,k3,iq loop indices

    ! initialize from k1 = 5
    do k3 = 0, N3
        do k2 = 0, N2
    !$omp simd
            do iq = 1, NQ
                c = spl%coeff(iq, 5, k2, k3, i1, i2, i3)
                coeff_23(iq, k2, k3) = c
                coeff_23_dx1(iq, k2, k3) = 5.0_dp*c
                coeff_23_dx1x1(iq, k2, k3) = 20.0_dp*c
            end do
        end do
    end do

    ! k1 = 4
    do k3 = 0, N3
        do k2 = 0, N2
    !$omp simd
            do iq = 1, NQ
                c = spl%coeff(iq, 4, k2, k3, i1, i2, i3)
                coeff_23(iq, k2, k3) = c + x1*coeff_23(iq, k2, k3)
                coeff_23_dx1(iq, k2, k3) = 4.0_dp*c + x1*coeff_23_dx1(iq, k2, k3)
                coeff_23_dx1x1(iq, k2, k3) = 12.0_dp*c + &
                                             x1*coeff_23_dx1x1(iq, k2, k3)
            end do
        end do
    end do

    ! k1 = 3
    do k3 = 0, N3
        do k2 = 0, N2
    !$omp simd
            do iq = 1, NQ
                c = spl%coeff(iq, 3, k2, k3, i1, i2, i3)
                coeff_23(iq, k2, k3) = c + x1*coeff_23(iq, k2, k3)
                coeff_23_dx1(iq, k2, k3) = 3.0_dp*c + x1*coeff_23_dx1(iq, k2, k3)
                coeff_23_dx1x1(iq, k2, k3) = 6.0_dp*c + &
                                             x1*coeff_23_dx1x1(iq, k2, k3)
            end do
        end do
    end do

    ! k1 = 2
    do k3 = 0, N3
        do k2 = 0, N2
    !$omp simd
            do iq = 1, NQ
                c = spl%coeff(iq, 2, k2, k3, i1, i2, i3)
                coeff_23(iq, k2, k3) = c + x1*coeff_23(iq, k2, k3)
                coeff_23_dx1(iq, k2, k3) = 2.0_dp*c + x1*coeff_23_dx1(iq, k2, k3)
                coeff_23_dx1x1(iq, k2, k3) = 2.0_dp*c + &
                                             x1*coeff_23_dx1x1(iq, k2, k3)
            end do
        end do
    end do

    ! k1 = 1
    do k3 = 0, N3
        do k2 = 0, N2
    !$omp simd
            do iq = 1, NQ
                c = spl%coeff(iq, 1, k2, k3, i1, i2, i3)
                coeff_23(iq, k2, k3) = c + x1*coeff_23(iq, k2, k3)
                coeff_23_dx1(iq, k2, k3) = c + x1*coeff_23_dx1(iq, k2, k3)
            end do
        end do
    end do

    ! k1 = 0
    do k3 = 0, N3
        do k2 = 0, N2
    !$omp simd
            do iq = 1, NQ
                c = spl%coeff(iq, 0, k2, k3, i1, i2, i3)
                coeff_23(iq, k2, k3) = c + x1*coeff_23(iq, k2, k3)
            end do
        end do
    end do
