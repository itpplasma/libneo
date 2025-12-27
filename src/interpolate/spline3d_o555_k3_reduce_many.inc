    ! Third reduction (k3): initialize from k3 = N3
    !$omp simd
    do iq = 1, NQ
        y_batch(iq) = coeff_3(iq, N3)
        dy_batch(1, iq) = coeff_3_dx1(iq, N3)
        dy_batch(2, iq) = coeff_3_dx2(iq, N3)
        dy_batch(3, iq) = real(N3, dp)*coeff_3(iq, N3)
        D2_11(iq) = coeff_3_dx1x1(iq, N3)
        D2_12(iq) = coeff_3_dx1x2(iq, N3)
        D2_13(iq) = real(N3, dp)*coeff_3_dx1(iq, N3)
#ifdef SPLINE3D_O555_FULL_D2
        D2_22(iq) = coeff_3_dx2x2(iq, N3)
        D2_23(iq) = real(N3, dp)*coeff_3_dx2(iq, N3)
        D2_33(iq) = real(N3*(N3 - 1), dp)*coeff_3(iq, N3)
#endif
    end do

    ! k3 = 4
    !$omp simd
    do iq = 1, NQ
        y_batch(iq) = coeff_3(iq, 4) + x3*y_batch(iq)
        dy_batch(1, iq) = coeff_3_dx1(iq, 4) + x3*dy_batch(1, iq)
        dy_batch(2, iq) = coeff_3_dx2(iq, 4) + x3*dy_batch(2, iq)
        dy_batch(3, iq) = 4.0_dp*coeff_3(iq, 4) + x3*dy_batch(3, iq)
        D2_11(iq) = coeff_3_dx1x1(iq, 4) + x3*D2_11(iq)
        D2_12(iq) = coeff_3_dx1x2(iq, 4) + x3*D2_12(iq)
        D2_13(iq) = 4.0_dp*coeff_3_dx1(iq, 4) + x3*D2_13(iq)
#ifdef SPLINE3D_O555_FULL_D2
        D2_22(iq) = coeff_3_dx2x2(iq, 4) + x3*D2_22(iq)
        D2_23(iq) = 4.0_dp*coeff_3_dx2(iq, 4) + x3*D2_23(iq)
        D2_33(iq) = 12.0_dp*coeff_3(iq, 4) + x3*D2_33(iq)
#endif
    end do

    ! k3 = 3
    !$omp simd
    do iq = 1, NQ
        y_batch(iq) = coeff_3(iq, 3) + x3*y_batch(iq)
        dy_batch(1, iq) = coeff_3_dx1(iq, 3) + x3*dy_batch(1, iq)
        dy_batch(2, iq) = coeff_3_dx2(iq, 3) + x3*dy_batch(2, iq)
        dy_batch(3, iq) = 3.0_dp*coeff_3(iq, 3) + x3*dy_batch(3, iq)
        D2_11(iq) = coeff_3_dx1x1(iq, 3) + x3*D2_11(iq)
        D2_12(iq) = coeff_3_dx1x2(iq, 3) + x3*D2_12(iq)
        D2_13(iq) = 3.0_dp*coeff_3_dx1(iq, 3) + x3*D2_13(iq)
#ifdef SPLINE3D_O555_FULL_D2
        D2_22(iq) = coeff_3_dx2x2(iq, 3) + x3*D2_22(iq)
        D2_23(iq) = 3.0_dp*coeff_3_dx2(iq, 3) + x3*D2_23(iq)
        D2_33(iq) = 6.0_dp*coeff_3(iq, 3) + x3*D2_33(iq)
#endif
    end do

    ! k3 = 2
    !$omp simd
    do iq = 1, NQ
        y_batch(iq) = coeff_3(iq, 2) + x3*y_batch(iq)
        dy_batch(1, iq) = coeff_3_dx1(iq, 2) + x3*dy_batch(1, iq)
        dy_batch(2, iq) = coeff_3_dx2(iq, 2) + x3*dy_batch(2, iq)
        dy_batch(3, iq) = 2.0_dp*coeff_3(iq, 2) + x3*dy_batch(3, iq)
        D2_11(iq) = coeff_3_dx1x1(iq, 2) + x3*D2_11(iq)
        D2_12(iq) = coeff_3_dx1x2(iq, 2) + x3*D2_12(iq)
        D2_13(iq) = 2.0_dp*coeff_3_dx1(iq, 2) + x3*D2_13(iq)
#ifdef SPLINE3D_O555_FULL_D2
        D2_22(iq) = coeff_3_dx2x2(iq, 2) + x3*D2_22(iq)
        D2_23(iq) = 2.0_dp*coeff_3_dx2(iq, 2) + x3*D2_23(iq)
        D2_33(iq) = 2.0_dp*coeff_3(iq, 2) + x3*D2_33(iq)
#endif
    end do

    ! k3 = 1
    !$omp simd
    do iq = 1, NQ
        y_batch(iq) = coeff_3(iq, 1) + x3*y_batch(iq)
        dy_batch(1, iq) = coeff_3_dx1(iq, 1) + x3*dy_batch(1, iq)
        dy_batch(2, iq) = coeff_3_dx2(iq, 1) + x3*dy_batch(2, iq)
        dy_batch(3, iq) = coeff_3(iq, 1) + x3*dy_batch(3, iq)
        D2_11(iq) = coeff_3_dx1x1(iq, 1) + x3*D2_11(iq)
        D2_12(iq) = coeff_3_dx1x2(iq, 1) + x3*D2_12(iq)
        D2_13(iq) = coeff_3_dx1(iq, 1) + x3*D2_13(iq)
#ifdef SPLINE3D_O555_FULL_D2
        D2_22(iq) = coeff_3_dx2x2(iq, 1) + x3*D2_22(iq)
        D2_23(iq) = coeff_3_dx2(iq, 1) + x3*D2_23(iq)
#endif
    end do

    ! k3 = 0
    !$omp simd
    do iq = 1, NQ
        y_batch(iq) = coeff_3(iq, 0) + x3*y_batch(iq)
        dy_batch(1, iq) = coeff_3_dx1(iq, 0) + x3*dy_batch(1, iq)
        dy_batch(2, iq) = coeff_3_dx2(iq, 0) + x3*dy_batch(2, iq)
        D2_11(iq) = coeff_3_dx1x1(iq, 0) + x3*D2_11(iq)
        D2_12(iq) = coeff_3_dx1x2(iq, 0) + x3*D2_12(iq)
#ifdef SPLINE3D_O555_FULL_D2
        D2_22(iq) = coeff_3_dx2x2(iq, 0) + x3*D2_22(iq)
#endif
    end do
