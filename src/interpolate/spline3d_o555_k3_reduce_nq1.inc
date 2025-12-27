      ! Third reduction (k3): initialize from k3 = N3
      y_batch(1) = coeff_3(N3)
      dy_batch(1, 1) = coeff_3_dx1(N3)
      dy_batch(2, 1) = coeff_3_dx2(N3)
      dy_batch(3, 1) = real(N3, dp)*coeff_3(N3)
      D2_11 = coeff_3_dx1x1(N3)
      D2_12 = coeff_3_dx1x2(N3)
      D2_13 = real(N3, dp)*coeff_3_dx1(N3)
#ifdef SPLINE3D_O555_FULL_D2
      D2_22 = coeff_3_dx2x2(N3)
      D2_23 = real(N3, dp)*coeff_3_dx2(N3)
      D2_33 = real(N3*(N3 - 1), dp)*coeff_3(N3)
#endif

      ! k3 = 4
      y_batch(1) = coeff_3(4) + x3*y_batch(1)
      dy_batch(1, 1) = coeff_3_dx1(4) + x3*dy_batch(1, 1)
      dy_batch(2, 1) = coeff_3_dx2(4) + x3*dy_batch(2, 1)
      dy_batch(3, 1) = 4.0_dp*coeff_3(4) + x3*dy_batch(3, 1)
      D2_11 = coeff_3_dx1x1(4) + x3*D2_11
      D2_12 = coeff_3_dx1x2(4) + x3*D2_12
      D2_13 = 4.0_dp*coeff_3_dx1(4) + x3*D2_13
#ifdef SPLINE3D_O555_FULL_D2
      D2_22 = coeff_3_dx2x2(4) + x3*D2_22
      D2_23 = 4.0_dp*coeff_3_dx2(4) + x3*D2_23
      D2_33 = 12.0_dp*coeff_3(4) + x3*D2_33
#endif

      ! k3 = 3
      y_batch(1) = coeff_3(3) + x3*y_batch(1)
      dy_batch(1, 1) = coeff_3_dx1(3) + x3*dy_batch(1, 1)
      dy_batch(2, 1) = coeff_3_dx2(3) + x3*dy_batch(2, 1)
      dy_batch(3, 1) = 3.0_dp*coeff_3(3) + x3*dy_batch(3, 1)
      D2_11 = coeff_3_dx1x1(3) + x3*D2_11
      D2_12 = coeff_3_dx1x2(3) + x3*D2_12
      D2_13 = 3.0_dp*coeff_3_dx1(3) + x3*D2_13
#ifdef SPLINE3D_O555_FULL_D2
      D2_22 = coeff_3_dx2x2(3) + x3*D2_22
      D2_23 = 3.0_dp*coeff_3_dx2(3) + x3*D2_23
      D2_33 = 6.0_dp*coeff_3(3) + x3*D2_33
#endif

      ! k3 = 2
      y_batch(1) = coeff_3(2) + x3*y_batch(1)
      dy_batch(1, 1) = coeff_3_dx1(2) + x3*dy_batch(1, 1)
      dy_batch(2, 1) = coeff_3_dx2(2) + x3*dy_batch(2, 1)
      dy_batch(3, 1) = 2.0_dp*coeff_3(2) + x3*dy_batch(3, 1)
      D2_11 = coeff_3_dx1x1(2) + x3*D2_11
      D2_12 = coeff_3_dx1x2(2) + x3*D2_12
      D2_13 = 2.0_dp*coeff_3_dx1(2) + x3*D2_13
#ifdef SPLINE3D_O555_FULL_D2
      D2_22 = coeff_3_dx2x2(2) + x3*D2_22
      D2_23 = 2.0_dp*coeff_3_dx2(2) + x3*D2_23
      D2_33 = 2.0_dp*coeff_3(2) + x3*D2_33
#endif

      ! k3 = 1
      y_batch(1) = coeff_3(1) + x3*y_batch(1)
      dy_batch(1, 1) = coeff_3_dx1(1) + x3*dy_batch(1, 1)
      dy_batch(2, 1) = coeff_3_dx2(1) + x3*dy_batch(2, 1)
      dy_batch(3, 1) = coeff_3(1) + x3*dy_batch(3, 1)
      D2_11 = coeff_3_dx1x1(1) + x3*D2_11
      D2_12 = coeff_3_dx1x2(1) + x3*D2_12
      D2_13 = coeff_3_dx1(1) + x3*D2_13
#ifdef SPLINE3D_O555_FULL_D2
      D2_22 = coeff_3_dx2x2(1) + x3*D2_22
      D2_23 = coeff_3_dx2(1) + x3*D2_23
#endif

      ! k3 = 0
      y_batch(1) = coeff_3(0) + x3*y_batch(1)
      dy_batch(1, 1) = coeff_3_dx1(0) + x3*dy_batch(1, 1)
      dy_batch(2, 1) = coeff_3_dx2(0) + x3*dy_batch(2, 1)
      D2_11 = coeff_3_dx1x1(0) + x3*D2_11
      D2_12 = coeff_3_dx1x2(0) + x3*D2_12
#ifdef SPLINE3D_O555_FULL_D2
      D2_22 = coeff_3_dx2x2(0) + x3*D2_22
#endif
