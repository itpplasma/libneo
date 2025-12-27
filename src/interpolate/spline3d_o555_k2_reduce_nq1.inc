      ! Second reduction (k2): initialize from k2 = N2
      do k3 = 0, N3
         coeff_3(k3) = coeff_23(N2, k3)
         coeff_3_dx1(k3) = coeff_23_dx1(N2, k3)
         coeff_3_dx1x1(k3) = coeff_23_dx1x1(N2, k3)
         coeff_3_dx2(k3) = real(N2, dp)*coeff_3(k3)
         coeff_3_dx1x2(k3) = real(N2, dp)*coeff_3_dx1(k3)
#ifdef SPLINE3D_O555_FULL_D2
         coeff_3_dx2x2(k3) = real(N2*(N2 - 1), dp)*coeff_3(k3)
#endif
      end do

      ! k2 = 4
      do k3 = 0, N3
         coeff_3(k3) = coeff_23(4, k3) + x2*coeff_3(k3)
         coeff_3_dx1(k3) = coeff_23_dx1(4, k3) + x2*coeff_3_dx1(k3)
         coeff_3_dx1x1(k3) = coeff_23_dx1x1(4, k3) + x2*coeff_3_dx1x1(k3)
         coeff_3_dx2(k3) = 4.0_dp*coeff_23(4, k3) + x2*coeff_3_dx2(k3)
         coeff_3_dx1x2(k3) = 4.0_dp*coeff_23_dx1(4, k3) + x2*coeff_3_dx1x2(k3)
#ifdef SPLINE3D_O555_FULL_D2
         coeff_3_dx2x2(k3) = 12.0_dp*coeff_23(4, k3) + x2*coeff_3_dx2x2(k3)
#endif
      end do

      ! k2 = 3
      do k3 = 0, N3
         coeff_3(k3) = coeff_23(3, k3) + x2*coeff_3(k3)
         coeff_3_dx1(k3) = coeff_23_dx1(3, k3) + x2*coeff_3_dx1(k3)
         coeff_3_dx1x1(k3) = coeff_23_dx1x1(3, k3) + x2*coeff_3_dx1x1(k3)
         coeff_3_dx2(k3) = 3.0_dp*coeff_23(3, k3) + x2*coeff_3_dx2(k3)
         coeff_3_dx1x2(k3) = 3.0_dp*coeff_23_dx1(3, k3) + x2*coeff_3_dx1x2(k3)
#ifdef SPLINE3D_O555_FULL_D2
         coeff_3_dx2x2(k3) = 6.0_dp*coeff_23(3, k3) + x2*coeff_3_dx2x2(k3)
#endif
      end do

      ! k2 = 2
      do k3 = 0, N3
         coeff_3(k3) = coeff_23(2, k3) + x2*coeff_3(k3)
         coeff_3_dx1(k3) = coeff_23_dx1(2, k3) + x2*coeff_3_dx1(k3)
         coeff_3_dx1x1(k3) = coeff_23_dx1x1(2, k3) + x2*coeff_3_dx1x1(k3)
         coeff_3_dx2(k3) = 2.0_dp*coeff_23(2, k3) + x2*coeff_3_dx2(k3)
         coeff_3_dx1x2(k3) = 2.0_dp*coeff_23_dx1(2, k3) + x2*coeff_3_dx1x2(k3)
#ifdef SPLINE3D_O555_FULL_D2
         coeff_3_dx2x2(k3) = 2.0_dp*coeff_23(2, k3) + x2*coeff_3_dx2x2(k3)
#endif
      end do

      ! k2 = 1
      do k3 = 0, N3
         coeff_3(k3) = coeff_23(1, k3) + x2*coeff_3(k3)
         coeff_3_dx1(k3) = coeff_23_dx1(1, k3) + x2*coeff_3_dx1(k3)
         coeff_3_dx1x1(k3) = coeff_23_dx1x1(1, k3) + x2*coeff_3_dx1x1(k3)
         coeff_3_dx2(k3) = coeff_23(1, k3) + x2*coeff_3_dx2(k3)
         coeff_3_dx1x2(k3) = coeff_23_dx1(1, k3) + x2*coeff_3_dx1x2(k3)
      end do

      ! k2 = 0
      do k3 = 0, N3
         coeff_3(k3) = coeff_23(0, k3) + x2*coeff_3(k3)
         coeff_3_dx1(k3) = coeff_23_dx1(0, k3) + x2*coeff_3_dx1(k3)
         coeff_3_dx1x1(k3) = coeff_23_dx1x1(0, k3) + x2*coeff_3_dx1x1(k3)
      end do
