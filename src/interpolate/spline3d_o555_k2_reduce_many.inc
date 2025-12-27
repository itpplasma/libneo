      ! Second reduction (k2): initialize from k2 = N2
      do k3 = 0, N3
!$omp simd
         do iq = 1, NQ
            coeff_3(iq, k3) = coeff_23(iq, N2, k3)
            coeff_3_dx1(iq, k3) = coeff_23_dx1(iq, N2, k3)
            coeff_3_dx1x1(iq, k3) = coeff_23_dx1x1(iq, N2, k3)
            coeff_3_dx2(iq, k3) = real(N2, dp)*coeff_3(iq, k3)
            coeff_3_dx1x2(iq, k3) = real(N2, dp)*coeff_3_dx1(iq, k3)
#ifdef SPLINE3D_O555_FULL_D2
            coeff_3_dx2x2(iq, k3) = real(N2*(N2 - 1), dp)*coeff_3(iq, k3)
#endif
         end do
      end do

      ! k2 = 4
      do k3 = 0, N3
!$omp simd
         do iq = 1, NQ
            coeff_3(iq, k3) = coeff_23(iq, 4, k3) + x2*coeff_3(iq, k3)
            coeff_3_dx1(iq, k3) = coeff_23_dx1(iq, 4, k3) + x2*coeff_3_dx1(iq, k3)
            coeff_3_dx1x1(iq, k3) = coeff_23_dx1x1(iq, 4, k3) + &
                                    x2*coeff_3_dx1x1(iq, k3)
            coeff_3_dx2(iq, k3) = 4.0_dp*coeff_23(iq, 4, k3) + &
                                  x2*coeff_3_dx2(iq, k3)
            coeff_3_dx1x2(iq, k3) = 4.0_dp*coeff_23_dx1(iq, 4, k3) + &
                                    x2*coeff_3_dx1x2(iq, k3)
#ifdef SPLINE3D_O555_FULL_D2
            coeff_3_dx2x2(iq, k3) = 12.0_dp*coeff_23(iq, 4, k3) + &
                                    x2*coeff_3_dx2x2(iq, k3)
#endif
         end do
      end do

      ! k2 = 3
      do k3 = 0, N3
!$omp simd
         do iq = 1, NQ
            coeff_3(iq, k3) = coeff_23(iq, 3, k3) + x2*coeff_3(iq, k3)
            coeff_3_dx1(iq, k3) = coeff_23_dx1(iq, 3, k3) + x2*coeff_3_dx1(iq, k3)
            coeff_3_dx1x1(iq, k3) = coeff_23_dx1x1(iq, 3, k3) + &
                                    x2*coeff_3_dx1x1(iq, k3)
            coeff_3_dx2(iq, k3) = 3.0_dp*coeff_23(iq, 3, k3) + &
                                  x2*coeff_3_dx2(iq, k3)
            coeff_3_dx1x2(iq, k3) = 3.0_dp*coeff_23_dx1(iq, 3, k3) + &
                                    x2*coeff_3_dx1x2(iq, k3)
#ifdef SPLINE3D_O555_FULL_D2
            coeff_3_dx2x2(iq, k3) = 6.0_dp*coeff_23(iq, 3, k3) + &
                                    x2*coeff_3_dx2x2(iq, k3)
#endif
         end do
      end do

      ! k2 = 2
      do k3 = 0, N3
!$omp simd
         do iq = 1, NQ
            coeff_3(iq, k3) = coeff_23(iq, 2, k3) + x2*coeff_3(iq, k3)
            coeff_3_dx1(iq, k3) = coeff_23_dx1(iq, 2, k3) + x2*coeff_3_dx1(iq, k3)
            coeff_3_dx1x1(iq, k3) = coeff_23_dx1x1(iq, 2, k3) + &
                                    x2*coeff_3_dx1x1(iq, k3)
            coeff_3_dx2(iq, k3) = 2.0_dp*coeff_23(iq, 2, k3) + &
                                  x2*coeff_3_dx2(iq, k3)
            coeff_3_dx1x2(iq, k3) = 2.0_dp*coeff_23_dx1(iq, 2, k3) + &
                                    x2*coeff_3_dx1x2(iq, k3)
#ifdef SPLINE3D_O555_FULL_D2
            coeff_3_dx2x2(iq, k3) = 2.0_dp*coeff_23(iq, 2, k3) + &
                                    x2*coeff_3_dx2x2(iq, k3)
#endif
         end do
      end do

      ! k2 = 1
      do k3 = 0, N3
!$omp simd
         do iq = 1, NQ
            coeff_3(iq, k3) = coeff_23(iq, 1, k3) + x2*coeff_3(iq, k3)
            coeff_3_dx1(iq, k3) = coeff_23_dx1(iq, 1, k3) + x2*coeff_3_dx1(iq, k3)
            coeff_3_dx1x1(iq, k3) = coeff_23_dx1x1(iq, 1, k3) + &
                                    x2*coeff_3_dx1x1(iq, k3)
            coeff_3_dx2(iq, k3) = coeff_23(iq, 1, k3) + x2*coeff_3_dx2(iq, k3)
            coeff_3_dx1x2(iq, k3) = coeff_23_dx1(iq, 1, k3) + &
                                    x2*coeff_3_dx1x2(iq, k3)
         end do
      end do

      ! k2 = 0
      do k3 = 0, N3
!$omp simd
         do iq = 1, NQ
            coeff_3(iq, k3) = coeff_23(iq, 0, k3) + x2*coeff_3(iq, k3)
            coeff_3_dx1(iq, k3) = coeff_23_dx1(iq, 0, k3) + x2*coeff_3_dx1(iq, k3)
            coeff_3_dx1x1(iq, k3) = coeff_23_dx1x1(iq, 0, k3) + &
                                    x2*coeff_3_dx1x1(iq, k3)
         end do
      end do
