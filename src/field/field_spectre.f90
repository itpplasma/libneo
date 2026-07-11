module field_spectre
    !> SPECTRE/SPEC magnetic field evaluation in the stacked-rho chart (#370).
    !>
    !> x = (rho_g, theta, zeta). Per point: lvol = min(int(rho_g) + 1, Mvol),
    !> local s = 2*(rho_g - lvol + 1) - 1 in [-1, 1], so ds/drho_g = 2 is the
    !> chain factor on every radial derivative. spectre_basis returns A and its
    !> derivatives in LOCAL s; the chart rescales only the radial coordinate, so
    !> covariant angle components A_theta, A_zeta and the poloidal/toroidal
    !> contravariant field are chart-invariant, while B^rho picks up 1/2.
    !>
    !> In the A_s = 0 gauge, the metric-independent 2-form F = dA gives
    !>     sqrtg B^rho   = d_theta A_zeta - d_zeta A_theta   (= F_tz, chart-free)
    !>     sqrtg B^theta = -d_rho  A_zeta  = -2 dAzt/ds
    !>     sqrtg B^zeta  =  d_rho  A_theta =  2 dAth/ds
    !> with sqrtg the chart Jacobian (2*|jac_local|, see #370). Hence B^theta and
    !> B^zeta match the local-s field, and B^rho equals the local-s B^s halved.
    !> Bmod**2 = B^i B^j g_ij; hcov_i = g_ij B^j / Bmod; Acov = (0, A_th, A_zt).
    !> SI units (Tesla, meter); consumers convert.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_base, only: magnetic_field_t
    use spectre_reader, only: spectre_data_t, load_spectre
    use spectre_basis, only: spectre_vecpot_t, eval_spectre_vector_potential
    use libneo_coordinates, only: make_spectre_coordinate_system, &
                                  spectre_coordinate_system_t

    implicit none
    private

    public :: spectre_field_t, create_spectre_field

    real(dp), parameter :: ds_drho = 2.0_dp

    type, extends(magnetic_field_t) :: spectre_field_t
        type(spectre_data_t) :: data
    contains
        procedure :: evaluate => spectre_evaluate
        procedure :: evaluate_der => spectre_evaluate_der
    end type spectre_field_t

contains

    subroutine create_spectre_field(field, filename, ierr)
        !> Load a SPECTRE equilibrium and build its stacked-rho coordinate chart.
        type(spectre_field_t), intent(out) :: field
        character(len=*), intent(in) :: filename
        integer, intent(out) :: ierr

        call load_spectre(filename, field%data, ierr)
        if (ierr /= 0) return
        call make_spectre_coordinate_system(field%coords, field%data)
    end subroutine create_spectre_field

    subroutine spectre_evaluate(self, x, Acov, hcov, Bmod, sqgBctr)
        class(spectre_field_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: Acov(3), hcov(3), Bmod
        real(dp), intent(out), optional :: sqgBctr(3)

        type(spectre_vecpot_t) :: av
        integer :: lvol
        real(dp) :: s, g(3, 3), ginv(3, 3), sqrtg
        real(dp) :: Bctr(3), Bcov(3), sqgB(3)

        call locate_volume(self%data%Mvol, x(1), lvol, s)
        call eval_spectre_vector_potential(self%data, lvol, s, x(2), x(3), av)
        call self%coords%metric_tensor(x, g, ginv, sqrtg)

        sqgB(1) = av%dAzt(2) - av%dAth(3)
        sqgB(2) = -ds_drho*av%dAzt(1)
        sqgB(3) = ds_drho*av%dAth(1)
        Bctr = sqgB/sqrtg

        Bcov = matmul(g, Bctr)
        Bmod = sqrt(dot_product(Bctr, Bcov))

        Acov = [0.0_dp, av%Ath, av%Azt]
        hcov = Bcov/Bmod

        if (present(sqgBctr)) sqgBctr = sqgB
    end subroutine spectre_evaluate

    subroutine spectre_evaluate_der(self, x, hcov, Bmod, sqgBctr, sqrtg, dhcov, dBmod)
        !> Field quantities and the analytic x-derivatives of covariant h and |B|:
        !> dhcov(i,j) = d hcov_i/dx_j, dBmod(j) = d|B|/dx_j (SI). sqgBctr = sqrtg B^i
        !> and the chart Jacobian sqrtg are returned so a drift caller needs no
        !> separate metric_tensor call. The 2-form sqrtg B^i is metric-free; |B| and
        !> covariant h carry the metric, so their derivatives combine the vector-
        !> potential Hessian with the analytic metric derivative dg.
        class(spectre_field_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: hcov(3), Bmod, sqgBctr(3), sqrtg
        real(dp), intent(out) :: dhcov(3, 3), dBmod(3)

        type(spectre_vecpot_t) :: av
        integer :: lvol, i, j, k
        real(dp) :: s
        real(dp) :: g(3, 3), ginv(3, 3), dg(3, 3, 3), dsqrtg(3)
        real(dp) :: sqgB(3), dsqgB(3, 3), Bctr(3), dBctr(3, 3)
        real(dp) :: Bcov(3), dBcov(3, 3), B2, dB2(3)

        call locate_volume(self%data%Mvol, x(1), lvol, s)
        call eval_spectre_vector_potential(self%data, lvol, s, x(2), x(3), av)

        select type (cs => self%coords)
        type is (spectre_coordinate_system_t)
            call cs%metric_tensor_der(x, g, ginv, sqrtg, dg, dsqrtg)
        class default
            error stop 'spectre_evaluate_der: coords not spectre_coordinate_system_t'
        end select

        ! Metric-free 2-form sqrtg B^i (see spectre_evaluate); d_rho = ds_drho d_s.
        sqgB(1) = av%dAzt(2) - av%dAth(3)
        sqgB(2) = -ds_drho*av%dAzt(1)
        sqgB(3) = ds_drho*av%dAth(1)

        ! d sqgB(i)/dx_j from the A Hessian (slots ss,st,sz,tt,tz,zz); the radial
        ! index carries ds_drho, so radial-radial entries carry ds_drho**2.
        dsqgB(1, 1) = ds_drho*(av%d2Azt(2) - av%d2Ath(3))
        dsqgB(1, 2) = av%d2Azt(4) - av%d2Ath(5)
        dsqgB(1, 3) = av%d2Azt(5) - av%d2Ath(6)
        dsqgB(2, 1) = -ds_drho*ds_drho*av%d2Azt(1)
        dsqgB(2, 2) = -ds_drho*av%d2Azt(2)
        dsqgB(2, 3) = -ds_drho*av%d2Azt(3)
        dsqgB(3, 1) = ds_drho*ds_drho*av%d2Ath(1)
        dsqgB(3, 2) = ds_drho*av%d2Ath(2)
        dsqgB(3, 3) = ds_drho*av%d2Ath(3)

        Bctr = sqgB/sqrtg
        do j = 1, 3
            dBctr(:, j) = (dsqgB(:, j) - Bctr*dsqrtg(j))/sqrtg
        end do

        Bcov = matmul(g, Bctr)
        do j = 1, 3
            do i = 1, 3
                dBcov(i, j) = 0.0_dp
                do k = 1, 3
                    dBcov(i, j) = dBcov(i, j) + dg(i, k, j)*Bctr(k) &
                                  + g(i, k)*dBctr(k, j)
                end do
            end do
        end do

        B2 = dot_product(Bctr, Bcov)
        Bmod = sqrt(B2)
        do j = 1, 3
            dB2(j) = dot_product(dBctr(:, j), Bcov) + dot_product(Bctr, dBcov(:, j))
        end do
        dBmod = dB2/(2.0_dp*Bmod)

        hcov = Bcov/Bmod
        do j = 1, 3
            dhcov(:, j) = (dBcov(:, j) - hcov*dBmod(j))/Bmod
        end do

        sqgBctr = sqgB
    end subroutine spectre_evaluate_der

    pure subroutine locate_volume(Mvol, rho_g, lvol, s)
        !> Stacked chart: rho_g in [0, Mvol] -> volume index and local s.
        integer, intent(in) :: Mvol
        real(dp), intent(in) :: rho_g
        integer, intent(out) :: lvol
        real(dp), intent(out) :: s

        lvol = min(int(rho_g) + 1, Mvol)
        lvol = max(lvol, 1)
        s = ds_drho*(rho_g - real(lvol - 1, dp)) - 1.0_dp
    end subroutine locate_volume

end module field_spectre
