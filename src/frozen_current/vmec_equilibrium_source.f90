!> VMEC equilibrium source of the frozen-current field.
!>
!> Uses libneo's VMEC reader (spline_vmec_sub, global state set by
!> spline_vmec_data), in its units: cm and G.  B is VMEC's field
!> B = grad A_theta x grad vartheta + grad A_phi x grad varphi with
!> vartheta = theta + lambda; its components and the Jacobian come from the
!> generated leaf vmec_source_field.
module neo_vmec_equilibrium_source
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neo_equilibrium_source, only: equilibrium_source_t, source_points_t, &
                                      source_resolution_t, flux_quadrature
    use neo_frozen_current_leaves, only: vmec_source_field, vmec_surface_normal
    implicit none
    private

    public :: vmec_equilibrium_source_t, vmec_equilibrium_source_init

    type, extends(equilibrium_source_t) :: vmec_equilibrium_source_t
        real(dp) :: orientation = 1.0_dp !! sign of the Jacobian
    contains
        procedure :: sample_volume => vmec_sample_volume
        procedure :: sample_surface => vmec_sample_surface
        procedure :: boundary_curve => vmec_boundary_curve
    end type vmec_equilibrium_source_t

contains

    !> Read and spline the wout file (netCDF) through libneo's VMEC reader.
    subroutine vmec_equilibrium_source_init(self, wout_file)
        use new_vmec_stuff_mod, only: netcdffile, nper
        use spline_vmec_sub, only: spline_vmec_data
        type(vmec_equilibrium_source_t), intent(out) :: self
        character(*), intent(in) :: wout_file
        real(dp) :: sqg, b(3)

        netcdffile = wout_file
        call spline_vmec_data
        self%nfp = nper
        call vmec_point(0.5_dp, 0.3_dp, 0.1_dp, sqg, b)
        self%orientation = sign(1.0_dp, sqg)
    end subroutine vmec_equilibrium_source_init

    subroutine vmec_sample_volume(self, res, pts)
        class(vmec_equilibrium_source_t), intent(in) :: self
        type(source_resolution_t), intent(in) :: res
        type(source_points_t), intent(out) :: pts

        real(dp), allocatable :: s(:), ws(:), theta(:), zeta(:)
        real(dp) :: dtheta, dzeta, r, z, sqg, b(3)
        integer :: i, j, k, m, n

        call flux_quadrature(res, self%nfp, s, ws, theta, zeta, dtheta, dzeta)
        n = res%ns*res%ntheta*res%nzeta
        pts%nfp = self%nfp
        allocate (pts%x(3, n), pts%b(3, n), pts%sqg(n), pts%w(n))
        !$omp parallel do collapse(3) default(shared) private(i, j, k, m, r, z, sqg, b)
        do k = 1, res%nzeta
            do j = 1, res%ntheta
                do i = 1, res%ns
                    m = i + res%ns*((j - 1) + res%ntheta*(k - 1))
                    call vmec_point(s(i), theta(j), zeta(k), sqg, b, r, z)
                    pts%x(:, m) = [r, zeta(k), z]
                    pts%b(:, m) = b
                    pts%sqg(m) = sqg
                    pts%w(m) = abs(sqg)*ws(i)*dtheta*dzeta
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine vmec_sample_volume

    subroutine vmec_sample_surface(self, res, pts)
        class(vmec_equilibrium_source_t), intent(in) :: self
        type(source_resolution_t), intent(in) :: res
        type(source_points_t), intent(out) :: pts

        real(dp), allocatable :: s(:), ws(:), theta(:), zeta(:)
        real(dp) :: dtheta, dzeta, r, z, sqg, b(3), a(3)
        integer :: j, k, m, n

        call flux_quadrature(res, self%nfp, s, ws, theta, zeta, dtheta, dzeta)
        n = res%ntheta*res%nzeta
        pts%nfp = self%nfp
        allocate (pts%x(3, n), pts%b(3, n), pts%sqg(n), pts%w(n), pts%area(3, n))
        !$omp parallel do collapse(2) default(shared) private(j, k, m, r, z, sqg, b, a)
        do k = 1, res%nzeta
            do j = 1, res%ntheta
                m = j + res%ntheta*(k - 1)
                call vmec_point(1.0_dp, theta(j), zeta(k), sqg, b, r, z, a)
                pts%x(:, m) = [r, zeta(k), z]
                pts%b(:, m) = b
                pts%sqg(m) = sqg
                pts%w(m) = dtheta*dzeta
                pts%area(:, m) = self%orientation*a*dtheta*dzeta
            end do
        end do
        !$omp end parallel do
    end subroutine vmec_sample_surface

    subroutine vmec_boundary_curve(self, phi, theta, r, z)
        class(vmec_equilibrium_source_t), intent(in) :: self
        real(dp), intent(in) :: phi, theta(:)
        real(dp), intent(out) :: r(:), z(:)
        real(dp) :: sqg, b(3)
        integer :: j

        do j = 1, size(theta)
            call vmec_point(1.0_dp, theta(j), phi, sqg, b, r(j), z(j))
        end do
    end subroutine vmec_boundary_curve

    !> Field (B_R, B_phi, B_Z), Jacobian, position and the area vector
    !> e_theta x e_zeta of the surface s = const at one VMEC point.
    subroutine vmec_point(s, theta, zeta, sqg, b, r, z, area)
        use spline_vmec_sub, only: splint_vmec_data
        real(dp), intent(in) :: s, theta, zeta
        real(dp), intent(out) :: sqg, b(3)
        real(dp), intent(out), optional :: r, z, area(3)

        real(dp) :: a_phi, a_theta, da_phi, da_theta, iota, rr, zz, lam
        real(dp) :: rs, rt, rp, zs, zt, zp, ls, lt, lp

        call splint_vmec_data(s, theta, zeta, a_phi, a_theta, da_phi, da_theta, &
                              iota, rr, zz, lam, rs, rt, rp, zs, zt, zp, ls, lt, lp)
        call vmec_source_field(rr, rs, rt, rp, zs, zt, zp, lt, lp, da_theta, &
                               da_phi, sqg, b(1), b(2), b(3))
        if (present(r)) r = rr
        if (present(z)) z = zz
        if (present(area)) call vmec_surface_normal(rr, rt, rp, zt, zp, area(1), &
                                                    area(2), area(3))
    end subroutine vmec_point

end module neo_vmec_equilibrium_source
