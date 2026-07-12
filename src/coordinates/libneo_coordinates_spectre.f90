module libneo_coordinates_spectre
    !> SPECTRE/SPEC stacked-volume coordinate chart.
    !>
    !> One chart covers all Mvol volumes through a stacked radial coordinate
    !> rho_g in [0, Mvol]. Per point: lvol = min(int(rho_g) + 1, Mvol),
    !> local s = 2*(rho_g - lvol + 1) - 1 in [-1, 1], so ds/drho_g = 2 enters
    !> every radial derivative as a chain factor.
    !>
    !> Geometry is the SPEC interface Fourier blend (get_position_vetor_x in
    !> SPECTRE coordinate_mod.F90): interior volumes blend interfaces lvol-1 and
    !> lvol linearly in s; the axis volume (lvol = 1) uses the sbar**m power law
    !> (sbar**2 for m = 0, Igeometry = 3) with sbar = (s + 1)/2. Cylindrical
    !> embedding is phi = +zeta: x = (R cos zeta, R sin zeta, Z), matching
    !> ca/01_coordinates.wl and the reference Jacobian generator.
    !>
    !> Jacobian sign/scale relation to the SPECTRE reference `jac` column
    !> (get_jacobian, evaluated at local s): the reference is the local-s
    !> Jacobian jac = R (R_t Z_s - R_s Z_t) = det(dx/d(s,theta,zeta)). This
    !> chart's covariant basis is dx/d(rho_g,theta,zeta), so
    !>     det(e_cov) = (ds/drho_g) * jac = 2 * jac,
    !>     sqrtg = |det(e_cov)| = 2 * |jac|  (no sign flip: same orientation).
    !> The factor 2 is the deliberate chart convention, not an error; it is what
    !> keeps sqrtg consistent with finite differences of evaluate_cart, which is
    !> parametrized by rho_g. Geometry is C0 across interfaces; metric radial
    !> derivatives jump there and that is physical.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_coordinates_base, only: coordinate_system_t
    use spectre_reader, only: spectre_data_t
    use cylindrical_cartesian, only: cyl_to_cart
    use math_constants, only: TWOPI
    implicit none
    private

    public :: spectre_coordinate_system_t, make_spectre_coordinate_system

    type, extends(coordinate_system_t) :: spectre_coordinate_system_t
        integer :: Mvol = 0
        integer :: mn = 0
        integer, allocatable :: im(:), in(:)
        real(dp), allocatable :: Rbc(:, :), Rbs(:, :), Zbc(:, :), Zbs(:, :)
    contains
        procedure :: evaluate_cart => spectre_evaluate_cart
        procedure :: evaluate_cyl => spectre_evaluate_cyl
        procedure :: covariant_basis => spectre_covariant_basis
        procedure :: metric_tensor => spectre_metric_tensor
        procedure :: metric_tensor_der => spectre_metric_tensor_der
        procedure :: from_cyl => spectre_from_cyl
        procedure :: from_cyl_warm => spectre_from_cyl_warm
    end type spectre_coordinate_system_t

contains

    subroutine make_spectre_coordinate_system(cs, data)
        class(coordinate_system_t), allocatable, intent(out) :: cs
        type(spectre_data_t), intent(in) :: data

        type(spectre_coordinate_system_t) :: tmp

        tmp%Mvol = data%Mvol
        tmp%mn = data%mn
        allocate (tmp%im(data%mn), tmp%in(data%mn))
        allocate (tmp%Rbc(data%mn, 0:data%Mvol))
        allocate (tmp%Rbs(data%mn, 0:data%Mvol))
        allocate (tmp%Zbc(data%mn, 0:data%Mvol))
        allocate (tmp%Zbs(data%mn, 0:data%Mvol))
        tmp%im = data%im
        tmp%in = data%in
        tmp%Rbc = data%Rbc
        tmp%Rbs = data%Rbs
        tmp%Zbc = data%Zbc
        tmp%Zbs = data%Zbs
        allocate (cs, source=tmp)
    end subroutine make_spectre_coordinate_system

    pure subroutine spectre_locate(Mvol, rho_g, lvol, s)
        integer, intent(in) :: Mvol
        real(dp), intent(in) :: rho_g
        integer, intent(out) :: lvol
        real(dp), intent(out) :: s

        lvol = min(int(rho_g) + 1, Mvol)
        lvol = max(lvol, 1)
        s = 2.0_dp*(rho_g - real(lvol - 1, dp)) - 1.0_dp
    end subroutine spectre_locate

    pure subroutine spectre_geom(self, lvol, s, theta, zeta, R, Z, &
                                 dRds, dRdt, dRdz, dZds, dZdt, dZdz)
        !> R, Z and their derivatives with respect to local (s, theta, zeta).
        class(spectre_coordinate_system_t), intent(in) :: self
        integer, intent(in) :: lvol
        real(dp), intent(in) :: s, theta, zeta
        real(dp), intent(out) :: R, Z, dRds, dRdt, dRdz, dZds, dZdt, dZdz

        integer :: ii, m, n
        real(dp) :: arg, c, sn, sbar, fj, dfj, alss, blss
        real(dp) :: rc, rs, zc, zs, drc, drs, dzc, dzs

        R = 0.0_dp; Z = 0.0_dp
        dRds = 0.0_dp; dRdt = 0.0_dp; dRdz = 0.0_dp
        dZds = 0.0_dp; dZdt = 0.0_dp; dZdz = 0.0_dp

        do ii = 1, self%mn
            m = self%im(ii)
            n = self%in(ii)
            arg = m*theta - n*zeta
            c = cos(arg)
            sn = sin(arg)

            if (lvol == 1) then
                sbar = 0.5_dp*(s + 1.0_dp)
                if (m == 0) then
                    fj = sbar*sbar
                    dfj = sbar
                else
                    fj = sbar**m
                    dfj = 0.5_dp*real(m, dp)*sbar**(m - 1)
                end if
                rc = self%Rbc(ii, 0) + (self%Rbc(ii, 1) - self%Rbc(ii, 0))*fj
                rs = self%Rbs(ii, 0) + (self%Rbs(ii, 1) - self%Rbs(ii, 0))*fj
                zc = self%Zbc(ii, 0) + (self%Zbc(ii, 1) - self%Zbc(ii, 0))*fj
                zs = self%Zbs(ii, 0) + (self%Zbs(ii, 1) - self%Zbs(ii, 0))*fj
                drc = (self%Rbc(ii, 1) - self%Rbc(ii, 0))*dfj
                drs = (self%Rbs(ii, 1) - self%Rbs(ii, 0))*dfj
                dzc = (self%Zbc(ii, 1) - self%Zbc(ii, 0))*dfj
                dzs = (self%Zbs(ii, 1) - self%Zbs(ii, 0))*dfj
            else
                alss = 0.5_dp*(1.0_dp - s)
                blss = 0.5_dp*(1.0_dp + s)
                rc = alss*self%Rbc(ii, lvol - 1) + blss*self%Rbc(ii, lvol)
                rs = alss*self%Rbs(ii, lvol - 1) + blss*self%Rbs(ii, lvol)
                zc = alss*self%Zbc(ii, lvol - 1) + blss*self%Zbc(ii, lvol)
                zs = alss*self%Zbs(ii, lvol - 1) + blss*self%Zbs(ii, lvol)
                drc = 0.5_dp*(self%Rbc(ii, lvol) - self%Rbc(ii, lvol - 1))
                drs = 0.5_dp*(self%Rbs(ii, lvol) - self%Rbs(ii, lvol - 1))
                dzc = 0.5_dp*(self%Zbc(ii, lvol) - self%Zbc(ii, lvol - 1))
                dzs = 0.5_dp*(self%Zbs(ii, lvol) - self%Zbs(ii, lvol - 1))
            end if

            R = R + rc*c + rs*sn
            Z = Z + zc*c + zs*sn
            dRds = dRds + drc*c + drs*sn
            dZds = dZds + dzc*c + dzs*sn
            dRdt = dRdt + real(m, dp)*(-rc*sn + rs*c)
            dZdt = dZdt + real(m, dp)*(-zc*sn + zs*c)
            dRdz = dRdz + real(n, dp)*(rc*sn - rs*c)
            dZdz = dZdz + real(n, dp)*(zc*sn - zs*c)
        end do
    end subroutine spectre_geom

    pure subroutine spectre_geom_d2(self, lvol, s, theta, zeta, R, Z, &
                                    dRds, dRdt, dRdz, dZds, dZdt, dZdz, d2R, d2Z)
        !> R, Z, their first derivatives and their second derivatives with respect
        !> to local (s, theta, zeta), second-derivative slot order (ss, st, sz, tt,
        !> tz, zz). Kept separate from spectre_geom so the existing first-derivative
        !> geometry path stays byte-for-byte unchanged. The radial blend is linear
        !> in s in interior volumes (zero second radial derivative); the axis volume
        !> uses the sbar**m power law.
        class(spectre_coordinate_system_t), intent(in) :: self
        integer, intent(in) :: lvol
        real(dp), intent(in) :: s, theta, zeta
        real(dp), intent(out) :: R, Z, dRds, dRdt, dRdz, dZds, dZdt, dZdz
        real(dp), intent(out) :: d2R(6), d2Z(6)

        integer :: ii, m, n
        real(dp) :: arg, c, sn, sbar, fj, dfj, d2fj, alss, blss, rm, rn
        real(dp) :: rc, rs, zc, zs, drc, drs, dzc, dzs
        real(dp) :: d2rc, d2rs, d2zc, d2zs, valr, valz

        R = 0.0_dp; Z = 0.0_dp
        dRds = 0.0_dp; dRdt = 0.0_dp; dRdz = 0.0_dp
        dZds = 0.0_dp; dZdt = 0.0_dp; dZdz = 0.0_dp
        d2R = 0.0_dp; d2Z = 0.0_dp
        d2rc = 0.0_dp; d2rs = 0.0_dp; d2zc = 0.0_dp; d2zs = 0.0_dp

        do ii = 1, self%mn
            m = self%im(ii)
            n = self%in(ii)
            rm = real(m, dp)
            rn = real(n, dp)
            arg = m*theta - n*zeta
            c = cos(arg)
            sn = sin(arg)

            if (lvol == 1) then
                sbar = 0.5_dp*(s + 1.0_dp)
                if (m == 0) then
                    fj = sbar*sbar; dfj = sbar; d2fj = 0.5_dp
                else if (m == 1) then
                    fj = sbar; dfj = 0.5_dp; d2fj = 0.0_dp
                else
                    fj = sbar**m
                    dfj = 0.5_dp*rm*sbar**(m - 1)
                    d2fj = 0.25_dp*real(m*(m - 1), dp)*sbar**(m - 2)
                end if
                rc = self%Rbc(ii, 0) + (self%Rbc(ii, 1) - self%Rbc(ii, 0))*fj
                rs = self%Rbs(ii, 0) + (self%Rbs(ii, 1) - self%Rbs(ii, 0))*fj
                zc = self%Zbc(ii, 0) + (self%Zbc(ii, 1) - self%Zbc(ii, 0))*fj
                zs = self%Zbs(ii, 0) + (self%Zbs(ii, 1) - self%Zbs(ii, 0))*fj
                drc = (self%Rbc(ii, 1) - self%Rbc(ii, 0))*dfj
                drs = (self%Rbs(ii, 1) - self%Rbs(ii, 0))*dfj
                dzc = (self%Zbc(ii, 1) - self%Zbc(ii, 0))*dfj
                dzs = (self%Zbs(ii, 1) - self%Zbs(ii, 0))*dfj
                d2rc = (self%Rbc(ii, 1) - self%Rbc(ii, 0))*d2fj
                d2rs = (self%Rbs(ii, 1) - self%Rbs(ii, 0))*d2fj
                d2zc = (self%Zbc(ii, 1) - self%Zbc(ii, 0))*d2fj
                d2zs = (self%Zbs(ii, 1) - self%Zbs(ii, 0))*d2fj
            else
                alss = 0.5_dp*(1.0_dp - s)
                blss = 0.5_dp*(1.0_dp + s)
                rc = alss*self%Rbc(ii, lvol - 1) + blss*self%Rbc(ii, lvol)
                rs = alss*self%Rbs(ii, lvol - 1) + blss*self%Rbs(ii, lvol)
                zc = alss*self%Zbc(ii, lvol - 1) + blss*self%Zbc(ii, lvol)
                zs = alss*self%Zbs(ii, lvol - 1) + blss*self%Zbs(ii, lvol)
                drc = 0.5_dp*(self%Rbc(ii, lvol) - self%Rbc(ii, lvol - 1))
                drs = 0.5_dp*(self%Rbs(ii, lvol) - self%Rbs(ii, lvol - 1))
                dzc = 0.5_dp*(self%Zbc(ii, lvol) - self%Zbc(ii, lvol - 1))
                dzs = 0.5_dp*(self%Zbs(ii, lvol) - self%Zbs(ii, lvol - 1))
                d2rc = 0.0_dp; d2rs = 0.0_dp; d2zc = 0.0_dp; d2zs = 0.0_dp
            end if

            R = R + rc*c + rs*sn
            Z = Z + zc*c + zs*sn
            dRds = dRds + drc*c + drs*sn
            dZds = dZds + dzc*c + dzs*sn
            dRdt = dRdt + rm*(-rc*sn + rs*c)
            dZdt = dZdt + rm*(-zc*sn + zs*c)
            dRdz = dRdz + rn*(rc*sn - rs*c)
            dZdz = dZdz + rn*(zc*sn - zs*c)

            valr = rc*c + rs*sn
            valz = zc*c + zs*sn
            d2R(1) = d2R(1) + d2rc*c + d2rs*sn
            d2R(2) = d2R(2) + rm*(-drc*sn + drs*c)
            d2R(3) = d2R(3) + rn*(drc*sn - drs*c)
            d2R(4) = d2R(4) - rm*rm*valr
            d2R(5) = d2R(5) + rm*rn*valr
            d2R(6) = d2R(6) - rn*rn*valr
            d2Z(1) = d2Z(1) + d2zc*c + d2zs*sn
            d2Z(2) = d2Z(2) + rm*(-dzc*sn + dzs*c)
            d2Z(3) = d2Z(3) + rn*(dzc*sn - dzs*c)
            d2Z(4) = d2Z(4) - rm*rm*valz
            d2Z(5) = d2Z(5) + rm*rn*valz
            d2Z(6) = d2Z(6) - rn*rn*valz
        end do
    end subroutine spectre_geom_d2

    subroutine spectre_evaluate_cyl(self, u, x)
        class(spectre_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: x(3)

        integer :: lvol
        real(dp) :: s, R, Z
        real(dp) :: dRds, dRdt, dRdz, dZds, dZdt, dZdz

        call spectre_locate(self%Mvol, u(1), lvol, s)
        call spectre_geom(self, lvol, s, u(2), u(3), R, Z, &
                          dRds, dRdt, dRdz, dZds, dZdt, dZdz)
        x(1) = R
        x(2) = u(3)
        x(3) = Z
    end subroutine spectre_evaluate_cyl

    subroutine spectre_evaluate_cart(self, u, x)
        class(spectre_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: x(3)

        real(dp) :: xcyl(3)

        call self%evaluate_cyl(u, xcyl)
        call cyl_to_cart(xcyl, x)
    end subroutine spectre_evaluate_cart

    subroutine spectre_covariant_basis(self, u, e_cov)
        class(spectre_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: e_cov(3, 3)

        integer :: lvol
        real(dp) :: s, R, Z, cz, sz
        real(dp) :: dRds, dRdt, dRdz, dZds, dZdt, dZdz

        call spectre_locate(self%Mvol, u(1), lvol, s)
        call spectre_geom(self, lvol, s, u(2), u(3), R, Z, &
                          dRds, dRdt, dRdz, dZds, dZdt, dZdz)

        cz = cos(u(3))
        sz = sin(u(3))

        e_cov(1, 1) = 2.0_dp*dRds*cz
        e_cov(2, 1) = 2.0_dp*dRds*sz
        e_cov(3, 1) = 2.0_dp*dZds

        e_cov(1, 2) = dRdt*cz
        e_cov(2, 2) = dRdt*sz
        e_cov(3, 2) = dZdt

        e_cov(1, 3) = dRdz*cz - R*sz
        e_cov(2, 3) = dRdz*sz + R*cz
        e_cov(3, 3) = dZdz
    end subroutine spectre_covariant_basis

    subroutine spectre_metric_tensor(self, u, g, ginv, sqrtg)
        class(spectre_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: g(3, 3), ginv(3, 3), sqrtg

        real(dp) :: e_cov(3, 3), det
        integer :: i, j

        call self%covariant_basis(u, e_cov)

        do i = 1, 3
            do j = 1, 3
                g(i, j) = dot_product(e_cov(:, i), e_cov(:, j))
            end do
        end do

        det = g(1, 1)*(g(2, 2)*g(3, 3) - g(2, 3)*g(3, 2)) &
              - g(1, 2)*(g(2, 1)*g(3, 3) - g(2, 3)*g(3, 1)) &
              + g(1, 3)*(g(2, 1)*g(3, 2) - g(2, 2)*g(3, 1))

        sqrtg = sqrt(abs(det))

        ginv(1, 1) = (g(2, 2)*g(3, 3) - g(2, 3)*g(3, 2))/det
        ginv(1, 2) = (g(1, 3)*g(3, 2) - g(1, 2)*g(3, 3))/det
        ginv(1, 3) = (g(1, 2)*g(2, 3) - g(1, 3)*g(2, 2))/det
        ginv(2, 1) = (g(2, 3)*g(3, 1) - g(2, 1)*g(3, 3))/det
        ginv(2, 2) = (g(1, 1)*g(3, 3) - g(1, 3)*g(3, 1))/det
        ginv(2, 3) = (g(1, 3)*g(2, 1) - g(1, 1)*g(2, 3))/det
        ginv(3, 1) = (g(2, 1)*g(3, 2) - g(2, 2)*g(3, 1))/det
        ginv(3, 2) = (g(1, 2)*g(3, 1) - g(1, 1)*g(3, 2))/det
        ginv(3, 3) = (g(1, 1)*g(2, 2) - g(1, 2)*g(2, 1))/det
    end subroutine spectre_metric_tensor

    pure subroutine metric_from_ecov(e_cov, g, ginv, sqrtg)
        !> Covariant metric g_ij = e_i . e_j, its inverse and Jacobian sqrt|det g|
        !> for the analytic-derivative path. spectre_metric_tensor keeps its own
        !> inline copy so the symplectic construction flowing through it stays
        !> byte-for-byte unchanged.
        real(dp), intent(in) :: e_cov(3, 3)
        real(dp), intent(out) :: g(3, 3), ginv(3, 3), sqrtg

        real(dp) :: det
        integer :: i, j

        do i = 1, 3
            do j = 1, 3
                g(i, j) = dot_product(e_cov(:, i), e_cov(:, j))
            end do
        end do

        det = g(1, 1)*(g(2, 2)*g(3, 3) - g(2, 3)*g(3, 2)) &
              - g(1, 2)*(g(2, 1)*g(3, 3) - g(2, 3)*g(3, 1)) &
              + g(1, 3)*(g(2, 1)*g(3, 2) - g(2, 2)*g(3, 1))

        sqrtg = sqrt(abs(det))

        ginv(1, 1) = (g(2, 2)*g(3, 3) - g(2, 3)*g(3, 2))/det
        ginv(1, 2) = (g(1, 3)*g(3, 2) - g(1, 2)*g(3, 3))/det
        ginv(1, 3) = (g(1, 2)*g(2, 3) - g(1, 3)*g(2, 2))/det
        ginv(2, 1) = (g(2, 3)*g(3, 1) - g(2, 1)*g(3, 3))/det
        ginv(2, 2) = (g(1, 1)*g(3, 3) - g(1, 3)*g(3, 1))/det
        ginv(2, 3) = (g(1, 3)*g(2, 1) - g(1, 1)*g(2, 3))/det
        ginv(3, 1) = (g(2, 1)*g(3, 2) - g(2, 2)*g(3, 1))/det
        ginv(3, 2) = (g(1, 2)*g(3, 1) - g(1, 1)*g(3, 2))/det
        ginv(3, 3) = (g(1, 1)*g(2, 2) - g(1, 2)*g(2, 1))/det
    end subroutine metric_from_ecov

    subroutine spectre_metric_tensor_der(self, u, g, ginv, sqrtg, dg, dsqrtg)
        !> Metric, inverse, Jacobian and their analytic first derivatives:
        !> dg(i,k,j) = d g_ik/du_j and dsqrtg(j) = d sqrtg/du_j. The chart chain
        !> ds/drho_g = 2 scales radial (j = 1) geometry derivatives; dsqrtg uses
        !> Jacobi's formula d sqrtg = (sqrtg/2) g^ik d g_ik.
        class(spectre_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: g(3, 3), ginv(3, 3), sqrtg
        real(dp), intent(out) :: dg(3, 3, 3), dsqrtg(3)

        integer :: lvol, i, j, k, a
        real(dp) :: s, R, Z, cz, sz
        real(dp) :: dRds, dRdt, dRdz, dZds, dZdt, dZdz
        real(dp) :: d2R(6), d2Z(6)
        real(dp) :: Ru(3), Zu(3), Ruu(3, 3), Zuu(3, 3), d3(3)
        real(dp) :: e_cov(3, 3), de_cov(3, 3, 3), tr

        call spectre_locate(self%Mvol, u(1), lvol, s)
        call spectre_geom_d2(self, lvol, s, u(2), u(3), R, Z, &
                             dRds, dRdt, dRdz, dZds, dZdt, dZdz, d2R, d2Z)

        cz = cos(u(3))
        sz = sin(u(3))

        ! Chart chain d/drho_g = 2 d/ds: radial first/second derivatives carry 2
        ! and 4, radial-angle mixed second derivatives carry 2.
        Ru = [2.0_dp*dRds, dRdt, dRdz]
        Zu = [2.0_dp*dZds, dZdt, dZdz]
        Ruu(1, 1) = 4.0_dp*d2R(1); Ruu(1, 2) = 2.0_dp*d2R(2); Ruu(1, 3) = 2.0_dp*d2R(3)
        Ruu(2, 2) = d2R(4); Ruu(2, 3) = d2R(5); Ruu(3, 3) = d2R(6)
        Ruu(2, 1) = Ruu(1, 2); Ruu(3, 1) = Ruu(1, 3); Ruu(3, 2) = Ruu(2, 3)
        Zuu(1, 1) = 4.0_dp*d2Z(1); Zuu(1, 2) = 2.0_dp*d2Z(2); Zuu(1, 3) = 2.0_dp*d2Z(3)
        Zuu(2, 2) = d2Z(4); Zuu(2, 3) = d2Z(5); Zuu(3, 3) = d2Z(6)
        Zuu(2, 1) = Zuu(1, 2); Zuu(3, 1) = Zuu(1, 3); Zuu(3, 2) = Zuu(2, 3)

        ! Cartesian embedding X = (R cos zeta, R sin zeta, Z); e_cov(a,i) = dX_a/du_i,
        ! de_cov(a,i,j) = d2X_a/du_i du_j. Only the zeta index (3) sees cos/sin.
        d3 = [0.0_dp, 0.0_dp, 1.0_dp]
        do i = 1, 3
            e_cov(1, i) = Ru(i)*cz - R*sz*d3(i)
            e_cov(2, i) = Ru(i)*sz + R*cz*d3(i)
            e_cov(3, i) = Zu(i)
        end do
        do j = 1, 3
            do i = 1, 3
                de_cov(1, i, j) = Ruu(i, j)*cz - sz*(Ru(i)*d3(j) + Ru(j)*d3(i)) &
                                  - R*cz*d3(i)*d3(j)
                de_cov(2, i, j) = Ruu(i, j)*sz + cz*(Ru(i)*d3(j) + Ru(j)*d3(i)) &
                                  - R*sz*d3(i)*d3(j)
                de_cov(3, i, j) = Zuu(i, j)
            end do
        end do

        call metric_from_ecov(e_cov, g, ginv, sqrtg)

        do j = 1, 3
            do k = 1, 3
                do i = 1, 3
                    dg(i, k, j) = 0.0_dp
                    do a = 1, 3
                        dg(i, k, j) = dg(i, k, j) + de_cov(a, i, j)*e_cov(a, k) &
                                      + e_cov(a, i)*de_cov(a, k, j)
                    end do
                end do
            end do
        end do

        do j = 1, 3
            tr = 0.0_dp
            do k = 1, 3
                do i = 1, 3
                    tr = tr + ginv(i, k)*dg(i, k, j)
                end do
            end do
            dsqrtg(j) = 0.5_dp*sqrtg*tr
        end do
    end subroutine spectre_metric_tensor_der

    subroutine spectre_from_cyl(self, xcyl, u, ierr)
        !> Inverse map: for each volume, seed theta by a coarse poloidal scan
        !> (the SPEC poloidal angle need not track the geometric angle), run a
        !> 2D Newton in (s, theta) at fixed zeta = phi, and accept the volume
        !> whose converged s lands in [-1, 1].
        class(spectre_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: xcyl(3)
        real(dp), intent(out) :: u(3)
        integer, intent(out) :: ierr

        real(dp), parameter :: tol_res = 1.0e-11_dp
        real(dp), parameter :: tol_det = 1.0e-14_dp
        real(dp), parameter :: s_slack = 1.0e-9_dp

        real(dp) :: zeta, theta0, s, theta
        integer :: lvol
        logical :: converged

        ierr = 1
        u = 0.0_dp
        zeta = xcyl(2)

        do lvol = 1, self%Mvol
            theta0 = spectre_seed_theta(self, lvol, xcyl(1), xcyl(3), zeta)
            call newton_volume(self, lvol, xcyl(1), xcyl(3), zeta, 0.0_dp, theta0, &
                               tol_res, tol_det, s, theta, converged)
            if (converged) then
                if (s >= -1.0_dp - s_slack) then
                    if (s <= 1.0_dp + s_slack) then
                        u(1) = real(lvol - 1, dp) + 0.5_dp*(s + 1.0_dp)
                        u(1) = min(max(u(1), 0.0_dp), real(self%Mvol, dp))
                        u(2) = modulo(theta, TWOPI)
                        u(3) = zeta
                        ierr = 0
                        return
                    end if
                end if
            end if
        end do
    end subroutine spectre_from_cyl

    subroutine spectre_from_cyl_warm(self, xcyl, u, lvol, ierr)
        class(spectre_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: xcyl(3)
        real(dp), intent(inout) :: u(3)
        integer, intent(in) :: lvol
        integer, intent(out) :: ierr

        real(dp), parameter :: tol_res = 1.0e-11_dp
        real(dp), parameter :: tol_det = 1.0e-14_dp
        real(dp), parameter :: s_slack = 1.0e-9_dp

        real(dp) :: s0, s, theta
        logical :: converged

        ierr = 1
        if (lvol < 1 .or. lvol > self%Mvol) return

        s0 = 2.0_dp*(u(1) - real(lvol - 1, dp)) - 1.0_dp
        call newton_volume(self, lvol, xcyl(1), xcyl(3), xcyl(2), s0, u(2), &
                           tol_res, tol_det, s, theta, converged)
        if (.not. converged) return
        if (s < -1.0_dp - s_slack .or. s > 1.0_dp + s_slack) return

        if (abs(s + 1.0_dp) <= s_slack) s = -1.0_dp
        if (abs(s - 1.0_dp) <= s_slack) s = 1.0_dp
        u(1) = real(lvol - 1, dp) + 0.5_dp*(s + 1.0_dp)
        u(2) = modulo(theta, TWOPI)
        u(3) = xcyl(2)
        ierr = 0
    end subroutine spectre_from_cyl_warm

    pure function spectre_seed_theta(self, lvol, R_tgt, Z_tgt, zeta) result(theta0)
        class(spectre_coordinate_system_t), intent(in) :: self
        integer, intent(in) :: lvol
        real(dp), intent(in) :: R_tgt, Z_tgt, zeta
        real(dp) :: theta0

        integer, parameter :: nscan = 16
        integer :: i
        real(dp) :: t, R, Z, res, best
        real(dp) :: dRds, dRdt, dRdz, dZds, dZdt, dZdz

        best = huge(1.0_dp)
        theta0 = 0.0_dp
        do i = 0, nscan - 1
            t = TWOPI*real(i, dp)/real(nscan, dp)
            call spectre_geom(self, lvol, 0.0_dp, t, zeta, R, Z, &
                              dRds, dRdt, dRdz, dZds, dZdt, dZdz)
            res = (R - R_tgt)**2 + (Z - Z_tgt)**2
            if (res < best) then
                best = res
                theta0 = t
            end if
        end do
    end function spectre_seed_theta

    subroutine newton_volume(self, lvol, R_tgt, Z_tgt, zeta, s0, theta0, &
                             tol_res, tol_det, s, theta, converged)
        class(spectre_coordinate_system_t), intent(in) :: self
        integer, intent(in) :: lvol
        real(dp), intent(in) :: R_tgt, Z_tgt, zeta, s0, theta0, tol_res, tol_det
        real(dp), intent(out) :: s, theta
        logical, intent(out) :: converged

        integer, parameter :: max_iter = 60
        real(dp) :: R, Z, dRds, dRdt, dRdz, dZds, dZdt, dZdz
        real(dp) :: f1, f2, det, ds, dt, res, res_try, alpha
        real(dp) :: s_try, theta_try
        integer :: iter, k

        converged = .false.
        s = min(1.0_dp, max(-1.0_dp, s0))
        theta = theta0

        do iter = 1, max_iter
            call spectre_geom(self, lvol, s, theta, zeta, R, Z, &
                              dRds, dRdt, dRdz, dZds, dZdt, dZdz)
            f1 = R - R_tgt
            f2 = Z - Z_tgt
            res = sqrt(f1*f1 + f2*f2)
            if (res < tol_res) then
                converged = .true.
                return
            end if

            det = dRds*dZdt - dRdt*dZds
            if (abs(det) < tol_det) return

            ds = (-f1*dZdt + f2*dRdt)/det
            dt = (-dRds*f2 + dZds*f1)/det

            alpha = 1.0_dp
            do k = 1, 12
                s_try = min(1.0_dp, max(-1.0_dp, s + alpha*ds))
                theta_try = theta + alpha*dt
                call spectre_geom(self, lvol, s_try, theta_try, zeta, R, Z, &
                                  dRds, dRdt, dRdz, dZds, dZdt, dZdz)
                res_try = sqrt((R - R_tgt)**2 + (Z - Z_tgt)**2)
                if (res_try < res) then
                    s = s_try
                    theta = theta_try
                    exit
                end if
                alpha = 0.5_dp*alpha
            end do
            if (k > 12) return
        end do
    end subroutine newton_volume

end module libneo_coordinates_spectre
