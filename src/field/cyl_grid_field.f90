!> Array-backed magnetic field on a uniform cylindrical (R, Z, phi) grid.
!>
!> Coordinates are x = (R, phi, Z) and B = (B_R, B_phi, B_Z) in the local
!> cylindrical basis, the convention of the field-line tools (poincare).
!> The grid covers one field period 0 <= phi < 2 pi / nfp with nphi planes
!> (phi_k = (k-1) * 2 pi / (nfp * nphi), the plane at the period end is the
!> first plane again).  The field is periodic in phi.  Interpolation is
!> 6-point (fifth order) Lagrange in each direction.  Grid nodes can be marked
!> invalid (outside the computational domain); a point whose interpolation
!> stencil touches an invalid node is "outside": compute_bfield returns B = 0
!> and contains_point is .false.
module neo_cyl_grid_field
    use neo_field_base, only: field_t
    implicit none
    private
    public :: cyl_grid_field_t

    integer, parameter :: dp = kind(1.0d0)
    integer, parameter :: np = 6

    type, extends(field_t) :: cyl_grid_field_t
        integer :: nr = 0, nz = 0, nphi = 0, nfp = 1
        real(dp) :: r_min = 0.0_dp, z_min = 0.0_dp
        real(dp) :: dr = 1.0_dp, dz = 1.0_dp, dphi = 1.0_dp
        !> b(1:3, ir, iz, k) = (B_R, B_phi, B_Z)
        real(dp), allocatable :: b(:, :, :, :)
        !> .false. for nodes outside the computational domain
        logical, allocatable :: valid(:, :, :)
    contains
        procedure :: cyl_grid_field_init
        procedure :: compute_afield
        procedure :: compute_bfield
        procedure :: compute_abfield
        procedure :: compute_bfield_jacobian
        procedure :: contains_point
    end type cyl_grid_field_t

contains

    !> r and z are uniform node coordinates, b(3, nr, nz, nphi) the field.
    !> Without valid, nodes with B = 0 exactly are treated as outside.
    subroutine cyl_grid_field_init(self, r, z, nfp, b, valid)
        class(cyl_grid_field_t), intent(out) :: self
        real(dp), intent(in) :: r(:), z(:)
        integer, intent(in) :: nfp
        real(dp), intent(in) :: b(:, :, :, :)
        logical, intent(in), optional :: valid(:, :, :)

        if (size(b, 1) /= 3 .or. size(b, 2) /= size(r) .or. &
            size(b, 3) /= size(z)) error stop 'cyl_grid_field: shape mismatch'
        if (size(r) < np .or. size(z) < np .or. size(b, 4) < np) &
            error stop 'cyl_grid_field: grid too small for 6-point stencil'
        if (nfp < 1) error stop 'cyl_grid_field: nfp must be positive'
        self%nr = size(r)
        self%nz = size(z)
        self%nphi = size(b, 4)
        self%nfp = nfp
        self%r_min = r(1)
        self%z_min = z(1)
        self%dr = (r(self%nr) - r(1))/real(self%nr - 1, dp)
        self%dz = (z(self%nz) - z(1))/real(self%nz - 1, dp)
        self%dphi = 8.0_dp*atan(1.0_dp)/real(nfp*self%nphi, dp)
        allocate (self%b, source=b)
        allocate (self%valid(self%nr, self%nz, self%nphi))
        if (present(valid)) then
            self%valid = valid
        else
            self%valid = any(b /= 0.0_dp, dim=1)
        end if
    end subroutine cyl_grid_field_init

    subroutine compute_afield(self, x, A)
        class(cyl_grid_field_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: A(3)
        A = 0.0_dp
        error stop 'cyl_grid_field: no vector potential available'
    end subroutine compute_afield

    subroutine compute_abfield(self, x, A, B)
        class(cyl_grid_field_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: A(3), B(3)
        A = 0.0_dp
        B = 0.0_dp
        error stop 'cyl_grid_field: no vector potential available'
    end subroutine compute_abfield

    !> Fifth-order Lagrange weights for nodes 0..5 at position t, and the
    !> derivative weights with respect to t.
    pure subroutine lag6(t, w, dw)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: w(np), dw(np)
        integer :: i, j, k
        real(dp) :: num, dnum, den

        do i = 1, np
            den = 1.0_dp
            num = 1.0_dp
            dnum = 0.0_dp
            do j = 1, np
                if (j == i) cycle
                den = den*real(i - j, dp)
                num = num*(t - real(j - 1, dp))
            end do
            do k = 1, np
                if (k == i) cycle
                num = 1.0_dp
                do j = 1, np
                    if (j == i .or. j == k) cycle
                    num = num*(t - real(j - 1, dp))
                end do
                dnum = dnum + num
            end do
            dw(i) = dnum/den
            num = 1.0_dp
            do j = 1, np
                if (j == i) cycle
                num = num*(t - real(j - 1, dp))
            end do
            w(i) = num/den
        end do
    end subroutine lag6

    !> Stencil: first node indices (1-based) in R and Z, phi node list, and
    !> weights.  inside is .false. when x leaves the grid in R or Z.
    subroutine stencil(self, x, ir0, iz0, ik, wr, dwr, wz, dwz, wp, dwp, inside)
        class(cyl_grid_field_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        integer, intent(out) :: ir0, iz0, ik(np)
        real(dp), intent(out) :: wr(np), dwr(np), wz(np), dwz(np)
        real(dp), intent(out) :: wp(np), dwp(np)
        logical, intent(out) :: inside
        real(dp) :: u, v, s, period, tp
        integer :: i0, j

        inside = .false.
        ir0 = 1
        iz0 = 1
        ik = 1
        wr = 0.0_dp
        dwr = 0.0_dp
        wz = 0.0_dp
        dwz = 0.0_dp
        wp = 0.0_dp
        dwp = 0.0_dp
        u = (x(1) - self%r_min)/self%dr
        v = (x(3) - self%z_min)/self%dz
        if (u < 0.0_dp .or. v < 0.0_dp) return
        if (u > real(self%nr - 1, dp) .or. v > real(self%nz - 1, dp)) return
        ir0 = max(0, min(self%nr - np, int(u) - 2))
        iz0 = max(0, min(self%nz - np, int(v) - 2))
        call lag6(u - real(ir0, dp), wr, dwr)
        call lag6(v - real(iz0, dp), wz, dwz)
        ir0 = ir0 + 1
        iz0 = iz0 + 1
        period = self%dphi*real(self%nphi, dp)
        s = x(2) - period*floor(x(2)/period)
        s = s/self%dphi
        i0 = int(s)
        tp = s - real(i0, dp)
        do j = 1, np
            ik(j) = modulo(i0 - 3 + j, self%nphi) + 1
        end do
        call lag6(tp + 2.0_dp, wp, dwp)
        inside = .true.
    end subroutine stencil

    logical function stencil_valid(self, ir0, iz0, ik)
        class(cyl_grid_field_t), intent(in) :: self
        integer, intent(in) :: ir0, iz0, ik(np)
        integer :: j
        stencil_valid = .true.
        do j = 1, np
            if (.not. all(self%valid(ir0:ir0 + np - 1, iz0:iz0 + np - 1, ik(j)))) then
                stencil_valid = .false.
                return
            end if
        end do
    end function stencil_valid

    !> True when the whole interpolation stencil lies in the valid domain.
    logical function contains_point(self, x)
        class(cyl_grid_field_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        integer :: ir0, iz0, ik(np)
        real(dp) :: wr(np), dwr(np), wz(np), dwz(np), wp(np), dwp(np)
        logical :: inside

        call stencil(self, x, ir0, iz0, ik, wr, dwr, wz, dwz, wp, dwp, inside)
        contains_point = .false.
        if (.not. inside) return
        contains_point = stencil_valid(self, ir0, iz0, ik)
    end function contains_point

    subroutine compute_bfield(self, x, B)
        class(cyl_grid_field_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: B(3)
        real(dp) :: dB(3, 3)
        call self%compute_bfield_jacobian(x, B, dB)
    end subroutine compute_bfield

    !> B and dB(i, j) = d B_i / d x_j with x = (R, phi, Z): the phi derivative
    !> is with respect to the angle (not divided by R) and acts on the
    !> components in the local basis.  Returns zeros outside the domain.
    subroutine compute_bfield_jacobian(self, x, B, dB)
        class(cyl_grid_field_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: B(3), dB(3, 3)
        integer :: ir0, iz0, ik(np), i, j, k, c
        real(dp) :: wr(np), dwr(np), wz(np), dwz(np), wp(np), dwp(np)
        real(dp) :: wzp, dwzp, wzpz, f(3)
        logical :: inside

        B = 0.0_dp
        dB = 0.0_dp
        call stencil(self, x, ir0, iz0, ik, wr, dwr, wz, dwz, wp, dwp, inside)
        if (.not. inside) return
        if (.not. stencil_valid(self, ir0, iz0, ik)) return
        do k = 1, np
            do j = 1, np
                wzp = wz(j)*wp(k)
                dwzp = dwz(j)*wp(k)
                wzpz = wz(j)*dwp(k)
                do i = 1, np
                    f = self%b(:, ir0 + i - 1, iz0 + j - 1, ik(k))
                    do c = 1, 3
                        B(c) = B(c) + f(c)*wr(i)*wzp
                        dB(c, 1) = dB(c, 1) + f(c)*dwr(i)*wzp
                        dB(c, 2) = dB(c, 2) + f(c)*wr(i)*wzpz
                        dB(c, 3) = dB(c, 3) + f(c)*wr(i)*dwzp
                    end do
                end do
            end do
        end do
        dB(:, 1) = dB(:, 1)/self%dr
        dB(:, 2) = dB(:, 2)/self%dphi
        dB(:, 3) = dB(:, 3)/self%dz
    end subroutine compute_bfield_jacobian

end module neo_cyl_grid_field
