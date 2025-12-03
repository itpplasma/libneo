module gframe_boundary
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: gframe_boundary_t
    public :: make_circular_boundary
    public :: boundary_from_modes
    public :: allocate_boundary

    type :: gframe_boundary_t
        integer :: ntheta = 0
        integer :: nzeta = 0
        real(dp), allocatable :: theta(:)
        real(dp), allocatable :: zeta(:)
        real(dp), allocatable :: Rb(:, :)
        real(dp), allocatable :: Zb(:, :)
    end type gframe_boundary_t

contains

    subroutine make_circular_boundary(boundary, ntheta, nzeta, r0, a)
        type(gframe_boundary_t), intent(out) :: boundary
        integer, intent(in) :: ntheta, nzeta
        real(dp), intent(in) :: r0, a

        integer :: jt, jz
        real(dp) :: theta_val, zeta_val

        if (ntheta < 2 .or. nzeta < 2) then
            error stop "make_circular_boundary: require ntheta,nzeta >= 2"
        end if

        call allocate_boundary(boundary, ntheta, nzeta)

        do jz = 1, nzeta
            zeta_val = two_pi() * real(jz - 1, dp) / real(nzeta - 1, dp)
            boundary%zeta(jz) = zeta_val
            do jt = 1, ntheta
                theta_val = two_pi() * real(jt - 1, dp) / real(ntheta - 1, dp)
                boundary%theta(jt) = theta_val
                boundary%Rb(jt, jz) = r0 + a * cos(theta_val)
                boundary%Zb(jt, jz) = a * sin(theta_val)
            end do
        end do
    end subroutine make_circular_boundary

    subroutine boundary_from_modes(m, n, rcos, rsin, zcos, zsin, nfp, ntheta, &
        nzeta, boundary)
        integer, intent(in) :: m(:), n(:)
        real(dp), intent(in) :: rcos(:), rsin(:), zcos(:), zsin(:)
        integer, intent(in) :: nfp, ntheta, nzeta
        type(gframe_boundary_t), intent(out) :: boundary

        integer :: nmodes
        integer :: jt, jz, idx
        real(dp) :: theta_val, zeta_val
        real(dp) :: phase

        nmodes = size(m)

        if (ntheta < 2 .or. nzeta < 2) then
            error stop "boundary_from_modes: require ntheta,nzeta >= 2"
        end if

        if (nmodes /= size(n) .or. nmodes /= size(rcos) .or. nmodes /= size(rsin) &
            .or. nmodes /= size(zcos) .or. nmodes /= size(zsin)) then
            error stop "boundary_from_modes: mode arrays must have equal length"
        end if

        call allocate_boundary(boundary, ntheta, nzeta)

        do jz = 1, nzeta
            zeta_val = two_pi() * real(jz - 1, dp) / real(nzeta - 1, dp)
            boundary%zeta(jz) = zeta_val
            do jt = 1, ntheta
                theta_val = two_pi() * real(jt - 1, dp) / real(ntheta - 1, dp)
                boundary%theta(jt) = theta_val
                boundary%Rb(jt, jz) = 0.0_dp
                boundary%Zb(jt, jz) = 0.0_dp
                do idx = 1, nmodes
                    phase = real(m(idx), dp) * theta_val - real(n(idx), dp) &
                        * real(nfp, dp) * zeta_val
                    boundary%Rb(jt, jz) = boundary%Rb(jt, jz) &
                        + rcos(idx) * cos(phase) + rsin(idx) * sin(phase)
                    boundary%Zb(jt, jz) = boundary%Zb(jt, jz) &
                        + zcos(idx) * cos(phase) + zsin(idx) * sin(phase)
                end do
            end do
        end do
    end subroutine boundary_from_modes

    subroutine allocate_boundary(boundary, ntheta, nzeta)
        type(gframe_boundary_t), intent(inout) :: boundary
        integer, intent(in) :: ntheta, nzeta

        if (allocated(boundary%theta)) deallocate(boundary%theta)
        if (allocated(boundary%zeta)) deallocate(boundary%zeta)
        if (allocated(boundary%Rb)) deallocate(boundary%Rb)
        if (allocated(boundary%Zb)) deallocate(boundary%Zb)

        allocate(boundary%theta(ntheta))
        allocate(boundary%zeta(nzeta))
        allocate(boundary%Rb(ntheta, nzeta))
        allocate(boundary%Zb(ntheta, nzeta))

        boundary%ntheta = ntheta
        boundary%nzeta = nzeta
    end subroutine allocate_boundary

    pure real(dp) function two_pi()
        two_pi = 2.0_dp * acos(-1.0_dp)
    end function two_pi

end module gframe_boundary
