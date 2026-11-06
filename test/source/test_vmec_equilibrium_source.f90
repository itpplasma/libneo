!> VMEC equilibrium source of the frozen-current field on the
!> LandremanPaul2021 QA reference equilibrium: B^varphi = B_phi/R and |B|
!> match vmec_field's symmetry-flux components and metric, the volume
!> quadrature integrates the wout volume, the boundary area vectors are
!> outward and close the volume by the divergence theorem, and B is
!> tangential to the boundary.
program test_vmec_equilibrium_source
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use netcdf
    use neo_equilibrium_source, only: source_points_t, source_resolution_t
    use neo_vmec_equilibrium_source, only: vmec_equilibrium_source_t, &
                                           vmec_equilibrium_source_init
    use spline_vmec_sub, only: vmec_field
    implicit none

    type(vmec_equilibrium_source_t) :: src
    type(source_resolution_t) :: res
    type(source_points_t) :: vol, surf
    real(dp) :: volume, volume_wout, closure, bn, bmax, bctr_phi, bmod
    integer :: k

    call vmec_equilibrium_source_init(src, "wout.nc")
    call check(src%nfp == 2, "nfp from the wout file")

    res%ns = 24
    res%ntheta = 48
    res%nzeta = 24
    call src%sample_volume(res, vol)
    call src%sample_surface(res, surf)

    do k = 1, size(vol%w), 97
        call vmec_reference(vol, k, bctr_phi, bmod)
        call check(abs(vol%b(2, k)/vol%x(1, k) - bctr_phi) <= &
                   1.0e-12_dp*abs(bctr_phi), "B^varphi agrees with vmec_field")
        call check(abs(norm2(vol%b(:, k)) - bmod) <= 1.0e-12_dp*bmod, &
                   "|B| agrees with vmec_field")
    end do

    volume = src%nfp*sum(vol%w)
    volume_wout = wout_volume("wout.nc")*1.0e6_dp
    print *, "volume relative error", volume/volume_wout - 1.0_dp
    call check(abs(volume/volume_wout - 1.0_dp) < 1.0e-4_dp, "volume")

    closure = src%nfp*sum(surf%x(1, :)*surf%area(1, :) + &
                          surf%x(3, :)*surf%area(3, :))/3.0_dp
    print *, "divergence-theorem volume relative error", closure/volume - 1.0_dp
    call check(abs(closure/volume - 1.0_dp) < 1.0e-6_dp, "outward area vectors")

    bmax = 0.0_dp
    do k = 1, size(surf%w)
        bmod = norm2(surf%b(:, k))*norm2(surf%area(:, k))
        bn = dot_product(surf%b(:, k), surf%area(:, k))
        bmax = max(bmax, abs(bn)/bmod)
    end do
    call check(bmax < 1.0e-12_dp, "B tangential to the boundary")
    print *, "PASS"

contains

    subroutine vmec_reference(pts, k, bctr_phi, bmod)
        type(source_points_t), intent(in) :: pts
        integer, intent(in) :: k
        real(dp), intent(out) :: bctr_phi, bmod
        real(dp) :: s, theta, a_theta, a_phi, da_theta, da_phi, iota, sqg, lam
        real(dp) :: dl_ds, dl_dt, dl_dp, bctr_t, bcov_s, bcov_t, bcov_p
        integer :: i, j

        i = mod(k - 1, res%ns) + 1
        j = mod((k - 1)/res%ns, res%ntheta) + 1
        s = gauss_node(i)
        theta = (j - 1)*8.0_dp*atan(1.0_dp)/res%ntheta
        call vmec_field(s, theta, pts%x(2, k), a_theta, a_phi, da_theta, da_phi, &
                        iota, sqg, lam, dl_ds, dl_dt, dl_dp, bctr_t, bctr_phi, &
                        bcov_s, bcov_t, bcov_p)
        bmod = sqrt(bctr_t*bcov_t + bctr_phi*bcov_p)
    end subroutine vmec_reference

    real(dp) function gauss_node(i) result(s)
        use fortnum_quadrature, only: gauss_legendre_ab
        integer, intent(in) :: i
        real(dp) :: x(res%ns), w(res%ns)

        call gauss_legendre_ab(res%ns, 0.0_dp, 1.0_dp, x, w)
        s = x(i)
    end function gauss_node

    real(dp) function wout_volume(path) result(v)
        character(*), intent(in) :: path
        integer :: ncid, vid

        call check(nf90_open(path, NF90_NOWRITE, ncid) == NF90_NOERR, "open wout")
        call check(nf90_inq_varid(ncid, "volume_p", vid) == NF90_NOERR, "volume_p")
        call check(nf90_get_var(ncid, vid, v) == NF90_NOERR, "volume_p")
        call check(nf90_close(ncid) == NF90_NOERR, "close wout")
    end function wout_volume

    subroutine check(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what

        if (.not. ok) then
            print *, "FAIL: ", what
            error stop 1
        end if
    end subroutine check

end program test_vmec_equilibrium_source
