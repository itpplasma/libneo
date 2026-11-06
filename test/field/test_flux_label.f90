!> Continued flux label on a circular torus boundary of minor radius b:
!> outside, s = (1 + d/b)^2 with d the distance to the circle; inside, the
!> given label is kept and clipped to 1.  The BEAMS3D export writes it.
module test_flux_label_circle
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neo_equilibrium_source, only: equilibrium_source_t, source_points_t, &
                                      source_resolution_t
    implicit none

    type, extends(equilibrium_source_t) :: circle_t
        real(dp) :: r0 = 300.0_dp, b = 80.0_dp
    contains
        procedure :: sample_volume => no_volume
        procedure :: sample_surface => no_surface
        procedure :: boundary_curve => circle_curve
    end type circle_t

contains

    subroutine circle_curve(self, phi, theta, r, z)
        class(circle_t), intent(in) :: self
        real(dp), intent(in) :: phi, theta(:)
        real(dp), intent(out) :: r(:), z(:)

        r = self%r0 + self%b*cos(theta) + 0.0_dp*phi
        z = self%b*sin(theta)
    end subroutine circle_curve

    subroutine no_volume(self, res, pts)
        class(circle_t), intent(in) :: self
        type(source_resolution_t), intent(in) :: res
        type(source_points_t), intent(out) :: pts
        error stop "not used"
    end subroutine no_volume

    subroutine no_surface(self, res, pts)
        class(circle_t), intent(in) :: self
        type(source_resolution_t), intent(in) :: res
        type(source_points_t), intent(out) :: pts
        error stop "not used"
    end subroutine no_surface

end module test_flux_label_circle

program test_flux_label
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use test_flux_label_circle, only: circle_t
    use hdf5_tools, only: hid_t, h5_init, h5_deinit, h5_create, h5_open, h5_close, &
                          h5_add, h5_get
    use neo_biotsavart, only: coils_init
    use neo_biotsavart_field, only: biotsavart_field_t
    use neo_beams3d_field_grid, only: write_beams3d_field
    use neo_flux_label, only: boundary_section_t, boundary_section, &
                              continued_flux_label
    implicit none

    type(circle_t) :: src
    type(boundary_section_t) :: sec
    type(biotsavart_field_t) :: loop
    real(dp) :: s, d, raxis(3), phiaxis(2), zaxis(1), sarr(3, 2, 1), sout(3, 2, 1)
    integer(hid_t) :: h5id
    integer :: i

    call boundary_section(src, 0.3_dp, 1024, sec)
    print *, "effective minor radius", sec%a
    call check(abs(sec%a/src%b - 1.0_dp) < 1.0e-5_dp, "a = sqrt(A/pi)")
    do i = 1, 4
        d = 10.0_dp*i
        s = continued_flux_label(sec, src%r0 + src%b + d, 0.0_dp, 0.0_dp)
        call check(abs(s/(1.0_dp + d/src%b)**2 - 1.0_dp) < 1.0e-4_dp, &
                   "outside: (1 + d/a)^2")
        s = continued_flux_label(sec, src%r0, -(src%b + d), 0.0_dp)
        call check(abs(s/(1.0_dp + d/src%b)**2 - 1.0_dp) < 1.0e-4_dp, &
                   "outside below")
    end do
    s = continued_flux_label(sec, src%r0 + 0.5_dp*src%b, 0.0_dp, 0.25_dp)
    call check(s == 0.25_dp, "inside: label kept")
    s = continued_flux_label(sec, src%r0 + 0.99_dp*src%b, 0.0_dp, 1.5_dp)
    call check(s == 1.0_dp, "inside: label at most 1")
    s = continued_flux_label(sec, src%r0 + 1.001_dp*src%b, 0.0_dp, 0.0_dp)
    call check(abs(s - 1.0_dp) < 3.0e-3_dp, "continuous at the boundary")

    ! the BEAMS3D export: R = 3.4, 4.0, 4.6 m on Z = 0, the first inside
    raxis = [3.4_dp, 4.0_dp, 4.6_dp]
    phiaxis = [0.0_dp, 0.5_dp]
    zaxis = [0.0_dp]
    sarr = 0.3_dp
    call h5_init()
    call h5_create('flux_label_template.h5', h5id)
    call h5_add(h5id, 'nr', 3)
    call h5_add(h5id, 'nphi', 2)
    call h5_add(h5id, 'nz', 1)
    call h5_add(h5id, 'raxis', raxis, lbound(raxis), ubound(raxis))
    call h5_add(h5id, 'phiaxis', phiaxis, lbound(phiaxis), ubound(phiaxis))
    call h5_add(h5id, 'zaxis', zaxis, lbound(zaxis), ubound(zaxis))
    call h5_add(h5id, 'S_ARR', sarr, lbound(sarr), ubound(sarr))
    call h5_add(h5id, 'U_ARR', sarr, lbound(sarr), ubound(sarr))
    call h5_close(h5id)
    call h5_deinit()
    call coils_init([200.0_dp, 0.0_dp, -200.0_dp], [0.0_dp, 200.0_dp, 0.0_dp], &
                    [500.0_dp, 500.0_dp, 500.0_dp], [1.0_dp, 1.0_dp, 0.0_dp], &
                    loop%coils)
    call write_beams3d_field('flux_label_template.h5', 'flux_label_out.h5', loop, src)
    call h5_init()
    call h5_open('flux_label_out.h5', h5id)
    call h5_get(h5id, 'S_ARR', sout)
    call h5_close(h5id)
    call h5_deinit()
    call check(all(sout(1, :, 1) == 0.3_dp), "export keeps the label inside")
    do i = 2, 3
        d = raxis(i)*100.0_dp - src%r0 - src%b
        call check(all(abs(sout(i, :, 1)/(1.0_dp + d/src%b)**2 - 1.0_dp) < 1.0e-4_dp), &
                   "export continues the label outside")
    end do
    print *, "PASS"

contains

    subroutine check(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what

        if (.not. ok) then
            print *, "FAIL: ", what
            error stop 1
        end if
    end subroutine check

end program test_flux_label
