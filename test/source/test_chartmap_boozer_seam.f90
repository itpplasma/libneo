program test_chartmap_boozer_seam
    ! Boozer-angle chart of the nfp = 3 stellarator fixture (libneo#304):
    ! one field period later the chart gives the rotated point, the covariant
    ! basis is continuous across the period seam, and the position error in
    ! the cells next to the seam is no larger than in the interior and falls
    ! with the grid like it.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use math_constants, only: TWOPI
    use libneo_coordinates_base, only: coordinate_system_t
    use libneo_coordinates_chartmap, only: make_chartmap_coordinate_system
    use chartmap_gate_fixtures, only: write_stellarator_boozer_chartmap, &
                                      stellarator_chart_position
    implicit none
    real(dp), parameter :: period = TWOPI/3
    class(coordinate_system_t), allocatable :: coarse, fine
    real(dp) :: u(3), x1(3), x2(3), rot(3, 3), e1(3, 3), e2(3, 3), worst
    real(dp) :: seam_c, inner_c, seam_f, inner_f
    integer :: i

    call write_stellarator_boozer_chartmap('seam_coarse.nc', 9, 16, 12, 1)
    call write_stellarator_boozer_chartmap('seam_fine.nc', 17, 32, 24, 1)
    call make_chartmap_coordinate_system(coarse, 'seam_coarse.nc')
    call make_chartmap_coordinate_system(fine, 'seam_fine.nc')

    rot = reshape([cos(period), sin(period), 0.0_dp, -sin(period), cos(period), &
                   0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 3])
    worst = 0
    do i = 1, 50
        u = [0.02_dp*i, 0.37_dp*i, 0.11_dp*i]
        call coarse%evaluate_cart(u, x1)
        call coarse%evaluate_cart(u + [0.0_dp, 0.0_dp, period], x2)
        worst = max(worst, maxval(abs(x2 - matmul(rot, x1))))
    end do
    print '(a,es10.3)', 'one period later minus rotated point (cm):', worst
    if (worst > 1.0e-10_dp) error stop 'field-period rotation'

    worst = 0
    do i = 1, 20
        u = [0.05_dp*i, 0.31_dp*i, period]
        call coarse%covariant_basis(u - [0.0_dp, 0.0_dp, 1.0e-9_dp], e1)
        call coarse%covariant_basis(u + [0.0_dp, 0.0_dp, 1.0e-9_dp], e2)
        worst = max(worst, maxval(abs(e2 - e1))/maxval(abs(e1)))
    end do
    print '(a,es10.3)', 'basis jump across the seam (relative):', worst
    if (worst > 1.0e-6_dp) error stop 'basis continuity'

    call cell_errors(coarse, 12, seam_c, inner_c)
    call cell_errors(fine, 24, seam_f, inner_f)
    print '(a,4es10.3)', 'error seam/inner coarse, seam/inner fine (cm):', &
        seam_c, inner_c, seam_f, inner_f
    if (seam_c > 2*inner_c .or. seam_f > 2*inner_f) error stop 'seam cell error'
    if (.not. seam_c > 8*seam_f) error stop 'seam convergence'
contains
    ! Largest position error at cell-centre zeta in the seam cells (first and
    ! last of the period) and in all other cells.
    subroutine cell_errors(chart, nz, seam, inner)
        class(coordinate_system_t), intent(in) :: chart
        integer, intent(in) :: nz
        real(dp), intent(out) :: seam, inner
        real(dp) :: v(3), xc(3), xa(3), e
        integer :: k, j, r
        seam = 0; inner = 0
        do k = 0, nz - 1
            do j = 0, 15
                do r = 1, 4
                    v = [0.2_dp*r, TWOPI*(j + 0.3_dp)/16, period*(k + 0.5_dp)/nz]
                    call chart%evaluate_cart(v, xc)
                    call stellarator_chart_position(v(1), v(2), v(3), 1, xa)
                    e = maxval(abs(xc - xa))
                    if (k == 0 .or. k == nz - 1) then
                        seam = max(seam, e)
                    else
                        inner = max(inner, e)
                    end if
                end do
            end do
        end do
    end subroutine cell_errors
end program test_chartmap_boozer_seam
