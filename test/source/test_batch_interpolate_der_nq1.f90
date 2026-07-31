program test_batch_interpolate_der_nq1
    ! Guards the single-quantity specialisation of the 3D first-derivative
    ! batch-spline evaluator (evaluate_batch_splines_3d_der_core_nq1).
    !
    ! The specialisation exists for speed: the general kernel indexes coefficient
    ! arrays whose leading dimension is MAX_QUANTITIES, so with num_quantities = 1
    ! every sweep is strided and the !$omp simd loops run a single iteration.
    ! Correctness therefore has to be pinned down by something other than the
    ! kernel itself.
    !
    ! Two independent oracles are used, neither of which is a recording of this
    ! code's own output:
    !
    !   1. Dispatch equivalence, bitwise. The same underlying data is splined
    !      twice: once as a 1-quantity spline (which takes the nq1 kernel) and
    !      once as the first quantity of a 2-quantity spline (which takes the
    !      general kernel). The nq1 kernel performs the same Horner recurrences
    !      in the same order, so the two must agree to the last bit, not merely
    !      to a tolerance. A tolerance check here would hide exactly the
    !      reassociation bugs that break SIMPLE's byte-exact golden records.
    !
    !   2. Cross-kernel agreement. Value and first derivatives from
    !      evaluate_batch_splines_3d_der must match those from
    !      evaluate_batch_splines_3d_der2, which reaches them through a
    !      completely separate kernel family, and the value must match
    !      evaluate_batch_splines_3d.

    use, intrinsic :: iso_fortran_env, only: dp => real64, int64, error_unit
    use batch_interpolate, only: BatchSplineData3D, construct_batch_splines_3d, &
                                 destroy_batch_splines_3d, &
                                 evaluate_batch_splines_3d, &
                                 evaluate_batch_splines_3d_der, &
                                 evaluate_batch_splines_3d_der2
    implicit none

    integer :: nfail

    nfail = 0
    call check_dispatch_equivalence_bitwise(nfail)
    call check_against_der2(nfail)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "test_batch_interpolate_der_nq1: all tests passed"

contains

    real(dp) function sample_field(x1, x2, x3, iq) result(v)
        real(dp), intent(in) :: x1, x2, x3
        integer, intent(in) :: iq
        v = cos(x1)*cos(x2)*cos(x3) + 0.25_dp*sin(2.0_dp*x1)*cos(x3) &
            + 0.1_dp*real(iq - 1, dp)*sin(x2)
    end function sample_field

    subroutine build_grid(ngrid, x_min, x_max, nq, y_grid)
        integer, intent(in) :: ngrid(3), nq
        real(dp), intent(in) :: x_min(3), x_max(3)
        real(dp), allocatable, intent(out) :: y_grid(:, :, :, :)

        integer :: i1, i2, i3, iq
        real(dp) :: p1, p2, p3

        allocate (y_grid(ngrid(1), ngrid(2), ngrid(3), nq))
        do iq = 1, nq
            do i3 = 1, ngrid(3)
                p3 = x_min(3) + real(i3 - 1, dp)*(x_max(3) - x_min(3)) / &
                     real(ngrid(3) - 1, dp)
                do i2 = 1, ngrid(2)
                    p2 = x_min(2) + real(i2 - 1, dp)*(x_max(2) - x_min(2)) / &
                         real(ngrid(2) - 1, dp)
                    do i1 = 1, ngrid(1)
                        p1 = x_min(1) + real(i1 - 1, dp)*(x_max(1) - x_min(1)) / &
                             real(ngrid(1) - 1, dp)
                        y_grid(i1, i2, i3, iq) = sample_field(p1, p2, p3, iq)
                    end do
                end do
            end do
        end do
    end subroutine build_grid

    ! Deterministic sample points spanning the periodic wrap in every direction,
    ! so the modulo branch of the kernel is exercised as well as the interior.
    subroutine eval_point(ip, npts, x_min, x_max, x)
        integer, intent(in) :: ip, npts
        real(dp), intent(in) :: x_min(3), x_max(3)
        real(dp), intent(out) :: x(3)

        real(dp) :: u
        integer :: j

        do j = 1, 3
            u = mod(real(ip, dp)*(0.6180339887498949_dp + 0.11_dp*real(j, dp)), 1.0_dp)
            ! Deliberately overshoot the domain on both sides.
            x(j) = x_min(j) - (x_max(j) - x_min(j)) + &
                   3.0_dp*(x_max(j) - x_min(j))*u
        end do
        associate (unused => npts); end associate
    end subroutine eval_point

    subroutine check_dispatch_equivalence_bitwise(nfail)
        integer, intent(inout) :: nfail

        integer, parameter :: order(3) = [5, 5, 5]
        integer, parameter :: ngrid(3) = [14, 12, 11]
        integer, parameter :: npts = 500
        logical, parameter :: periodic(3) = [.true., .true., .true.]
        real(dp), parameter :: x_min(3) = [0.0_dp, 0.0_dp, 0.0_dp]
        real(dp), parameter :: x_max(3) = [1.7_dp, 2.3_dp, 1.1_dp]

        type(BatchSplineData3D) :: spl1, spl2
        real(dp), allocatable :: y1(:, :, :, :), y2(:, :, :, :)
        real(dp) :: x(3)
        real(dp) :: v1(1), dv1(3, 1)
        real(dp) :: v2(2), dv2(3, 2)
        integer :: ip, nbad

        call build_grid(ngrid, x_min, x_max, 1, y1)
        call build_grid(ngrid, x_min, x_max, 2, y2)
        ! Quantity 1 of the 2-quantity grid is the same data as the 1-quantity
        ! grid (sample_field only varies with iq through a term that vanishes
        ! for iq = 1), so the splines agree on that quantity by construction.
        call construct_batch_splines_3d(x_min, x_max, y1, order, periodic, spl1)
        call construct_batch_splines_3d(x_min, x_max, y2, order, periodic, spl2)

        nbad = 0
        do ip = 1, npts
            call eval_point(ip, npts, x_min, x_max, x)
            call evaluate_batch_splines_3d_der(spl1, x, v1, dv1)
            call evaluate_batch_splines_3d_der(spl2, x, v2, dv2)
            if (v1(1) /= v2(1)) nbad = nbad + 1
            if (dv1(1, 1) /= dv2(1, 1)) nbad = nbad + 1
            if (dv1(2, 1) /= dv2(2, 1)) nbad = nbad + 1
            if (dv1(3, 1) /= dv2(3, 1)) nbad = nbad + 1
        end do

        if (nbad > 0) then
            write (error_unit, "(a,i0,a)") &
                "  nq1 kernel differs from general kernel in ", nbad, &
                " of 2000 compared values (must be bitwise identical)"
            nfail = nfail + 1
        else
            write (*, "(a)") "  dispatch equivalence (nq1 vs general): bitwise same"
        end if

        call destroy_batch_splines_3d(spl1)
        call destroy_batch_splines_3d(spl2)
    end subroutine check_dispatch_equivalence_bitwise

    subroutine check_against_der2(nfail)
        integer, intent(inout) :: nfail

        integer, parameter :: order(3) = [5, 5, 5]
        integer, parameter :: ngrid(3) = [14, 12, 11]
        integer, parameter :: npts = 500
        logical, parameter :: periodic(3) = [.true., .true., .true.]
        real(dp), parameter :: x_min(3) = [0.0_dp, 0.0_dp, 0.0_dp]
        real(dp), parameter :: x_max(3) = [1.7_dp, 2.3_dp, 1.1_dp]
        real(dp), parameter :: tol = 1.0e-12_dp

        type(BatchSplineData3D) :: spl
        real(dp), allocatable :: y(:, :, :, :)
        real(dp) :: x(3)
        real(dp) :: v(1), dv(3, 1)
        real(dp) :: v2(1), dv2(3, 1), d2v2(6, 1)
        real(dp) :: vplain(1)
        real(dp) :: worst_val, worst_der, scal
        integer :: ip, j

        call build_grid(ngrid, x_min, x_max, 1, y)
        call construct_batch_splines_3d(x_min, x_max, y, order, periodic, spl)

        worst_val = 0.0_dp
        worst_der = 0.0_dp
        do ip = 1, npts
            call eval_point(ip, npts, x_min, x_max, x)
            call evaluate_batch_splines_3d_der(spl, x, v, dv)
            call evaluate_batch_splines_3d_der2(spl, x, v2, dv2, d2v2)
            call evaluate_batch_splines_3d(spl, x, vplain)

            worst_val = max(worst_val, abs(v(1) - v2(1)))
            worst_val = max(worst_val, abs(v(1) - vplain(1)))
            do j = 1, 3
                scal = max(1.0_dp, abs(dv2(j, 1)))
                worst_der = max(worst_der, abs(dv(j, 1) - dv2(j, 1))/scal)
            end do
        end do

        write (*, "(a,es10.3,a,es10.3)") "  vs der2/plain: worst value diff ", &
            worst_val, "   worst rel derivative diff ", worst_der

        if (.not. (worst_val < tol)) then
            write (error_unit, "(a)") "  value disagrees with der2 or plain evaluator"
            nfail = nfail + 1
        end if
        if (.not. (worst_der < tol)) then
            write (error_unit, "(a)") "  first derivatives disagree with der2"
            nfail = nfail + 1
        end if

        call destroy_batch_splines_3d(spl)
    end subroutine check_against_der2

end program test_batch_interpolate_der_nq1
