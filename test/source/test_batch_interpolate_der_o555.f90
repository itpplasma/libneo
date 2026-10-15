program test_batch_interpolate_der_o555
    ! Guards the quintic multi-quantity specialisation of the 3D first-derivative
    ! batch-spline evaluator (evaluate_batch_splines_3d_der_core_o555).
    !
    ! The specialisation exists for speed: SIMPLE's SPECTRE full orbit spends
    ! most of its time in quintic four-quantity first-derivative evaluations.
    ! Two oracles, neither a recording of this kernel's output:
    !
    !   1. Agreement with the general reference kernel, which performs the same
    !      Horner recurrences in the same order for any degree: bitwise when the
    !      test is run with the argument "bitwise" (LIBNEO_DETERMINISTIC_FP
    !      builds, where FMA contraction is off), else to 1e-13 relative.
    !   2. Value and first derivatives agree with evaluate_batch_splines_3d_der2,
    !      a separate kernel family, and the value with evaluate_batch_splines_3d.
    !
    ! Cases cover 4 and MAX_QUANTITIES quantities, all-periodic axes and the
    ! non-periodic first axis used by SPECTRE's per-volume radial coordinate.

    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use batch_interpolate, only: BatchSplineData3D, construct_batch_splines_3d, &
                                 destroy_batch_splines_3d, &
                                 evaluate_batch_splines_3d, &
                                 evaluate_batch_splines_3d_der, &
                                 evaluate_batch_splines_3d_der_general, &
                                 evaluate_batch_splines_3d_der2
    implicit none

    integer :: nfail
    logical :: bitwise
    character(len=16) :: arg

    call get_command_argument(1, arg)
    bitwise = trim(arg) == "bitwise"
    nfail = 0
    call check_case(4, [.true., .true., .true.], nfail)
    call check_case(4, [.false., .true., .true.], nfail)
    call check_case(8, [.false., .true., .true.], nfail)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "test_batch_interpolate_der_o555: all tests passed"

contains

    real(dp) function sample_field(x1, x2, x3, iq) result(v)
        real(dp), intent(in) :: x1, x2, x3
        integer, intent(in) :: iq
        v = cos(x1 + 0.1_dp*iq)*cos(x2)*cos(x3) + 0.25_dp*sin(2.0_dp*x1)*cos(x3) &
            + 0.1_dp*real(iq - 1, dp)*sin(x2 + x3)
    end function sample_field

    subroutine check_case(nq, periodic, nfail)
        integer, intent(in) :: nq
        logical, intent(in) :: periodic(3)
        integer, intent(inout) :: nfail

        integer, parameter :: order(3) = [5, 5, 5]
        integer, parameter :: ngrid(3) = [14, 12, 11]
        integer, parameter :: npts = 500
        real(dp), parameter :: x_min(3) = [0.0_dp, 0.0_dp, 0.0_dp]
        real(dp), parameter :: x_max(3) = [1.7_dp, 2.3_dp, 1.1_dp]
        real(dp), parameter :: tol = 1.0e-12_dp

        type(BatchSplineData3D) :: spl
        real(dp), allocatable :: y(:, :, :, :)
        real(dp) :: x(3), p(3), u, scal, worst_val, worst_der, worst_ref
        real(dp) :: v(nq), dv(3, nq), vr(nq), dvr(3, nq)
        real(dp) :: v2(nq), dv2(3, nq), d2v2(6, nq), vplain(nq)
        integer :: i1, i2, i3, iq, ip, j, nbad

        allocate (y(ngrid(1), ngrid(2), ngrid(3), nq))
        do iq = 1, nq
            do i3 = 1, ngrid(3)
                do i2 = 1, ngrid(2)
                    do i1 = 1, ngrid(1)
                        p = x_min + real([i1, i2, i3] - 1, dp)*(x_max - x_min) &
                            /real(ngrid - 1, dp)
                        y(i1, i2, i3, iq) = sample_field(p(1), p(2), p(3), iq)
                    end do
                end do
            end do
        end do
        call construct_batch_splines_3d(x_min, x_max, y, order, periodic, spl)

        nbad = 0
        worst_ref = 0.0_dp
        worst_val = 0.0_dp
        worst_der = 0.0_dp
        do ip = 1, npts
            do j = 1, 3
                u = mod(real(ip, dp)*(0.6180339887498949_dp + 0.11_dp*real(j, dp)), &
                        1.0_dp)
                if (periodic(j)) then
                    x(j) = x_min(j) - (x_max(j) - x_min(j)) + &
                           3.0_dp*(x_max(j) - x_min(j))*u
                else
                    x(j) = x_min(j) + (x_max(j) - x_min(j))*u
                end if
            end do
            call evaluate_batch_splines_3d_der(spl, x, v, dv)
            call evaluate_batch_splines_3d_der_general(spl, x, vr, dvr)
            nbad = nbad + count(v /= vr) + count(dv /= dvr)
            worst_ref = max(worst_ref, maxval(abs(v - vr)/max(1.0_dp, abs(vr))), &
                            maxval(abs(dv - dvr)/max(1.0_dp, abs(dvr))))

            call evaluate_batch_splines_3d_der2(spl, x, v2, dv2, d2v2)
            call evaluate_batch_splines_3d(spl, x, vplain)
            worst_val = max(worst_val, maxval(abs(v - v2)), maxval(abs(v - vplain)))
            do iq = 1, nq
                do j = 1, 3
                    scal = max(1.0_dp, abs(dv2(j, iq)))
                    worst_der = max(worst_der, abs(dv(j, iq) - dv2(j, iq))/scal)
                end do
            end do
        end do

        write (*, "(a,i0,a,3l2,a,i0,a,es10.3,a,es10.3)") "  nq=", nq, &
            " periodic=", periodic, ": mismatches vs reference ", nbad, &
            ", worst value diff ", worst_val, ", worst rel derivative diff ", &
            worst_der
        if ((bitwise .and. nbad > 0) .or. .not. (worst_ref < tol)) then
            write (error_unit, "(a,es10.3)") &
                "  differs from the general reference kernel, worst rel ", worst_ref
            nfail = nfail + 1
        end if
        if (.not. (worst_val < tol .and. worst_der < tol)) then
            write (error_unit, "(a)") "  disagrees with der2 or plain evaluator"
            nfail = nfail + 1
        end if

        call destroy_batch_splines_3d(spl)
    end subroutine check_case

end program test_batch_interpolate_der_o555
