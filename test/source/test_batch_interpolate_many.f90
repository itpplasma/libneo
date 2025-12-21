program test_batch_interpolate_many
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use batch_interpolate, only: BatchSplineData1D, BatchSplineData2D, BatchSplineData3D
    use batch_interpolate, only: construct_batch_splines_1d, &
                                 construct_batch_splines_1d_lines, &
                                 construct_batch_splines_2d, &
                                 construct_batch_splines_2d_lines, &
                                 construct_batch_splines_3d
    use batch_interpolate, only: construct_batch_splines_3d_lines
    use batch_interpolate, only: construct_batch_splines_1d_resident, &
                                 construct_batch_splines_1d_resident_device, &
                                 construct_batch_splines_2d_resident, &
                                 construct_batch_splines_2d_resident_device, &
                                 construct_batch_splines_3d_resident, &
                                 construct_batch_splines_3d_resident_device
    use batch_interpolate, only: destroy_batch_splines_1d, destroy_batch_splines_2d, &
                                 destroy_batch_splines_3d
    use batch_interpolate, only: evaluate_batch_splines_1d, evaluate_batch_splines_2d, &
                                 evaluate_batch_splines_3d
    use batch_interpolate, only: evaluate_batch_splines_1d_many, &
                                 evaluate_batch_splines_2d_many, &
                                 evaluate_batch_splines_3d_many
    use batch_interpolate, only: evaluate_batch_splines_1d_many_resident, &
                                 evaluate_batch_splines_2d_many_resident, &
                                 evaluate_batch_splines_3d_many_resident
    implicit none

    type :: lcg_t
        integer(int64) :: state
    end type lcg_t

    call test_1d()
    call test_2d()
    call test_3d()

contains

    pure subroutine lcg_init(rng, seed)
        type(lcg_t), intent(inout) :: rng
        integer(int64), intent(in) :: seed
        rng%state = seed
    end subroutine lcg_init

    real(dp) function lcg_u01(rng) result(r)
        type(lcg_t), intent(inout) :: rng
        integer(int64), parameter :: a = 6364136223846793005_int64
        integer(int64), parameter :: c = 1442695040888963407_int64
        integer(int64) :: x
        rng%state = a*rng%state + c
        x = ishft(rng%state, -11)
        r = real(iand(x, int(z'1FFFFFFFFFFFFF', int64)), dp) / &
            real(int(z'20000000000000', int64), dp)
    end function lcg_u01

    subroutine assert_close(name, a, b, tol)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: a(:)
        real(dp), intent(in) :: b(:)
        real(dp), intent(in) :: tol
        real(dp) :: diff

        diff = maxval(abs(a - b))
        if (diff > tol) then
            write (*, *) "FAIL", trim(name), "diff", diff, "tol", tol
            error stop 1
        end if
    end subroutine assert_close

    subroutine test_1d()
        integer, parameter :: order = 5
        integer, parameter :: order3 = 3
        integer, parameter :: ngrid = 64
        integer, parameter :: nq = 7
        integer, parameter :: npts = 2000
        logical, parameter :: periodic = .true.

        real(dp), parameter :: x_min = 1.23d0
        real(dp), parameter :: x_max = x_min + 4.0d0

        type(BatchSplineData1D) :: spl
        type(BatchSplineData1D) :: spl3
        type(BatchSplineData1D) :: spl_lines
        type(BatchSplineData1D) :: spl3_lines
        real(dp), allocatable :: x_grid(:)
        real(dp), allocatable :: y_grid(:, :)
        real(dp), allocatable :: x_eval(:)
        real(dp), allocatable :: y_many(:, :)
        real(dp), allocatable :: y_many3(:, :)
        real(dp), allocatable :: y_many_res(:, :)
        real(dp), allocatable :: y_many_lines(:, :)
        real(dp), allocatable :: y_one(:)

        type(lcg_t) :: rng
        integer :: ip, iq
        real(dp) :: period
        real(dp), parameter :: tol = 1.0d-12

        allocate (x_grid(ngrid))
        allocate (y_grid(ngrid, nq))
        do ip = 1, ngrid
            x_grid(ip) = x_min + (x_max - x_min)*real(ip - 1, dp)/real(ngrid - 1, dp)
        end do
        do iq = 1, nq
            do ip = 1, ngrid
                y_grid(ip, iq) = cos(x_grid(ip) + 0.1d0*real(iq - 1, dp))
            end do
        end do

        call construct_batch_splines_1d(x_min, x_max, y_grid, order, periodic, spl)

        allocate (x_eval(npts))
        allocate (y_many(nq, npts))
        allocate (y_one(nq))

        call lcg_init(rng, 12345_int64)
        period = (x_max - x_min)
        do ip = 1, npts
            x_eval(ip) = x_min - 2.0d0*period + 5.0d0*period*lcg_u01(rng)
        end do

        call evaluate_batch_splines_1d_many(spl, x_eval, y_many)
        do ip = 1, npts
            call evaluate_batch_splines_1d(spl, x_eval(ip), y_one)
            call assert_close("1d", y_one, y_many(:, ip), tol)
        end do

        call construct_batch_splines_1d_lines(x_min, x_max, y_grid, order, periodic, spl_lines)
        allocate (y_many_lines(nq, npts))
        call evaluate_batch_splines_1d_many(spl_lines, x_eval, y_many_lines)
        call assert_close("1d_lines", reshape(y_many_lines, [nq*npts]), &
                          reshape(y_many, [nq*npts]), tol)
        call destroy_batch_splines_1d(spl_lines)
        deallocate (y_many_lines)

        call destroy_batch_splines_1d(spl)
        call construct_batch_splines_1d_resident(x_min, x_max, y_grid, &
                                                 order, periodic, spl)
        allocate (y_many_res(nq, npts))
        y_many_res = 0.0d0

        !$acc enter data copyin(x_eval) create(y_many_res)
        call evaluate_batch_splines_1d_many_resident(spl, x_eval, y_many_res)
        !$acc update self(y_many_res)
        !$acc exit data delete(y_many_res, x_eval)

        call assert_close("1d_resident", reshape(y_many_res, [nq*npts]), &
                          reshape(y_many, [nq*npts]), tol)

        call destroy_batch_splines_1d(spl)
        call construct_batch_splines_1d_resident_device(x_min, x_max, y_grid, &
                                                        order, periodic, spl, &
                                                        update_host=.true.)
        y_many_res = 0.0d0
        call evaluate_batch_splines_1d_many(spl, x_eval, y_many_res)
        call assert_close("1d_device_host_coeff", reshape(y_many_res, [nq*npts]), &
                          reshape(y_many, [nq*npts]), tol)

        call destroy_batch_splines_1d(spl)
        call construct_batch_splines_1d_resident_device(x_min, x_max, y_grid, &
                                                        order, periodic, spl, &
                                                        update_host=.false.)
        y_many_res = 0.0d0
        !$acc enter data copyin(x_eval) create(y_many_res)
        call evaluate_batch_splines_1d_many_resident(spl, x_eval, y_many_res)
        !$acc update self(y_many_res)
        !$acc exit data delete(y_many_res, x_eval)
        call assert_close("1d_device_resident_only", reshape(y_many_res, [nq*npts]), &
                          reshape(y_many, [nq*npts]), tol)

        call destroy_batch_splines_1d(spl)

        call construct_batch_splines_1d(x_min, x_max, y_grid, order3, periodic, spl3)
        allocate (y_many3(nq, npts))
        call evaluate_batch_splines_1d_many(spl3, x_eval, y_many3)

        call construct_batch_splines_1d_lines(x_min, x_max, y_grid, order3, periodic, &
                                              spl3_lines)
        allocate (y_many_lines(nq, npts))
        call evaluate_batch_splines_1d_many(spl3_lines, x_eval, y_many_lines)
        call assert_close("1d_lines_order3", reshape(y_many_lines, [nq*npts]), &
                          reshape(y_many3, [nq*npts]), tol)
        call destroy_batch_splines_1d(spl3_lines)
        deallocate (y_many_lines)

        call destroy_batch_splines_1d(spl3)
        call construct_batch_splines_1d_resident_device(x_min, x_max, y_grid, order3, &
                                                        periodic, spl3, update_host=.false.)
        y_many_res = 0.0d0
        !$acc enter data copyin(x_eval) create(y_many_res)
        call evaluate_batch_splines_1d_many_resident(spl3, x_eval, y_many_res)
        !$acc update self(y_many_res)
        !$acc exit data delete(y_many_res, x_eval)
        call assert_close("1d_device_resident_only_order3", reshape(y_many_res, [nq*npts]), &
                          reshape(y_many3, [nq*npts]), tol)

        call destroy_batch_splines_1d(spl3)
    end subroutine test_1d

    subroutine test_2d()
        integer, parameter :: order(2) = [5, 5]
        integer, parameter :: ngrid(2) = [32, 24]
        integer, parameter :: nq = 5
        integer, parameter :: npts = 2000
        logical, parameter :: periodic(2) = [.true., .false.]

        real(dp), parameter :: x_min(2) = [1.23d0, -0.7d0]
        real(dp), parameter :: x_max(2) = [x_min(1) + 3.0d0, x_min(2) + 2.0d0]

        type(BatchSplineData2D) :: spl
        type(BatchSplineData2D) :: spl_lines
        real(dp), allocatable :: y_grid(:, :, :)
        real(dp), allocatable :: x_eval(:, :)
        real(dp), allocatable :: y_many(:, :)
        real(dp), allocatable :: y_many_res(:, :)
        real(dp), allocatable :: y_many_lines(:, :)
        real(dp), allocatable :: y_one(:)

        type(lcg_t) :: rng
        integer :: i1, i2, iq, ip
        real(dp), parameter :: tol = 1.0d-12

        allocate (y_grid(ngrid(1), ngrid(2), nq))
        do iq = 1, nq
            do i2 = 1, ngrid(2)
                do i1 = 1, ngrid(1)
                    y_grid(i1, i2, iq) = &
                        cos(x_min(1) + real(i1 - 1, dp)*(x_max(1) - x_min(1)) / &
                            real(ngrid(1) - 1, dp)) * &
                        cos(x_min(2) + real(i2 - 1, dp)*(x_max(2) - x_min(2)) / &
                            real(ngrid(2) - 1, dp)) * &
                        (1.0d0 + 0.05d0*real(iq - 1, dp))
                end do
            end do
        end do

        call construct_batch_splines_2d(x_min, x_max, y_grid, order, periodic, spl)

        allocate (x_eval(2, npts))
        allocate (y_many(nq, npts))
        allocate (y_one(nq))

        call lcg_init(rng, 67890_int64)
        do ip = 1, npts
            x_eval(1, ip) = x_min(1) - 2.0d0*(x_max(1) - x_min(1)) + &
                            5.0d0*(x_max(1) - x_min(1))*lcg_u01(rng)
            x_eval(2, ip) = x_min(2) - 1.0d0*(x_max(2) - x_min(2)) + &
                            3.0d0*(x_max(2) - x_min(2))*lcg_u01(rng)
        end do

        call evaluate_batch_splines_2d_many(spl, x_eval, y_many)
        do ip = 1, npts
            call evaluate_batch_splines_2d(spl, x_eval(:, ip), y_one)
            call assert_close("2d", y_one, y_many(:, ip), tol)
        end do

        call construct_batch_splines_2d_lines(x_min, x_max, y_grid, order, periodic, spl_lines)
        allocate (y_many_lines(nq, npts))
        call evaluate_batch_splines_2d_many(spl_lines, x_eval, y_many_lines)
        call assert_close("2d_lines", reshape(y_many_lines, [nq*npts]), &
                          reshape(y_many, [nq*npts]), tol)
        call destroy_batch_splines_2d(spl_lines)
        deallocate (y_many_lines)

        call destroy_batch_splines_2d(spl)
        call construct_batch_splines_2d_resident(x_min, x_max, y_grid, &
                                                 order, periodic, spl)
        allocate (y_many_res(nq, npts))
        y_many_res = 0.0d0

        !$acc enter data copyin(x_eval) create(y_many_res)
        call evaluate_batch_splines_2d_many_resident(spl, x_eval, y_many_res)
        !$acc update self(y_many_res)
        !$acc exit data delete(y_many_res, x_eval)

        call assert_close("2d_resident", reshape(y_many_res, [nq*npts]), &
                          reshape(y_many, [nq*npts]), tol)

        call destroy_batch_splines_2d(spl)
        call construct_batch_splines_2d_resident_device(x_min, x_max, y_grid, &
                                                        order, periodic, spl, &
                                                        update_host=.true.)
        y_many_res = 0.0d0
        call evaluate_batch_splines_2d_many(spl, x_eval, y_many_res)
        call assert_close("2d_device_host_coeff", reshape(y_many_res, [nq*npts]), &
                          reshape(y_many, [nq*npts]), tol)

        call destroy_batch_splines_2d(spl)
        call construct_batch_splines_2d_resident_device(x_min, x_max, y_grid, &
                                                        order, periodic, spl, &
                                                        update_host=.false.)
        y_many_res = 0.0d0
        !$acc enter data copyin(x_eval) create(y_many_res)
        call evaluate_batch_splines_2d_many_resident(spl, x_eval, y_many_res)
        !$acc update self(y_many_res)
        !$acc exit data delete(y_many_res, x_eval)
        call assert_close("2d_device_resident_only", reshape(y_many_res, [nq*npts]), &
                          reshape(y_many, [nq*npts]), tol)

        call destroy_batch_splines_2d(spl)
    end subroutine test_2d

    subroutine test_3d()
        integer, parameter :: order(3) = [5, 3, 3]
        integer, parameter :: order_dev(3) = [5, 5, 5]
        integer, parameter :: ngrid(3) = [16, 12, 10]
        integer, parameter :: nq = 4
        integer, parameter :: npts = 800
        logical, parameter :: periodic(3) = [.true., .true., .true.]

        real(dp), parameter :: x_min(3) = [1.23d0, -0.7d0, 0.5d0]
        real(dp), parameter :: x_max(3) = [x_min(1) + 2.0d0, x_min(2) + 1.5d0, &
                                           x_min(3) + 1.0d0]

        type(BatchSplineData3D) :: spl
        type(BatchSplineData3D) :: spl_lines
        real(dp), allocatable :: y_grid(:, :, :, :)
        real(dp), allocatable :: x_eval(:, :)
        real(dp), allocatable :: y_many(:, :)
        real(dp), allocatable :: y_many_res(:, :)
        real(dp), allocatable :: y_many_lines(:, :)
        real(dp), allocatable :: y_one(:)

        type(lcg_t) :: rng
        integer :: i1, i2, i3, iq, ip
        real(dp), parameter :: tol = 1.0d-12

        allocate (y_grid(ngrid(1), ngrid(2), ngrid(3), nq))
        do iq = 1, nq
            do i3 = 1, ngrid(3)
                do i2 = 1, ngrid(2)
                    do i1 = 1, ngrid(1)
                        y_grid(i1, i2, i3, iq) = &
                            cos(x_min(1) + real(i1 - 1, dp)*(x_max(1) - x_min(1)) / &
                                real(ngrid(1) - 1, dp)) * &
                            cos(x_min(2) + real(i2 - 1, dp)*(x_max(2) - x_min(2)) / &
                                real(ngrid(2) - 1, dp)) * &
                            cos(x_min(3) + real(i3 - 1, dp)*(x_max(3) - x_min(3)) / &
                                real(ngrid(3) - 1, dp)) * &
                            (1.0d0 + 0.03d0*real(iq - 1, dp))
                    end do
                end do
            end do
        end do

        call construct_batch_splines_3d(x_min, x_max, y_grid, order, periodic, spl)

        allocate (x_eval(3, npts))
        allocate (y_many(nq, npts))
        allocate (y_one(nq))

        call lcg_init(rng, 13579_int64)
        do ip = 1, npts
            x_eval(1, ip) = x_min(1) - 2.0d0*(x_max(1) - x_min(1)) + &
                            5.0d0*(x_max(1) - x_min(1))*lcg_u01(rng)
            x_eval(2, ip) = x_min(2) - 2.0d0*(x_max(2) - x_min(2)) + &
                            5.0d0*(x_max(2) - x_min(2))*lcg_u01(rng)
            x_eval(3, ip) = x_min(3) - 2.0d0*(x_max(3) - x_min(3)) + &
                            5.0d0*(x_max(3) - x_min(3))*lcg_u01(rng)
        end do

        call evaluate_batch_splines_3d_many(spl, x_eval, y_many)
        do ip = 1, npts
            call evaluate_batch_splines_3d(spl, x_eval(:, ip), y_one)
            call assert_close("3d", y_one, y_many(:, ip), tol)
        end do

        call construct_batch_splines_3d_lines(x_min, x_max, y_grid, order, periodic, spl_lines)
        allocate (y_many_lines(nq, npts))
        call evaluate_batch_splines_3d_many(spl_lines, x_eval, y_many_lines)
        call assert_close("3d_lines", reshape(y_many_lines, [nq*npts]), &
                          reshape(y_many, [nq*npts]), tol)
        call destroy_batch_splines_3d(spl_lines)
        deallocate (y_many_lines)

        call destroy_batch_splines_3d(spl)
        call construct_batch_splines_3d_resident(x_min, x_max, y_grid, &
                                                 order, periodic, spl)
        allocate (y_many_res(nq, npts))
        y_many_res = 0.0d0

        !$acc enter data copyin(x_eval) create(y_many_res)
        call evaluate_batch_splines_3d_many_resident(spl, x_eval, y_many_res)
        !$acc update self(y_many_res)
        !$acc exit data delete(y_many_res, x_eval)

        call assert_close("3d_resident", reshape(y_many_res, [nq*npts]), &
                          reshape(y_many, [nq*npts]), tol)

        call destroy_batch_splines_3d(spl)

        y_many_res = 0.0d0
        call construct_batch_splines_3d(x_min, x_max, y_grid, order_dev, periodic, spl)
        call evaluate_batch_splines_3d_many(spl, x_eval, y_many_res)
        call destroy_batch_splines_3d(spl)

        call construct_batch_splines_3d_resident_device(x_min, x_max, y_grid, &
                                                        order_dev, periodic, spl, &
                                                        update_host=.true.)
        y_many = 0.0d0
        call evaluate_batch_splines_3d_many(spl, x_eval, y_many)
        call assert_close("3d_device_host_coeff", reshape(y_many, [nq*npts]), &
                          reshape(y_many_res, [nq*npts]), tol)
        call destroy_batch_splines_3d(spl)

        call construct_batch_splines_3d_resident_device(x_min, x_max, y_grid, &
                                                        order_dev, periodic, spl, &
                                                        update_host=.false.)
        y_many = 0.0d0
        !$acc enter data copyin(x_eval) create(y_many)
        call evaluate_batch_splines_3d_many_resident(spl, x_eval, y_many)
        !$acc update self(y_many)
        !$acc exit data delete(y_many, x_eval)
        call assert_close("3d_device_resident_only", reshape(y_many, [nq*npts]), &
                          reshape(y_many_res, [nq*npts]), tol)
        call destroy_batch_splines_3d(spl)
    end subroutine test_3d

end program test_batch_interpolate_many
