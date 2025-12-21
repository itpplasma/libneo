module util_timer
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: iso_c_binding, only: c_int, c_long
    implicit none
    private

    public :: wall_time

    integer(c_int), parameter :: CLOCK_MONOTONIC = 1

    type, bind(c) :: timespec_t
        integer(c_long) :: tv_sec
        integer(c_long) :: tv_nsec
    end type timespec_t

    interface
        integer(c_int) function clock_gettime(clk_id, tp) bind(c, name="clock_gettime")
            import :: c_int, timespec_t
            integer(c_int), value :: clk_id
            type(timespec_t) :: tp
        end function clock_gettime
    end interface

contains

    real(dp) function wall_time() result(t)
        type(timespec_t) :: ts
        integer(c_int) :: ierr

        ierr = clock_gettime(CLOCK_MONOTONIC, ts)
        if (ierr /= 0) then
            error stop "util_timer: clock_gettime failed"
        end if

        t = real(int(ts%tv_sec, int64), dp) + real(int(ts%tv_nsec, int64), dp)*1.0d-9
    end function wall_time

end module util_timer

