module util
    use iso_fortran_env, only: dp => real64
    
    implicit none

contains

    function linspace(lo, hi, cnt, excl_lo, excl_hi)
        real(dp), intent(in) :: lo, hi
        integer, intent(in) :: cnt
        integer, intent(in), optional :: excl_lo, excl_hi
        real(dp) :: linspace(cnt)
        real(dp) :: step
        integer :: k, omit_lo, omit_hi

        omit_lo = 0
        if (present(excl_lo)) omit_lo = excl_lo
        omit_hi = 0
        if (present(excl_hi)) omit_hi = excl_hi
        step = (hi - lo) / dble(cnt - 1 + omit_lo + omit_hi)
        linspace = lo + [(k * step, k = omit_lo, cnt - 1 + omit_lo)]
        if (omit_hi == 0) linspace(cnt) = hi
    end function linspace

end module util
