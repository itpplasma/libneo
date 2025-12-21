module util_lcg_rng
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    implicit none
    private

    public :: lcg_state_t
    public :: lcg_init
    public :: lcg_uniform_0_1

    type :: lcg_state_t
        integer(int64) :: state
    end type lcg_state_t

contains

    pure subroutine lcg_init(rng, seed)
        type(lcg_state_t), intent(inout) :: rng
        integer(int64), intent(in) :: seed

        if (seed == 0_int64) then
            rng%state = 1_int64
        else
            rng%state = seed
        end if
    end subroutine lcg_init

    real(dp) function lcg_uniform_0_1(rng) result(u)
        type(lcg_state_t), intent(inout) :: rng
        integer(int64) :: x
        integer(int64), parameter :: a = 6364136223846793005_int64
        integer(int64), parameter :: c = 1442695040888963407_int64

        x = a*rng%state + c
        rng%state = x
        u = real(ior(ishft(x, -11), 1_int64), dp)*(1.0d0/real(huge(1_int64), dp))
    end function lcg_uniform_0_1

end module util_lcg_rng
