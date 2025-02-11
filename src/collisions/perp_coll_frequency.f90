module libneo_collisions

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_kinds, only : real_kind, complex_kind
    implicit none


    contains

    subroutine calc_coulomb_log()

        implicit none

    end subroutine

    subroutine calc_perp_coll_freq(vel)

        implicit none
        real(kind=real_kind), intent(in) :: vel

        print *, vel

    end subroutine

    subroutine calc_perp_coll_freq_slow_limit(vel)

        implicit none
        real(kind=real_kind), intent(in) :: vel

        print *, vel

    end subroutine

    subroutine calc_perp_coll_freq_fast_limit(vel)

        implicit none
        real(kind=real_kind), intent(in) :: vel

        print *, vel

    end subroutine

end module