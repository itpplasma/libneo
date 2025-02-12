program test_perp_coll_frequency

    use libneo_kinds, only : real_kind
    use libneo_collisions, only: calc_perp_coll_freq_slow_limit_ee, calc_perp_coll_freq_fast_limit_ee

    implicit none

    print *, "works!"
    !call test_calc_coulomb_log
    !call test_calc_perp_coll_freq
    !call test_calc_perp_coll_freq_slow_limit
    !call test_calc_perp_coll_freq_fast_limit
    call test_calc_perp_coll_freq_fast_limit
    !call test_calc_perp_coll_freq_slow_limit

    contains

    subroutine test_calc_perp_coll_freq_fast_limit

        implicit none
        real(kind=real_kind) :: vel, dens_e, coulomb_log
        real(kind=real_kind) :: coll_freq = 0.0d0

        vel = 1.0e8
        dens_e = 1.0e13
        coulomb_log = 15.0

        call calc_perp_coll_freq_fast_limit_ee(vel, dens_e, coulomb_log, coll_freq)
        print *, "coll_freq = ", coll_freq

    end subroutine

end program