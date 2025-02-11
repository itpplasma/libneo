program test_perp_coll_frequency

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use util_for_test
    use libneo_collisions, only: calc_perp_coll_freq, calc_perp_coll_freq_slow_limit, &
    calc_perp_coll_freq_fast_limit, calc_coulomb_log

    implicit none

    print *, "works!"
    !call test_calc_coulomb_log
    !call test_calc_perp_coll_freq
    !call test_calc_perp_coll_freq_slow_limit
    !call test_calc_perp_coll_freq_fast_limit

end program