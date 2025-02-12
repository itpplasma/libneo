program test_perp_coll_frequency

    use libneo_kinds, only : real_kind
    use libneo_collisions, only: calc_perp_coll_freq_slow_limit_ee, calc_perp_coll_freq_fast_limit_ee, &
        calc_coulomb_log

    implicit none

    call test_calc_coulomb_log
    !call test_calc_perp_coll_freq
    !call test_calc_perp_coll_freq_slow_limit
    !call test_calc_perp_coll_freq_fast_limit
    call test_calc_perp_coll_freq_fast_limit
    !call test_calc_perp_coll_freq_slow_limit

    contains

    subroutine test_calc_coulomb_log

        use math_constants, only: m_e, m_p, m_D
        implicit none

        character(len=2) :: interaction_type
        real(kind=real_kind) :: temp_a, dens_a, temp_b, dens_b
        real(kind=real_kind) :: mass_a, mass_b
        integer :: charge_num_a, charge_num_b
        real(kind=real_kind) :: coulomb_log

        ! electron-electron interaction
        temp_a = 1.0d3
        dens_a = 1.0d13
        mass_a = m_e
        charge_num_a = 0
        temp_b = 0d0
        dens_b = 0d0
        mass_b = 0d0
        charge_num_b = 0

        interaction_type = "ee"
        call calc_coulomb_log(interaction_type, temp_a, dens_a, mass_a, charge_num_a, &
        temp_b, dens_b,  mass_b, charge_num_b, coulomb_log)
        print *, "interaction_type ", interaction_type, ": ", "coulomb_log = ", coulomb_log

        ! electron-ion interaction
        temp_b = 1d3
        dens_b = dens_a
        mass_b = m_D
        charge_num_b = 1

        interaction_type = "ei"
        call calc_coulomb_log(interaction_type, temp_a, dens_a, mass_a, charge_num_a, &
        temp_b, dens_b,  mass_b, charge_num_b, coulomb_log)
        print *, "interaction_type ", interaction_type, ": ", "coulomb_log = ", coulomb_log

        ! ion-ion interaction
        mass_a = m_p
        charge_num_a = 1

        interaction_type = "ii"
        call calc_coulomb_log(interaction_type, temp_a, dens_a, mass_a, charge_num_a, &
        temp_b, dens_b,  mass_b, charge_num_b, coulomb_log)
        print *, "interaction_type ", interaction_type, ": ", "coulomb_log = ", coulomb_log
        

    end subroutine

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