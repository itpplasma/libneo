program test_perp_coll_frequency

    use libneo_kinds, only : real_kind
    use libneo_collisions, only: calc_perp_coll_freq_slow_limit_ee, calc_perp_coll_freq_fast_limit_ee, &
        calc_coulomb_log
    use util_for_test, only: print_test, print_ok, print_fail

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

        call print_test("test_calc_coulomb_log")

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
        if (coulomb_log - 15.940948099344785 > 1e-6) then ! number is for the parameters above
            call print_fail
            stop "coulomb_log != 15.940948099344785"
        else
            call print_ok
        end if

        ! electron-ion interaction (here, the formula is conditional, c.f. NRL formulary)
        temp_b = 1.0d3
        dens_b = dens_a
        mass_b = m_D
        charge_num_b = 1

        interaction_type = "ei"

        ! case 1:
        temp_a = 9.0d0
        call calc_coulomb_log(interaction_type, temp_a, dens_a, mass_a, charge_num_a, &
        temp_b, dens_b,  mass_b, charge_num_b, coulomb_log)
        print *, "interaction_type ", interaction_type, ", case 1: ", "coulomb_log = ", coulomb_log
        if (coulomb_log - 15.9409521745208 > 1e-6) then ! number is for the parameters above
            call print_fail
            stop "coulomb_log != 15.9409521745208"
        else
            call print_ok
        end if

        ! case 2:
        temp_a = 1.0d3
        call calc_coulomb_log(interaction_type, temp_a, dens_a, mass_a, charge_num_a, &
        temp_b, dens_b,  mass_b, charge_num_b, coulomb_log)
        print *, "interaction_type ", interaction_type, ", case 2: ", "coulomb_log = ", coulomb_log
        if (coulomb_log - 15.9409521745208 > 1e-6) then ! number is for the parameters above
            call print_fail
            stop "coulomb_log != 15.9409521745208"
        else
            call print_ok
        end if

        ! case 3
        temp_a = 0.1d0
        
        call calc_coulomb_log(interaction_type, temp_a, dens_a, mass_a, charge_num_a, &
        temp_b, dens_b,  mass_b, charge_num_b, coulomb_log)
        print *, "interaction_type ", interaction_type, ", case 3: ", "coulomb_log = ", coulomb_log
        if (coulomb_log - 10.7021790072942 > 1e-6) then ! number is for the parameters above
            call print_fail
            stop "coulomb_log != 10.7021790072942"
        else
            call print_ok
        end if

        ! ion-ion interaction
        temp_a = 1.0d3
        mass_a = m_p
        charge_num_a = 1
        mass_a = m_D
        mass_b = m_p

        interaction_type = "ii"
        call calc_coulomb_log(interaction_type, temp_a, dens_a, mass_a, charge_num_a, &
        temp_b, dens_b,  mass_b, charge_num_b, coulomb_log)
        print *, "interaction_type ", interaction_type, ": ", "coulomb_log = ", coulomb_log
        if (coulomb_log - 18.0482562237319 > 1e-6) then ! number is for the parameters above
            call print_fail
            stop "coulomb_log != 18.0482562237319"
        else
            call print_ok
        end if

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