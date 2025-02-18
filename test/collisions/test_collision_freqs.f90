program test_collision_freqs

    use libneo_kinds, only : real_kind
    use libneo_collisions, only: calc_perp_coll_freq_slow_limit_ee, calc_perp_coll_freq_fast_limit_ee, &
        calc_coulomb_log, calc_perp_coll_freq
    use util_for_test, only: print_test, print_ok, print_fail

    implicit none

    call test_calc_coulomb_log
    call test_calc_perp_coll_freq
    !call test_calc_perp_coll_freq_slow_limit
    call test_calc_perp_coll_freq_fast_limit
    call test_fill_species_arr_coulomb_log

    contains

    subroutine test_calc_coulomb_log

        use math_constants, only: m_e, m_p, m_D
        use libneo_species, only: species_t

        implicit none

        type(species_t) :: species_a, species_b

        character(len=2) :: interaction_type
        real(kind=real_kind) :: coulomb_log

        call print_test("test_calc_coulomb_log")

        ! electron-electron interaction
        species_a%temp = 1.0d3
        species_a%dens = 1.0d13
        species_a%mass = m_e
        species_a%charge_num = 0
        species_b%temp = 0d0
        species_b%dens = 0d0
        species_b%mass = 0d0
        species_b%charge_num = 0

        interaction_type = "ee"
        call calc_coulomb_log(interaction_type, species_a, species_b, coulomb_log)
        print *, "interaction_type ", interaction_type, ": ", "coulomb_log = ", coulomb_log
        if (abs(coulomb_log - 15.940948099344785) > 1e-6) then ! number is for the parameters above
            call print_fail
            stop "coulomb_log != 15.940948099344785"
        else
            call print_ok
        end if

        ! electron-ion interaction (here, the formula is conditional, c.f. NRL formulary)
        species_b%temp = 1.0d3
        species_b%dens = species_a%dens
        species_b%mass = m_D
        species_b%charge_num = 1

        interaction_type = "ei"

        ! case 1:
        species_a%temp = 9.0d0
        call calc_coulomb_log(interaction_type, species_a, species_b, coulomb_log)
        print *, "interaction_type ", interaction_type, ", case 1: ", "coulomb_log = ", coulomb_log
        if (abs(coulomb_log - 11.329033761543) > 1e-6) then ! number is for the parameters above
            call print_fail
            stop "coulomb_log != 11.329033761543"
        else
            call print_ok
        end if

        ! case 2:
        species_a%temp = 1.0d3
        call calc_coulomb_log(interaction_type, species_a, species_b, coulomb_log)
        print *, "interaction_type ", interaction_type, ", case 2: ", "coulomb_log = ", coulomb_log
        if (abs(coulomb_log - 15.9409521745208) > 1e-6) then ! number is for the parameters above
            call print_fail
            stop "coulomb_log != 15.9409521745208"
        else
            call print_ok
        end if

        ! case 3
        species_a%temp = 0.1d0
        
        call calc_coulomb_log(interaction_type, species_a, species_b, coulomb_log)
        print *, "interaction_type ", interaction_type, ", case 3: ", "coulomb_log = ", coulomb_log
        if (abs(coulomb_log - 10.7021790072942) > 1e-6) then ! number is for the parameters above
            call print_fail
            stop "coulomb_log != 10.7021790072942"
        else
            call print_ok
        end if

        ! ion-ion interaction
        species_a%temp = 1.0d3
        species_a%mass = m_p
        species_a%charge_num = 1
        species_a%mass = m_D
        species_b%mass = m_p

        interaction_type = "ii"
        call calc_coulomb_log(interaction_type, species_a, species_b, coulomb_log)
        print *, "interaction_type ", interaction_type, ": ", "coulomb_log = ", coulomb_log
        if (abs(coulomb_log - 18.0482562237319) > 1e-6) then ! number is for the parameters above
            call print_fail
            stop "coulomb_log != 18.0482562237319"
        else
            call print_ok
        end if

    end subroutine

    subroutine test_fill_species_arr_coulomb_log

        use libneo_species, only: species_t, init_deuterium_plasma
        use libneo_collisions, only: fill_species_arr_coulomb_log

        implicit none

        integer, parameter :: num_species = 2
        type(species_t) :: species_arr(num_species)
        integer :: i, j
        real(kind=real_kind), dimension(2,2) :: coulomb_log_comparison

        coulomb_log_comparison = reshape([15.940948099344892, 15.940952174520840, 15.940952174520840, 18.048256223731936], [2, 2])

        call print_test("test_fill_species_arr_coulomb_log")

        call init_deuterium_plasma(1d3, 1d3, 1d13, species_arr)

        call fill_species_arr_coulomb_log(num_species, species_arr)

        do i=1, num_species
            do j=1, num_species
                print *, "coulomb_log(", i, ",", j, ") = ", species_arr(i)%coulomb_log(j), &
                    " compared to ", coulomb_log_comparison(i, j)
                if (abs(species_arr(i)%coulomb_log(j) - coulomb_log_comparison(i,j)) > 1e-6) then ! number is for the parameters above
                    call print_fail
                    stop "coulomb_log not correct"
                else
                    call print_ok
                end if
            end do
        end do

    end subroutine

    subroutine test_calc_perp_coll_freq

        use math_constants, only: m_e, m_D
        use libneo_species, only: species_t

        implicit none

        type(species_t) :: species_a, species_b
        character(len=2) :: interaction_type
        real(kind=real_kind) :: coulomb_log, coll_freq, vel
        
        call print_test("test_calc_perp_coll_freq")

        ! electron-electron interaction
        species_a%temp = 1.0d3
        species_a%dens = 1.0d13
        species_a%mass = m_e
        species_a%charge_num = -1
        species_b%temp = 1.0d3
        species_b%dens = species_a%dens
        species_b%mass = m_e
        species_b%charge_num = -1

        vel = 1.0d9

        ! electron electron interaction
        interaction_type = "ee"
        call calc_coulomb_log(interaction_type, species_a, species_b, coulomb_log)
        call calc_perp_coll_freq(vel, species_a, species_b, coulomb_log, coll_freq)

        print *, "interaction_type ", interaction_type, ": ", "coll_freq = ", coll_freq
        if (abs(coll_freq - 97549.2133117624) > 1d0) then ! number is for the parameters above
            call print_fail
            stop "coll_freq != 97549.2133117624"
        else
            call print_ok
        end if

        ! electron ion interaction (case 2 of coulomb log)
        interaction_type = "ei"
        species_b%mass = m_D
        species_b%charge_num = 1

        call calc_coulomb_log(interaction_type, species_a, species_b, coulomb_log)
        call calc_perp_coll_freq(vel, species_a, species_b, coulomb_log, coll_freq)

        print *, "interaction_type ", interaction_type, ": ", "coll_freq = ", coll_freq
        if (abs(coll_freq - 256856.553227032) > 1d0) then ! number is for the parameters above
            call print_fail
            stop "coll_freq != 256856.553227032"
        else
            call print_ok
        end if

        ! ion electron interaction (case 2 of coulomb log, same value as for ei)
        interaction_type = "ie"
        species_a%mass = m_D
        species_b%mass = m_e
        species_a%charge_num = 1
        species_b%charge_num = -1

        call calc_perp_coll_freq(vel, species_a, species_b, coulomb_log, coll_freq)

        print *, "interaction_type ", interaction_type, ": ", "coll_freq = ", coll_freq
        if (abs(coll_freq - 0.00724064837092233) > 1d0) then ! number is for the parameters above
            call print_fail
            stop "coll_freq != 0.00724064837092233"
        else
            call print_ok
        end if

        ! ion-ion interaction
        interaction_type = "ii"
        species_a%mass = m_D
        species_b%mass = m_D
        species_a%charge_num = 1
        species_b%charge_num = 1

        call calc_coulomb_log(interaction_type, species_a, species_b, coulomb_log)
        call calc_perp_coll_freq(vel, species_a, species_b, coulomb_log, coll_freq)

        print *, "interaction_type ", interaction_type, ": ", "coll_freq = ", coll_freq
        if (abs(coll_freq - 0.0215856541453319) > 1d-6) then ! number is for the parameters above
            call print_fail
            stop "coll_freq != 0.0215856541453319"
        else
            call print_ok
        end if

    end subroutine

    subroutine test_calc_perp_coll_freq_fast_limit

        implicit none
        real(kind=real_kind) :: vel, dens_e, coulomb_log
        real(kind=real_kind) :: coll_freq

        vel = 1.0d8
        dens_e = 1.0d13
        coulomb_log = 15.94094809934479d0

        call print_test("test_calc_perp_coll_freq_fast_limit")

        call calc_perp_coll_freq_fast_limit_ee(vel, dens_e, coulomb_log, coll_freq)
        print *, "Fast limit ee: coll_freq = ", coll_freq
        if (abs(coll_freq - 256089317.82351d0) > 10.0d0) then
            call print_fail
            print *, "Difference = " , abs(coll_freq - 256089317.82351d0)
            stop "coll_freq != 256089317.82351"
        else
            call print_ok
        end if

    end subroutine

end program