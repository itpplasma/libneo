program test_species
    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    logical :: all_passed

    all_passed = .true.

    call test_init_electron_species(all_passed)
    call test_init_deuterium_species(all_passed)
    call test_init_deuterium_plasma(all_passed)
    call test_electron_physical_constants(all_passed)
    call test_deuterium_physical_constants(all_passed)
    call test_species_type_assignments(all_passed)

    if (.not. all_passed) then
        error stop "One or more species tests failed"
    end if

contains

    subroutine test_init_electron_species(passed)
        use libneo_species, only: species_t, init_electron_species
        use math_constants, only: m_e
        logical, intent(inout) :: passed

        type(species_t) :: electron
        real(dp), parameter :: temp = 1000.0d0, dens = 1.0d14

        call init_electron_species(temp, dens, electron)

        if (trim(electron%name) /= "electron") then
            write(*,*) "FAIL: electron name. Got:", trim(electron%name)
            passed = .false.
            return
        end if

        if (electron%typ /= "e") then
            write(*,*) "FAIL: electron type. Got:", electron%typ
            passed = .false.
            return
        end if

        if (electron%typ_int /= 1) then
            write(*,*) "FAIL: electron typ_int. Got:", electron%typ_int, " expected 1"
            passed = .false.
            return
        end if

        if (abs(electron%temp - temp) > 1.0d-12) then
            write(*,*) "FAIL: electron temp. Got:", electron%temp
            passed = .false.
            return
        end if

        if (abs(electron%dens - dens) > 1.0d-12) then
            write(*,*) "FAIL: electron dens. Got:", electron%dens
            passed = .false.
            return
        end if

        if (abs(electron%mass - m_e) > 1.0d-40) then
            write(*,*) "FAIL: electron mass. Got:", electron%mass, " expected:", m_e
            passed = .false.
            return
        end if

        if (electron%charge_num /= -1) then
            write(*,*) "FAIL: electron charge_num. Got:", electron%charge_num
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_init_electron_species"
    end subroutine test_init_electron_species

    subroutine test_init_deuterium_species(passed)
        use libneo_species, only: species_t, init_deuterium_species
        use math_constants, only: m_D
        logical, intent(inout) :: passed

        type(species_t) :: deuterium
        real(dp), parameter :: temp = 1500.0d0, dens = 5.0d13

        call init_deuterium_species(temp, dens, deuterium)

        if (trim(deuterium%name) /= "deuterium") then
            write(*,*) "FAIL: deuterium name. Got:", trim(deuterium%name)
            passed = .false.
            return
        end if

        if (deuterium%typ /= "i") then
            write(*,*) "FAIL: deuterium type. Got:", deuterium%typ
            passed = .false.
            return
        end if

        if (deuterium%typ_int /= 2) then
            write(*,*) "FAIL: deuterium typ_int. Got:", deuterium%typ_int, " expected 2"
            passed = .false.
            return
        end if

        if (abs(deuterium%temp - temp) > 1.0d-12) then
            write(*,*) "FAIL: deuterium temp. Got:", deuterium%temp
            passed = .false.
            return
        end if

        if (abs(deuterium%dens - dens) > 1.0d-12) then
            write(*,*) "FAIL: deuterium dens. Got:", deuterium%dens
            passed = .false.
            return
        end if

        if (abs(deuterium%mass - m_D) > 1.0d-35) then
            write(*,*) "FAIL: deuterium mass. Got:", deuterium%mass, " expected:", m_D
            passed = .false.
            return
        end if

        if (deuterium%charge_num /= 1) then
            write(*,*) "FAIL: deuterium charge_num. Got:", deuterium%charge_num
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_init_deuterium_species"
    end subroutine test_init_deuterium_species

    subroutine test_init_deuterium_plasma(passed)
        use libneo_species, only: species_t, init_deuterium_plasma
        logical, intent(inout) :: passed

        type(species_t) :: plasma(2)
        real(dp), parameter :: temp_e = 3000.0d0, temp_D = 1500.0d0, dens = 1.0d14

        call init_deuterium_plasma(temp_e, temp_D, dens, plasma)

        if (plasma(1)%typ /= "e") then
            write(*,*) "FAIL: plasma(1) should be electron. Got:", plasma(1)%typ
            passed = .false.
            return
        end if

        if (plasma(2)%typ /= "i") then
            write(*,*) "FAIL: plasma(2) should be ion. Got:", plasma(2)%typ
            passed = .false.
            return
        end if

        if (abs(plasma(1)%temp - temp_e) > 1.0d-12) then
            write(*,*) "FAIL: electron temp in plasma. Got:", plasma(1)%temp
            passed = .false.
            return
        end if

        if (abs(plasma(2)%temp - temp_D) > 1.0d-12) then
            write(*,*) "FAIL: deuterium temp in plasma. Got:", plasma(2)%temp
            passed = .false.
            return
        end if

        if (abs(plasma(1)%dens - dens) > 1.0d-12) then
            write(*,*) "FAIL: electron dens in plasma. Got:", plasma(1)%dens
            passed = .false.
            return
        end if

        if (abs(plasma(2)%dens - dens) > 1.0d-12) then
            write(*,*) "FAIL: deuterium dens in plasma. Got:", plasma(2)%dens
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_init_deuterium_plasma"
    end subroutine test_init_deuterium_plasma

    subroutine test_electron_physical_constants(passed)
        use math_constants, only: m_e
        logical, intent(inout) :: passed

        real(dp), parameter :: m_e_nist = 9.1093837139d-28

        if (abs(m_e - m_e_nist) / m_e_nist > 1.0d-6) then
            write(*,*) "FAIL: electron mass differs from NIST value"
            write(*,*) "  Library:", m_e
            write(*,*) "  NIST:   ", m_e_nist
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_electron_physical_constants"
    end subroutine test_electron_physical_constants

    subroutine test_deuterium_physical_constants(passed)
        use math_constants, only: m_D, m_p
        logical, intent(inout) :: passed

        real(dp), parameter :: m_D_nist = 3.343583719d-24
        real(dp), parameter :: m_p_nist = 1.672621637d-24

        if (abs(m_D - m_D_nist) / m_D_nist > 1.0d-6) then
            write(*,*) "FAIL: deuterium mass differs from NIST value"
            write(*,*) "  Library:", m_D
            write(*,*) "  NIST:   ", m_D_nist
            passed = .false.
            return
        end if

        if (abs(m_p - m_p_nist) / m_p_nist > 1.0d-6) then
            write(*,*) "FAIL: proton mass differs from NIST value"
            write(*,*) "  Library:", m_p
            write(*,*) "  NIST:   ", m_p_nist
            passed = .false.
            return
        end if

        if (m_D <= m_p) then
            write(*,*) "FAIL: deuterium mass should be greater than proton mass"
            passed = .false.
            return
        end if

        if (m_D >= 2.0d0 * m_p) then
            write(*,*) "FAIL: deuterium mass should be less than 2x proton mass"
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_deuterium_physical_constants"
    end subroutine test_deuterium_physical_constants

    subroutine test_species_type_assignments(passed)
        use libneo_species, only: species_t, init_electron_species, init_deuterium_species
        logical, intent(inout) :: passed

        type(species_t) :: electron, deuterium
        integer :: ee_type, ei_type, ii_type

        call init_electron_species(1000.0d0, 1.0d13, electron)
        call init_deuterium_species(1000.0d0, 1.0d13, deuterium)

        ee_type = electron%typ_int + electron%typ_int
        ei_type = electron%typ_int + deuterium%typ_int
        ii_type = deuterium%typ_int + deuterium%typ_int

        if (ee_type /= 2) then
            write(*,*) "FAIL: e-e interaction type should be 2. Got:", ee_type
            passed = .false.
            return
        end if

        if (ei_type /= 3) then
            write(*,*) "FAIL: e-i interaction type should be 3. Got:", ei_type
            passed = .false.
            return
        end if

        if (ii_type /= 4) then
            write(*,*) "FAIL: i-i interaction type should be 4. Got:", ii_type
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_species_type_assignments"
    end subroutine test_species_type_assignments

end program test_species
