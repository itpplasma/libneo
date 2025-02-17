program test_transport

    !use libneo_kinds, only: real_kind
    use libneo_transport, only: D_one_over_nu_11, init_gauss_laguerre_integration

    implicit none

    call test_D_one_over_nu_11

    contains

    subroutine test_D_one_over_nu_11
        
        implicit none

        call init_gauss_laguerre_integration(5.0d0/2.0d0)
        call D_one_over_nu_11

    end subroutine

end program