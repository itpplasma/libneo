program test_jorek_model303_field
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use jorek_model303_field, only: evaluate_jorek_model303_a, &
        evaluate_jorek_model303_b, evaluate_jorek_model303_at
    use jorek_restart, only: jorek_restart_t
    use util_for_test, only: print_test, print_ok, print_fail

    implicit none

    real(dp), parameter :: tol = 3.0e-14_dp
    integer :: nfail

    nfail = 0
    call test_reduced_mhd_field
    call test_reduced_mhd_potential
    call test_global_field_query
    call test_rejected_model
    call test_rejected_singular_geometry
    if (nfail > 0) error stop

contains

    subroutine test_reduced_mhd_field
        type(jorek_restart_t) :: data
        real(dp) :: b_r_z_phi(3), r, z
        integer :: ierr, nfail_before

        call print_test('model 303 field follows the reduced-MHD flux convention')
        nfail_before = nfail
        call make_model303_element(data)
        call evaluate_jorek_model303_b(data, 1, 0.25_dp, 0.6_dp, 0.4_dp, &
            b_r_z_phi, ierr)
        r = 10.25_dp
        z = 0.6_dp
        call check_int('ierr', ierr, 0)
        call check_real('B_R', b_r_z_phi(1), 1.0_dp)
        call check_real('B_Z', b_r_z_phi(2), -z/r)
        call check_real('B_phi', b_r_z_phi(3), data%F0/r)
        call report(nfail_before)
    end subroutine test_reduced_mhd_field

    subroutine test_reduced_mhd_potential
        type(jorek_restart_t) :: data
        real(dp) :: a_r_phi_z(3), r, z
        integer :: ierr, nfail_before

        call print_test('model 303 potential has the reduced-MHD curl')
        nfail_before = nfail
        call make_model303_element(data)
        call evaluate_jorek_model303_a(data, 1, 0.25_dp, 0.6_dp, 0.4_dp, &
            a_r_phi_z, ierr)
        r = 10.25_dp
        z = 0.6_dp
        call check_int('ierr', ierr, 0)
        call check_real('A_R', a_r_phi_z(1), 0.0_dp)
        call check_real('A_phi', a_r_phi_z(2), -r*z)
        call check_real('A_Z', a_r_phi_z(3), -data%F0*log(r))
        call report(nfail_before)
    end subroutine test_reduced_mhd_potential

    subroutine test_global_field_query
        type(jorek_restart_t) :: data
        real(dp) :: a_r_phi_z(3), b_r_z_phi(3), st(2), r, z
        integer :: element, ierr, nfail_before

        call print_test('global model 303 query returns compatible A and B')
        nfail_before = nfail
        call make_model303_element(data)
        r = 10.25_dp
        z = 0.6_dp
        call evaluate_jorek_model303_at(data, [r, z], 0.4_dp, a_r_phi_z, &
            b_r_z_phi, element, st, ierr)
        call check_int('ierr', ierr, 0)
        call check_int('element', element, 1)
        call check_real('s', st(1), 0.25_dp)
        call check_real('t', st(2), 0.6_dp)
        call check_real('A_phi', a_r_phi_z(2), -r*z)
        call check_real('B_R', b_r_z_phi(1), 1.0_dp)
        call check_real('B_Z', b_r_z_phi(2), -z/r)
        call check_real('B_phi', b_r_z_phi(3), data%F0/r)
        call report(nfail_before)
    end subroutine test_global_field_query

    subroutine test_rejected_model
        type(jorek_restart_t) :: data
        real(dp) :: b_r_z_phi(3)
        integer :: ierr, nfail_before

        call print_test('model 303 field rejects another JOREK model')
        nfail_before = nfail
        call make_model303_element(data)
        data%jorek_model = 199
        call evaluate_jorek_model303_b(data, 1, 0.5_dp, 0.5_dp, 0.0_dp, &
            b_r_z_phi, ierr)
        call check_int('ierr', ierr, 7)
        call report(nfail_before)
    end subroutine test_rejected_model

    subroutine test_rejected_singular_geometry
        type(jorek_restart_t) :: data
        real(dp) :: b_r_z_phi(3)
        integer :: ierr, nfail_before

        call print_test('model 303 field rejects a singular element map')
        nfail_before = nfail
        call make_model303_element(data)
        data%x = 0.0_dp
        call evaluate_jorek_model303_b(data, 1, 0.5_dp, 0.5_dp, 0.0_dp, &
            b_r_z_phi, ierr)
        call check_int('ierr', ierr, 8)
        call report(nfail_before)
    end subroutine test_rejected_singular_geometry

    subroutine make_model303_element(data)
        type(jorek_restart_t), intent(out) :: data

        real(dp), parameter :: r_node(4) = [10.0_dp, 11.0_dp, 11.0_dp, 10.0_dp]
        real(dp), parameter :: z_node(4) = [0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp]
        integer :: node

        data%jorek_model = 303
        data%n_order = 3
        data%n_degrees = 4
        data%n_coord_tor = 1
        data%n_dim = 2
        data%n_tor = 1
        data%n_period = 1
        data%n_var = 7
        data%n_nodes = 4
        data%n_elements = 1
        data%n_vertex_max = 4
        data%F0 = 20.0_dp
        allocate (data%x(4, 1, 4, 2), data%values(4, 1, 4, 7))
        allocate (data%vertex(1, 4), data%size(1, 4, 4))
        data%x = 0.0_dp
        data%values = 0.0_dp
        data%vertex(1, :) = [1, 2, 3, 4]
        data%size = 1.0_dp
        do node = 1, 4
            data%x(node, 1, 1, 1) = r_node(node)
            data%x(node, 1, 1, 2) = z_node(node)
            data%x(node, 1, 2, 1) = 1.0_dp/3.0_dp
            data%x(node, 1, 3, 2) = 1.0_dp/3.0_dp
            data%values(node, 1, 1, 1) = r_node(node)*z_node(node)
            data%values(node, 1, 2, 1) = z_node(node)/3.0_dp
            data%values(node, 1, 3, 1) = r_node(node)/3.0_dp
            data%values(node, 1, 4, 1) = 1.0_dp/9.0_dp
        end do
    end subroutine make_model303_element

    subroutine check_real(name, actual, expected)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: actual, expected

        if (abs(actual - expected) > tol*max(1.0_dp, abs(expected))) then
            call fail(name//' mismatch')
        end if
    end subroutine check_real

    subroutine check_int(name, actual, expected)
        character(len=*), intent(in) :: name
        integer, intent(in) :: actual, expected

        if (actual /= expected) call fail(name//' mismatch')
    end subroutine check_int

    subroutine fail(message)
        character(len=*), intent(in) :: message

        print *, '    ', message
        nfail = nfail + 1
    end subroutine fail

    subroutine report(nfail_before)
        integer, intent(in) :: nfail_before

        if (nfail > nfail_before) then
            call print_fail
        else
            call print_ok
        end if
    end subroutine report

end program test_jorek_model303_field
