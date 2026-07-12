program test_jorek_bezier_geometry
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use jorek_bezier, only: cubic_jorek_basis, evaluate_jorek_geometry
    use jorek_restart, only: jorek_restart_t
    use util_for_test, only: print_test, print_ok, print_fail

    implicit none

    real(dp), parameter :: tol = 2.0e-14_dp
    integer :: nfail

    nfail = 0
    call test_jorek_basis_oracle
    call test_affine_geometry
    call test_shared_edge
    call test_rejected_inputs
    if (nfail > 0) error stop

contains

    subroutine test_jorek_basis_oracle
        real(dp), parameter :: basis_ref(4, 4) = reshape([ &
            0.193536_dp, 0.022464_dp, 0.081536_dp, 0.702464_dp, &
            0.082944_dp, -0.020736_dp, -0.075264_dp, 0.301056_dp, &
            0.169344_dp, 0.019656_dp, -0.045864_dp, -0.395136_dp, &
            0.072576_dp, -0.018144_dp, 0.042336_dp, -0.169344_dp], [4, 4])
        real(dp), parameter :: basis_s_ref(4, 4) = reshape([ &
            -0.20736_dp, 0.20736_dp, 0.75264_dp, -0.75264_dp, &
            0.20736_dp, -0.18144_dp, -0.65856_dp, 0.75264_dp, &
            -0.18144_dp, 0.18144_dp, -0.42336_dp, 0.42336_dp, &
            0.18144_dp, -0.15876_dp, 0.37044_dp, -0.42336_dp], [4, 4])
        real(dp), parameter :: basis_t_ref(4, 4) = reshape([ &
            -1.12896_dp, -0.13104_dp, 0.13104_dp, 1.12896_dp, &
            -0.48384_dp, 0.12096_dp, -0.12096_dp, 0.48384_dp, &
            -0.88704_dp, -0.10296_dp, 0.02184_dp, 0.18816_dp, &
            -0.38016_dp, 0.09504_dp, -0.02016_dp, 0.08064_dp], [4, 4])

        real(dp) :: basis(4, 4), basis_s(4, 4), basis_t(4, 4)
        integer :: ierr, nfail_before

        call print_test('cubic basis matches JOREK oracle at s=0.2, t=0.7')
        nfail_before = nfail
        call cubic_jorek_basis(0.2_dp, 0.7_dp, basis, basis_s, basis_t, ierr)
        call check_int('ierr', ierr, 0)
        call check_array('basis', basis, basis_ref)
        call check_array('basis_s', basis_s, basis_s_ref)
        call check_array('basis_t', basis_t, basis_t_ref)
        call report(nfail_before)
    end subroutine test_jorek_basis_oracle

    subroutine test_affine_geometry
        type(jorek_restart_t) :: data
        real(dp) :: rz(2), rz_st(2, 2)
        integer :: ierr, nfail_before

        call print_test('geometry reproduces an affine element and Jacobian')
        nfail_before = nfail
        call make_two_element_strip(data)
        call evaluate_jorek_geometry(data, 1, 0.25_dp, 0.6_dp, rz, rz_st, ierr)
        call check_int('ierr', ierr, 0)
        call check_real('R', rz(1), 10.25_dp)
        call check_real('Z', rz(2), 0.6_dp)
        call check_real('R_s', rz_st(1, 1), 1.0_dp)
        call check_real('Z_s', rz_st(2, 1), 0.0_dp)
        call check_real('R_t', rz_st(1, 2), 0.0_dp)
        call check_real('Z_t', rz_st(2, 2), 1.0_dp)
        call report(nfail_before)
    end subroutine test_affine_geometry

    subroutine test_shared_edge
        type(jorek_restart_t) :: data
        real(dp) :: left(2), right(2), left_st(2, 2), right_st(2, 2)
        integer :: ierr, nfail_before

        call print_test('adjacent elements agree on their shared edge')
        nfail_before = nfail
        call make_two_element_strip(data)
        call evaluate_jorek_geometry(data, 1, 1.0_dp, 0.37_dp, &
            left, left_st, ierr)
        call check_int('left ierr', ierr, 0)
        call evaluate_jorek_geometry(data, 2, 0.0_dp, 0.37_dp, &
            right, right_st, ierr)
        call check_int('right ierr', ierr, 0)
        call check_vector('position', left, right)
        call check_vector('s tangent', left_st(:, 1), right_st(:, 1))
        call check_vector('t tangent', left_st(:, 2), right_st(:, 2))
        call report(nfail_before)
    end subroutine test_shared_edge

    subroutine test_rejected_inputs
        type(jorek_restart_t) :: data
        real(dp) :: rz(2), rz_st(2, 2)
        integer :: ierr, nfail_before

        call print_test('geometry rejects invalid element, coordinates, and order')
        nfail_before = nfail
        call make_two_element_strip(data)
        call evaluate_jorek_geometry(data, 0, 0.5_dp, 0.5_dp, rz, rz_st, ierr)
        call check_int('element ierr', ierr, 1)
        call evaluate_jorek_geometry(data, 1, -0.1_dp, 0.5_dp, rz, rz_st, ierr)
        call check_int('coordinate ierr', ierr, 2)
        data%n_order = 5
        call evaluate_jorek_geometry(data, 1, 0.5_dp, 0.5_dp, rz, rz_st, ierr)
        call check_int('order ierr', ierr, 3)
        call report(nfail_before)
    end subroutine test_rejected_inputs

    subroutine make_two_element_strip(data)
        type(jorek_restart_t), intent(out) :: data

        real(dp), parameter :: r_node(6) = [10.0_dp, 11.0_dp, 11.0_dp, &
            10.0_dp, 12.0_dp, 12.0_dp]
        real(dp), parameter :: z_node(6) = [0.0_dp, 0.0_dp, 1.0_dp, &
            1.0_dp, 0.0_dp, 1.0_dp]
        integer :: i

        data%n_order = 3
        data%n_degrees = 4
        data%n_coord_tor = 1
        data%n_dim = 2
        data%n_nodes = 6
        data%n_elements = 2
        data%n_vertex_max = 4
        allocate (data%x(6, 1, 4, 2), data%vertex(2, 4), data%size(2, 4, 4))
        data%x = 0.0_dp
        data%size = 1.0_dp
        data%vertex(1, :) = [1, 2, 3, 4]
        data%vertex(2, :) = [2, 5, 6, 3]
        do i = 1, 6
            data%x(i, 1, 1, 1) = r_node(i)
            data%x(i, 1, 1, 2) = z_node(i)
            data%x(i, 1, 2, 1) = 1.0_dp/3.0_dp
            data%x(i, 1, 3, 2) = 1.0_dp/3.0_dp
        end do
    end subroutine make_two_element_strip

    subroutine check_array(name, actual, expected)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: actual(:, :), expected(:, :)

        if (maxval(abs(actual - expected)) > tol) call fail(name//' mismatch')
    end subroutine check_array

    subroutine check_vector(name, actual, expected)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: actual(:), expected(:)

        if (maxval(abs(actual - expected)) > tol) call fail(name//' mismatch')
    end subroutine check_vector

    subroutine check_real(name, actual, expected)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: actual, expected

        if (abs(actual - expected) > tol) call fail(name//' mismatch')
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

end program test_jorek_bezier_geometry
