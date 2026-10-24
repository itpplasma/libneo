module libneo_coordinate_jets
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private
    public :: rotate_z_vector, rotate_z_vector_der, rotate_z_vector_der2
    public :: cylindrical_field_direction, cylindrical_field_direction_der
    public :: pullback_covector, pullback_covector_der, scale_basis_der
contains
    include '../generated/coordinate_jets_leaves.inc'

    pure subroutine cylindrical_field_direction(position, field, direction)
        real(dp), intent(in) :: position(3), field(3)
        real(dp), intent(out) :: direction(3)
        if (sum(position(1:2)**2) <= 0.0_dp) error stop 'cylindrical direction: R=0'
        if (sum(field**2) <= 0.0_dp) error stop 'cylindrical direction: B=0'
        call cylindrical_direction_value(position(1), position(2), &
                                          field(1), field(2), field(3), &
                                          direction(1), direction(2), direction(3))
    end subroutine

    pure subroutine cylindrical_field_direction_der(position, dposition, field, &
                                                      dfield, direction, ddirection)
        real(dp), intent(in) :: position(3), dposition(3, 3), field(3), dfield(3, 3)
        real(dp), intent(out) :: direction(3), ddirection(3, 3)
        integer :: j
        if (sum(position(1:2)**2) <= 0.0_dp) error stop 'cylindrical direction: R=0'
        if (sum(field**2) <= 0.0_dp) error stop 'cylindrical direction: B=0'
        do j = 1, 3
            call cylindrical_direction_der(position(1), position(2), &
                       field(1), field(2), field(3), dposition(1, j), dposition(2, j), &
                       dfield(1, j), dfield(2, j), dfield(3, j), &
                       direction(1), direction(2), direction(3), &
                       ddirection(1, j), ddirection(2, j), ddirection(3, j))
        end do
    end subroutine

    pure subroutine rotate_z_vector(angle, vector, rotated)
        real(dp), intent(in) :: angle, vector(3)
        real(dp), intent(out) :: rotated(3)
        call rotate_pair_value(angle, vector(1), vector(2), rotated(1), rotated(2))
        rotated(3) = vector(3)
    end subroutine

    pure subroutine rotate_z_vector_der(angle, vector, derivative, rotated, drotated)
        real(dp), intent(in) :: angle, vector(3), derivative(3, 3)
        real(dp), intent(out) :: rotated(3), drotated(3, 3)
        real(dp) :: angle_derivative(3)
        integer :: j
        angle_derivative = [0.0_dp, 0.0_dp, 1.0_dp]
        do j = 1, 3
            call rotate_pair_der(angle, vector(1), vector(2), angle_derivative(j), &
                           derivative(1, j), derivative(2, j), rotated(1), rotated(2), &
                                 drotated(1, j), drotated(2, j))
        end do
        rotated(3) = vector(3)
        drotated(3, :) = derivative(3, :)
    end subroutine

    pure subroutine rotate_z_vector_der2(angle, vector, derivative, second, &
                                         rotated, drotated, d2rotated)
        real(dp), intent(in) :: angle, vector(3), derivative(3, 3), second(3, 6)
        real(dp), intent(out) :: rotated(3), drotated(3, 3), d2rotated(3, 3, 3)
       integer, parameter :: packed(3, 3) = reshape([1, 2, 3, 2, 4, 5, 3, 5, 6], [3, 3])
        real(dp) :: angle_derivative(3), first_i(2), first_j(2)
        integer :: i, j, k
        angle_derivative = [0.0_dp, 0.0_dp, 1.0_dp]
        do j = 1, 3
            do i = 1, j
                k = packed(i, j)
                call rotate_pair_der2(angle, vector(1), vector(2), &
                              angle_derivative(i), derivative(1, i), derivative(2, i), &
                              angle_derivative(j), derivative(1, j), derivative(2, j), &
                           0.0_dp, second(1, k), second(2, k), rotated(1), rotated(2), &
                                      first_i(1), first_i(2), first_j(1), first_j(2), &
                                      d2rotated(1, i, j), d2rotated(2, i, j))
                drotated(1:2, i) = first_i
                drotated(1:2, j) = first_j
                d2rotated(3, i, j) = second(3, k)
                d2rotated(:, j, i) = d2rotated(:, i, j)
            end do
        end do
        rotated(3) = vector(3)
        drotated(3, :) = derivative(3, :)
    end subroutine

    pure subroutine pullback_covector(vector, basis, covector)
        real(dp), intent(in) :: vector(3), basis(3, 3)
        real(dp), intent(out) :: covector(3)
        integer :: i
        do i = 1, 3
            call pullback_value(basis(1, i), basis(2, i), basis(3, i), &
                                vector(1), vector(2), vector(3), covector(i))
        end do
    end subroutine

    pure subroutine pullback_covector_der(vector, derivative, basis, dbasis, &
                                          covector, dcovector)
        real(dp), intent(in) :: vector(3), derivative(3, 3)
        real(dp), intent(in) :: basis(3, 3), dbasis(3, 3, 3)
        real(dp), intent(out) :: covector(3), dcovector(3, 3)
        integer :: i, j
        do j = 1, 3
            do i = 1, 3
                ! Output convention is dcovector(derivative index,component).
                call pullback_der(basis(1, i), basis(2, i), basis(3, i), &
                                  vector(1), vector(2), vector(3), &
                                  dbasis(1, i, j), dbasis(2, i, j), dbasis(3, i, j), &
                                 derivative(1, j), derivative(2, j), derivative(3, j), &
                                  covector(i), dcovector(j, i))
            end do
        end do
    end subroutine

    pure subroutine scale_basis_der(jacobian, second, basis, dbasis, scaled, dscaled)
        real(dp), intent(in) :: jacobian(3), second(3), basis(3, 3), dbasis(3, 3, 3)
        real(dp), intent(out) :: scaled(3, 3), dscaled(3, 3, 3)
        real(dp) :: diagonal_second
        integer :: component, i, j
        do j = 1, 3
            do i = 1, 3
                diagonal_second = 0.0_dp
                if (i == j) diagonal_second = second(i)
                do component = 1, 3
                 call diagonal_basis_der(basis(component, i), dbasis(component, i, j), &
                                            jacobian(i), jacobian(j), diagonal_second, &
                                         scaled(component, i), dscaled(component, i, j))
                end do
            end do
        end do
    end subroutine
end module libneo_coordinate_jets
