module jorek_field_values
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use jorek_bezier, only: cubic_jorek_basis
    use jorek_restart, only: jorek_restart_t

    implicit none
    private

    public :: evaluate_jorek_variable

contains

    pure subroutine evaluate_jorek_variable(data, element, variable, s, t, phi, &
            value, derivative, ierr)
        type(jorek_restart_t), intent(in) :: data
        integer, intent(in) :: element, variable
        real(dp), intent(in) :: s, t, phi
        real(dp), intent(out) :: value, derivative(3)
        integer, intent(out) :: ierr

        real(dp) :: basis(4, 4), basis_s(4, 4), basis_t(4, 4)
        real(dp), allocatable :: weights(:), weights_phi(:)
        real(dp) :: coefficient, coefficient_phi, scale
        integer :: degree, node, vertex

        value = 0.0_dp
        derivative = 0.0_dp
        ierr = 0
        if (element < 1 .or. element > data%n_elements) then
            ierr = 1
            return
        end if
        call cubic_jorek_basis(s, t, basis, basis_s, basis_t, ierr)
        if (ierr /= 0) return
        if (.not. supported_value_layout(data)) then
            ierr = 3
            return
        end if
        if (variable < 1 .or. variable > data%n_var) then
            ierr = 5
            return
        end if
        if (data%n_tor < 1 .or. mod(data%n_tor, 2) /= 1 &
            .or. data%n_period < 1) then
            ierr = 6
            return
        end if
        if (size(data%values, 2) < data%n_tor) then
            ierr = 6
            return
        end if

        allocate (weights(data%n_tor), weights_phi(data%n_tor))
        call toroidal_weights(data%n_tor, data%n_period, phi, weights, &
            weights_phi)
        do vertex = 1, 4
            node = data%vertex(element, vertex)
            if (node < 1 .or. node > data%n_nodes) then
                ierr = 4
                value = 0.0_dp
                derivative = 0.0_dp
                return
            end if
            do degree = 1, 4
                scale = data%size(element, vertex, degree)
                coefficient = dot_product(data%values(node, :, degree, variable), &
                    weights)*scale
                coefficient_phi = &
                    dot_product(data%values(node, :, degree, variable), &
                    weights_phi)*scale
                value = value + coefficient*basis(vertex, degree)
                derivative(1) = derivative(1) &
                    + coefficient*basis_s(vertex, degree)
                derivative(2) = derivative(2) &
                    + coefficient*basis_t(vertex, degree)
                derivative(3) = derivative(3) &
                    + coefficient_phi*basis(vertex, degree)
            end do
        end do
    end subroutine evaluate_jorek_variable

    pure subroutine toroidal_weights(n_tor, n_period, phi, weights, weights_phi)
        integer, intent(in) :: n_tor, n_period
        real(dp), intent(in) :: phi
        real(dp), intent(out) :: weights(n_tor), weights_phi(n_tor)

        real(dp) :: angle, mode
        integer :: harmonic

        weights = 0.0_dp
        weights_phi = 0.0_dp
        weights(1) = 1.0_dp
        do harmonic = 1, (n_tor - 1)/2
            mode = real(n_period*harmonic, dp)
            angle = mode*phi
            weights(2*harmonic) = cos(angle)
            weights(2*harmonic + 1) = sin(angle)
            weights_phi(2*harmonic) = -mode*sin(angle)
            weights_phi(2*harmonic + 1) = mode*cos(angle)
        end do
    end subroutine toroidal_weights

    pure logical function supported_value_layout(data)
        type(jorek_restart_t), intent(in) :: data

        supported_value_layout = .false.
        if (data%n_order /= 3) return
        if (data%n_degrees /= 4) return
        if (data%n_vertex_max /= 4) return
        if (data%n_var < 1) return
        if (.not. allocated(data%values)) return
        if (.not. allocated(data%vertex)) return
        if (.not. allocated(data%size)) return
        if (size(data%values, 1) < data%n_nodes) return
        if (size(data%values, 3) < 4) return
        if (size(data%values, 4) < data%n_var) return
        if (any(shape(data%vertex) < [data%n_elements, 4])) return
        if (any(shape(data%size) < [data%n_elements, 4, 4])) return
        supported_value_layout = .true.
    end function supported_value_layout

end module jorek_field_values
