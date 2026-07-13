module jorek_bezier
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use jorek_restart, only: jorek_restart_t

    implicit none
    private

    public :: cubic_jorek_basis, evaluate_jorek_geometry, locate_jorek_element

contains

    pure subroutine cubic_jorek_basis(s, t, basis, basis_s, basis_t, ierr)
        real(dp), intent(in) :: s, t
        real(dp), intent(out) :: basis(4, 4), basis_s(4, 4), basis_t(4, 4)
        integer, intent(out) :: ierr

        integer, parameter :: s_end(4) = [1, 2, 2, 1]
        integer, parameter :: t_end(4) = [1, 1, 2, 2]
        real(dp) :: bs(2, 2), bt(2, 2), dbs(2, 2), dbt(2, 2)
        integer :: vertex

        basis = 0.0_dp
        basis_s = 0.0_dp
        basis_t = 0.0_dp
        ierr = 0
        if (s < 0.0_dp .or. s > 1.0_dp .or. &
            t < 0.0_dp .or. t > 1.0_dp) then
            ierr = 2
            return
        end if

        call cubic_endpoint_basis(s, bs, dbs)
        call cubic_endpoint_basis(t, bt, dbt)
        do vertex = 1, 4
            basis(vertex, 1) = bs(s_end(vertex), 1)*bt(t_end(vertex), 1)
            basis(vertex, 2) = bs(s_end(vertex), 2)*bt(t_end(vertex), 1)
            basis(vertex, 3) = bs(s_end(vertex), 1)*bt(t_end(vertex), 2)
            basis(vertex, 4) = bs(s_end(vertex), 2)*bt(t_end(vertex), 2)

            basis_s(vertex, 1) = dbs(s_end(vertex), 1)*bt(t_end(vertex), 1)
            basis_s(vertex, 2) = dbs(s_end(vertex), 2)*bt(t_end(vertex), 1)
            basis_s(vertex, 3) = dbs(s_end(vertex), 1)*bt(t_end(vertex), 2)
            basis_s(vertex, 4) = dbs(s_end(vertex), 2)*bt(t_end(vertex), 2)

            basis_t(vertex, 1) = bs(s_end(vertex), 1)*dbt(t_end(vertex), 1)
            basis_t(vertex, 2) = bs(s_end(vertex), 2)*dbt(t_end(vertex), 1)
            basis_t(vertex, 3) = bs(s_end(vertex), 1)*dbt(t_end(vertex), 2)
            basis_t(vertex, 4) = bs(s_end(vertex), 2)*dbt(t_end(vertex), 2)
        end do
    end subroutine cubic_jorek_basis

    pure subroutine evaluate_jorek_geometry(data, element, s, t, rz, rz_st, ierr)
        type(jorek_restart_t), intent(in) :: data
        integer, intent(in) :: element
        real(dp), intent(in) :: s, t
        real(dp), intent(out) :: rz(2), rz_st(2, 2)
        integer, intent(out) :: ierr

        real(dp) :: basis(4, 4), basis_s(4, 4), basis_t(4, 4), coefficient
        integer :: coordinate, degree, node, vertex

        rz = 0.0_dp
        rz_st = 0.0_dp
        ierr = 0
        if (element < 1 .or. element > data%n_elements) then
            ierr = 1
            return
        end if
        call cubic_jorek_basis(s, t, basis, basis_s, basis_t, ierr)
        if (ierr /= 0) return
        if (.not. supported_geometry_layout(data)) then
            ierr = 3
            return
        end if

        do vertex = 1, 4
            node = data%vertex(element, vertex)
            if (node < 1 .or. node > data%n_nodes) then
                ierr = 4
                rz = 0.0_dp
                rz_st = 0.0_dp
                return
            end if
            do coordinate = 1, 2
                do degree = 1, 4
                    coefficient = data%x(node, 1, degree, coordinate) &
                        *data%size(element, vertex, degree)
                    rz(coordinate) = rz(coordinate) &
                        + coefficient*basis(vertex, degree)
                    rz_st(coordinate, 1) = rz_st(coordinate, 1) &
                        + coefficient*basis_s(vertex, degree)
                    rz_st(coordinate, 2) = rz_st(coordinate, 2) &
                        + coefficient*basis_t(vertex, degree)
                end do
            end do
        end do
    end subroutine evaluate_jorek_geometry

    pure subroutine locate_jorek_element(data, target_rz, element, s, t, ierr)
        type(jorek_restart_t), intent(in) :: data
        real(dp), intent(in) :: target_rz(2)
        integer, intent(out) :: element
        real(dp), intent(out) :: s, t
        integer, intent(out) :: ierr

        real(dp), parameter :: seeds(2, 9) = reshape([ &
            0.5_dp, 0.5_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
            1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.5_dp, 0.0_dp, &
            1.0_dp, 0.5_dp, 0.5_dp, 1.0_dp, 0.0_dp, 0.5_dp], [2, 9])
        integer, parameter :: max_iterations = 30

        real(dp) :: delta(2), determinant, determinant_scale, residual(2)
        real(dp) :: rz(2), rz_st(2, 2), scale, tolerance, trial_s, trial_t
        integer :: candidate, iteration, seed, geometry_ierr

        element = 0
        s = 0.0_dp
        t = 0.0_dp
        ierr = 1
        if (data%n_elements < 1) then
            ierr = 2
            return
        end if
        call evaluate_jorek_geometry(data, 1, 0.5_dp, 0.5_dp, rz, rz_st, &
            geometry_ierr)
        if (geometry_ierr /= 0) then
            ierr = 2
            return
        end if

        do candidate = 1, data%n_elements
            do seed = 1, size(seeds, 2)
                trial_s = seeds(1, seed)
                trial_t = seeds(2, seed)
                do iteration = 1, max_iterations
                    call evaluate_jorek_geometry(data, candidate, trial_s, &
                        trial_t, rz, rz_st, geometry_ierr)
                    if (geometry_ierr /= 0) exit
                    residual = target_rz - rz
                    scale = max(1.0_dp, maxval(abs(target_rz)), maxval(abs(rz)))
                    tolerance = 1024.0_dp*epsilon(1.0_dp)**0.75_dp*scale
                    if (maxval(abs(residual)) <= tolerance) then
                        element = candidate
                        s = trial_s
                        t = trial_t
                        ierr = 0
                        return
                    end if
                    determinant = rz_st(1, 1)*rz_st(2, 2) &
                        - rz_st(1, 2)*rz_st(2, 1)
                    determinant_scale = maxval(abs(rz_st(:, 1))) &
                        *maxval(abs(rz_st(:, 2)))
                    if (abs(determinant) <= 64.0_dp*epsilon(1.0_dp) &
                        *determinant_scale) exit
                    delta(1) = (residual(1)*rz_st(2, 2) &
                        - residual(2)*rz_st(1, 2))/determinant
                    delta(2) = (-residual(1)*rz_st(2, 1) &
                        + residual(2)*rz_st(1, 1))/determinant
                    trial_s = min(1.0_dp, max(0.0_dp, trial_s + delta(1)))
                    trial_t = min(1.0_dp, max(0.0_dp, trial_t + delta(2)))
                end do
            end do
        end do
    end subroutine locate_jorek_element

    pure subroutine cubic_endpoint_basis(u, basis, derivative)
        real(dp), intent(in) :: u
        real(dp), intent(out) :: basis(2, 2), derivative(2, 2)

        basis(1, 1) = 2.0_dp*u**3 - 3.0_dp*u**2 + 1.0_dp
        basis(1, 2) = 3.0_dp*(u**3 - 2.0_dp*u**2 + u)
        basis(2, 1) = -2.0_dp*u**3 + 3.0_dp*u**2
        basis(2, 2) = 3.0_dp*(u**3 - u**2)

        derivative(1, 1) = 6.0_dp*u**2 - 6.0_dp*u
        derivative(1, 2) = 9.0_dp*u**2 - 12.0_dp*u + 3.0_dp
        derivative(2, 1) = -6.0_dp*u**2 + 6.0_dp*u
        derivative(2, 2) = 9.0_dp*u**2 - 6.0_dp*u
    end subroutine cubic_endpoint_basis

    pure logical function supported_geometry_layout(data)
        type(jorek_restart_t), intent(in) :: data

        supported_geometry_layout = .false.
        if (data%n_order /= 3) return
        if (data%n_degrees /= 4) return
        if (data%n_vertex_max /= 4) return
        if (data%n_coord_tor < 1) return
        if (data%n_dim < 2) return
        if (.not. allocated(data%x)) return
        if (.not. allocated(data%vertex)) return
        if (.not. allocated(data%size)) return
        if (any(shape(data%x) < [data%n_nodes, 1, 4, 2])) return
        if (any(shape(data%vertex) < [data%n_elements, 4])) return
        if (any(shape(data%size) < [data%n_elements, 4, 4])) return
        supported_geometry_layout = .true.
    end function supported_geometry_layout

end module jorek_bezier
