module jorek_bezier
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use jorek_restart, only: jorek_restart_t

    implicit none
    private

    type :: jorek_locator_t
        integer :: n_elements = 0
        integer :: n_bins(2) = 0
        real(dp) :: grid_lower(2) = 0.0_dp
        real(dp) :: bin_width(2) = 0.0_dp
        logical :: has_axis = .false.
        real(dp) :: axis_rz(2) = 0.0_dp
        real(dp), allocatable :: bounds(:, :, :)
        integer, allocatable :: bin_offsets(:)
        integer, allocatable :: bin_elements(:)
        integer, allocatable :: axis_elements(:)
        integer, allocatable :: axis_sides(:)
    end type jorek_locator_t

    public :: jorek_locator_t, cubic_jorek_basis, evaluate_jorek_geometry
    public :: build_jorek_locator, locate_jorek_element
    public :: locate_jorek_element_indexed
    public :: jorek_locator_candidate_count
    public :: is_jorek_axis_target

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

        integer :: candidate, candidate_ierr

        element = 0
        s = 0.0_dp
        t = 0.0_dp
        ierr = 1
        if (data%n_elements < 1) then
            ierr = 2
            return
        end if
        do candidate = 1, data%n_elements
            call invert_jorek_element(data, candidate, target_rz, s, t, &
                candidate_ierr)
            if (candidate_ierr == 0) then
                element = candidate
                ierr = 0
                return
            end if
            if (candidate_ierr == 2) then
                ierr = 2
                return
            end if
        end do
    end subroutine locate_jorek_element

    subroutine build_jorek_locator(data, locator, ierr)
        type(jorek_restart_t), intent(in) :: data
        type(jorek_locator_t), intent(out) :: locator
        integer, intent(out) :: ierr

        integer :: element

        ierr = 0
        if (.not. supported_geometry_layout(data) .or. data%n_elements < 1) then
            ierr = 2
            return
        end if
        locator%n_elements = data%n_elements
        allocate (locator%bounds(2, 2, data%n_elements))
        locator%bounds(:, 1, :) = 0.0_dp
        locator%bounds(:, 2, :) = 0.0_dp
        do element = 1, data%n_elements
            call build_element_bounds(data, element, &
                locator%bounds(:, :, element), ierr)
            if (ierr /= 0) return
        end do
        call build_axis_metadata(data, locator, ierr)
        if (ierr /= 0) return
        call build_locator_bins(locator, ierr)
    end subroutine build_jorek_locator

    subroutine build_axis_metadata(data, locator, ierr)
        type(jorek_restart_t), intent(in) :: data
        type(jorek_locator_t), intent(inout) :: locator
        integer, intent(out) :: ierr

        integer, allocatable :: elements(:), sides(:)
        real(dp) :: edge_rz(2), reference_rz(2), scale, tolerance
        integer :: count, element, side, edge_ierr

        ierr = 0
        count = 0
        allocate (elements(data%n_elements), sides(data%n_elements))
        do element = 1, data%n_elements
            call collapsed_element_side(data, element, side, edge_rz, edge_ierr)
            if (edge_ierr /= 0) then
                ierr = 2
                return
            end if
            if (side == 0) cycle
            if (count == 0) then
                reference_rz = edge_rz
            else
                scale = max(1.0_dp, maxval(abs(reference_rz)), &
                    maxval(abs(edge_rz)))
                tolerance = 64.0_dp*epsilon(1.0_dp)*scale
                if (maxval(abs(edge_rz - reference_rz)) > tolerance) then
                    ierr = 2
                    return
                end if
            end if
            count = count + 1
            elements(count) = element
            sides(count) = side
        end do
        if (count == 0) return
        locator%has_axis = .true.
        locator%axis_rz = reference_rz
        allocate (locator%axis_elements(count), locator%axis_sides(count))
        locator%axis_elements = elements(:count)
        locator%axis_sides = sides(:count)
    end subroutine build_axis_metadata

    pure subroutine collapsed_element_side(data, element, side, edge_rz, ierr)
        type(jorek_restart_t), intent(in) :: data
        integer, intent(in) :: element
        integer, intent(out) :: side
        real(dp), intent(out) :: edge_rz(2)
        integer, intent(out) :: ierr

        real(dp), parameter :: corners(2, 4) = reshape([ &
            0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
            1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp], [2, 4])
        integer, parameter :: edge_vertices(2, 4) = reshape([ &
            1, 4, 2, 3, 1, 2, 4, 3], [2, 4])
        real(dp) :: corner_rz(2, 4), rz_st(2, 2), scale, tolerance
        integer :: corner, vertex_1, vertex_2

        side = 0
        edge_rz = 0.0_dp
        ierr = 0
        do corner = 1, 4
            call evaluate_jorek_geometry(data, element, corners(1, corner), &
                corners(2, corner), corner_rz(:, corner), rz_st, ierr)
            if (ierr /= 0) return
        end do
        scale = max(1.0_dp, maxval(abs(corner_rz)))
        tolerance = 64.0_dp*epsilon(1.0_dp)*scale
        do side = 1, 4
            vertex_1 = edge_vertices(1, side)
            vertex_2 = edge_vertices(2, side)
            if (maxval(abs(corner_rz(:, vertex_1) &
                    - corner_rz(:, vertex_2))) <= tolerance) then
                edge_rz = 0.5_dp*(corner_rz(:, vertex_1) &
                    + corner_rz(:, vertex_2))
                return
            end if
        end do
        side = 0
    end subroutine collapsed_element_side

    pure logical function is_jorek_axis_target(locator, target_rz)
        type(jorek_locator_t), intent(in) :: locator
        real(dp), intent(in) :: target_rz(2)

        real(dp) :: scale, tolerance

        is_jorek_axis_target = .false.
        if (.not. locator%has_axis) return
        if (.not. allocated(locator%axis_elements)) return
        scale = max(1.0_dp, maxval(abs(locator%axis_rz)))
        tolerance = 64.0_dp*epsilon(1.0_dp)*scale
        is_jorek_axis_target = &
            maxval(abs(target_rz - locator%axis_rz)) <= tolerance
    end function is_jorek_axis_target

    pure subroutine build_element_bounds(data, element, bounds, ierr)
        type(jorek_restart_t), intent(in) :: data
        integer, intent(in) :: element
        real(dp), intent(out) :: bounds(2, 2)
        integer, intent(out) :: ierr

        real(dp) :: coefficient, controls(2, 4, 4), s_weights(4), t_weights(4)
        integer :: coordinate, degree, i, j, node, vertex

        controls = 0.0_dp
        ierr = 0
        do vertex = 1, 4
            node = data%vertex(element, vertex)
            if (node < 1 .or. node > data%n_nodes) then
                ierr = 2
                return
            end if
            do degree = 1, 4
                call jorek_bernstein_weights(vertex, degree, s_weights, &
                    t_weights)
                do coordinate = 1, 2
                    coefficient = data%x(node, 1, degree, coordinate) &
                        *data%size(element, vertex, degree)
                    do j = 1, 4
                        do i = 1, 4
                            controls(coordinate, i, j) = &
                                controls(coordinate, i, j) &
                                + coefficient*s_weights(i)*t_weights(j)
                        end do
                    end do
                end do
            end do
        end do
        do coordinate = 1, 2
            bounds(coordinate, 1) = minval(controls(coordinate, :, :))
            bounds(coordinate, 2) = maxval(controls(coordinate, :, :))
        end do
    end subroutine build_element_bounds

    pure subroutine jorek_bernstein_weights(vertex, degree, s_weights, &
            t_weights)
        integer, intent(in) :: vertex, degree
        real(dp), intent(out) :: s_weights(4), t_weights(4)

        integer, parameter :: s_end(4) = [1, 2, 2, 1]
        integer, parameter :: t_end(4) = [1, 1, 2, 2]
        integer, parameter :: s_kind(4) = [1, 2, 1, 2]
        integer, parameter :: t_kind(4) = [1, 1, 2, 2]

        call endpoint_bernstein_weights(s_end(vertex), s_kind(degree), &
            s_weights)
        call endpoint_bernstein_weights(t_end(vertex), t_kind(degree), &
            t_weights)
    end subroutine jorek_bernstein_weights

    pure subroutine endpoint_bernstein_weights(endpoint, kind, weights)
        integer, intent(in) :: endpoint, kind
        real(dp), intent(out) :: weights(4)

        weights = 0.0_dp
        if (endpoint == 1 .and. kind == 1) then
            weights(1:2) = 1.0_dp
        else if (endpoint == 1) then
            weights(2) = 1.0_dp
        else if (kind == 1) then
            weights(3:4) = 1.0_dp
        else
            weights(3) = -1.0_dp
        end if
    end subroutine endpoint_bernstein_weights

    pure subroutine locate_jorek_element_indexed(data, locator, target_rz, &
            element, s, t, ierr)
        type(jorek_restart_t), intent(in) :: data
        type(jorek_locator_t), intent(in) :: locator
        real(dp), intent(in) :: target_rz(2)
        integer, intent(out) :: element
        real(dp), intent(out) :: s, t
        integer, intent(out) :: ierr

        real(dp) :: scale, tolerance
        integer :: bin, candidate, candidate_ierr, index

        element = 0
        s = 0.0_dp
        t = 0.0_dp
        ierr = 1
        if (.not. valid_locator(locator, data%n_elements)) then
            ierr = 2
            return
        end if
        bin = locator_bin(locator, target_rz)
        if (bin == 0) return
        do index = locator%bin_offsets(bin), locator%bin_offsets(bin + 1) - 1
            candidate = locator%bin_elements(index)
            scale = max(1.0_dp, maxval(abs(target_rz)), &
                maxval(abs(locator%bounds(:, :, candidate))))
            tolerance = 64.0_dp*epsilon(1.0_dp)*scale
            if (any(target_rz < locator%bounds(:, 1, candidate) - tolerance) &
                .or. any(target_rz > locator%bounds(:, 2, candidate) &
                + tolerance)) cycle
            call invert_jorek_element(data, candidate, target_rz, s, t, &
                candidate_ierr)
            if (candidate_ierr == 0) then
                element = candidate
                ierr = 0
                return
            end if
            if (candidate_ierr == 2) then
                ierr = 2
                return
            end if
        end do
    end subroutine locate_jorek_element_indexed

    pure integer function jorek_locator_candidate_count(locator, target_rz)
        type(jorek_locator_t), intent(in) :: locator
        real(dp), intent(in) :: target_rz(2)

        integer :: bin

        jorek_locator_candidate_count = 0
        if (.not. valid_locator(locator, locator%n_elements)) return
        bin = locator_bin(locator, target_rz)
        if (bin == 0) return
        jorek_locator_candidate_count = locator%bin_offsets(bin + 1) &
            - locator%bin_offsets(bin)
    end function jorek_locator_candidate_count

    subroutine build_locator_bins(locator, ierr)
        type(jorek_locator_t), intent(inout) :: locator
        integer, intent(out) :: ierr

        integer, allocatable :: counts(:), cursor(:)
        integer :: bin, element, first(2), last(2), n_total

        ierr = 0
        call set_locator_grid(locator)
        n_total = product(locator%n_bins)
        allocate(counts(n_total), cursor(n_total))
        counts = 0
        do element = 1, locator%n_elements
            call element_bin_range(locator, element, first, last)
            call add_element_bins(locator%n_bins, first, last, counts)
        end do
        allocate(locator%bin_offsets(n_total + 1))
        locator%bin_offsets(1) = 1
        do bin = 1, n_total
            locator%bin_offsets(bin + 1) = locator%bin_offsets(bin) + counts(bin)
        end do
        allocate(locator%bin_elements(locator%bin_offsets(n_total + 1) - 1))
        cursor = locator%bin_offsets(:n_total)
        do element = 1, locator%n_elements
            call element_bin_range(locator, element, first, last)
            call store_element_bins(locator%n_bins, first, last, element, &
                cursor, locator%bin_elements)
        end do
    end subroutine build_locator_bins

    subroutine set_locator_grid(locator)
        type(jorek_locator_t), intent(inout) :: locator

        real(dp) :: extent(2), ratio
        integer :: coordinate, element
        real(dp) :: lower(2), upper(2)

        locator%grid_lower = huge(1.0_dp)
        upper = -huge(1.0_dp)
        do element = 1, locator%n_elements
            call expanded_bounds(locator, element, lower, extent)
            locator%grid_lower = min(locator%grid_lower, lower)
            upper = max(upper, extent)
        end do
        extent = upper - locator%grid_lower
        if (all(extent > 0.0_dp)) then
            ratio = extent(1)/extent(2)
            locator%n_bins(1) = max(1, min(locator%n_elements, &
                nint(sqrt(real(locator%n_elements, dp)*ratio))))
            locator%n_bins(2) = max(1, ceiling(real(locator%n_elements, dp) &
                /locator%n_bins(1)))
        else
            locator%n_bins = 1
            coordinate = maxloc(extent, 1)
            locator%n_bins(coordinate) = locator%n_elements
        end if
        do coordinate = 1, 2
            locator%bin_width(coordinate) = extent(coordinate) &
                /locator%n_bins(coordinate)
            if (locator%bin_width(coordinate) <= 0.0_dp) &
                locator%bin_width(coordinate) = 1.0_dp
        end do
    end subroutine set_locator_grid

    pure subroutine expanded_bounds(locator, element, lower, upper)
        type(jorek_locator_t), intent(in) :: locator
        integer, intent(in) :: element
        real(dp), intent(out) :: lower(2), upper(2)

        real(dp) :: scale, tolerance

        scale = max(1.0_dp, maxval(abs(locator%bounds(:, :, element))))
        tolerance = 64.0_dp*epsilon(1.0_dp)*scale
        lower = locator%bounds(:, 1, element) - tolerance
        upper = locator%bounds(:, 2, element) + tolerance
    end subroutine expanded_bounds

    pure subroutine element_bin_range(locator, element, first, last)
        type(jorek_locator_t), intent(in) :: locator
        integer, intent(in) :: element
        integer, intent(out) :: first(2), last(2)

        real(dp) :: lower(2), upper(2)
        integer :: coordinate

        call expanded_bounds(locator, element, lower, upper)
        do coordinate = 1, 2
            first(coordinate) = coordinate_bin(locator, coordinate, &
                lower(coordinate))
            last(coordinate) = coordinate_bin(locator, coordinate, &
                upper(coordinate))
        end do
    end subroutine element_bin_range

    pure integer function coordinate_bin(locator, coordinate, value)
        type(jorek_locator_t), intent(in) :: locator
        integer, intent(in) :: coordinate
        real(dp), intent(in) :: value

        coordinate_bin = int((value - locator%grid_lower(coordinate)) &
            /locator%bin_width(coordinate)) + 1
        coordinate_bin = min(locator%n_bins(coordinate), &
            max(1, coordinate_bin))
    end function coordinate_bin

    pure integer function locator_bin(locator, target_rz)
        type(jorek_locator_t), intent(in) :: locator
        real(dp), intent(in) :: target_rz(2)

        real(dp) :: upper(2)
        integer :: indices(2)

        upper = locator%grid_lower + locator%bin_width*locator%n_bins
        locator_bin = 0
        if (any(target_rz < locator%grid_lower) .or. any(target_rz > upper)) return
        indices(1) = coordinate_bin(locator, 1, target_rz(1))
        indices(2) = coordinate_bin(locator, 2, target_rz(2))
        locator_bin = (indices(2) - 1)*locator%n_bins(1) + indices(1)
    end function locator_bin

    pure subroutine add_element_bins(n_bins, first, last, counts)
        integer, intent(in) :: n_bins(2), first(2), last(2)
        integer, intent(inout) :: counts(:)

        integer :: bin, i, j

        do j = first(2), last(2)
            do i = first(1), last(1)
                bin = (j - 1)*n_bins(1) + i
                counts(bin) = counts(bin) + 1
            end do
        end do
    end subroutine add_element_bins

    pure subroutine store_element_bins(n_bins, first, last, element, cursor, &
            elements)
        integer, intent(in) :: n_bins(2), first(2), last(2), element
        integer, intent(inout) :: cursor(:), elements(:)

        integer :: bin, i, j

        do j = first(2), last(2)
            do i = first(1), last(1)
                bin = (j - 1)*n_bins(1) + i
                elements(cursor(bin)) = element
                cursor(bin) = cursor(bin) + 1
            end do
        end do
    end subroutine store_element_bins

    pure logical function valid_locator(locator, n_elements)
        type(jorek_locator_t), intent(in) :: locator
        integer, intent(in) :: n_elements

        integer :: n_total

        valid_locator = .false.
        if (n_elements < 1 .or. locator%n_elements /= n_elements) return
        if (any(locator%n_bins <= 0)) return
        if (.not. allocated(locator%bounds) &
            .or. .not. allocated(locator%bin_offsets) &
            .or. .not. allocated(locator%bin_elements)) return
        if (any(shape(locator%bounds) /= [2, 2, n_elements])) return
        n_total = product(locator%n_bins)
        if (size(locator%bin_offsets) /= n_total + 1) return
        if (size(locator%bin_elements) &
            /= locator%bin_offsets(n_total + 1) - 1) return
        if (locator%has_axis) then
            if (.not. allocated(locator%axis_elements) &
                .or. .not. allocated(locator%axis_sides)) return
            if (size(locator%axis_elements) < 1 &
                .or. size(locator%axis_sides) &
                /= size(locator%axis_elements)) return
        end if
        valid_locator = .true.
    end function valid_locator

    pure subroutine invert_jorek_element(data, element, target_rz, s, t, ierr)
        type(jorek_restart_t), intent(in) :: data
        integer, intent(in) :: element
        real(dp), intent(in) :: target_rz(2)
        real(dp), intent(out) :: s, t
        integer, intent(out) :: ierr

        real(dp), parameter :: seeds(2, 9) = reshape([ &
            0.5_dp, 0.5_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
            1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.5_dp, 0.0_dp, &
            1.0_dp, 0.5_dp, 0.5_dp, 1.0_dp, 0.0_dp, 0.5_dp], [2, 9])
        integer, parameter :: max_iterations = 30

        real(dp) :: delta(2), determinant, determinant_scale, residual(2)
        real(dp) :: rz(2), rz_st(2, 2), scale, tolerance
        integer :: geometry_ierr, iteration, seed

        s = 0.0_dp
        t = 0.0_dp
        ierr = 1
        do seed = 1, size(seeds, 2)
            s = seeds(1, seed)
            t = seeds(2, seed)
            do iteration = 1, max_iterations
                call evaluate_jorek_geometry(data, element, s, t, rz, rz_st, &
                    geometry_ierr)
                if (geometry_ierr /= 0) then
                    ierr = 2
                    return
                end if
                residual = target_rz - rz
                scale = max(1.0_dp, maxval(abs(target_rz)), maxval(abs(rz)))
                tolerance = 1024.0_dp*epsilon(1.0_dp)**0.75_dp*scale
                if (maxval(abs(residual)) <= tolerance) then
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
                s = min(1.0_dp, max(0.0_dp, s + delta(1)))
                t = min(1.0_dp, max(0.0_dp, t + delta(2)))
            end do
        end do
    end subroutine invert_jorek_element

    pure subroutine cubic_endpoint_basis(u, basis, derivative)
        real(dp), intent(in) :: u
        real(dp), intent(out) :: basis(2, 2), derivative(2, 2)

        basis(1, 1) = 2.0_dp*u**3 - 3.0_dp*u**2 + 1.0_dp
        basis(1, 2) = 3.0_dp*(u**3 - 2.0_dp*u**2 + u)
        basis(2, 1) = -2.0_dp*u**3 + 3.0_dp*u**2
        basis(2, 2) = 3.0_dp*(u**2 - u**3)

        derivative(1, 1) = 6.0_dp*u**2 - 6.0_dp*u
        derivative(1, 2) = 9.0_dp*u**2 - 12.0_dp*u + 3.0_dp
        derivative(2, 1) = -6.0_dp*u**2 + 6.0_dp*u
        derivative(2, 2) = 6.0_dp*u - 9.0_dp*u**2
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
