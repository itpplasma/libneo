module geoflux_coordinates

    use, intrinsic :: iso_fortran_env, only : dp => real64
    use cylindrical_cartesian, only : cyl_to_cart, cart_to_cyl
    use math_constants, only : pi
    use geqdsk_tools, only : geqdsk_t, geqdsk_read, geqdsk_standardise, geqdsk_deinit
    use binsrc_sub, only : binsrc
    use interpolate, only : SplineData1D, construct_splines_1d, &
        destroy_splines_1d, evaluate_splines_1d

    implicit none

    real(dp), parameter :: tol_s = 1.0d-10
    real(dp), parameter :: tol_root = 1.0d-8
    integer, parameter :: max_bracket_iter = 200
    integer, parameter :: max_bisect_iter = 60
    integer, parameter :: default_ns_cache = 65
    integer, parameter :: default_ntheta_cache = 128

    abstract interface
        function psi_evaluator_i(R, Z) result(psi)
            import :: dp
            real(dp), intent(in) :: R, Z
            real(dp) :: psi
        end function psi_evaluator_i
    end interface

    type :: geoflux_context_t
        type(geqdsk_t) :: geqdsk
        procedure(psi_evaluator_i), pointer, nopass :: psi_eval => null()
        logical :: use_geqdsk = .true.
        real(dp) :: psi_axis = 0.0_dp
        real(dp) :: psi_sep = 0.0_dp
        real(dp) :: psi_tor_edge = 0.0_dp
        real(dp) :: R_axis = 0.0_dp
        real(dp) :: Z_axis = 0.0_dp
        real(dp) :: R_min = 0.0_dp
        real(dp) :: R_max = 0.0_dp
        real(dp) :: Z_min = 0.0_dp
        real(dp) :: Z_max = 0.0_dp
        real(dp) :: ray_step = 0.0_dp
        real(dp) :: max_radius = 0.0_dp
        logical :: has_toroidal = .false.
        logical :: initialised = .false.
        real(dp), allocatable :: psi_grid(:)
        real(dp), allocatable :: s_grid(:)
        integer :: ns_cache = 0
        integer :: ntheta_cache = 0
        real(dp), allocatable :: s_nodes(:)
        real(dp), allocatable :: theta_nodes(:)
        real(dp), allocatable :: R_cache(:,:)
        real(dp), allocatable :: Z_cache(:,:)
        type(SplineData1D) :: psi_of_s_spline
        type(SplineData1D) :: s_of_psi_spline
        logical :: psi_of_s_ready = .false.
        logical :: s_of_psi_ready = .false.
        logical :: cache_built = .false.
    end type geoflux_context_t

    type(geoflux_context_t), save :: ctx

    abstract interface
        subroutine transform_i(xfrom, xto, dxto_dxfrom)
            import :: dp
            real(dp), intent(in) :: xfrom(3)
            real(dp), intent(out) :: xto(3)
            real(dp), intent(out), optional :: dxto_dxfrom(3,3)
        end subroutine transform_i
    end interface

contains

    subroutine init_geoflux_coordinates(filename, ns_cache_in, ntheta_cache_in)
        character(len=*), intent(in) :: filename
        integer, intent(in), optional :: ns_cache_in, ntheta_cache_in
        integer :: ns_val, ntheta_val

        call cleanup_context()

        ctx%use_geqdsk = .true.
        ctx%psi_eval => null()

        call geqdsk_read(ctx%geqdsk, filename)
        call geqdsk_standardise(ctx%geqdsk)

        ctx%psi_axis = ctx%geqdsk%simag
        ctx%psi_sep = ctx%geqdsk%sibry
        ctx%R_axis = ctx%geqdsk%rmaxis
        ctx%Z_axis = ctx%geqdsk%zmaxis

        ctx%R_min = minval(ctx%geqdsk%R_eqd)
        ctx%R_max = maxval(ctx%geqdsk%R_eqd)
        ctx%Z_min = minval(ctx%geqdsk%Z_eqd)
        ctx%Z_max = maxval(ctx%geqdsk%Z_eqd)

        ctx%max_radius = sqrt((ctx%R_max - ctx%R_axis)**2 + (ctx%Z_max - ctx%Z_axis)**2)
        ctx%max_radius = ctx%max_radius + 1.0_dp
        ctx%ray_step = max(ctx%max_radius/100.0_dp, 1.0d-3)

        call build_radial_mapping()
        call build_radial_splines()

        ns_val = default_ns_cache
        ntheta_val = default_ntheta_cache
        if (present(ns_cache_in)) ns_val = max(3, ns_cache_in)
        if (present(ntheta_cache_in)) ntheta_val = max(4, ntheta_cache_in)
        ctx%ns_cache = ns_val
        ctx%ntheta_cache = ntheta_val

        call build_flux_surface_cache()

        ctx%initialised = .true.
    end subroutine init_geoflux_coordinates

    subroutine initialize_analytical_geoflux(psi_evaluator, R0, epsilon, kappa, delta, &
                                            psi_axis, psi_sep, &
                                            ns_cache_in, ntheta_cache_in, npsi_grid_in)
        procedure(psi_evaluator_i) :: psi_evaluator
        real(dp), intent(in) :: R0, epsilon, kappa, delta
        real(dp), intent(in) :: psi_axis, psi_sep
        integer, intent(in), optional :: ns_cache_in, ntheta_cache_in, npsi_grid_in
        integer :: ns_val, ntheta_val, npsi_val
        integer :: i
        real(dp) :: a_minor

        call cleanup_context()

        ctx%use_geqdsk = .false.
        ctx%psi_eval => psi_evaluator

        a_minor = R0 * epsilon
        ctx%R_axis = R0
        ctx%Z_axis = 0.0_dp
        ctx%psi_axis = psi_axis
        ctx%psi_sep = psi_sep

        ctx%R_min = R0 - 1.5_dp * a_minor * (1.0_dp + abs(delta))
        ctx%R_max = R0 + 1.5_dp * a_minor * (1.0_dp + abs(delta))
        ctx%Z_min = -1.5_dp * a_minor * kappa
        ctx%Z_max = 1.5_dp * a_minor * kappa

        ctx%max_radius = sqrt((ctx%R_max - ctx%R_axis)**2 + (ctx%Z_max - ctx%Z_axis)**2)
        ctx%max_radius = ctx%max_radius + 1.0_dp
        ctx%ray_step = max(ctx%max_radius/100.0_dp, 1.0d-3)

        npsi_val = 129
        if (present(npsi_grid_in)) npsi_val = max(10, npsi_grid_in)

        allocate(ctx%psi_grid(npsi_val))
        allocate(ctx%s_grid(npsi_val))

        do i = 1, npsi_val
            ctx%s_grid(i) = real(i - 1, dp) / real(npsi_val - 1, dp)
            ctx%psi_grid(i) = ctx%psi_axis + ctx%s_grid(i) * (ctx%psi_sep - ctx%psi_axis)
        end do

        ctx%s_grid(1) = 0.0_dp
        ctx%s_grid(npsi_val) = 1.0_dp
        ctx%has_toroidal = .false.

        call build_radial_splines()

        ns_val = default_ns_cache
        ntheta_val = default_ntheta_cache
        if (present(ns_cache_in)) ns_val = max(3, ns_cache_in)
        if (present(ntheta_cache_in)) ntheta_val = max(4, ntheta_cache_in)
        ctx%ns_cache = ns_val
        ctx%ntheta_cache = ntheta_val

        call build_flux_surface_cache()

        ctx%initialised = .true.
    end subroutine initialize_analytical_geoflux

    subroutine geoflux_to_cyl(xfrom, xto, dxto_dxfrom)
        real(dp), intent(in) :: xfrom(3)
        real(dp), intent(out) :: xto(3)
        real(dp), intent(out), optional :: dxto_dxfrom(3,3)

        real(dp) :: s_val, theta_val, phi_val
        real(dp) :: R_val, Z_val

        call ensure_initialised()

        s_val = clamp01(xfrom(1))
        theta_val = xfrom(2)
        phi_val = xfrom(3)

        if (ctx%cache_built) then
            call interpolate_cached_surface(s_val, theta_val, R_val, Z_val)
        else
            call locate_flux_surface(s_val, theta_val, R_val, Z_val)
        end if

        xto(1) = R_val
        xto(2) = phi_val
        xto(3) = Z_val

        if (present(dxto_dxfrom)) then
            call assign_geoflux_to_cyl_jacobian(s_val, theta_val, phi_val, R_val, Z_val, dxto_dxfrom)
        end if
    end subroutine geoflux_to_cyl

    subroutine geoflux_to_cart(xfrom, xto, dxto_dxfrom)
        real(dp), intent(in) :: xfrom(3)
        real(dp), intent(out) :: xto(3)
        real(dp), intent(out), optional :: dxto_dxfrom(3,3)

        real(dp) :: xcyl(3)
        real(dp) :: dxcyl_dxgeo(3,3)
        real(dp) :: dxcart_dxcyl(3,3)

        if (present(dxto_dxfrom)) then
            call geoflux_to_cyl(xfrom, xcyl, dxcyl_dxgeo)
            call cyl_to_cart(xcyl, xto, dxcart_dxcyl)
            dxto_dxfrom = matmul(dxcart_dxcyl, dxcyl_dxgeo)
        else
            call geoflux_to_cyl(xfrom, xcyl)
            call cyl_to_cart(xcyl, xto)
        end if
    end subroutine geoflux_to_cart

    subroutine cyl_to_geoflux(xfrom, xto, dxto_dxfrom)
        real(dp), intent(in) :: xfrom(3)
        real(dp), intent(out) :: xto(3)
        real(dp), intent(out), optional :: dxto_dxfrom(3,3)

        real(dp) :: R_val, phi_in, Z_val
        real(dp) :: s_val, theta_val

        call ensure_initialised()

        R_val = clamp(xfrom(1), ctx%R_min, ctx%R_max)
        phi_in = xfrom(2)
        Z_val = clamp(xfrom(3), ctx%Z_min, ctx%Z_max)

        s_val = s_from_position(R_val, Z_val)
        theta_val = atan2(Z_val - ctx%Z_axis, R_val - ctx%R_axis)

        xto(1) = s_val
        xto(2) = theta_val
        xto(3) = phi_in

        if (present(dxto_dxfrom)) then
            call assign_cyl_to_geoflux_jacobian(R_val, phi_in, Z_val, dxto_dxfrom)
        end if
    end subroutine cyl_to_geoflux

    function get_transform(from, to) result(func)
        procedure(transform_i), pointer :: func
        character(*), intent(in) :: from, to

        func => null()

        select case (trim(from))
        case ('geoflux')
            select case (trim(to))
            case ('cyl')
                func => geoflux_to_cyl
            case ('cart')
                func => geoflux_to_cart
            case default
                call report_unknown(from, to)
            end select
        case ('cyl')
            select case (trim(to))
            case ('geoflux')
                func => cyl_to_geoflux
            case ('cart')
                func => cyl_to_cart
            case default
                call report_unknown(from, to)
            end select
        case ('cart')
            if (trim(to) == 'cyl') then
                func => cart_to_cyl
            else
                call report_unknown(from, to)
            end if
        case default
            call report_unknown(from, to)
        end select
    end function get_transform

    subroutine report_unknown(from, to)
        character(*), intent(in) :: from
        character(*), intent(in), optional :: to

        if (present(to)) then
            write(*,*) 'Unknown transform from ', from, ' to ', to
        else
            write(*,*) 'Unknown transform from ', from
        end if
        error stop
    end subroutine report_unknown

    subroutine cleanup_context()
        if (ctx%initialised) then
            call geqdsk_deinit(ctx%geqdsk)
        end if
        if (ctx%psi_of_s_ready) then
            call destroy_splines_1d(ctx%psi_of_s_spline)
            ctx%psi_of_s_ready = .false.
        end if
        if (ctx%s_of_psi_ready) then
            call destroy_splines_1d(ctx%s_of_psi_spline)
            ctx%s_of_psi_ready = .false.
        end if
        if (allocated(ctx%psi_grid)) deallocate(ctx%psi_grid)
        if (allocated(ctx%s_grid)) deallocate(ctx%s_grid)
        if (allocated(ctx%s_nodes)) deallocate(ctx%s_nodes)
        if (allocated(ctx%theta_nodes)) deallocate(ctx%theta_nodes)
        if (allocated(ctx%R_cache)) deallocate(ctx%R_cache)
        if (allocated(ctx%Z_cache)) deallocate(ctx%Z_cache)
        ctx%initialised = .false.
        ctx%has_toroidal = .false.
        ctx%cache_built = .false.
    end subroutine cleanup_context

    subroutine ensure_initialised()
        if (.not. ctx%initialised) then
            write(*,*) 'geoflux_coordinates: call init_geoflux_coordinates first.'
            error stop
        end if
    end subroutine ensure_initialised

    subroutine build_radial_mapping()
        integer :: npsi
        real(dp), allocatable :: tor_flux(:)
        real(dp) :: denom, scale

        npsi = size(ctx%geqdsk%psi_eqd)
        if (npsi < 2) then
            write(*,*) 'geoflux_coordinates: insufficient psi grid.'
            error stop
        end if

        allocate(ctx%psi_grid(npsi))
        allocate(ctx%s_grid(npsi))
        allocate(tor_flux(npsi))

        ctx%psi_grid = ctx%geqdsk%psi_eqd
        tor_flux = 0.0_dp

        ctx%psi_tor_edge = 0.0_dp

        if (allocated(ctx%geqdsk%qpsi)) then
            call integrate_toroidal_flux(tor_flux)
        end if

        ctx%psi_tor_edge = tor_flux(npsi)

        if (abs(ctx%psi_tor_edge) > 0.0_dp) then
            ctx%has_toroidal = .true.
            scale = 1.0_dp / abs(ctx%psi_tor_edge)
            ctx%s_grid = tor_flux * scale
            if (ctx%psi_tor_edge < 0.0_dp) ctx%s_grid = -ctx%s_grid
        else
            ctx%has_toroidal = .false.
        end if

        if (.not. ctx%has_toroidal) then
            denom = ctx%psi_sep - ctx%psi_axis
            if (abs(denom) < 1.0d-12) then
                ctx%s_grid = 0.0_dp
            else
                ctx%s_grid = (ctx%psi_grid - ctx%psi_axis) / denom
            end if
        end if

        ctx%s_grid(1) = 0.0_dp
        ctx%s_grid(npsi) = 1.0_dp

        deallocate(tor_flux)
    end subroutine build_radial_mapping

    subroutine build_radial_splines()
        integer :: npsi, i
        real(dp), allocatable :: psi_samples(:)
        real(dp) :: s_uniform
        real(dp) :: psi_min, psi_max
        integer, parameter :: spline_order = 5

        npsi = size(ctx%psi_grid)
        if (npsi < 2) then
            ctx%psi_of_s_ready = .false.
            ctx%s_of_psi_ready = .false.
            return
        end if

        if (ctx%psi_of_s_ready) then
            call destroy_splines_1d(ctx%psi_of_s_spline)
            ctx%psi_of_s_ready = .false.
        end if
        if (ctx%s_of_psi_ready) then
            call destroy_splines_1d(ctx%s_of_psi_spline)
            ctx%s_of_psi_ready = .false.
        end if

        allocate(psi_samples(npsi))
        do i = 1, npsi
            s_uniform = real(i - 1, dp) / real(npsi - 1, dp)
            psi_samples(i) = linear_interp_monotonic(ctx%s_grid, ctx%psi_grid, s_uniform)
        end do

        call construct_splines_1d(0.0_dp, 1.0_dp, psi_samples, spline_order, &
            .false., ctx%psi_of_s_spline)
        ctx%psi_of_s_ready = .true.
        deallocate(psi_samples)

        psi_min = minval(ctx%psi_grid)
        psi_max = maxval(ctx%psi_grid)
        call construct_splines_1d(psi_min, psi_max, ctx%s_grid, spline_order, &
            .false., ctx%s_of_psi_spline)
        ctx%s_of_psi_ready = .true.
    end subroutine build_radial_splines

    subroutine build_flux_surface_cache()
        integer :: is, itheta
        real(dp) :: s_val, theta_val
        real(dp) :: R_tmp, Z_tmp
        real(dp) :: two_pi, dtheta

        if (ctx%ns_cache < 3 .or. ctx%ntheta_cache < 4) then
            ctx%cache_built = .false.
            return
        end if

        if (allocated(ctx%s_nodes)) deallocate(ctx%s_nodes)
        if (allocated(ctx%theta_nodes)) deallocate(ctx%theta_nodes)
        if (allocated(ctx%R_cache)) deallocate(ctx%R_cache)
        if (allocated(ctx%Z_cache)) deallocate(ctx%Z_cache)

        allocate(ctx%s_nodes(ctx%ns_cache))
        allocate(ctx%theta_nodes(ctx%ntheta_cache))
        allocate(ctx%R_cache(ctx%ns_cache, ctx%ntheta_cache))
        allocate(ctx%Z_cache(ctx%ns_cache, ctx%ntheta_cache))

        do is = 1, ctx%ns_cache
            ctx%s_nodes(is) = real(is - 1, dp) / &
                real(ctx%ns_cache - 1, dp)
        end do

        two_pi = 2.0_dp * pi
        dtheta = two_pi / real(ctx%ntheta_cache, dp)
        do itheta = 1, ctx%ntheta_cache
            ctx%theta_nodes(itheta) = -pi + real(itheta - 1, dp) * dtheta
        end do

        do itheta = 1, ctx%ntheta_cache
            ctx%R_cache(1, itheta) = ctx%R_axis
            ctx%Z_cache(1, itheta) = ctx%Z_axis
        end do

        do is = 2, ctx%ns_cache
            s_val = ctx%s_nodes(is)
            do itheta = 1, ctx%ntheta_cache
                theta_val = ctx%theta_nodes(itheta)
                call locate_flux_surface(s_val, theta_val, R_tmp, Z_tmp)
                ctx%R_cache(is, itheta) = R_tmp
                ctx%Z_cache(is, itheta) = Z_tmp
            end do
        end do

        ctx%cache_built = .true.
    end subroutine build_flux_surface_cache

    subroutine integrate_toroidal_flux(tor_flux)
        real(dp), intent(inout) :: tor_flux(:)

        integer :: i, npsi
        real(dp) :: dpsi, q_lo, q_hi

        npsi = size(tor_flux)
        tor_flux(1) = 0.0_dp

        do i = 2, npsi
            dpsi = ctx%psi_grid(i) - ctx%psi_grid(i-1)
            q_lo = ctx%geqdsk%qpsi(i-1)
            q_hi = ctx%geqdsk%qpsi(i)
            tor_flux(i) = tor_flux(i-1) + 0.5_dp * (q_lo + q_hi) * dpsi / (2.0_dp * pi)
        end do
    end subroutine integrate_toroidal_flux

    subroutine locate_flux_surface(s_val, theta_val, R_val, Z_val)
        real(dp), intent(in) :: s_val, theta_val
        real(dp), intent(out) :: R_val, Z_val

        real(dp) :: target_psi, target_norm
        real(dp) :: r_low, r_high, flux_low, flux_high, r_mid, flux_mid
        integer :: iter

        if (s_val <= tol_s) then
            R_val = ctx%R_axis
            Z_val = ctx%Z_axis
            return
        end if

        target_psi = psi_from_s(s_val)
        target_norm = normalized_poloidal(target_psi)

        r_low = 0.0_dp
        flux_low = 0.0_dp
        r_high = ctx%ray_step

        flux_high = flux_along_ray(r_high, theta_val)

        do iter = 1, max_bracket_iter
            if (flux_high >= target_norm) exit
            r_high = r_high + ctx%ray_step
            if (r_high >= ctx%max_radius) then
                r_high = ctx%max_radius
                exit
            end if
            flux_high = flux_along_ray(r_high, theta_val)
        end do

        if (flux_high < target_norm) then
            r_high = ctx%max_radius
            flux_high = flux_along_ray(r_high, theta_val)
            if (flux_high < target_norm) then
                R_val = clamp(ctx%R_axis + r_high * cos(theta_val), ctx%R_min, ctx%R_max)
                Z_val = clamp(ctx%Z_axis + r_high * sin(theta_val), ctx%Z_min, ctx%Z_max)
                return
            end if
        end if

        do iter = 1, max_bisect_iter
            r_mid = 0.5_dp * (r_low + r_high)
            flux_mid = flux_along_ray(r_mid, theta_val)

            if (abs(flux_mid - target_norm) <= tol_root) exit

            if (flux_mid >= target_norm) then
                r_high = r_mid
                flux_high = flux_mid
            else
                r_low = r_mid
                flux_low = flux_mid
            end if
        end do

        R_val = clamp(ctx%R_axis + r_high * cos(theta_val), ctx%R_min, ctx%R_max)
        Z_val = clamp(ctx%Z_axis + r_high * sin(theta_val), ctx%Z_min, ctx%Z_max)
    end subroutine locate_flux_surface

    function flux_along_ray(r_val, theta_val) result(s_norm)
        real(dp), intent(in) :: r_val, theta_val
        real(dp) :: s_norm
        real(dp) :: R_pt, Z_pt

        R_pt = ctx%R_axis + r_val * cos(theta_val)
        Z_pt = ctx%Z_axis + r_val * sin(theta_val)
        s_norm = normalized_poloidal(psi_from_position(R_pt, Z_pt))
    end function flux_along_ray

    function psi_from_s(s_val) result(psi_val)
        real(dp), intent(in) :: s_val
        real(dp) :: psi_val

        if (ctx%psi_of_s_ready) then
            call evaluate_splines_1d(ctx%psi_of_s_spline, clamp01(s_val), psi_val)
        else
            psi_val = linear_interp_monotonic(ctx%s_grid, ctx%psi_grid, &
                clamp01(s_val))
        end if
    end function psi_from_s

    function s_from_position(R_val, Z_val) result(s_val)
        real(dp), intent(in) :: R_val, Z_val
        real(dp) :: s_val
        real(dp) :: psi_val

        psi_val = psi_from_position(R_val, Z_val)
        if (ctx%s_of_psi_ready) then
            call evaluate_splines_1d(ctx%s_of_psi_spline, &
                clamp(psi_val, minval(ctx%psi_grid), maxval(ctx%psi_grid)), s_val)
        else
            s_val = linear_interp_monotonic(ctx%psi_grid, ctx%s_grid, psi_val)
        end if
        s_val = clamp01(s_val)
    end function s_from_position

    function psi_from_position(R_val, Z_val) result(psi_val)
        real(dp), intent(in) :: R_val, Z_val
        real(dp) :: psi_val
        integer :: i_hi, j_hi, i_lo, j_lo
        real(dp) :: t_R, t_Z
        real(dp) :: R_clamped, Z_clamped

        if (ctx%use_geqdsk) then
            R_clamped = clamp(R_val, ctx%R_min, ctx%R_max)
            Z_clamped = clamp(Z_val, ctx%Z_min, ctx%Z_max)

            call binsrc(ctx%geqdsk%R_eqd, 1, size(ctx%geqdsk%R_eqd), R_clamped, i_hi)
            call binsrc(ctx%geqdsk%Z_eqd, 1, size(ctx%geqdsk%Z_eqd), Z_clamped, j_hi)

            i_lo = max(1, min(i_hi - 1, size(ctx%geqdsk%R_eqd) - 1))
            i_hi = i_lo + 1
            j_lo = max(1, min(j_hi - 1, size(ctx%geqdsk%Z_eqd) - 1))
            j_hi = j_lo + 1

            t_R = (R_clamped - ctx%geqdsk%R_eqd(i_lo)) / &
                  max(ctx%geqdsk%R_eqd(i_hi) - ctx%geqdsk%R_eqd(i_lo), 1.0d-12)
            t_Z = (Z_clamped - ctx%geqdsk%Z_eqd(j_lo)) / &
                  max(ctx%geqdsk%Z_eqd(j_hi) - ctx%geqdsk%Z_eqd(j_lo), 1.0d-12)

            psi_val = bilinear(ctx%geqdsk%psirz, i_lo, j_lo, t_R, t_Z)
        else
            if (.not. associated(ctx%psi_eval)) then
                error stop 'geoflux_coordinates: psi evaluator not set'
            end if
            psi_val = ctx%psi_eval(R_val, Z_val)
        end if
    end function psi_from_position

    function bilinear(psirz, i_lo, j_lo, t_R, t_Z) result(value)
        real(dp), intent(in) :: psirz(:,:)
        integer, intent(in) :: i_lo, j_lo
        real(dp), intent(in) :: t_R, t_Z
        real(dp) :: value
        real(dp) :: f00, f10, f01, f11

        f00 = psirz(i_lo    , j_lo)
        f10 = psirz(i_lo + 1, j_lo)
        f01 = psirz(i_lo    , j_lo + 1)
        f11 = psirz(i_lo + 1, j_lo + 1)

        value = (1.0_dp - t_R) * (1.0_dp - t_Z) * f00 + &
                t_R * (1.0_dp - t_Z) * f10 + &
                (1.0_dp - t_R) * t_Z * f01 + &
                t_R * t_Z * f11
    end function bilinear

    function normalized_poloidal(psi_val) result(s_norm)
        real(dp), intent(in) :: psi_val
        real(dp) :: s_norm
        real(dp) :: denom

        denom = ctx%psi_sep - ctx%psi_axis
        if (abs(denom) < 1.0d-12) then
            s_norm = 0.0_dp
        else
            s_norm = (psi_val - ctx%psi_axis) / denom
        end if
        s_norm = clamp01(s_norm)
    end function normalized_poloidal

    subroutine interpolate_cached_surface(s_val, theta_val, R_val, Z_val)
        real(dp), intent(in) :: s_val, theta_val
        real(dp), intent(out) :: R_val, Z_val
        real(dp) :: s_use, theta_use, s_pos, theta_pos
        real(dp) :: ws, wt, dtheta, two_pi
        integer :: i_lo, i_hi, j_lo, j_hi

        if (.not. ctx%cache_built) then
            call locate_flux_surface(s_val, theta_val, R_val, Z_val)
            return
        end if

        s_use = clamp01(s_val)
        if (s_use <= tol_s) then
            R_val = ctx%R_axis
            Z_val = ctx%Z_axis
            return
        end if

        s_pos = s_use * real(ctx%ns_cache - 1, dp)
        i_lo = int(floor(s_pos)) + 1
        if (i_lo >= ctx%ns_cache) then
            i_lo = ctx%ns_cache - 1
            i_hi = ctx%ns_cache
            ws = 1.0_dp
        else
            i_hi = i_lo + 1
            ws = s_pos - real(i_lo - 1, dp)
        end if

        two_pi = 2.0_dp * pi
        dtheta = two_pi / real(ctx%ntheta_cache, dp)
        theta_use = wrap_theta(theta_val)
        theta_pos = (theta_use + pi) / dtheta
        j_lo = int(floor(theta_pos)) + 1
        wt = theta_pos - real(j_lo - 1, dp)
        if (j_lo > ctx%ntheta_cache) then
            j_lo = 1
            wt = theta_pos - floor(theta_pos)
        end if
        j_hi = j_lo + 1
        if (j_hi > ctx%ntheta_cache) j_hi = 1

        R_val = (1.0_dp - ws) * (1.0_dp - wt) * ctx%R_cache(i_lo, j_lo) &
            + ws * (1.0_dp - wt) * ctx%R_cache(i_hi, j_lo) &
            + (1.0_dp - ws) * wt * ctx%R_cache(i_lo, j_hi) &
            + ws * wt * ctx%R_cache(i_hi, j_hi)

        Z_val = (1.0_dp - ws) * (1.0_dp - wt) * ctx%Z_cache(i_lo, j_lo) &
            + ws * (1.0_dp - wt) * ctx%Z_cache(i_hi, j_lo) &
            + (1.0_dp - ws) * wt * ctx%Z_cache(i_lo, j_hi) &
            + ws * wt * ctx%Z_cache(i_hi, j_hi)
    end subroutine interpolate_cached_surface

    pure function wrap_theta(theta) result(theta_wrapped)
        real(dp), intent(in) :: theta
        real(dp) :: theta_wrapped
        real(dp) :: two_pi

        two_pi = 2.0_dp * pi
        theta_wrapped = modulo(theta + pi, two_pi) - pi
    end function wrap_theta

    pure function clamp(x, xmin, xmax) result(val)
        real(dp), intent(in) :: x, xmin, xmax
        real(dp) :: val

        val = min(max(x, xmin), xmax)
    end function clamp

    pure function clamp01(x) result(val)
        real(dp), intent(in) :: x
        real(dp) :: val

        val = min(max(x, 0.0_dp), 1.0_dp)
    end function clamp01

    function linear_interp_monotonic(x, y, xval) result(yval)
        real(dp), intent(in) :: x(:), y(:), xval
        real(dp) :: yval
        integer :: n, idx
        real(dp) :: x1, x2, w
        real(dp) :: x_clamped

        n = size(x)
        if (n /= size(y)) then
            error stop 'geoflux_coordinates: inconsistent interpolation arrays.'
        end if

        if (n == 0) then
            yval = 0.0_dp
            return
        end if

        if (n == 1) then
            yval = y(1)
            return
        end if

        do idx = 2, n
            if (x(idx) < x(idx - 1)) then
                error stop 'geoflux_coordinates: interpolation abscissae must be non-decreasing.'
            end if
        end do

        x_clamped = clamp(xval, x(1), x(n))

        if (x_clamped <= x(1)) then
            yval = y(1)
            return
        else if (x_clamped >= x(n)) then
            yval = y(n)
            return
        end if

        call binsrc(x, 1, n, x_clamped, idx)
        idx = max(2, min(n, idx))

        x1 = x(idx - 1)
        x2 = x(idx)

        if (abs(x2 - x1) < 1.0d-14) then
            yval = 0.5_dp * (y(idx - 1) + y(idx))
        else
            w = (x_clamped - x1) / (x2 - x1)
            yval = (1.0_dp - w) * y(idx - 1) + w * y(idx)
        end if
    end function linear_interp_monotonic

    subroutine assign_geoflux_to_cyl_jacobian(s_val, theta_val, phi_val, R_val, Z_val, jac)
        real(dp), intent(in) :: s_val, theta_val, phi_val, R_val, Z_val
        real(dp), intent(out) :: jac(3,3)
        real(dp) :: ds, dt
        real(dp) :: xp(3), xm(3)
        real(dp) :: rp(3), rm(3)

        jac = 0.0_dp

        ds = max(1.0d-5, 1.0d-3 * max(1.0_dp, s_val))
        dt = 1.0d-4

        xp = [clamp01(s_val + ds), theta_val, phi_val]
        xm = [clamp01(s_val - ds), theta_val, phi_val]
        call geoflux_to_cyl_internal(xp, rp)
        call geoflux_to_cyl_internal(xm, rm)
        jac(1,1) = (rp(1) - rm(1)) / (2.0_dp * ds)
        jac(3,1) = (rp(3) - rm(3)) / (2.0_dp * ds)
        jac(2,1) = 0.0_dp

        xp = [s_val, theta_val + dt, phi_val]
        xm = [s_val, theta_val - dt, phi_val]
        call geoflux_to_cyl_internal(xp, rp)
        call geoflux_to_cyl_internal(xm, rm)
        jac(1,2) = (rp(1) - rm(1)) / (2.0_dp * dt)
        jac(3,2) = (rp(3) - rm(3)) / (2.0_dp * dt)
        jac(2,2) = 0.0_dp

        jac(1,3) = 0.0_dp
        jac(2,3) = 1.0_dp
        jac(3,3) = 0.0_dp
    end subroutine assign_geoflux_to_cyl_jacobian

    subroutine geoflux_to_cyl_internal(xfrom, xto)
        real(dp), intent(in) :: xfrom(3)
        real(dp), intent(out) :: xto(3)

        real(dp) :: R_val, Z_val
        if (ctx%cache_built) then
            call interpolate_cached_surface(xfrom(1), xfrom(2), R_val, Z_val)
        else
            call locate_flux_surface(xfrom(1), xfrom(2), R_val, Z_val)
        end if
        xto(1) = R_val
        xto(2) = xfrom(3)
        xto(3) = Z_val
    end subroutine geoflux_to_cyl_internal

    subroutine assign_cyl_to_geoflux_jacobian(R_val, phi_in, Z_val, jac)
        real(dp), intent(in) :: R_val, phi_in, Z_val
        real(dp), intent(out) :: jac(3,3)
        real(dp) :: dR, dZ
        real(dp) :: xp(3), xm(3)
        real(dp) :: sp(3), sm(3)

        jac = 0.0_dp

        dR = max(1.0d-5, 1.0d-3 * max(1.0_dp, abs(R_val)))
        dZ = max(1.0d-5, 1.0d-3 * max(1.0_dp, abs(Z_val)))

        xp = [clamp(R_val + dR, ctx%R_min, ctx%R_max), phi_in, Z_val]
        xm = [clamp(R_val - dR, ctx%R_min, ctx%R_max), phi_in, Z_val]
        call cyl_to_geoflux_internal(xp, sp)
        call cyl_to_geoflux_internal(xm, sm)
        jac(1,1) = (sp(1) - sm(1)) / (2.0_dp * dR)
        jac(2,1) = (sp(2) - sm(2)) / (2.0_dp * dR)
        jac(3,1) = 0.0_dp

        xp = [R_val, phi_in, clamp(Z_val + dZ, ctx%Z_min, ctx%Z_max)]
        xm = [R_val, phi_in, clamp(Z_val - dZ, ctx%Z_min, ctx%Z_max)]
        call cyl_to_geoflux_internal(xp, sp)
        call cyl_to_geoflux_internal(xm, sm)
        jac(1,3) = (sp(1) - sm(1)) / (2.0_dp * dZ)
        jac(2,3) = (sp(2) - sm(2)) / (2.0_dp * dZ)
        jac(3,3) = 0.0_dp

        jac(3,2) = 1.0_dp
        jac(2,2) = 0.0_dp
        jac(1,2) = 0.0_dp
    end subroutine assign_cyl_to_geoflux_jacobian

    subroutine cyl_to_geoflux_internal(xfrom, xto)
        real(dp), intent(in) :: xfrom(3)
        real(dp), intent(out) :: xto(3)

        real(dp) :: s_val, theta_val
        s_val = s_from_position(xfrom(1), xfrom(3))
        theta_val = atan2(xfrom(3) - ctx%Z_axis, xfrom(1) - ctx%R_axis)
        xto(1) = s_val
        xto(2) = theta_val
        xto(3) = xfrom(2)
    end subroutine cyl_to_geoflux_internal

    subroutine geoflux_get_axis(R_axis, Z_axis)
        real(dp), intent(out) :: R_axis, Z_axis

        call ensure_initialised()

        R_axis = ctx%R_axis
        Z_axis = ctx%Z_axis
    end subroutine geoflux_get_axis

end module geoflux_coordinates
