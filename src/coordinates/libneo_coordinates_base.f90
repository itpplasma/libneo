module libneo_coordinates_base
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_coordinate_conventions
    implicit none
    private

    public :: coordinate_system_t
    public :: UNKNOWN, CYL, VMEC, BOOZER
    public :: RHO_TOR, RHO_POL, PSI_TOR_NORM, PSI_POL_NORM

    public :: chartmap_netcdf_spec, chartmap_from_cyl_ierr_spec
    public :: chartmap_from_cyl_ok, chartmap_from_cyl_err_max_iter
    public :: chartmap_from_cyl_err_singular, chartmap_from_cyl_err_out_of_bounds
    public :: chartmap_from_cyl_err_invalid
    public :: refcoords_file_unknown, refcoords_file_chartmap, refcoords_file_vmec_wout

    character(len=*), parameter :: chartmap_netcdf_spec = &
        "Chartmap NetCDF conventions:"//new_line('a')// &
        "- Dimensions: rho, theta, zeta"//new_line('a')// &
        "- Variables (required):"//new_line('a')// &
        "  rho(rho), theta(theta)"//new_line('a')// &
        "  zeta(zeta)"//new_line('a')// &
        "  x(zeta,theta,rho)"//new_line('a')// &
        "  y(zeta,theta,rho)"//new_line('a')// &
        "  z(zeta,theta,rho)"//new_line('a')// &
        "- Optional: num_field_periods (integer >= 1)"//new_line('a')// &
        "- Required: zeta_convention (global attribute)"//new_line('a')// &
        "  allowed: cyl, vmec"//new_line('a')// &
        "- Ranges:"//new_line('a')// &
        "  rho in [0,1]"//new_line('a')// &
        "  theta in [0,2pi)"//new_line('a')// &
        "  zeta in [0,2pi/num_field_periods)"//new_line('a')// &
        "periodic dims exclude endpoint"//new_line('a')// &
        "- Storage order:"//new_line('a')// &
        "  file dims (zeta,theta,rho)"//new_line('a')// &
        "Fortran reads x(rho,theta,zeta)"//new_line('a')// &
        "- Units: x,y,z in cm"

    integer, parameter :: chartmap_from_cyl_ok = 0
    integer, parameter :: chartmap_from_cyl_err_max_iter = 1
    integer, parameter :: chartmap_from_cyl_err_singular = 2
    integer, parameter :: chartmap_from_cyl_err_out_of_bounds = 3
    integer, parameter :: chartmap_from_cyl_err_invalid = 4

    integer, parameter :: refcoords_file_unknown = 0
    integer, parameter :: refcoords_file_chartmap = 1
    integer, parameter :: refcoords_file_vmec_wout = 2

    character(len=*), parameter :: chartmap_from_cyl_ierr_spec = &
        "chartmap_from_cyl ierr codes:"//new_line('a')// &
        "0 ok"//new_line('a')// &
        "1 max iterations or no progress"//new_line('a')// &
        "2 singular normal equations"//new_line('a')// &
        "3 step out of bounds for rho"//new_line('a')// &
        "4 invalid mapping slice"

    type, abstract :: coordinate_system_t
    contains
        procedure(evaluate_cart_if), deferred :: evaluate_cart
        procedure(evaluate_cyl_if), deferred :: evaluate_cyl
        procedure(covariant_basis_if), deferred :: covariant_basis
        procedure(metric_tensor_if), deferred :: metric_tensor
        procedure(from_cyl_if), deferred :: from_cyl
        procedure :: cov_to_cart => coordinate_system_cov_to_cart
        procedure :: ctr_to_cart => coordinate_system_ctr_to_cart
    end type coordinate_system_t

    abstract interface
        subroutine evaluate_cart_if(self, u, x)
            import :: coordinate_system_t, dp
            class(coordinate_system_t), intent(in) :: self
            real(dp), intent(in) :: u(3)
            real(dp), intent(out) :: x(3)
        end subroutine

        subroutine evaluate_cyl_if(self, u, x)
            import :: coordinate_system_t, dp
            class(coordinate_system_t), intent(in) :: self
            real(dp), intent(in) :: u(3)
            real(dp), intent(out) :: x(3)
        end subroutine

        subroutine covariant_basis_if(self, u, e_cov)
            import :: coordinate_system_t, dp
            class(coordinate_system_t), intent(in) :: self
            real(dp), intent(in) :: u(3)
            real(dp), intent(out) :: e_cov(3, 3)
        end subroutine

        subroutine metric_tensor_if(self, u, g, ginv, sqrtg)
            import :: coordinate_system_t, dp
            class(coordinate_system_t), intent(in) :: self
            real(dp), intent(in) :: u(3)
            real(dp), intent(out) :: g(3, 3), ginv(3, 3), sqrtg
        end subroutine

        subroutine from_cyl_if(self, xcyl, u, ierr)
            import :: coordinate_system_t, dp
            class(coordinate_system_t), intent(in) :: self
            real(dp), intent(in) :: xcyl(3)
            real(dp), intent(out) :: u(3)
            integer, intent(out) :: ierr
        end subroutine
    end interface

contains

    subroutine coordinate_system_cov_to_cart(self, u, v_cov, v_cart)
        class(coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(in) :: v_cov(3)
        real(dp), intent(out) :: v_cart(3)

        real(dp) :: e_cov(3, 3), g(3, 3), ginv(3, 3), sqrtg
        real(dp) :: v_ctr(3)

        call self%metric_tensor(u, g, ginv, sqrtg)
        v_ctr = matmul(ginv, v_cov)
        call self%covariant_basis(u, e_cov)
        v_cart = matmul(e_cov, v_ctr)
    end subroutine coordinate_system_cov_to_cart

    subroutine coordinate_system_ctr_to_cart(self, u, v_ctr, v_cart)
        class(coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(in) :: v_ctr(3)
        real(dp), intent(out) :: v_cart(3)

        real(dp) :: e_cov(3, 3)

        call self%covariant_basis(u, e_cov)
        v_cart = matmul(e_cov, v_ctr)
    end subroutine coordinate_system_ctr_to_cart

end module libneo_coordinates_base
