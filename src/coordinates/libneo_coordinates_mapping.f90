module libneo_coordinates_mapping
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_coordinates, only: vmec_coordinate_system_t, &
                                  chartmap_coordinate_system_t, &
                                  CYL, VMEC
    implicit none
    private

    public :: map_vmec_u_to_chartmap_u

contains

    subroutine map_vmec_u_to_chartmap_u(vmec_cs, ccs, u_vmec, u_chart, ierr, message)
        type(vmec_coordinate_system_t), intent(in) :: vmec_cs
        type(chartmap_coordinate_system_t), intent(in) :: ccs
        real(dp), intent(in) :: u_vmec(3)
        real(dp), intent(out) :: u_chart(3)
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: message

        real(dp) :: xcyl(3)
        real(dp) :: xcart(3)

        ierr = 0
        message = ""
        u_chart = 0.0_dp

        call vmec_cs%evaluate_cyl(u_vmec, xcyl)

        if (ccs%zeta_convention == CYL .or. ccs%zeta_convention == VMEC) then
            call ccs%from_cyl(xcyl, u_chart, ierr)
            if (ierr /= 0) then
                write (message, '(a,i0)') "chartmap from_cyl failed ierr=", ierr
                return
            end if
            return
        end if

        call cyl_to_cart(xcyl, xcart)
        call ccs%from_cart(xcart, u_chart, ierr)
        if (ierr /= 0) then
            write (message, '(a,i0)') "chartmap from_cart failed ierr=", ierr
        end if
    end subroutine map_vmec_u_to_chartmap_u

    pure subroutine cyl_to_cart(xcyl, xcart)
        real(dp), intent(in) :: xcyl(3)
        real(dp), intent(out) :: xcart(3)

        xcart(1) = xcyl(1)*cos(xcyl(2))
        xcart(2) = xcyl(1)*sin(xcyl(2))
        xcart(3) = xcyl(3)
    end subroutine cyl_to_cart

end module libneo_coordinates_mapping
