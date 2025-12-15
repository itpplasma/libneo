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

        ierr = 1
        message = "unsupported zeta_convention for vmec mapping"
    end subroutine map_vmec_u_to_chartmap_u

end module libneo_coordinates_mapping
