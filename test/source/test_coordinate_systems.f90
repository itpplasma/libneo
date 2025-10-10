program test_coordinate_systems
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_coordinates
    implicit none

    integer :: nerrors
    logical :: all_passed

    nerrors = 0
    all_passed = .true.

    print *, "Testing coordinate_system_t abstraction..."
    print *, ""

    call test_coordinate_interface(nerrors)
    print *, ""

    if (nerrors > 0) then
        print *, "FAILED: ", nerrors, " error(s) detected"
        all_passed = .false.
    else
        print *, "All coordinate system interface tests passed!"
    end if

    if (.not. all_passed) then
        error stop 1
    end if

contains

    subroutine test_coordinate_interface(nerrors)
        integer, intent(inout) :: nerrors
        class(coordinate_system_t), allocatable :: cs_vmec, cs_geo

        print *, "Testing factory functions..."

        call make_vmec_coordinate_system(cs_vmec)
        if (.not. allocated(cs_vmec)) then
            print *, "  FAIL: make_vmec_coordinate_system did not allocate cs"
            nerrors = nerrors + 1
        else
            print *, "  PASS: make_vmec_coordinate_system allocated"
        end if

        call make_geoflux_coordinate_system(cs_geo)
        if (.not. allocated(cs_geo)) then
            print *, "  FAIL: make_geoflux_coordinate_system did not allocate cs"
            nerrors = nerrors + 1
        else
            print *, "  PASS: make_geoflux_coordinate_system allocated"
        end if

        if (allocated(cs_vmec)) deallocate(cs_vmec)
        if (allocated(cs_geo)) deallocate(cs_geo)

        print *, ""
        print *, "Coordinate system interface verified successfully."
        print *, "Note: Full integration tests with real data are available through:"
        print *, "  - test_vmec_modules (VMEC coordinate transforms)"
        print *, "  - test_geoflux (Geoflux coordinate transforms)"

    end subroutine test_coordinate_interface

end program test_coordinate_systems
