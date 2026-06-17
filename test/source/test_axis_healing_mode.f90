program test_axis_healing_mode
    !> resolve_axis_healing_mode maps the axis_healing / axis_healing_boundary
    !> selectors and the deprecated boolean switches to a (mode, boundary) pair.
    !> Explicit selectors win and are case-insensitive; when axis_healing is empty
    !> the legacy booleans are mapped for back-compat.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use new_vmec_stuff_mod, only: axis_healing, axis_healing_boundary, &
                                  old_axis_healing, old_axis_healing_boundary, &
                                  axis_healing_power_law, axis_healing_polyfit
    use spline_vmec_sub, only: resolve_axis_healing_mode
    implicit none

    character(len=16) :: mode, boundary

    ! Defaults (everything unset) reproduce the legacy/fixed behavior.
    call reset()
    call resolve_axis_healing_mode(mode, boundary)
    call expect(mode, 'legacy', boundary, 'fixed', 'defaults')

    ! Deprecated booleans map to the right mode.
    call reset(); axis_healing_power_law = .True.
    call resolve_axis_healing_mode(mode, boundary)
    call expect(mode, 'powerlaw', boundary, 'fixed', 'legacy bool powerlaw')

    call reset(); axis_healing_polyfit = .True.
    call resolve_axis_healing_mode(mode, boundary)
    call expect(mode, 'polyfit', boundary, 'fixed', 'legacy bool polyfit')

    call reset(); old_axis_healing_boundary = .False.
    call resolve_axis_healing_mode(mode, boundary)
    call expect(mode, 'legacy', boundary, 'adaptive', 'legacy bool adaptive boundary')

    ! Explicit selector overrides the booleans and is case-insensitive.
    call reset(); axis_healing = 'polyfit'; axis_healing_power_law = .True.
    call resolve_axis_healing_mode(mode, boundary)
    call expect(mode, 'polyfit', boundary, 'fixed', 'explicit overrides bool')

    call reset(); axis_healing = 'POWERLAW'
    call resolve_axis_healing_mode(mode, boundary)
    call expect(mode, 'powerlaw', boundary, 'fixed', 'case-insensitive')

    call reset(); axis_healing = 'legacy'; axis_healing_boundary = 'Adaptive'
    call resolve_axis_healing_mode(mode, boundary)
    call expect(mode, 'legacy', boundary, 'adaptive', 'explicit boundary')

    print *, 'axis healing mode resolution: all checks passed'

contains

    subroutine reset()
        axis_healing = ''
        axis_healing_boundary = ''
        old_axis_healing = .True.
        old_axis_healing_boundary = .True.
        axis_healing_power_law = .False.
        axis_healing_polyfit = .False.
    end subroutine reset

    subroutine expect(mode, mode_exp, boundary, boundary_exp, label)
        character(len=*), intent(in) :: mode, mode_exp, boundary, boundary_exp, label
        if (trim(mode) /= mode_exp .or. trim(boundary) /= boundary_exp) then
            print *, 'FAIL [', label, ']: got mode=', trim(mode), ' boundary=', trim(boundary), &
                ' expected mode=', mode_exp, ' boundary=', boundary_exp
            error stop 1
        end if
    end subroutine expect

end program test_axis_healing_mode
