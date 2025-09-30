program test_vmec_modules
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use vmec_coordinates
    implicit none

    procedure(transform_i), pointer :: fptr
    real(dp) :: cyl(3)

    fptr => get_transform('cyl', 'cart')
    call cyl_to_cart([1.0_dp, 0.0_dp, 0.0_dp], cyl)
    if (abs(cyl(1) - 1.0_dp) > 1.0e-12_dp) then
        error stop 'cyl_to_cart sanity check failed'
    end if

    print *, 'vmec modules available'
end program test_vmec_modules
