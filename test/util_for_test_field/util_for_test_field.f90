module util_for_test_field
use, intrinsic :: iso_fortran_env, only: dp => real64

implicit none


contains


function compute_cartesian_curla(field, x, tol) result(curla)
    use neo_field_base, only: field_t

    class(field_t), intent(in) :: field
    real(dp), intent(in) :: x(3)
    real(dp), intent(in) :: tol
    real(dp) :: curla(3)

    real(dp) :: x_temp(3), A_temp1(3), A_temp2(3)
    real(dp) :: dA_dx(3,3), delta
    integer :: i

    delta = sqrt(tol)
    x_temp = x
    dA_dx = 0.0_dp
    do i = 1, 3
        x_temp = x
        x_temp(i) = x(i) + delta
        call field%compute_afield(x_temp, A_temp1)
        x_temp(i) = x(i) - delta
        call field%compute_afield(x_temp, A_temp2)
        dA_dx(i,:) = (A_temp1 - A_temp2) / (2.0_dp * delta)
    end do

    curla(1) = dA_dx(2,3) - dA_dx(3,2)
    curla(2) = dA_dx(3,1) - dA_dx(1,3)
    curla(3) = dA_dx(1,2) - dA_dx(2,1)
end function compute_cartesian_curla


function compute_cylindrical_curla(field, x, tol) result(curla)

    use neo_field_base, only: field_t

    class(field_t), intent(in) :: field
    real(dp), intent(in) :: x(3)
    real(dp), intent(in) :: tol
    real(dp) :: curla(3)

    real(dp) :: x_temp(3), A_temp1(3), A_temp2(3)
    real(dp) :: dA_dx(3,3), A(3), delta
    integer :: i

    delta = sqrt(tol)
    x_temp = x
    dA_dx = 0.0_dp
    do i = 1, 3
        x_temp = x
        x_temp(i) = x(i) + delta
        call field%compute_afield(x_temp, A_temp1)
        x_temp(i) = x(i) - delta
        call field%compute_afield(x_temp, A_temp2)
        dA_dx(i,:) = (A_temp1 - A_temp2) / (2.0_dp * delta)
    end do

    curla(1) = dA_dx(2,3)/x(1) - dA_dx(3,2)             ! (dA_Z/dphi)/R - dA_phi/dZ
    curla(2) = dA_dx(3,1) - dA_dx(1,3)                  ! dA_Z/dR - dA_R/dZ
    call field%compute_afield(x, A)
    curla(3) = dA_dx(1,2) + A(2)/x(1) - dA_dx(2,1)/x(1) ! (d(R*A_phi)/dR - dA_R/dphi)/R
end function compute_cylindrical_curla


function compute_cartesian_divb(field, x, tol) result(divb)
    use neo_field_base, only: field_t

    class(field_t), intent(in) :: field
    real(dp), intent(in) :: x(3)
    real(dp), intent(in) :: tol
    real(dp) :: divb

    real(dp) :: x_temp(3), B_temp1(3), B_temp2(3)
    real(dp) :: dB_dx(3), delta
    integer :: i

    delta = sqrt(tol)
    x_temp = x
    dB_dx = 0.0_dp
    do i = 1, 3
        x_temp = x
        x_temp(i) = x(i) + delta
        call field%compute_bfield(x_temp, B_temp1)
        x_temp(i) = x(i) - delta
        call field%compute_bfield(x_temp, B_temp2)
        dB_dx(i) = (B_temp1(i) - B_temp2(i)) / (2.0_dp * delta)
    end do

    divb = sum(dB_dx)
end function compute_cartesian_divb


end module util_for_test_field
