subroutine spl_per(ns, n, h, splcoe)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use spl_three_to_five_sub, only: spl_per_mod => spl_per
    implicit none
    integer, intent(in) :: ns, n
    real(dp), intent(in) :: h
    real(dp), intent(inout) :: splcoe(0:ns, n)

    call spl_per_mod(ns, n, h, splcoe)
end subroutine spl_per

subroutine spl_reg(ns, n, h, splcoe)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use spl_three_to_five_sub, only: spl_reg_mod => spl_reg
    implicit none
    integer, intent(in) :: ns, n
    real(dp), intent(in) :: h
    real(dp), intent(inout) :: splcoe(0:ns, n)

    call spl_reg_mod(ns, n, h, splcoe)
end subroutine spl_reg
