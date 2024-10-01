program test_poincare_plot
use, intrinsic :: iso_fortran_env, only: dp => real64
use util_for_test, only: print_test, print_ok, print_fail
implicit none

call test_poincare_plot_call

contains


subroutine test_poincare_plot_call
    use neo_poincare_plot, only: make_poincare_plot

    call print_test("test_poincare_plot_call")

    call print_ok
end subroutine test_poincare_plot_call

end program test_poincare_plot