!> Guiding-centre field quantities in a curvilinear chart from the Cartesian
!> field and its Jacobian, for tracers that tabulate them once per chart.
module neo_guiding_centre_field
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: guiding_centre_quantities

    !> Order of q in guiding_centre_quantities.
    integer, parameter, public :: GC_BMOD = 1, GC_SQRTG = 2, GC_BDER = 3, &
                                  GC_HCOV = 6, GC_HCTR = 9, GC_HCURL = 12, GC_COUNT = 14

contains

    !> e(i, j) = dx_i/du_j, b and db(i, j) = db_i/dx_j at the same point. Returns
    !> |B|, sqrt(g), d_j log|B|, h_j, h^j and (curl h)^j in q(1:14).
    pure subroutine guiding_centre_quantities(e, b, db, q)
        real(dp), intent(in) :: e(3, 3), b(3), db(3, 3)
        real(dp), intent(out) :: q(GC_COUNT)

        call guiding_centre_field(e(1, 1), e(1, 2), e(1, 3), e(2, 1), e(2, 2), &
                                  e(2, 3), e(3, 1), e(3, 2), e(3, 3), b(1), b(2), b(3), &
                                  db(1, 1), db(1, 2), db(1, 3), db(2, 1), db(2, 2), &
                                  db(2, 3), db(3, 1), db(3, 2), db(3, 3), q(1), q(2), &
                                  q(3), q(4), q(5), q(6), q(7), q(8), q(9), q(10), &
                                  q(11), q(12), q(13), q(14))
    end subroutine guiding_centre_quantities

    include '../generated/guiding_centre_field_leaves.inc'

end module neo_guiding_centre_field
