module libneo_lagrange_point
    ! Generated uniform-stencil Lagrange point kernels (gen/gen_lagrange_point).
    implicit none
    private
    public :: lagrange_point_3d_n4, lagrange_point_3d_n6, lagrange_value_3d_n6
contains
    include '../generated/lagrange_point_leaves.inc'
end module libneo_lagrange_point
