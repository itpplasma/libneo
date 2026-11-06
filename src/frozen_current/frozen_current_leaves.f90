!> Generated pointwise leaves of the frozen-current field.
module neo_frozen_current_leaves
    implicit none
    private
    public :: vmec_source_field, vmec_surface_normal, area_cross_field
    public :: cyl_to_cart, cart_to_cyl
contains
    include '../generated/frozen_current_source_leaves.inc'
end module neo_frozen_current_leaves
