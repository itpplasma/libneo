!> Save and load a solved frozen-current field (patch, per-mode V1
!> coefficients and coils), so consumers need not repeat the solve.
module neo_frozen_current_io
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use hdf5_tools, only: hid_t, h5_init, h5_deinit, h5_create, h5_open, h5_close, &
                          h5_add, h5_get
    use fortfem_meridian_bspline_n, only: meridian_bspline_build
    use neo_biotsavart, only: coils_init
    use neo_frozen_current_field, only: frozen_current_field_t, frozen_current_patch_t
    implicit none
    private

    public :: save_frozen_current_field, load_frozen_current_field

contains

    subroutine save_frozen_current_field(filename, patch, field)
        character(*), intent(in) :: filename
        type(frozen_current_patch_t), intent(in) :: patch
        type(frozen_current_field_t), intent(in) :: field

        integer(hid_t) :: h5id
        real(dp), allocatable :: x(:), y(:), z(:), current(:)
        integer :: npoints

        call h5_init()
        call h5_create(filename, h5id)
        call h5_add(h5id, 'r0', patch%r0)
        call h5_add(h5id, 'r1', patch%r1)
        call h5_add(h5id, 'z0', patch%z0)
        call h5_add(h5id, 'z1', patch%z1)
        call h5_add(h5id, 'nr', patch%nr)
        call h5_add(h5id, 'nz', patch%nz)
        call h5_add(h5id, 'degree', patch%degree)
        call h5_add(h5id, 'kmax', field%plasma%kmax)
        call h5_add(h5id, 'nfp', field%plasma%nfp)
        call h5_add(h5id, 'u', field%plasma%u, lbound(field%plasma%u), &
                    ubound(field%plasma%u))
        npoints = 0
        if (allocated(field%coils%x)) npoints = size(field%coils%x)
        call h5_add(h5id, 'coil_points', npoints)
        if (npoints > 0) then
            allocate (x, source=field%coils%x)
            allocate (y, source=field%coils%y)
            allocate (z, source=field%coils%z)
            allocate (current, source=field%coils%current)
            call h5_add(h5id, 'coil_x', x)
            call h5_add(h5id, 'coil_y', y)
            call h5_add(h5id, 'coil_z', z)
            call h5_add(h5id, 'coil_current', current)
        end if
        call h5_close(h5id)
        call h5_deinit()
    end subroutine save_frozen_current_field

    subroutine load_frozen_current_field(filename, field, status)
        character(*), intent(in) :: filename
        type(frozen_current_field_t), intent(out) :: field
        integer, intent(out) :: status

        type(frozen_current_patch_t) :: patch
        integer(hid_t) :: h5id
        real(dp), allocatable :: x(:), y(:), z(:), current(:)
        integer :: npoints

        call h5_init()
        call h5_open(filename, h5id)
        call h5_get(h5id, 'r0', patch%r0)
        call h5_get(h5id, 'r1', patch%r1)
        call h5_get(h5id, 'z0', patch%z0)
        call h5_get(h5id, 'z1', patch%z1)
        call h5_get(h5id, 'nr', patch%nr)
        call h5_get(h5id, 'nz', patch%nz)
        call h5_get(h5id, 'degree', patch%degree)
        call h5_get(h5id, 'kmax', field%plasma%kmax)
        call h5_get(h5id, 'nfp', field%plasma%nfp)
        call meridian_bspline_build(patch%r0, patch%r1, patch%z0, patch%z1, &
                                    patch%nr, patch%nz, patch%degree, &
                                    field%plasma%sp, status)
        if (status == 0) then
            allocate (field%plasma%u(field%plasma%sp%n1, 0:field%plasma%kmax))
            call h5_get(h5id, 'u', field%plasma%u)
            call h5_get(h5id, 'coil_points', npoints)
            if (npoints > 0) then
                allocate (x(npoints), y(npoints), z(npoints), current(npoints))
                call h5_get(h5id, 'coil_x', x)
                call h5_get(h5id, 'coil_y', y)
                call h5_get(h5id, 'coil_z', z)
                call h5_get(h5id, 'coil_current', current)
                call coils_init(x, y, z, current, field%coils)
            end if
        end if
        call h5_close(h5id)
        call h5_deinit()
    end subroutine load_frozen_current_field

end module neo_frozen_current_io
