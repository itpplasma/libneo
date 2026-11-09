module libneo_chartmap_field
    !> Magnetic field on the nodes of a chartmap file: physical cylindrical
    !> components of A (G cm) and B (G) as variables A_R, A_phi, A_Z, B_R,
    !> B_phi, B_Z on (rho, theta, zeta), marked by the global attribute
    !> magnetic_field = "cylindrical_A_B". Cylindrical components are periodic
    !> in theta and over one field period in zeta for any zeta convention, so
    !> a reader splines them on the chart grid itself and never resamples.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use netcdf, only: NF90_CLOBBER, NF90_NETCDF4, NF90_DOUBLE, NF90_INT, &
                      NF90_GLOBAL, NF90_NOERR, NF90_NOWRITE, nf90_close, nf90_create, &
                      nf90_def_dim, nf90_def_var, nf90_enddef, nf90_get_att, &
                      nf90_get_var, nf90_inq_dimid, nf90_inq_varid, &
                      nf90_inquire_dimension, nf90_open, &
                      nf90_put_att, nf90_put_var, nf90_strerror
    use neo_field_base, only: field_t
    use libneo_coordinates_base, only: coordinate_system_t
    use libneo_coordinates_chartmap, only: make_chartmap_coordinate_system
    use libneo_coordinates_validator, only: validate_chartmap_file
    use batch_interpolate, only: BatchSplineData3D, construct_batch_splines_3d, &
                                 evaluate_batch_splines_3d, destroy_batch_splines_3d
    implicit none
    private

    public :: chartmap_field_t, chartmap_has_field, read_chartmap_field, &
              sample_chartmap_field, write_chartmap_field

    character(len=*), parameter :: marker = 'cylindrical_A_B'
    character(len=5), parameter :: names(6) = &
                                   ['A_R  ', 'A_phi', 'A_Z  ', 'B_R  ', 'B_phi', 'B_Z  ']
    real(dp), parameter :: twopi = 6.283185307179586476925286766559_dp

    !> The chart grid of the file and the field at its nodes; a and b are
    !> (component R, phi, Z; rho; theta; zeta).
    type :: chartmap_field_t
        real(dp), allocatable :: rho(:), theta(:), zeta(:)
        integer :: num_field_periods = 1
        real(dp), allocatable :: a(:, :, :, :), b(:, :, :, :)
        !> Largest relative error of the order-5 spline of B at cell centres,
        !> and of A relative to the largest |A| (negative if not measured).
        real(dp) :: spline_error_b = -1.0_dp, spline_error_a = -1.0_dp
        character(len=1024) :: source = ''
    end type chartmap_field_t

contains

    !> True if filename is a valid chartmap carrying the field variables.
    logical function chartmap_has_field(filename) result(has)
        character(len=*), intent(in) :: filename
        integer :: ncid, i, varid, ierr
        character(len=256) :: message, value

        has = .false.
        call validate_chartmap_file(filename, ierr, message)
        if (ierr /= 0) return
        if (nf90_open(trim(filename), NF90_NOWRITE, ncid) /= NF90_NOERR) return
        value = ''
        if (nf90_get_att(ncid, NF90_GLOBAL, 'magnetic_field', value) == NF90_NOERR) then
            has = trim(value) == marker
            do i = 1, 6
                has = has .and. nf90_inq_varid(ncid, trim(names(i)), varid) == NF90_NOERR
            end do
        end if
        if (nf90_close(ncid) /= NF90_NOERR) has = .false.
    end function chartmap_has_field

    !> Field at the nodes of chartmap_file, from the Cartesian field. The
    !> spline error is measured at all cell centres against field itself.
    subroutine sample_chartmap_field(field, chartmap_file, source, f)
        class(field_t), intent(in) :: field
        character(len=*), intent(in) :: chartmap_file, source
        type(chartmap_field_t), intent(out) :: f

        real(dp), allocatable :: x(:, :, :), y(:, :, :), z(:, :, :)
        integer :: nr, nt, nz, i, j, k

        call read_grid(chartmap_file, f, x, y, z)
        nr = size(f%rho); nt = size(f%theta); nz = size(f%zeta)
        allocate (f%a(3, nr, nt, nz), f%b(3, nr, nt, nz))
        !$omp parallel do collapse(2) private(i)
        do k = 1, nz
            do j = 1, nt
                do i = 1, nr
                    call cylindrical_ab(field, [x(i, j, k), y(i, j, k), z(i, j, k)], &
                                        f%a(:, i, j, k), f%b(:, i, j, k))
                end do
            end do
        end do
        !$omp end parallel do
        f%source = source
        call measure_spline_error(field, chartmap_file, f)
    end subroutine sample_chartmap_field

    subroutine cylindrical_ab(field, xc, a, b)
        class(field_t), intent(in) :: field
        real(dp), intent(in) :: xc(3)
        real(dp), intent(out) :: a(3), b(3)
        real(dp) :: ac(3), bc(3), c, s, r

        call field%compute_abfield(xc, ac, bc)
        r = hypot(xc(1), xc(2))
        c = xc(1)/r; s = xc(2)/r
        a = [c*ac(1) + s*ac(2), -s*ac(1) + c*ac(2), ac(3)]
        b = [c*bc(1) + s*bc(2), -s*bc(1) + c*bc(2), bc(3)]
    end subroutine cylindrical_ab

    subroutine measure_spline_error(field, chartmap_file, f)
        class(field_t), intent(in) :: field
        character(len=*), intent(in) :: chartmap_file
        type(chartmap_field_t), intent(inout) :: f

        class(coordinate_system_t), allocatable :: chart
        type(BatchSplineData3D) :: spl
        real(dp) :: u(3), xc(3), a(3), b(3), y(6), amax, eb, ea, period
        integer :: nr, nt, nz, i, j, k

        nr = size(f%rho); nt = size(f%theta); nz = size(f%zeta)
        period = twopi/f%num_field_periods
        call field_spline(f, spl)
        call make_chartmap_coordinate_system(chart, chartmap_file)
        amax = maxval(norm2(f%a, 1))
        eb = 0; ea = 0
        !$omp parallel do collapse(2) private(i, u, xc, a, b, y) &
        !$omp reduction(max: eb, ea)
        do k = 1, nz
            do j = 1, nt
                do i = 1, nr - 1
                    u = [0.5_dp*(f%rho(i) + f%rho(i + 1)), &
                         f%theta(j) + 0.5_dp*twopi/nt, f%zeta(k) + 0.5_dp*period/nz]
                    call chart%evaluate_cart(u, xc)
                    call cylindrical_ab(field, xc, a, b)
                    call evaluate_batch_splines_3d(spl, u, y)
                    eb = max(eb, norm2(y(4:6) - b)/norm2(b))
                    ea = max(ea, norm2(y(1:3) - a)/amax)
                end do
            end do
        end do
        !$omp end parallel do
        f%spline_error_b = eb
        f%spline_error_a = ea
        call destroy_batch_splines_3d(spl)
    end subroutine measure_spline_error

    !> Order-5 spline of the nodal A and B on the chart grid, with the
    !> periodic theta and zeta endpoints appended.
    subroutine field_spline(f, spl)
        type(chartmap_field_t), intent(in) :: f
        type(BatchSplineData3D), intent(out) :: spl
        real(dp), allocatable :: yb(:, :, :, :)
        integer :: nr, nt, nz, q
        real(dp) :: drho

        nr = size(f%rho); nt = size(f%theta); nz = size(f%zeta)
        drho = (f%rho(nr) - f%rho(1))/(nr - 1)
        if (maxval(abs(f%rho(2:) - f%rho(:nr - 1) - drho)) > 1.0e-10_dp*drho) &
            error stop 'chartmap_field: rho grid must be uniform for the spline'
        allocate (yb(nr, nt + 1, nz + 1, 6))
        do q = 1, 3
            yb(:, :nt, :nz, q) = f%a(q, :, :, :)
            yb(:, :nt, :nz, q + 3) = f%b(q, :, :, :)
        end do
        yb(:, nt + 1, :nz, :) = yb(:, 1, :nz, :)
        yb(:, :, nz + 1, :) = yb(:, :, 1, :)
        call construct_batch_splines_3d([f%rho(1), f%theta(1), f%zeta(1)], &
                                        [f%rho(nr), f%theta(1) + twopi, &
                                         f%zeta(1) + twopi/f%num_field_periods], &
                                        yb, [5, 5, 5], [.false., .true., .true.], spl)
    end subroutine field_spline

    subroutine read_grid(filename, f, x, y, z)
        character(len=*), intent(in) :: filename
        type(chartmap_field_t), intent(inout) :: f
        real(dp), allocatable, intent(out) :: x(:, :, :), y(:, :, :), z(:, :, :)
        integer :: ncid, ierr, varid, nr, nt, nz
        character(len=2048) :: message

        call validate_chartmap_file(filename, ierr, message)
        if (ierr /= 0) then
            print *, 'chartmap_field: invalid chartmap ', trim(filename), ': ', &
                trim(message)
            error stop
        end if
        call check(nf90_open(trim(filename), NF90_NOWRITE, ncid), filename)
        nr = dim_len(ncid, 'rho'); nt = dim_len(ncid, 'theta'); nz = dim_len(ncid, 'zeta')
        allocate (f%rho(nr), f%theta(nt), f%zeta(nz))
        allocate (x(nr, nt, nz), y(nr, nt, nz), z(nr, nt, nz))
        call get1(ncid, 'rho', f%rho); call get1(ncid, 'theta', f%theta)
        call get1(ncid, 'zeta', f%zeta)
        call get3(ncid, 'x', x); call get3(ncid, 'y', y); call get3(ncid, 'z', z)
        f%num_field_periods = 1
        if (nf90_inq_varid(ncid, 'num_field_periods', varid) == NF90_NOERR) &
            call check(nf90_get_var(ncid, varid, f%num_field_periods), filename)
        call check(nf90_close(ncid), filename)
    end subroutine read_grid

    !> Read the field of a chartmap file written by write_chartmap_field.
    subroutine read_chartmap_field(filename, f)
        character(len=*), intent(in) :: filename
        type(chartmap_field_t), intent(out) :: f
        real(dp), allocatable :: x(:, :, :), y(:, :, :), z(:, :, :), v(:, :, :)
        integer :: ncid, q

        if (.not. chartmap_has_field(filename)) then
            print *, 'read_chartmap_field: ', trim(filename), &
                ' is not a chartmap with magnetic_field = "', marker, '"'
            error stop
        end if
        call read_grid(filename, f, x, y, z)
        allocate (v, mold=x)
        allocate (f%a(3, size(x, 1), size(x, 2), size(x, 3)), mold=0.0_dp)
        allocate (f%b, mold=f%a)
        call check(nf90_open(trim(filename), NF90_NOWRITE, ncid), filename)
        do q = 1, 3
            call get3(ncid, trim(names(q)), v); f%a(q, :, :, :) = v
            call get3(ncid, trim(names(q + 3)), v); f%b(q, :, :, :) = v
        end do
        if (nf90_get_att(ncid, NF90_GLOBAL, 'field_spline_error_B', f%spline_error_b) &
            /= NF90_NOERR) f%spline_error_b = -1
        if (nf90_get_att(ncid, NF90_GLOBAL, 'field_spline_error_A', f%spline_error_a) &
            /= NF90_NOERR) f%spline_error_a = -1
        if (nf90_get_att(ncid, NF90_GLOBAL, 'field_source', f%source) /= NF90_NOERR) &
            f%source = ''
        call check(nf90_close(ncid), filename)
    end subroutine read_chartmap_field

    !> Write the chart of chartmap_file together with the field f to output.
    subroutine write_chartmap_field(chartmap_file, f, output)
        character(len=*), intent(in) :: chartmap_file, output
        type(chartmap_field_t), intent(in) :: f

        type(chartmap_field_t) :: g
        real(dp), allocatable :: x(:, :, :), y(:, :, :), z(:, :, :)
        integer :: ncid, dims(3), var(3), vxyz(3), vf(6), vnfp, q, src
        character(len=16) :: zconv, rconv

        call read_grid(chartmap_file, g, x, y, z)
        if (any(shape(x) /= shape(f%a(1, :, :, :)))) &
            error stop 'write_chartmap_field: field and chart grids differ'
        call check(nf90_open(trim(chartmap_file), NF90_NOWRITE, src), chartmap_file)
        zconv = ''; rconv = ''
        if (nf90_get_att(src, NF90_GLOBAL, 'zeta_convention', zconv) /= NF90_NOERR) &
            zconv = ''
        if (nf90_get_att(src, NF90_GLOBAL, 'rho_convention', rconv) /= NF90_NOERR) &
            rconv = ''
        call check(nf90_close(src), chartmap_file)

        call check(nf90_create(trim(output), ior(NF90_CLOBBER, NF90_NETCDF4), ncid), &
                   output)
        if (len_trim(zconv) > 0) &
            call check(nf90_put_att(ncid, NF90_GLOBAL, 'zeta_convention', trim(zconv)), &
                       output)
        if (len_trim(rconv) > 0) &
            call check(nf90_put_att(ncid, NF90_GLOBAL, 'rho_convention', trim(rconv)), &
                       output)
        call check(nf90_put_att(ncid, NF90_GLOBAL, 'magnetic_field', marker), output)
        call check(nf90_put_att(ncid, NF90_GLOBAL, 'field_source', &
                                trim(f%source)), output)
        call check(nf90_put_att(ncid, NF90_GLOBAL, 'field_spline_error_B', &
                                f%spline_error_b), output)
        call check(nf90_put_att(ncid, NF90_GLOBAL, 'field_spline_error_A', &
                                f%spline_error_a), output)
        call check(nf90_def_dim(ncid, 'rho', size(g%rho), dims(1)), output)
        call check(nf90_def_dim(ncid, 'theta', size(g%theta), dims(2)), output)
        call check(nf90_def_dim(ncid, 'zeta', size(g%zeta), dims(3)), output)
        call check(nf90_def_var(ncid, 'rho', NF90_DOUBLE, dims(1:1), var(1)), output)
        call check(nf90_def_var(ncid, 'theta', NF90_DOUBLE, dims(2:2), var(2)), output)
        call check(nf90_def_var(ncid, 'zeta', NF90_DOUBLE, dims(3:3), var(3)), output)
        call check(nf90_def_var(ncid, 'x', NF90_DOUBLE, dims, vxyz(1)), output)
        call check(nf90_def_var(ncid, 'y', NF90_DOUBLE, dims, vxyz(2)), output)
        call check(nf90_def_var(ncid, 'z', NF90_DOUBLE, dims, vxyz(3)), output)
        do q = 1, 3
            call check(nf90_put_att(ncid, vxyz(q), 'units', 'cm'), output)
        end do
        call check(nf90_def_var(ncid, 'num_field_periods', NF90_INT, vnfp), output)
        do q = 1, 6
            call check(nf90_def_var(ncid, trim(names(q)), NF90_DOUBLE, dims, vf(q)), &
                       output)
            call check(nf90_put_att(ncid, vf(q), 'units', &
                                    merge('G cm', 'G   ', q <= 3)), output)
        end do
        call check(nf90_enddef(ncid), output)
        call check(nf90_put_var(ncid, var(1), g%rho), output)
        call check(nf90_put_var(ncid, var(2), g%theta), output)
        call check(nf90_put_var(ncid, var(3), g%zeta), output)
        call check(nf90_put_var(ncid, vxyz(1), x), output)
        call check(nf90_put_var(ncid, vxyz(2), y), output)
        call check(nf90_put_var(ncid, vxyz(3), z), output)
        call check(nf90_put_var(ncid, vnfp, g%num_field_periods), output)
        do q = 1, 3
            call check(nf90_put_var(ncid, vf(q), f%a(q, :, :, :)), output)
            call check(nf90_put_var(ncid, vf(q + 3), f%b(q, :, :, :)), output)
        end do
        call check(nf90_close(ncid), output)
    end subroutine write_chartmap_field

    integer function dim_len(ncid, name) result(n)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: name
        integer :: dimid
        call check(nf90_inq_dimid(ncid, name, dimid), name)
        call check(nf90_inquire_dimension(ncid, dimid, len=n), name)
    end function dim_len

    subroutine get1(ncid, name, v)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: name
        real(dp), intent(out) :: v(:)
        integer :: varid
        call check(nf90_inq_varid(ncid, name, varid), name)
        call check(nf90_get_var(ncid, varid, v), name)
    end subroutine get1

    subroutine get3(ncid, name, v)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: name
        real(dp), intent(out) :: v(:, :, :)
        integer :: varid
        call check(nf90_inq_varid(ncid, name, varid), name)
        call check(nf90_get_var(ncid, varid, v), name)
    end subroutine get3

    subroutine check(status, what)
        integer, intent(in) :: status
        character(len=*), intent(in) :: what
        if (status /= NF90_NOERR) then
            print *, 'chartmap_field: ', trim(what), ': ', trim(nf90_strerror(status))
            error stop
        end if
    end subroutine check

end module libneo_chartmap_field
