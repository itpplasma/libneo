module bc_file
    !> Reader for Boozer-coordinate .bc ASCII files (Strumberger/NEO-2 format).
    !>
    !> Returns the raw file quantities without unit conversion: currents stay
    !> in the file's [A] columns and no derived major radius is computed.
    !> Consumers apply their own conventions (e.g. -mu0/(2 pi) factors).

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    private
    public :: bc_data_t, read_bc_file

    integer, parameter :: max_header_lines = 100

    type :: bc_data_t
        integer :: m0b, n0b, nsurf, nper, nmode
        real(dp) :: flux !< toroidal flux [Tm^2] from global header
        real(dp) :: a, R !< minor/major radius [m] from global header
        integer, allocatable :: m(:), n(:) !< mode numbers, n without nper
        real(dp), allocatable :: s(:) !< per-surface normalized toroidal flux
        real(dp), allocatable :: iota(:)
        real(dp), allocatable :: Jpol_over_nper(:), Itor(:) !< [A]
        real(dp), allocatable :: rmnc(:, :), zmns(:, :) !< (nsurf, nmode) [m]
        real(dp), allocatable :: vmns(:, :) !< (nsurf, nmode) [ ]
        real(dp), allocatable :: bmnc(:, :) !< (nsurf, nmode) [T]
    end type bc_data_t

contains

    subroutine read_bc_file(filename, d)
        character(len=*), intent(in) :: filename
        type(bc_data_t), intent(out) :: d

        integer :: iunit, ios, isurf, imode, m_int, n_int
        real(dp) :: surf_vals(6), row(4)
        character(len=512) :: line

        open (newunit=iunit, file=trim(filename), status='old', &
              action='read', iostat=ios)
        if (ios /= 0) then
            print *, "bc_file: cannot open file: ", trim(filename)
            error stop
        end if

        call next_matching_line(iunit, line, is_global_header, &
                                "global header (m0b n0b nsurf nper flux a R)")
        read (line, *) d%m0b, d%n0b, d%nsurf, d%nper, d%flux, d%a, d%R

        do isurf = 1, d%nsurf
            call next_matching_line(iunit, line, is_real_row_6, &
                                    "surface parameter line")
            read (line, *) surf_vals
            if (isurf == 1) call first_surface_pass(d, iunit, line)
            call store_surface_params(d, isurf, surf_vals)
            if (isurf == 1) cycle
            call next_matching_line(iunit, line, is_mode_row, "mode table row")
            do imode = 1, d%nmode
                if (imode > 1) then
                    read (iunit, '(A)', iostat=ios) line
                    if (ios /= 0) error stop "bc_file: truncated mode table"
                end if
                read (line, *, iostat=ios) m_int, n_int, row(1:4)
                if (ios /= 0) error stop "bc_file: malformed mode table row"
                if (d%m(imode) /= m_int .or. d%n(imode) /= n_int) then
                    error stop "bc_file: mode table differs between surfaces"
                end if
                d%rmnc(isurf, imode) = row(1)
                d%zmns(isurf, imode) = row(2)
                d%vmns(isurf, imode) = row(3)
                d%bmnc(isurf, imode) = row(4)
            end do
        end do
        close (iunit)
    end subroutine read_bc_file

    !> Read the first surface block to discover the mode table size, detect
    !! non-stellarator-symmetric files, and allocate all storage.
    subroutine first_surface_pass(d, iunit, line)
        type(bc_data_t), intent(inout) :: d
        integer, intent(in) :: iunit
        character(len=512), intent(inout) :: line

        integer :: ios, m_int, n_int, nmode
        integer, parameter :: max_modes = 100000
        real(dp) :: probe(10), row(4)
        real(dp), allocatable :: m_tmp(:), n_tmp(:), table(:, :)

        allocate (m_tmp(1024), n_tmp(1024), table(1024, 4))

        call next_matching_line(iunit, line, is_mode_row, "first mode table row")
        read (line, *, iostat=ios) probe
        if (ios == 0) then
            error stop "bc_file: non-stellarator-symmetric .bc files "// &
                "(10 columns) are not supported"
        end if

        nmode = 0
        do
            read (line, *, iostat=ios) m_int, n_int, row
            if (ios /= 0) exit
            nmode = nmode + 1
            if (nmode > max_modes) error stop "bc_file: mode table too large"
            if (nmode > size(m_tmp)) call grow(m_tmp, n_tmp, table)
            m_tmp(nmode) = real(m_int, dp)
            n_tmp(nmode) = real(n_int, dp)
            table(nmode, :) = row
            read (iunit, '(A)', iostat=ios) line
            if (ios /= 0) exit
        end do
        if (nmode == 0) error stop "bc_file: empty mode table"

        d%nmode = nmode
        allocate (d%m(nmode), d%n(nmode))
        allocate (d%s(d%nsurf), d%iota(d%nsurf))
        allocate (d%Jpol_over_nper(d%nsurf), d%Itor(d%nsurf))
        allocate (d%rmnc(d%nsurf, nmode), d%zmns(d%nsurf, nmode))
        allocate (d%vmns(d%nsurf, nmode), d%bmnc(d%nsurf, nmode))
        d%m = nint(m_tmp(1:nmode))
        d%n = nint(n_tmp(1:nmode))
        d%rmnc(1, :) = table(1:nmode, 1)
        d%zmns(1, :) = table(1:nmode, 2)
        d%vmns(1, :) = table(1:nmode, 3)
        d%bmnc(1, :) = table(1:nmode, 4)
    end subroutine first_surface_pass

    subroutine store_surface_params(d, isurf, surf_vals)
        type(bc_data_t), intent(inout) :: d
        integer, intent(in) :: isurf
        real(dp), intent(in) :: surf_vals(6)

        d%s(isurf) = surf_vals(1)
        d%iota(isurf) = surf_vals(2)
        d%Jpol_over_nper(isurf) = surf_vals(3)
        d%Itor(isurf) = surf_vals(4)
    end subroutine store_surface_params

    !> Advance to the next line for which matches() is true, skipping headers
    !! and comments. Aborts with what_expected in the message on EOF.
    subroutine next_matching_line(iunit, line, matches, what_expected)
        integer, intent(in) :: iunit
        character(len=512), intent(out) :: line
        interface
            logical function matches(line)
                character(len=512), intent(in) :: line
            end function matches
        end interface
        character(len=*), intent(in) :: what_expected

        integer :: ios, tries

        do tries = 1, max_header_lines
            read (iunit, '(A)', iostat=ios) line
            if (ios /= 0) then
                print *, "bc_file: unexpected end of file, expected ", &
                    what_expected
                error stop
            end if
            if (matches(line)) return
        end do
        print *, "bc_file: could not find ", what_expected
        error stop
    end subroutine next_matching_line

    logical function is_global_header(line)
        character(len=512), intent(in) :: line

        integer :: ios, ints(4)
        real(dp) :: reals(3)

        read (line, *, iostat=ios) ints, reals
        is_global_header = (ios == 0)
    end function is_global_header

    logical function is_real_row_6(line)
        character(len=512), intent(in) :: line

        integer :: ios
        real(dp) :: vals(6)

        read (line, *, iostat=ios) vals
        is_real_row_6 = (ios == 0)
    end function is_real_row_6

    logical function is_mode_row(line)
        character(len=512), intent(in) :: line

        integer :: ios, ints(2)
        real(dp) :: vals(4)

        read (line, *, iostat=ios) ints, vals
        is_mode_row = (ios == 0)
    end function is_mode_row

    subroutine grow(m_tmp, n_tmp, table)
        real(dp), allocatable, intent(inout) :: m_tmp(:), n_tmp(:), table(:, :)

        real(dp), allocatable :: m_new(:), n_new(:), table_new(:, :)
        integer :: n_old

        n_old = size(m_tmp)
        allocate (m_new(2*n_old), n_new(2*n_old), table_new(2*n_old, 4))
        m_new(1:n_old) = m_tmp
        n_new(1:n_old) = n_tmp
        table_new(1:n_old, :) = table
        call move_alloc(m_new, m_tmp)
        call move_alloc(n_new, n_tmp)
        call move_alloc(table_new, table)
    end subroutine grow

end module bc_file
