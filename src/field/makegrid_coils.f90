!> STELLOPT/MAKEGRID coils file reader with coil groups (circuits).
!>
!> Parsing follows python/libneo/coils.py (CoilsFile.from_file with
!> split_on_zero_current): '#' and '!' start comments, 'periods' and 'mirror'
!> are recorded but not applied, 'begin filament' opens and 'end filament' or
!> 'end' closes a block, points are 'x y z I', and every point with I exactly 0
!> terminates a coil without belonging to it. Coils with fewer than two points
!> are dropped. The terminating line 'x y z 0 group name' assigns the group id
!> (circuit index, MAKEGRID EXTCUR index) and its optional name. A coil without
!> such a terminator gets group 0. The fourth column is kept per point, so a
!> group field built from it is the field per unit EXTCUR (for W7-X style files
!> the fourth column is the winding count).
module neo_makegrid_coils
    use neo_biotsavart, only: dp, coils_t, coils_init

    implicit none

    private
    public :: makegrid_coils_t, read_makegrid_coils, makegrid_to_biotsavart

    integer, parameter :: name_len = 64
    integer, parameter :: line_len = 4096

    type :: makegrid_coils_t
        integer :: periods = 0
        character(len=name_len) :: mirror = ""
        integer :: n_coils = 0
        integer :: n_groups = 0
        integer, allocatable :: first(:), last(:), group(:)
        real(dp), allocatable :: x(:), y(:), z(:), current(:)
        character(len=name_len), allocatable :: group_name(:)
    end type makegrid_coils_t

    type :: point_buffer_t
        integer :: n = 0
        real(dp), allocatable :: v(:, :)
        integer, allocatable :: group(:)
        character(len=name_len), allocatable :: name(:)
    end type point_buffer_t

contains

    !> Reads a MAKEGRID coils file. Without ierr, a malformed file stops with
    !> an error message naming file and line.
    subroutine read_makegrid_coils(filename, mc, ierr)
        character(*), intent(in) :: filename
        type(makegrid_coils_t), intent(out) :: mc
        integer, intent(out), optional :: ierr

        type(point_buffer_t) :: block, coils_pts
        integer, allocatable :: first(:), last(:), group(:)
        character(len=name_len), allocatable :: names(:)
        character(len=line_len) :: raw, line, low
        integer :: unit, ios, lineno, stat, n_coils
        logical :: in_filament
        character(len=256) :: msg

        stat = 0
        msg = ""
        lineno = 0
        n_coils = 0
        in_filament = .false.
        allocate(first(0), last(0), group(0), names(0))
        coils_pts%n = 0
        allocate(coils_pts%v(4, 0), coils_pts%group(0), coils_pts%name(0))
        call reset_buffer(block)

        open(newunit=unit, file=filename, status="old", action="read", iostat=ios)
        if (ios /= 0) then
            call fail(ierr, "cannot open MAKEGRID coils file "//trim(filename))
            return
        end if
        do
            read(unit, '(a)', iostat=ios) raw
            if (ios /= 0) exit
            lineno = lineno + 1
            line = strip_comment(raw)
            if (len_trim(line) == 0) cycle
            low = to_lower(line)
            if (starts_with(low, "periods")) then
                call parse_periods(line, mc%periods)
            else if (starts_with(low, "mirror")) then
                mc%mirror = second_token(line)
            else if (starts_with(low, "begin filament")) then
                if (in_filament) then
                    if (block%n == 0) then
                        stat = 1
                        exit
                    end if
                    call split_block(block, coils_pts, first, last, group, names)
                end if
                in_filament = .true.
            else if (starts_with(low, "end filament") .or. low == "end") then
                if (.not. in_filament) cycle
                if (block%n == 0) then
                    stat = 1
                    exit
                end if
                call split_block(block, coils_pts, first, last, group, names)
                in_filament = .false.
            else if (in_filament) then
                call append_point(line, block, stat)
                if (stat /= 0) exit
            end if
        end do
        close(unit)

        if (stat == 0 .and. in_filament) then
            if (block%n == 0) then
                stat = 1
            else
                call split_block(block, coils_pts, first, last, group, names)
            end if
        end if
        if (stat /= 0) then
            write(msg, '(a,a,a,i0,a)') "malformed MAKEGRID coils file ", &
                trim(filename), " at line ", lineno, &
                ": empty filament or fewer than 4 numeric fields (x y z I)"
            call fail(ierr, trim(msg))
            return
        end if

        call finalize(mc, coils_pts, first, last, group, names)
        if (present(ierr)) ierr = 0
    end subroutine read_makegrid_coils

    !> Builds a Biot-Savart polygon from all coils or from one group. Each coil
    !> is closed through its first point and ends with a zero-current point, as
    !> CoilsFile.write_simple does, so segment i carries the current of point i.
    subroutine makegrid_to_biotsavart(mc, coils, group)
        type(makegrid_coils_t), intent(in) :: mc
        type(coils_t), intent(out) :: coils
        integer, intent(in), optional :: group

        real(dp), allocatable :: x(:), y(:), z(:), cur(:)
        integer :: ic, n, np, i0, i1, ntot

        ntot = 0
        do ic = 1, mc%n_coils
            if (selected(mc, ic, group)) ntot = ntot + closed_count(mc, ic)
        end do
        allocate(x(ntot), y(ntot), z(ntot), cur(ntot))
        n = 0
        do ic = 1, mc%n_coils
            if (.not. selected(mc, ic, group)) cycle
            np = closed_count(mc, ic)
            i0 = mc%first(ic)
            i1 = i0 + np - 2
            x(n + 1:n + np - 1) = mc%x(i0:i1)
            y(n + 1:n + np - 1) = mc%y(i0:i1)
            z(n + 1:n + np - 1) = mc%z(i0:i1)
            cur(n + 1:n + np - 1) = mc%current(i0:i1)
            x(n + np) = mc%x(i0)
            y(n + np) = mc%y(i0)
            z(n + np) = mc%z(i0)
            cur(n + np) = 0.0_dp
            n = n + np
        end do
        call coils_init(x, y, z, cur, coils)
    end subroutine makegrid_to_biotsavart

    pure logical function selected(mc, ic, group)
        type(makegrid_coils_t), intent(in) :: mc
        integer, intent(in) :: ic
        integer, intent(in), optional :: group

        selected = .true.
        if (present(group)) selected = mc%group(ic) == group
    end function selected

    pure integer function closed_count(mc, ic)
        type(makegrid_coils_t), intent(in) :: mc
        integer, intent(in) :: ic

        integer :: i0, i1

        i0 = mc%first(ic)
        i1 = mc%last(ic)
        closed_count = i1 - i0 + 2
        if (i1 == i0) return
        if (mc%x(i0) == mc%x(i1) .and. mc%y(i0) == mc%y(i1) .and. &
            mc%z(i0) == mc%z(i1)) closed_count = closed_count - 1
    end function closed_count

    subroutine reset_buffer(buf)
        type(point_buffer_t), intent(inout) :: buf

        buf%n = 0
        if (allocated(buf%v)) deallocate(buf%v, buf%group, buf%name)
        allocate(buf%v(4, 64), buf%group(64), buf%name(64))
    end subroutine reset_buffer

    subroutine push(buf, v, group, name)
        type(point_buffer_t), intent(inout) :: buf
        real(dp), intent(in) :: v(4)
        integer, intent(in) :: group
        character(*), intent(in) :: name

        real(dp), allocatable :: v_new(:, :)
        integer, allocatable :: g_new(:)
        character(len=name_len), allocatable :: n_new(:)
        integer :: cap

        cap = size(buf%group)
        if (buf%n == cap) then
            cap = max(64, 2*cap)
            allocate(v_new(4, cap), g_new(cap), n_new(cap))
            v_new(:, 1:buf%n) = buf%v(:, 1:buf%n)
            g_new(1:buf%n) = buf%group(1:buf%n)
            n_new(1:buf%n) = buf%name(1:buf%n)
            call move_alloc(v_new, buf%v)
            call move_alloc(g_new, buf%group)
            call move_alloc(n_new, buf%name)
        end if
        buf%n = buf%n + 1
        buf%v(:, buf%n) = v
        buf%group(buf%n) = group
        buf%name(buf%n) = name
    end subroutine push

    subroutine append_point(line, block, stat)
        character(*), intent(in) :: line
        type(point_buffer_t), intent(inout) :: block
        integer, intent(out) :: stat

        character(len=line_len) :: tok(6)
        integer :: ntok, k, group, ios
        real(dp) :: v(4)

        call split_tokens(line, tok, ntok)
        stat = 1
        if (ntok < 4) return
        do k = 1, 4
            read(tok(k), *, iostat=ios) v(k)
            if (ios /= 0) return
        end do
        group = 0
        if (ntok >= 5) then
            read(tok(5), *, iostat=ios) group
            if (ios /= 0) group = 0
        end if
        if (ntok >= 6) then
            call push(block, v, group, trim(tok(6)))
        else
            call push(block, v, group, "")
        end if
        stat = 0
    end subroutine append_point

    subroutine split_block(block, pts, first, last, group, names)
        type(point_buffer_t), intent(inout) :: block, pts
        integer, allocatable, intent(inout) :: first(:), last(:), group(:)
        character(len=name_len), allocatable, intent(inout) :: names(:)

        integer :: i, start

        if (all(block%v(4, 1:block%n) /= 0.0_dp)) then
            call add_coil(block, 1, block%n, 0, "", 1, pts, first, last, group, &
                names)
            call reset_buffer(block)
            return
        end if
        start = 1
        do i = 1, block%n
            if (block%v(4, i) /= 0.0_dp) cycle
            call add_coil(block, start, i - 1, block%group(i), block%name(i), 2, &
                pts, first, last, group, names)
            start = i + 1
        end do
        if (start <= block%n) call add_coil(block, start, block%n, 0, "", 2, pts, &
            first, last, group, names)
        call reset_buffer(block)
    end subroutine split_block

    subroutine add_coil(block, i0, i1, g, name, min_points, pts, first, last, &
            group, names)
        type(point_buffer_t), intent(in) :: block
        integer, intent(in) :: i0, i1, g, min_points
        character(*), intent(in) :: name
        type(point_buffer_t), intent(inout) :: pts
        integer, allocatable, intent(inout) :: first(:), last(:), group(:)
        character(len=name_len), allocatable, intent(inout) :: names(:)

        integer :: i
        character(len=name_len) :: nm

        if (i1 - i0 + 1 < min_points) return
        nm = name
        call grow_int(first, pts%n + 1)
        do i = i0, i1
            call push(pts, block%v(:, i), g, "")
        end do
        call grow_int(last, pts%n)
        call grow_int(group, g)
        call grow_name(names, nm)
    end subroutine add_coil

    subroutine grow_int(a, v)
        integer, allocatable, intent(inout) :: a(:)
        integer, intent(in) :: v

        integer, allocatable :: b(:)

        allocate(b(size(a) + 1))
        b(1:size(a)) = a
        b(size(b)) = v
        call move_alloc(b, a)
    end subroutine grow_int

    subroutine grow_name(a, v)
        character(len=name_len), allocatable, intent(inout) :: a(:)
        character(len=name_len), intent(in) :: v

        character(len=name_len), allocatable :: b(:)

        allocate(b(size(a) + 1))
        b(1:size(a)) = a
        b(size(b)) = v
        call move_alloc(b, a)
    end subroutine grow_name

    subroutine finalize(mc, pts, first, last, group, names)
        type(makegrid_coils_t), intent(inout) :: mc
        type(point_buffer_t), intent(in) :: pts
        integer, intent(in) :: first(:), last(:), group(:)
        character(len=name_len), intent(in) :: names(:)

        integer :: ic, g

        mc%n_coils = size(first)
        allocate(mc%first(mc%n_coils), mc%last(mc%n_coils), mc%group(mc%n_coils))
        allocate(mc%x(pts%n), mc%y(pts%n), mc%z(pts%n), mc%current(pts%n))
        mc%first = first
        mc%last = last
        mc%group = group
        mc%x = pts%v(1, 1:pts%n)
        mc%y = pts%v(2, 1:pts%n)
        mc%z = pts%v(3, 1:pts%n)
        mc%current = pts%v(4, 1:pts%n)
        mc%n_groups = 0
        if (mc%n_coils > 0) mc%n_groups = max(0, maxval(group))
        allocate(mc%group_name(mc%n_groups))
        mc%group_name = ""
        do ic = 1, mc%n_coils
            g = group(ic)
            if (g < 1) cycle
            if (len_trim(mc%group_name(g)) == 0) mc%group_name(g) = names(ic)
        end do
    end subroutine finalize

    subroutine parse_periods(line, periods)
        character(*), intent(in) :: line
        integer, intent(inout) :: periods

        character(len=line_len) :: tok(6)
        integer :: ntok, value, ios

        call split_tokens(line, tok, ntok)
        if (ntok < 2) return
        if (verify(trim(tok(2)), "+-0123456789") /= 0) return
        read(tok(2), *, iostat=ios) value
        if (ios == 0) periods = value
    end subroutine parse_periods

    function second_token(line) result(tok2)
        character(*), intent(in) :: line
        character(len=name_len) :: tok2

        character(len=line_len) :: tok(6)
        integer :: ntok

        call split_tokens(line, tok, ntok)
        tok2 = ""
        if (ntok >= 2) tok2 = tok(2)
    end function second_token

    pure subroutine split_tokens(line, tok, ntok)
        character(*), intent(in) :: line
        character(len=*), intent(out) :: tok(:)
        integer, intent(out) :: ntok

        integer :: i, j, n

        tok = ""
        ntok = 0
        n = len_trim(line)
        i = 1
        do while (i <= n)
            if (is_blank(line(i:i))) then
                i = i + 1
                cycle
            end if
            j = i
            do while (j < n)
                if (is_blank(line(j + 1:j + 1))) exit
                j = j + 1
            end do
            ntok = ntok + 1
            if (ntok <= size(tok)) tok(ntok) = line(i:j)
            i = j + 1
        end do
    end subroutine split_tokens

    pure logical function is_blank(c)
        character, intent(in) :: c

        is_blank = c == " " .or. c == achar(9) .or. c == achar(13)
    end function is_blank

    pure function strip_comment(raw) result(line)
        character(*), intent(in) :: raw
        character(len=line_len) :: line

        integer :: ic

        ic = scan(raw, "#!")
        if (ic > 0) then
            line = trim_blanks(raw(1:ic - 1))
        else
            line = trim_blanks(raw)
        end if
    end function strip_comment

    pure function trim_blanks(s) result(t)
        character(*), intent(in) :: s
        character(len=line_len) :: t

        integer :: i0, i1

        i0 = 1
        i1 = len(s)
        do while (i0 <= i1)
            if (.not. is_blank(s(i0:i0))) exit
            i0 = i0 + 1
        end do
        do while (i1 >= i0)
            if (.not. is_blank(s(i1:i1))) exit
            i1 = i1 - 1
        end do
        t = s(i0:i1)
    end function trim_blanks

    pure function to_lower(s) result(t)
        character(*), intent(in) :: s
        character(len=len(s)) :: t

        integer :: i, c

        t = s
        do i = 1, len(s)
            c = iachar(s(i:i))
            if (c >= iachar("A") .and. c <= iachar("Z")) t(i:i) = achar(c + 32)
        end do
    end function to_lower

    pure logical function starts_with(s, prefix)
        character(*), intent(in) :: s, prefix

        starts_with = .false.
        if (len(s) < len(prefix)) return
        starts_with = s(1:len(prefix)) == prefix
    end function starts_with

    subroutine fail(ierr, msg)
        integer, intent(out), optional :: ierr
        character(*), intent(in) :: msg

        if (present(ierr)) then
            ierr = 1
            write(*, '(a)') trim(msg)
            return
        end if
        error stop msg
    end subroutine fail

end module neo_makegrid_coils
