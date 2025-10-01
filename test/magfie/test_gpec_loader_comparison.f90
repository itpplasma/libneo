program test_gpec_loader_comparison
    use iso_fortran_env, only: dp => real64, error_unit
    use coil_tools, only: coil_t, coils_read_GPEC, coil_deinit, coils_append
    use neo_biotsavart, only: coils_t, load_coils_from_gpec_file, coils_deinit

    implicit none

    character(len=*), parameter :: coil_files(2) = [character(len=32) :: 'aug_bu.dat', 'aug_bl.dat']
    real(dp), parameter :: tol = 1.0e-12_dp
    character(len=512) :: data_dir
    integer :: argc

    argc = command_argument_count()
    if (argc >= 1) then
        call get_command_argument(1, data_dir)
    else
        data_dir = 'test/magfie/test_data'
    end if
    call ensure_trailing_slash(data_dir)

    call test_each_file(data_dir, coil_files, tol)
    call test_combined_files(data_dir, coil_files, tol)

    write(*, '(A)') 'GPEC loader comparison passed.'

contains

    subroutine test_each_file(base, files, tolerance)
        character(len=*), intent(in) :: base
        character(len=*), intent(in) :: files(:)
        real(dp), intent(in) :: tolerance
        type(coil_t), allocatable :: coils_ct(:)
        type(coils_t) :: coils_neo
        integer :: i

        do i = 1, size(files)
            call coils_read_GPEC(trim(base) // trim(files(i)), coils_ct)
            call load_coils_from_gpec_file(trim(base) // trim(files(i)), coils_neo)
            call compare_coil_sets(coils_ct, coils_neo, tolerance, 'file '//trim(files(i)))
            call cleanup_coiltools(coils_ct)
            call coils_deinit(coils_neo)
        end do
    end subroutine test_each_file

    subroutine test_combined_files(base, files, tolerance)
        character(len=*), intent(in) :: base
        character(len=*), intent(in) :: files(:)
        real(dp), intent(in) :: tolerance
        type(coil_t), allocatable :: combined_ct(:), temp_ct(:)
        type(coils_t) :: combined_neo, temp_neo
        integer :: i

        if (allocated(combined_ct)) call cleanup_coiltools(combined_ct)
        if (allocated(temp_ct)) call cleanup_coiltools(temp_ct)
        if (allocated(combined_neo%x)) call coils_deinit(combined_neo)
        if (allocated(temp_neo%x)) call coils_deinit(temp_neo)

        do i = 1, size(files)
            call coils_read_GPEC(trim(base) // trim(files(i)), temp_ct)
            if (i == 1) then
                call move_alloc(temp_ct, combined_ct)
            else
                call coils_append(combined_ct, temp_ct)
            end if

            call load_coils_from_gpec_file(trim(base) // trim(files(i)), temp_neo)
            if (.not. allocated(combined_neo%x)) then
                call move_alloc(temp_neo%x, combined_neo%x)
                call move_alloc(temp_neo%y, combined_neo%y)
                call move_alloc(temp_neo%z, combined_neo%z)
                call move_alloc(temp_neo%current, combined_neo%current)
            else
                call append_neo_coils(combined_neo, temp_neo)
            end if
            if (allocated(temp_neo%x)) call coils_deinit(temp_neo)
        end do

        call compare_coil_sets(combined_ct, combined_neo, tolerance, 'combined files')
        call cleanup_coiltools(combined_ct)
        call coils_deinit(combined_neo)
    end subroutine test_combined_files

    subroutine compare_coil_sets(coils_ct, coils_neo, tolerance, context)
        type(coil_t), allocatable, intent(in) :: coils_ct(:)
        type(coils_t), intent(in) :: coils_neo
        real(dp), intent(in) :: tolerance
        character(len=*), intent(in) :: context
        integer :: total_points_expected
        integer :: ic, kc, idx, closing_idx, nseg
        real(dp) :: diff_coord
        real(dp) :: xyz_neo(3)

        if (.not. allocated(coils_ct)) then
            write(error_unit, '(A)') 'coil_tools returned no coils for '//trim(context)
            error stop
        end if
        if (.not. allocated(coils_neo%x)) then
            write(error_unit, '(A)') 'neo loader returned no coils for '//trim(context)
            error stop
        end if

        total_points_expected = 0
        do ic = 1, size(coils_ct)
            total_points_expected = total_points_expected + coils_ct(ic)%nseg + 1
        end do

        if (size(coils_neo%x) /= total_points_expected) then
            write(error_unit, '(A)') 'Point count mismatch for '//trim(context)
            write(error_unit, '(A,I0)') 'neo count: ', size(coils_neo%x)
            write(error_unit, '(A,I0)') 'expected: ', total_points_expected
            error stop
        end if

        idx = 0
        do ic = 1, size(coils_ct)
            nseg = coils_ct(ic)%nseg
            do kc = 1, nseg
                if (abs(coils_neo%current(idx + kc)) < tolerance) then
                    write(error_unit, '(A,I0,A)') 'Zero current within coil segments at index ', idx + kc, &
                        ' for '//trim(context)
                    error stop
                end if
                if (abs(coils_neo%current(idx + kc) - real(coils_ct(ic)%nwind, dp)) > tolerance) then
                    write(error_unit, '(A,I0,A,F12.6,A,I0)') 'Current mismatch at index ', idx + kc, &
                        ' neo=', coils_neo%current(idx + kc), ' nwind=', coils_ct(ic)%nwind
                    error stop
                end if
                xyz_neo = [coils_neo%x(idx + kc), coils_neo%y(idx + kc), &
                    coils_neo%z(idx + kc)]
                diff_coord = maxval(abs(xyz_neo - coils_ct(ic)%XYZ(:, kc)))
                if (diff_coord > tolerance) then
                    write(error_unit, '(A,I0,A,1PE12.3)') 'Coordinate mismatch at index ', idx + kc, &
                        ' diff=', diff_coord
                    error stop
                end if
            end do
            closing_idx = idx + nseg + 1
            if (abs(coils_neo%current(closing_idx)) > tolerance) then
                write(error_unit, '(A,I0,A,1PE12.3)') 'Closing segment current not zero at index ', &
                    closing_idx, ' value=', coils_neo%current(closing_idx)
                error stop
            end if
            xyz_neo = [coils_neo%x(closing_idx), coils_neo%y(closing_idx), &
                coils_neo%z(closing_idx)]
            diff_coord = maxval(abs(xyz_neo - coils_ct(ic)%XYZ(:, 1)))
            if (diff_coord > tolerance) then
                write(error_unit, '(A,1PE12.3)') 'Closing point mismatch diff=', diff_coord
                error stop
            end if
            idx = idx + nseg + 1
        end do
    end subroutine compare_coil_sets

    subroutine cleanup_coiltools(coils)
        type(coil_t), allocatable, intent(inout) :: coils(:)
        integer :: i

        if (.not. allocated(coils)) return
        do i = 1, size(coils)
            call coil_deinit(coils(i))
        end do
        deallocate(coils)
    end subroutine cleanup_coiltools

    subroutine append_neo_coils(dst, src)
        type(coils_t), intent(inout) :: dst
        type(coils_t), intent(inout) :: src
        integer :: n_old, n_new
        real(dp), allocatable :: tmp(:)

        if (.not. allocated(src%x)) return
        if (.not. allocated(dst%x)) then
            call move_alloc(src%x, dst%x)
            call move_alloc(src%y, dst%y)
            call move_alloc(src%z, dst%z)
            call move_alloc(src%current, dst%current)
            return
        end if

        n_old = size(dst%x)
        n_new = n_old + size(src%x)

        allocate(tmp(n_new))
        tmp(1:n_old) = dst%x
        tmp(n_old + 1:n_new) = src%x
        call move_alloc(tmp, dst%x)

        allocate(tmp(n_new))
        tmp(1:n_old) = dst%y
        tmp(n_old + 1:n_new) = src%y
        call move_alloc(tmp, dst%y)

        allocate(tmp(n_new))
        tmp(1:n_old) = dst%z
        tmp(n_old + 1:n_new) = src%z
        call move_alloc(tmp, dst%z)

        allocate(tmp(n_new))
        tmp(1:n_old) = dst%current
        tmp(n_old + 1:n_new) = src%current
        call move_alloc(tmp, dst%current)
    end subroutine append_neo_coils

    subroutine ensure_trailing_slash(path)
        character(len=*), intent(inout) :: path
        integer :: l

        l = len_trim(path)
        if (l == 0) return
        if (path(l:l) /= '/' .and. path(l:l) /= '\\') then
            path = trim(path) // '/'
        end if
    end subroutine ensure_trailing_slash

end program test_gpec_loader_comparison
