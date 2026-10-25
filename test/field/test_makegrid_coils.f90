program test_makegrid_coils
    use util_for_test, only: print_test, print_ok, print_fail

    implicit none

    integer, parameter :: dp = kind(1.0d0)
    real(dp), parameter :: pi = 3.141592653589793238462643383279502884197d0
    character(*), parameter :: coils_file = "coils.makegrid_test"
    integer, parameter :: n_seg = 24

    call write_synthetic_coils_file
    call test_metadata_and_groups
    call test_group_superposition
    call test_group_point_list
    call test_rejects_short_point_line

contains

    subroutine write_synthetic_coils_file
        integer :: unit, k

        open(newunit=unit, file=coils_file, status="replace", action="write")
        write(unit, '(a)') "# synthetic MAKEGRID test file"
        write(unit, '(a)') "periods 3"
        write(unit, '(a)') "begin filament"
        write(unit, '(a)') "mirror NIL"
        do k = 0, 2
            call write_planar_tf_coil(unit, 2.0d0*pi*k/3.0d0, 7.0d0, 1, "TF")
        end do
        call write_ring(unit, 1.6d0, 0.5d0, -3.0d0, 2, "PF")
        call write_ring(unit, 1.6d0, -0.5d0, -3.0d0, 2, "PF")
        call write_ring(unit, 0.9d0, 0.0d0, 2.0d0, 3, "")
        write(unit, '(a)') "end"
        close(unit)
    end subroutine write_synthetic_coils_file

    subroutine write_planar_tf_coil(unit, phi, current, group, name)
        integer, intent(in) :: unit, group
        real(dp), intent(in) :: phi, current
        character(*), intent(in) :: name
        integer :: i
        real(dp) :: t, r, z

        do i = 0, n_seg
            t = 2.0d0*pi*i/n_seg
            r = 1.0d0 + 0.45d0*cos(t)
            z = 0.6d0*sin(t)
            call write_point(unit, r*cos(phi), r*sin(phi), z, current, i, group, name)
        end do
    end subroutine write_planar_tf_coil

    subroutine write_ring(unit, radius, z, current, group, name)
        integer, intent(in) :: unit, group
        real(dp), intent(in) :: radius, z, current
        character(*), intent(in) :: name
        integer :: i
        real(dp) :: t

        do i = 0, n_seg
            t = 2.0d0*pi*i/n_seg
            call write_point(unit, radius*cos(t), radius*sin(t), z, current, i, &
                group, name)
        end do
    end subroutine write_ring

    subroutine write_point(unit, x, y, z, current, i, group, name)
        integer, intent(in) :: unit, i, group
        real(dp), intent(in) :: x, y, z, current
        character(*), intent(in) :: name

        if (i < n_seg) then
            write(unit, '(4es24.15)') x, y, z, current
        else
            write(unit, '(4es24.15,i4,1x,a)') x, y, z, 0.0d0, group, name
        end if
    end subroutine write_point

    subroutine test_metadata_and_groups
        use neo_makegrid_coils, only: makegrid_coils_t, read_makegrid_coils

        type(makegrid_coils_t) :: mc

        call print_test("makegrid_metadata_and_groups")
        call read_makegrid_coils(coils_file, mc)
        call check(mc%periods == 3, "periods")
        call check(trim(mc%mirror) == "NIL", "mirror")
        call check(mc%n_coils == 6, "coil count")
        call check(mc%n_groups == 3, "group count")
        call check(all(mc%group == [1, 1, 1, 2, 2, 3]), "coil groups")
        call check(trim(mc%group_name(1)) == "TF", "group name 1")
        call check(trim(mc%group_name(2)) == "PF", "group name 2")
        call check(len_trim(mc%group_name(3)) == 0, "group name 3")
        call check(all(mc%last - mc%first + 1 == n_seg), "points per coil")
        call check(abs(mc%current(mc%first(4)) + 3.0d0) == 0.0d0, "coil current")
        call print_ok
    end subroutine test_metadata_and_groups

    subroutine test_group_superposition
        use neo_biotsavart, only: coils_t, compute_magnetic_field
        use neo_makegrid_coils, only: makegrid_coils_t, read_makegrid_coils, &
            makegrid_to_biotsavart

        type(makegrid_coils_t) :: mc
        type(coils_t) :: all_coils, group_coils(3)
        real(dp) :: x(3), b_all(3), b_sum(3)
        integer :: g, ip

        call print_test("makegrid_group_superposition")
        call read_makegrid_coils(coils_file, mc)
        call makegrid_to_biotsavart(mc, all_coils)
        do g = 1, 3
            call makegrid_to_biotsavart(mc, group_coils(g), group=g)
        end do
        do ip = 1, 5
            x = [0.3d0 + 0.21d0*ip, 0.17d0*ip - 0.4d0, 0.05d0*ip - 0.1d0]
            b_all = compute_magnetic_field(all_coils, x)
            b_sum = 0.0d0
            do g = 1, 3
                b_sum = b_sum + compute_magnetic_field(group_coils(g), x)
            end do
            call check(norm2(b_all - b_sum) <= 1.0d-14*norm2(b_all), &
                "superposition")
        end do
        call print_ok
    end subroutine test_group_superposition

    subroutine test_group_point_list
        use neo_biotsavart, only: coils_t, compute_magnetic_field, &
            load_coils_from_file
        use neo_makegrid_coils, only: makegrid_coils_t, read_makegrid_coils, &
            makegrid_to_biotsavart

        character(*), parameter :: list_file = "coils.makegrid_test_group2"
        type(makegrid_coils_t) :: mc
        type(coils_t) :: from_group, from_list
        real(dp) :: x(3), b_group(3), b_list(3)
        integer :: unit, i, ip, k
        real(dp) :: t

        call print_test("makegrid_group_vs_point_list")
        open(newunit=unit, file=list_file, status="replace", action="write")
        write(unit, *) 2*(n_seg + 1)
        do k = 1, 2
            do i = 0, n_seg
                t = 2.0d0*pi*mod(i, n_seg)/n_seg
                write(unit, '(4es24.15)') 1.6d0*cos(t), 1.6d0*sin(t), &
                    merge(0.5d0, -0.5d0, k == 1), merge(-3.0d0, 0.0d0, i < n_seg)
            end do
        end do
        close(unit)

        call read_makegrid_coils(coils_file, mc)
        call makegrid_to_biotsavart(mc, from_group, group=2)
        call load_coils_from_file(list_file, from_list)
        call check(size(from_group%x) == size(from_list%x), "point count")
        do ip = 1, 5
            x = [0.2d0*ip, 0.1d0 - 0.13d0*ip, 0.07d0*ip]
            b_group = compute_magnetic_field(from_group, x)
            b_list = compute_magnetic_field(from_list, x)
            call check(norm2(b_group - b_list) <= 1.0d-14*norm2(b_list), &
                "group field vs point list")
        end do
        call print_ok
    end subroutine test_group_point_list

    subroutine test_rejects_short_point_line
        use neo_makegrid_coils, only: makegrid_coils_t, read_makegrid_coils

        character(*), parameter :: bad_file = "coils.makegrid_test_bad"
        type(makegrid_coils_t) :: mc
        integer :: unit, ierr

        call print_test("makegrid_rejects_short_point_line")
        open(newunit=unit, file=bad_file, status="replace", action="write")
        write(unit, '(a)') "begin filament"
        write(unit, '(a)') "1.0 0.0 0.0"
        close(unit)
        call read_makegrid_coils(bad_file, mc, ierr)
        call check(ierr /= 0, "short point line rejected")
        call print_ok
    end subroutine test_rejects_short_point_line

    subroutine check(condition, what)
        logical, intent(in) :: condition
        character(*), intent(in) :: what

        if (.not. condition) then
            print *, "FAILED: ", what
            call print_fail
            error stop
        end if
    end subroutine check

end program test_makegrid_coils
