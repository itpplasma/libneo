!> Dumps a MAKEGRID coils file as parsed by neo_makegrid_coils, and the
!> Biot-Savart polygon of one group, for test/scripts/compare_makegrid_coils.py.
program makegrid_coils_oracle
    use neo_biotsavart, only: dp, coils_t, save_coils_to_file
    use neo_makegrid_coils, only: makegrid_coils_t, read_makegrid_coils, &
        makegrid_to_biotsavart

    implicit none

    character(len=4096) :: coils_file, dump_file, group_file, arg
    type(makegrid_coils_t) :: mc
    type(coils_t) :: group_coils
    integer :: unit, ic, i, group

    call get_command_argument(1, coils_file)
    call get_command_argument(2, dump_file)
    call get_command_argument(3, group_file)
    call get_command_argument(4, arg)
    read(arg, *) group

    call read_makegrid_coils(trim(coils_file), mc)
    open(newunit=unit, file=trim(dump_file), status="replace", action="write")
    write(unit, '(a,1x,i0)') "periods", mc%periods
    write(unit, '(a,1x,a)') "mirror", trim(mc%mirror)
    write(unit, '(a,1x,i0,1x,i0)') "coils", mc%n_coils, mc%n_groups
    do ic = 1, mc%n_coils
        write(unit, '(a,1x,i0,1x,i0,1x,a)') "coil", mc%last(ic) - mc%first(ic) + 1, &
            mc%group(ic), trim(group_label(mc, ic))
        do i = mc%first(ic), mc%last(ic)
            write(unit, '(4es25.17)') mc%x(i), mc%y(i), mc%z(i), mc%current(i)
        end do
    end do
    close(unit)

    call makegrid_to_biotsavart(mc, group_coils, group=group)
    call save_coils_to_file(trim(group_file), group_coils)

contains

    function group_label(mc, ic) result(label)
        type(makegrid_coils_t), intent(in) :: mc
        integer, intent(in) :: ic
        character(len=64) :: label

        label = "-"
        if (mc%group(ic) < 1) return
        if (len_trim(mc%group_name(mc%group(ic))) > 0) &
            label = mc%group_name(mc%group(ic))
    end function group_label

end program makegrid_coils_oracle
