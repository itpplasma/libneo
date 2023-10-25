program coil_convert

  use iso_fortran_env, only: dp => real64, error_unit
  use coil_tools, only: coil_t, coil_init, coil_deinit, coils_append, &
       process_fixed_number_of_args, check_number_of_args, &
       AUG_coils_write, AUG_coils_read, AUG_coils_write_Nemov, AUG_coils_read_Nemov, &
       AUG_coils_write_GPEC, AUG_coils_read_GPEC

  implicit none

  character(len = *), parameter :: incommensurable_fmt = &
       '("Cannot distribute ", i0, " coil(s) to ", i0, " output file(s).")'

  integer :: in_ncoil, out_ncoil, kc
  character(len = 1024) :: in_type, out_type
  character(len = :), dimension(:), allocatable :: in_files, out_files
  type(coil_t), dimension(:), allocatable :: coils, more_coils

  call check_number_of_args(1)
  call get_command_argument(1, in_type)
  call process_fixed_number_of_args(2, in_ncoil, in_files)
  call check_number_of_args(3 + in_ncoil)
  call get_command_argument(3 + in_ncoil, out_type)
  call process_fixed_number_of_args(4 + in_ncoil, out_ncoil, out_files)

  if (in_type == 'AUG') then
     allocate(coils(in_ncoil))
     do kc = 1, in_ncoil
        call AUG_coils_read(trim(in_files(kc)), coils(kc))
     end do
  else if (in_type == 'GPEC') then
     call AUG_coils_read_GPEC(trim(in_files(1)), coils)
     do kc = 2, in_ncoil
        call AUG_coils_read_GPEC(trim(in_files(kc)), more_coils)
        call coils_append(coils, more_coils)
     end do
  else if (in_type == 'Nemov') then
     call AUG_coils_read_Nemov(trim(in_files(1)), coils)
     do kc = 2, in_ncoil
        call AUG_coils_read_Nemov(trim(in_files(kc)), more_coils)
        call coils_append(coils, more_coils)
     end do
  else
     write (error_unit, '("Unknown input type ", a, ".")') trim(in_type)
     error stop
  end if
  ! not actually necessary for format AUG
  if (.not. allocated(coils)) then
     write (error_unit, '("Unexpected error: coils array not allocated.")')
     error stop
  end if
  in_ncoil = size(coils)

  if (out_type == 'AUG') then
     if (in_ncoil /= out_ncoil) then
        write (error_unit, incommensurable_fmt) in_ncoil, out_ncoil
        error stop
     end if
     do kc = 1, out_ncoil
        call AUG_coils_write(trim(out_files(kc)), coils(kc))
     end do
  else if (out_type == 'GPEC') then
     if (out_ncoil > 1) then
        if (in_ncoil /= out_ncoil) then
           write (error_unit, incommensurable_fmt) in_ncoil, out_ncoil
           error stop
        end if
        do kc = 1, out_ncoil
           call AUG_coils_write_GPEC(trim(out_files(kc)), coils(kc:kc))
        end do
     else
        call AUG_coils_write_GPEC(trim(out_files(1)), coils)
     end if
  else if (out_type == 'Nemov') then
     if (out_ncoil > 1) then
        if (in_ncoil /= out_ncoil) then
           write (error_unit, incommensurable_fmt) in_ncoil, out_ncoil
           error stop
        end if
        do kc = 1, out_ncoil
           call AUG_coils_write_Nemov(trim(out_files(kc)), coils(kc:kc))
        end do
     else
        call AUG_coils_write_Nemov(trim(out_files(1)), coils)
     end if
  else
     write (error_unit, '("Unknown output type ", a, ".")') trim(out_type)
     error stop
  end if

  if (allocated(coils)) then
     do kc = 1, size(coils)
        call coil_deinit(coils(kc))
     end do
     deallocate(coils)
  end if
  if (allocated(in_files)) deallocate(in_files)
  if (allocated(out_files)) deallocate(out_files)

end program coil_convert
