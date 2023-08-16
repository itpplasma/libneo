!> Program to test the geqdsk_tools module.
!>
!> Ideas for further tests:
!> - read test with 'defect' efit file (should fail)
!> - that geqdsk_classify after geqdsk_standardise gives cocos index 3
!> - that a second call of geqdsk_check_consistency does nothing
!> - that a second call of geqdsk_standardise does nothing
program test_geqdsk_tools
  implicit none

  logical :: sucess, sucess_all

  sucess_all = .true.

  call test_read_write(sucess)
  sucess_all = sucess_all .and. sucess

  call test_classify(sucess)
  sucess_all = sucess_all .and. sucess

  if (.not. sucess_all) then
    error stop "ERROR At least one test of geqdsk_tools failed."
  end if

contains

  !> Check if eqdsk file can be read, written, and read again.
  !>
  !> This check tries:
  !>   - to read an eqdsk file
  !>   - write it to file
  !>   - read the just written eqdsk file
  !>   - checks if some quantities are the same between read and reread
  !>     file.
  !>
  !> input:
  !> ------
  !> none
  !>
  !> output:
  !> -------
  !> sucess: logical, true if sucessfull, false if not.
  !>
  !> sideeffects:
  !> ------------
  !> Leaves written file behind, if check of values was not sucessful,
  !> for debugging purposes.
  !>
  !> limitations:
  !> ------------
  !> Checks only limited amount of quantities for equality.
  subroutine test_read_write(sucess)
    use geqdsk_tools, only : geqdsk_t ! variables/types
    use geqdsk_tools, only : geqdsk_read, geqdsk_write ! subroutines

    implicit none

    logical, intent(out) :: sucess

    type(geqdsk_t) :: geqdsk_in, geqdsk_out

    character(len=*), parameter :: infilename = '../tests/resources/input_efit_file.dat', &
                                 & outfilename = 'geqdsk_output'
    integer :: file_unit, ios

    sucess = .false.

    call geqdsk_read(geqdsk_in, infilename)
    call geqdsk_write(geqdsk_in, outfilename)
    call geqdsk_read(geqdsk_out, outfilename)

    ! Note: only integer values checked.
    if (geqdsk_in%nw == geqdsk_out%nw &
        & .and. geqdsk_in%nh == geqdsk_out%nh &
        & .and. geqdsk_in%nbbbs == geqdsk_out%nbbbs &
        & .and. geqdsk_in%limitr == geqdsk_in%limitr) then
      sucess = .true.

      open(newunit=file_unit, iostat=ios, file=outfilename, status='old')
      if (ios == 0) close(unit=file_unit, status='delete')
    else
      write(*,*) "ERROR values of eqdsk files do not match."
    end if
  end subroutine test_read_write

  !> Check classification of geqdsk file.
  !>
  !> A file with known cocos class is read and classified. Then it is
  !> checked if the cocos index is as expected.
  !>
  !> input:
  !> ------
  !> none
  !>
  !> output:
  !> -------
  !> sucess: logical, true if sucessfull, false if not.
  !>
  !> sideeffects:
  !> ------------
  !> none
  !>
  !> limitations:
  !> ------------
  !> Checks only single file, and thus only single configuration.
  subroutine test_classify(sucess)
    use geqdsk_tools, only : geqdsk_t ! variables/types
    use geqdsk_tools, only : geqdsk_read, geqdsk_classify ! subroutines

    implicit none

    logical, intent(out) :: sucess

    type(geqdsk_t) :: geqdsk

    character(len=*), parameter :: infilename = '../tests/resources/input_efit_file.dat'

    sucess = .false.

    call geqdsk_read(geqdsk, infilename)

    call geqdsk_classify(geqdsk)

    if (geqdsk%cocos%index == 15) then
      sucess = .true.
    else
      write(*,*) "ERROR cocos index ", geqdsk%cocos%index, " does not&
          & match the expected value 15!"
    end if

  end subroutine test_classify

end program test_geqdsk_tools
