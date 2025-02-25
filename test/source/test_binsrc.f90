!> Program to test the binsrc subroutine.
!>
!> Creates simple array (value = index) and checks that for values within
!> the minimum and maximum values, as well as for smaller and larger
!> values the expected index is returned.
!>
!> There is no test with incorrect input.
program test_binsrc

  use binsrc_sub, only : binsrc
  use libneo_kinds, only : dp

  implicit none

  integer, parameter :: min_index = 1, max_index = 20

  logical :: sucess

  integer :: found_index
  integer :: k

  real(dp) :: value_to_find
  real(dp), dimension(min_index:max_index) :: values

  sucess = .true.

  values = (/ (k, k=min_index, max_index) /)

  ! Test with value between minimum and maximum of array.
  value_to_find = 5.5
  call binsrc(values, min_index, max_index, value_to_find, found_index)
  if (found_index /= floor(value_to_find)+1) then
    sucess = .false.
    write(*,*) "Error found index, ", found_index, &
        & ", differs from expected value, ", floor(value_to_find)+1, &
        & ", when looking for value ", value_to_find
  end if

  ! Test with value below minimum value of array.
  value_to_find = min_index - 1.0
  call binsrc(values, min_index, max_index, value_to_find, found_index)
  if (found_index /= min_index+1) then
    sucess = .false.
    write(*,*) "Error found index, ", found_index, &
        & ", differs from expected value, ", min_index, &
        & ", when looking for value ", value_to_find
  end if

  ! Test with value larger than maximum value of array.
  value_to_find = max_index + 1.0
  call binsrc(values, min_index, max_index, value_to_find, found_index)
  if (found_index /= max_index) then
    sucess = .false.
    write(*,*) "Error found index, ", found_index, &
        & ", differs from expected value, ", max_index, &
        & ", when looking for value ", value_to_find
  end if

  if (.not. sucess) then
    error stop "ERROR At least one test of binsrc failed."
  end if
end program test_binsrc
