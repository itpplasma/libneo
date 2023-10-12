subroutine getran(irand,ur)
  !  Produces the random number with zero deviation and unit square
  !  deviation
  !
  !  Input parameters: irand - 0 for continious, 1 for +1 -1,
  !  Output parameters: ur   - random number

  integer irand
  real(KIND(1.0)) ur

  call random_number(ur)

  if (irand.eq.0) then

    !  continuos random number, constant is sqrt(12)
    ur=3.464102*(ur-.5)

  else

    !  discrete random number
    if(ur.gt..5) then
      ur=1.
    else
      ur=-1.
    end if
  end if

end subroutine getran
