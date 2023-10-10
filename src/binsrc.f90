! Finds the index  i  of the array of increasing numbers   p  with dimension  n 
! which satisfies   p(i-1) <  xi  <  p(i) . Uses binary search algorithm.
!
!
! Note: there is no check if p(nmin) < xi < p(nmax).
!   If xi <= p(nmin) then nmin+1 is returned.
!   If xi >= p(nmax) then nmax is returned.
subroutine binsrc(p, nmin, nmax, xi, i)

  implicit none

  integer, intent(in) :: nmin, nmax
  integer, intent(out) :: i
  double precision, intent(in) :: xi
  double precision, dimension(nmin:nmax), intent(in) :: p

  integer :: n, imin, imax, k

  imin = nmin
  imax = nmax
  n = nmax-nmin

  do k=1,n
    i = (imax-imin)/2 + imin
    if (p(i) > xi) then
      imax = i
    else
      imin = i
    end if
    if (imax == imin+1) exit
  end do

  i = imax

end
