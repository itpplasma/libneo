!> Matrix module, which can pack itself
module matrix_module
  use packable_module
  
  implicit none
  
  !> Class matrix
  type, extends(packable) :: matrix
     real, dimension(:,:), allocatable :: mat    !< data
   contains     
     procedure :: print  => print_matrix
     procedure :: pack   => pack_matrix
     procedure :: unpack => unpack_matrix
     procedure :: free   => free_matrix
     
     procedure :: alloc => alloc_matrix
  end type matrix
  
contains
  
  !> Wrapper for allocation
  subroutine alloc_matrix(this, a, b, c, d)
    class(matrix) :: this
    integer, intent(in) :: a, b, c, d
    
    allocate(this%mat(a:b, c:d))
  end subroutine alloc_matrix
  
  !> Free memory
  subroutine free_matrix(this)
    class(matrix) :: this
    
    if (allocated(this%mat)) deallocate(this%mat)
  end subroutine free_matrix
  
  subroutine unpack_matrix(this)
    class(matrix) :: this
    write (*,*) "Matrix: Unpack not yet implemented"
    
  end subroutine unpack_matrix
  
  subroutine pack_matrix(this)
    class(matrix) :: this
    write (*,*) "Matrix: Pack not yet implemented!"
    
  end subroutine pack_matrix
  
  subroutine print_matrix(this)
    class(matrix) :: this
    integer :: ub1, ub2, lb1, lb2
    
    lb1 = lbound(this%mat, 1)
    ub1 = ubound(this%mat, 1)
    
    lb2 = lbound(this%mat, 2)
    ub2 = ubound(this%mat, 2)

    write (*,*) this%uid, "Result dimensions: ", lb1, ub1, lb2, ub2

    ! Nice printout
    !do i = lb1, ub1
    !   do j = lb2, ub2
    !      write (*, '(E15.5)', advance='no') this%mat(i, j)
    !   end do
    !   write (*,*)
    !end do
    
  end subroutine print_matrix
  
end module matrix_module
