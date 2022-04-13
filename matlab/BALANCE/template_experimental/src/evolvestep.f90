  subroutine evolvestep(timstep,eps)
!
  use grid_mod, only : nbaleqs,neqset,iboutype,npoic,params,y,dery
  USE sparse_mod, ONLY : sparse_talk,sparse_solve_method,sparse_solve, &
       column_full2pointer,remap_rc,sparse_solver_test
  use matrix_mod
  use recstep_mod, only : timstep_arr
!
  implicit none
!
  external :: rhs_balance !, rhs_func, gslint
  integer :: ipoi,ieq,i,k,npoi,iopt,nz_sp,nz_sq,nrow,ncol
  double precision :: timstep,x1,x2,eps,time_start,time_factorization,time_solver
  INTEGER,          DIMENSION(:),   ALLOCATABLE :: ipcol,irow_sp,icol_sp
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: amat_sp,bvec_sp
!
  x1=0.d0
!  x2=timstep
!
  if(iboutype.eq.1) then
    npoi=npoic-1
  else
    npoi=npoic
  endif
!
  do ipoi=1,npoi
    do ieq=1,nbaleqs
      i=nbaleqs*(ipoi-1)+ieq
      y(i)=params(ieq,ipoi)
    enddo
  enddo
!
  call initialize_rhs(y,dery)
  call rhs_balance(x1,y,dery)
!
  nz_sp=nz+nsize
  nrow=nsize
  ncol=nsize
!
  allocate(ipcol(nsize))
  allocate(amat_sp(nz_sp),irow_sp(nz_sp),icol_sp(nz_sp),bvec_sp(nsize))
  irow_sp(1:nz)=irow
  icol_sp(1:nz)=icol
!  amat_sp(1:nz)=-timstep*amat
  amat_sp(1:nz)=-timstep_arr(irow)*amat
!
  k=nz
  do i=1,nsize
    k=k+1
    irow_sp(k)=i 
    icol_sp(k)=i 
    amat_sp(k)=1.d0
  enddo
!
!  bvec_sp=y+timstep*dery
  bvec_sp=y+timstep_arr*dery
!
  call  remap_rc(nz_sp,nz_sq,irow_sp,icol_sp,amat_sp)
!
  nz_sp=nz_sq
!
  CALL column_full2pointer(icol_sp(1:nz_sp),ipcol)
!
!  iopt=1
  iopt=0
!
  call cpu_time(time_start)
  CALL sparse_solve(nrow,ncol,nz_sp,irow_sp(1:nz_sp),ipcol,amat_sp(1:nz_sp),       &
                    bvec_sp,iopt)
!
  call cpu_time(time_factorization)
!  print *,'factorization completed ',time_factorization - time_start,' sec'
!
!  iopt=2
!
! Solution of inhomogeneus equation (account of sources):
!
!  CALL sparse_solve(nrow,ncol,nz_sp,irow_sp(1:nz_sp),ipcol,amat_sp(1:nz_sp),       &
!                    bvec_sp,iopt)
!
  y=bvec_sp
!  iopt=3
!
!  CALL sparse_solve(nrow,ncol,nz_sp,irow_sp(1:nz),ipcol,amat_sp(1:nz_sp),bvec_sp,iopt)
!
  deallocate(ipcol)
  deallocate(amat_sp,irow_sp,icol_sp,bvec_sp)
!
  do ipoi=1,npoi
    do ieq=1,nbaleqs
      i=nbaleqs*(ipoi-1)+ieq
      params(ieq,ipoi)=y(i)
    enddo
  enddo
!
  return
  end subroutine evolvestep
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine genstartsource
!
  use grid_mod, only : y,dery,dery_equisource &
                     , nbaleqs,neqset,iboutype,npoic,params 

  use control_mod, only: iwrite
  use matrix_mod
!
  implicit none
!
  integer :: ipoi,ieq,i,npoi,icount,k
  double precision :: x
!
  if(iboutype.eq.1) then
    npoi=npoic-1
  else
    npoi=npoic
  endif
!
  do ipoi=1,npoi
    do ieq=1,nbaleqs
      i=nbaleqs*(ipoi-1)+ieq
      y(i)=params(ieq,ipoi)
    enddo
  enddo
!  
  call initialize_rhs(y,dery)
!
  dery_equisource=0.d0
  call rhs_balance(x,y,dery)
!
  do k=1,nz
     dery_equisource(irow(k))=dery_equisource(irow(k))-amat(k)*y(icol(k))
  end do
!
  dery_equisource=dery_equisource-rhsvec

!print *, '----'
  open(666,file='equisource.dat')
  do ipoi=1,npoi
    do ieq=1,nbaleqs
      i=nbaleqs*(ipoi-1)+ieq
      params(ieq,ipoi)=y(i)
    enddo
    write(666,*) dery_equisource(4*(ipoi-1)+1:4*ipoi)
  enddo
  close (666)
!
!
if (iwrite .eq. 1) then
    write(666,*) dery_equisource
    close (666)
end if

return
end subroutine genstartsource
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
