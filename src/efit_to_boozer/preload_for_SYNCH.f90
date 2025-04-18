  program preload_for_SYNCH
!
  use efit_to_boozer_mod
!
  implicit none
!
  integer :: nstep,nsurfmax,i
!
  double precision, dimension(:,:), allocatable :: R_st,Z_st,bmod_st,sqgnorm_st,Gfunc
!
  open(1,file='preload_for_SYNCH.inp')
  read (1,*) nstep    !number of integration steps
  read (1,*) nlabel   !grid size over radial variabl
  read (1,*) ntheta   !grid size over poloidal angle
  read (1,*) nsurfmax !number of starting points between the
                      !magnetic axis and right box boundary
                      !when searching for the separatrix
  close(1)
!
  allocate(rbeg(nlabel),rsmall(nlabel),qsaf(nlabel),psi_pol(0:nlabel))
  allocate(psi_tor_vac(nlabel),psi_tor(0:nlabel),C_const(nlabel))
!
  allocate(R_spl(0:nspl,0:ntheta,nlabel),Z_spl(0:nspl,0:ntheta,nlabel),bmod_spl(0:nspl,0:ntheta,nlabel))
  allocate(sqgnorm_spl(0:nspl,0:ntheta,nlabel),Gfunc_spl(0:nspl,0:ntheta,nlabel))
!
  call flint_for_Boozer(nstep,nsurfmax,nlabel,ntheta,     &
                        rmn,rmx,zmn,zmx,raxis,zaxis,      &
                        rbeg,rsmall,qsaf,                 &
                        psi_pol(1:nlabel),psi_tor_vac,    &
                        psi_tor(1:nlabel),C_const,        &
                        R_spl(0,1:ntheta,:),              &
                        Z_spl(0,1:ntheta,:),              &
                        bmod_spl(0,1:ntheta,:),           &
                        sqgnorm_spl(0,1:ntheta,:),        &
                        Gfunc_spl(0,1:ntheta,:))
!
  call spline_magdata_in_symfluxcoord
!
  allocate(R_st(nlabel,ntheta),Z_st(nlabel,ntheta),bmod_st(nlabel,ntheta),sqgnorm_st(nlabel,ntheta))
  allocate(Gfunc(nlabel,ntheta))
!
  R_st=transpose(R_spl(0,1:ntheta,:))
  Z_st=transpose(Z_spl(0,1:ntheta,:))
  bmod_st=transpose(bmod_spl(0,1:ntheta,:))
  sqgnorm_st=transpose(sqgnorm_spl(0,1:ntheta,:))
!
  open(1,form='formatted',file='box_size_axis.dat')
  write (1,*) rmn,rmx, '<= rmn, rmx (cm)'
  write (1,*) zmn,zmx, '<= zmn, zmx (cm)'
  write (1,*) raxis,zaxis, '<= raxis, zaxis (cm)'
  close(1)
!
  open(1,form='formatted',file='flux_functions.dat')
  write (1,*) '# R_beg, r,  q, psi_pol, psi_tor_vac, psi_tor'
  do i=1,nlabel
    write (1,*) rbeg(i),rsmall(i),qsaf(i),psi_pol(i),psi_tor_vac(i),psi_tor(i)*psitor_max
  enddo
  close(1)
!
  open(1,form='formatted',file='twodim_functions.dat')
  write (1,*) nlabel, ntheta, '<= nlabel, ntheta'
  write (1,*) 'R(label,theta)'
  do i=1,nlabel
    write (1,*) R_st(i,:)
  enddo
  write (1,*) 'Z(label,theta)'
  do i=1,nlabel
    write (1,*) Z_st(i,:)
  enddo
  write (1,*) 'B(label,theta)'
  do i=1,nlabel
    write (1,*) bmod_st(i,:)
  enddo
  write (1,*) 'sqrtg_norm(label,theta)'
  do i=1,nlabel
    write (1,*) sqgnorm_st(i,:)
  enddo
  close(1)
!
  deallocate(rbeg,rsmall,qsaf,psi_pol,psi_tor_vac,psi_tor,C_const)
  deallocate(R_st,Z_st,bmod_st,sqgnorm_st,Gfunc)
!
  end program preload_for_SYNCH
