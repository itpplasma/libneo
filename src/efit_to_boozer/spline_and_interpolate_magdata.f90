!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine spline_magdata_in_symfluxcoord
!
! Splines magnetic data over theta
!
  use spl_three_to_five_sub
  use efit_to_boozer_mod
!
  implicit none
!
  integer :: i,nthetap1
  double precision, dimension(:,:), allocatable :: splcoe
!
!-----------------------------------------------------------------------
!
  allocate(splcoe(0:nspl,0:ntheta))
  h_theta = twopi/dfloat(ntheta)
  nthetap1=ntheta+1
!
  do i=1,nlabel
    splcoe(0,1:ntheta)=R_spl(0,1:ntheta,i)
    splcoe(0,0)=splcoe(0,ntheta)
!
    call spl_per(nspl,nthetap1,h_theta,splcoe)
!
    R_spl(:,:,i)=splcoe
!
    splcoe(0,1:ntheta)=Z_spl(0,1:ntheta,i)
    splcoe(0,0)=splcoe(0,ntheta)
!
    call spl_per(nspl,nthetap1,h_theta,splcoe)
!
    Z_spl(:,:,i)=splcoe
!
    splcoe(0,1:ntheta)=bmod_spl(0,1:ntheta,i)
    splcoe(0,0)=splcoe(0,ntheta)
!
    call spl_per(nspl,nthetap1,h_theta,splcoe)
!
    bmod_spl(:,:,i)=splcoe
!
    splcoe(0,1:ntheta)=sqgnorm_spl(0,1:ntheta,i)
    splcoe(0,0)=splcoe(0,ntheta)
!
    call spl_per(nspl,nthetap1,h_theta,splcoe)
!
    sqgnorm_spl(:,:,i)=splcoe
!
    splcoe(0,1:ntheta)=Gfunc_spl(0,1:ntheta,i)
    splcoe(0,0)=splcoe(0,ntheta)
!
    call spl_per(nspl,nthetap1,h_theta,splcoe)
!
    Gfunc_spl(:,:,i)=splcoe
  enddo
!
  deallocate(splcoe)
!
!-----------------------------------------------------------------------
!
  psi_pol(0)=0.d0
  psi_tor(0)=0.d0
  psipol_max=psi_pol(nlabel)
  psitor_max=psi_tor(nlabel)
  psi_pol=psi_pol/psipol_max
  psi_tor=psi_tor/psitor_max
!
  end subroutine spline_magdata_in_symfluxcoord
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine magdata_in_symfluxcoord_ext(inp_label,s,psi,theta,q,dq_ds,C_norm,dC_norm_ds, &
                                         sqrtg,bmod,dbmod_dtheta,R,dR_ds,dR_dtheta,       &
                                         Z,dZ_ds,dZ_dtheta,G,dG_ds,dG_dtheta)
!
! Computes safety factor, sqrt(g) of symmetry flux coordinates, module of B, cylindrical
! coordinates R and Z and transformation function G as functions of the flux surface label
! and poloidal angle of symmetry flux coordinates.
! Computes also derivative of module-B over the poloidal angle and derivatives of R, Z and transformation
! function G over the nomalized toroidal or poloidal flux and poloidal angle theta.
! Uses periodic spline interpolation over theta and Lagrange polynomial interpolation over flux srface label.
! Two flux surface labels can be used as an input: normalized toroidal flux s and dimensional poloidal
! flux psi, depending on the input switch inp_label. If inp_label=1 s is an input variable and psi=psi(s) is computed.
! If inp_label=2, psi is an input variable, and s=s(psi) is computed.
!
! Input (inout) arguments:
!                 inp_label - input switch: 1 for s and 2 for psi
!                 s         - normalized toroidal flux (inout)
!                 psi       - poloidal flux (inout)
!                 theta     - poloidal angle of symmetry flux coordinates
! Output arguments:
!                 q         - safety factor
!                 dq_ds     - derivative of safety factor over normalized toroidal or poloidal flux
!                 C_norm    - normalization constant of transformation C
!                 dC_norm_ds- derivative of C over normalized toroidal or poloidal flux
!                 sqrtg     - normalized metric determinant sqrt(g)
!                 bmod      - module of B
!                 dbmod_dt  - derivative of module of B over theta
!                 R         - radius R of cylindrical coordinates
!                 dR_ds     - derivative of R over normalized toroidal or poloidal flux
!                 dR_dtheta - derivative of R over polidal angle of symmetry flux coordinates theta
!                 Z         - height Z of cylindrical coordinates
!                 dZ_ds     - derivative of Z over normalized toroidal or poloidal flux
!                 dZ_dtheta - derivative of Z over polidal angle of symmetry flux coordinates theta
!                 G         - transformation function to Boozer coordinates
!                 dG_ds     - derivative of G over normalized toroidal or poloidal flux
!                 dG_dtheta - derivative of G over polidal angle of symmetry flux coordinates theta
!
  use plag_coeff_sub
  use binsrc_sub
  use efit_to_boozer_mod
!
  implicit none
!
  integer, intent(in) :: inp_label
  double precision, intent(inout) :: s,psi,theta
  double precision, intent(out) :: q,dq_ds,C_norm,dC_norm_ds,sqrtg,bmod,dbmod_dtheta, &
                      R,dR_ds,dR_dtheta,Z,dZ_ds,dZ_dtheta,G,dG_ds,dG_dtheta
  integer :: it,k,km1,ibeg,iend
  double precision :: dtheta
!
  if(inp_label.eq.1) then
!
    call binsrc(psi_tor(0:nlabel),0,nlabel,s,ibeg)
!
    ibeg=max(1,ibeg-nplag/2)
    iend=ibeg+nplag-1
    if(iend.gt.nlabel) then
      iend=nlabel
      ibeg=iend+1-nplag
    endif
!
    call plag_coeff(nplag,nder,s,psi_tor(ibeg:iend),coef)
!
    psi=sum(coef(0,:)*psi_pol(ibeg:iend))*psipol_max
  elseif(inp_label.eq.2) then
    s=psi/psipol_max
!
    call binsrc(psi_pol(0:nlabel),0,nlabel,s,ibeg)
!
    ibeg=max(1,ibeg-nplag/2)
    iend=ibeg+nplag-1
    if(iend.gt.nlabel) then
      iend=nlabel
      ibeg=iend+1-nplag
    endif
!
    call plag_coeff(nplag,nder,s,psi_pol(ibeg:iend),coef)
!
    s=sum(coef(0,:)*psi_tor(ibeg:iend))
  else
    print *,'unknown mode for inp_label =',inp_label
    return
  endif
!
  q=sum(coef(0,:)*qsaf(ibeg:iend))
  dq_ds=sum(coef(1,:)*qsaf(ibeg:iend))
  C_norm=sum(coef(0,:)*C_const(ibeg:iend))
  dC_norm_ds=sum(coef(1,:)*C_const(ibeg:iend))
!
  dtheta=modulo(theta,twopi)/h_theta
  it=max(0,min(ntheta-1,int(dtheta)))
  dtheta=(dtheta-dfloat(it))*h_theta
!
  sqrtg_lag=sqgnorm_spl(nspl,it,ibeg:iend)
  bmod_lag=bmod_spl(nspl,it,ibeg:iend)
  R_lag=R_spl(nspl,it,ibeg:iend)
  Z_lag=Z_spl(nspl,it,ibeg:iend)
  G_lag=Gfunc_spl(nspl,it,ibeg:iend)
!
  dbmod_dt_lag=0.d0
  dR_dt_lag=0.d0
  dZ_dt_lag=0.d0
  dG_dt_lag=0.d0
!
  do k=nspl,1,-1
    km1=k-1
    sqrtg_lag=sqrtg_lag*dtheta+sqgnorm_spl(km1,it,ibeg:iend)
    bmod_lag=bmod_lag*dtheta+bmod_spl(km1,it,ibeg:iend)
    R_lag=R_lag*dtheta+R_spl(km1,it,ibeg:iend)
    Z_lag=Z_lag*dtheta+Z_spl(km1,it,ibeg:iend)
    G_lag=G_lag*dtheta+Gfunc_spl(km1,it,ibeg:iend)
!
    dbmod_dt_lag=dbmod_dt_lag*dtheta+bmod_spl(k,it,ibeg:iend)*dfloat(k)
    dR_dt_lag=dR_dt_lag*dtheta+R_spl(k,it,ibeg:iend)*dfloat(k)
    dZ_dt_lag=dZ_dt_lag*dtheta+Z_spl(k,it,ibeg:iend)*dfloat(k)
    dG_dt_lag=dG_dt_lag*dtheta+Gfunc_spl(k,it,ibeg:iend)*dfloat(k)
  enddo
!
  sqrtg=sum(coef(0,:)*sqrtg_lag)
  bmod=sum(coef(0,:)*bmod_lag)
  R=sum(coef(0,:)*R_lag)
  Z=sum(coef(0,:)*Z_lag)
  G=sum(coef(0,:)*G_lag)
!
  dR_ds=sum(coef(1,:)*R_lag)
  dZ_ds=sum(coef(1,:)*Z_lag)
  dG_ds=sum(coef(1,:)*G_lag)
!
  dbmod_dtheta=sum(coef(0,:)*dbmod_dt_lag)
  dR_dtheta=sum(coef(0,:)*dR_dt_lag)
  dZ_dtheta=sum(coef(0,:)*dZ_dt_lag)
  dG_dtheta=sum(coef(0,:)*dG_dt_lag)
!
  end subroutine magdata_in_symfluxcoord_ext
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
