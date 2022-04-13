
subroutine rhs_balance(x,y,dy)

use grid_mod, only :   nbaleqs,neqset,iboutype,npoic,npoib &
                     , params,ddr_params,Sc,Sb,deriv_coef  &
                     , ipbeg,ipend,rb,params_b,reint_coef  &
                     , fluxes_dif,fluxes_con,rc,dot_params &
                     , dae11,dae12,dae22,dai11,dai12,dai22 &
                     , dni22,visca,gpp_av,dery_equisource  &
                     , dqle11,dqle12,dqle21,dqle22         &
                     , dqli11,dqli12,dqli21,dqli22         &
                     , sqg_bthet_overc,Ercov,polforce,qlheat_e,qlheat_i &
                     , params_lin,Ercov_lin,ddr_params_nl,fluxes_con_nl &
                     , init_params, params_b_lin
use baseparam_mod, only : Z_i,e_charge,am,p_mass,c,btor
use control_mod, only : iwrite
use wave_code_data, only : q,Vth

use matrix_mod, only : isw_rhs,nz,nsize,irow,icol,amat,rhsvec

implicit none

integer :: ipoi,ieq,i,npoi,ibeg,iend,nshift,ibegb,iendb,ibegtot,iendtot,k,iprobe
double precision :: x,A_noE_1e,A_noE_2e,A_noE_1i,A_noE_2i,convel
double precision :: A_noE_1e_nl,A_noE_2e_nl,A_noE_1i_nl,A_noE_2i_nl
double precision :: gamma_e,gamma_i,dfluxvphi,Q_e,Q_i,A_1e,A_1i
double precision :: gamma_e_nl,Q_e_nl,Q_i_nl,A_1e_nl,A_1i_nl
double precision :: gamma_ql_e,gamma_ql_i
double precision, dimension(neqset) :: y,dy,y_lin
!
  if(iboutype.eq.1) then
    npoi=npoic-1
  else
    npoi=npoic
  endif
!
! equilibrium parameters:
!
y_lin = 0.0d0
params_lin = 0.0d0
ddr_params=0.0d0
Ercov_lin=0.0d0

do ipoi=1,npoi
    do ieq=1,nbaleqs
      i=nbaleqs*(ipoi-1)+ieq
      params(ieq,ipoi)=y(i)
    end do
end do

!
! Interpolation:
!
  do ipoi=1,npoib
    do ieq=1,nbaleqs
! radial derivatives of equilibrium parameters at cell boundaries:
      ddr_params_nl(ieq,ipoi) &
         =sum(params(ieq,ipbeg(ipoi):ipend(ipoi))*deriv_coef(:,ipoi))
! equilibrium parameters at cell boundaries:
      params_b(ieq,ipoi)   &
         =sum(params(ieq,ipbeg(ipoi):ipend(ipoi))*reint_coef(:,ipoi))
    enddo
  enddo
!
! Compute radial electric field:
!
!   Ercov=sqg_bthet_overc*params_b(2,:)                                 & !OLD
   Ercov=sqg_bthet_overc*(params_b(2,:)-Vth*q/rb)                           &
       +(params_b(4,:)*ddr_params_nl(1,:)/params_b(1,:)+ddr_params_nl(4,:)) &
       /(Z_i*e_charge)


! Compute diffusion coefficient matrices:
!
  call calc_dequi
!
  do ipoi=1,npoib
!
! Thermodynamic forces for zero radial electric field:
!
    A_noE_1e_nl=ddr_params_nl(1,ipoi)/params_b(1,ipoi)  &
            -1.5d0*ddr_params_nl(3,ipoi)/params_b(3,ipoi)
    A_noE_2e_nl=ddr_params_nl(3,ipoi)/params_b(3,ipoi)
    A_noE_1i_nl=ddr_params_nl(1,ipoi)/params_b(1,ipoi)  &
            -1.5d0*ddr_params_nl(4,ipoi)/params_b(4,ipoi)
    A_noE_2i_nl=ddr_params_nl(4,ipoi)/params_b(4,ipoi)
!
! Thermodynamic forces for finite radial electric field:
!
    A_1e_nl=A_noE_1e_nl+Ercov(ipoi)*e_charge/params_b(3,ipoi)
!ERROR    A_1i_nl=A_noE_1e_nl-Ercov(ipoi)*e_charge*Z_i/params_b(4,ipoi)
    A_1i_nl=A_noE_1i_nl-Ercov(ipoi)*e_charge*Z_i/params_b(4,ipoi) !<-FIXED

!
! particle flux densities:
!
    gamma_e_nl=-(dae11(ipoi)*A_noE_1e_nl+dae12(ipoi)*A_noE_2e_nl      &
              +  dqle11(ipoi)*A_1e_nl+dqle12(ipoi)*A_noE_2e_nl        &
                 )*params_b(1,ipoi)

! total particle flux:
    fluxes_con_nl(1,ipoi)=(Sb(ipoi)*gamma_e_nl - &
     (-Sb(ipoi)*ddr_params_nl(1,ipoi)*(dae11(ipoi)              &
      +dqle11(ipoi)*(1.d0+params_b(4,ipoi)/params_b(3,ipoi)/Z_i))))/params_b(1,ipoi)
!
! toroidal moment flux density divided by mass:
! total toroidal moment flux:
    fluxes_con_nl(2,ipoi)=0.d0
!
! electron heat flux density:
!
!colli    Q_e_nl=-(dae12(ipoi)*A_noE_1e_nl+dqle12(ipoi)*A_1e_nl    &
    Q_e_nl=-(dae12(ipoi)*A_noE_1e_nl+dqle21(ipoi)*A_1e_nl    &
       +  (dae22(ipoi)+dqle22(ipoi))*A_noE_2e_nl)            &
       * params_b(1,ipoi)*params_b(3,ipoi)
!
! ion heat flux density:
!
!colli    Q_i_nl=-(dai12(ipoi)*A_noE_1i_nl+dqli12(ipoi)*A_1i_nl    &
    Q_i_nl=-(dai12(ipoi)*A_noE_1i_nl+dqli21(ipoi)*A_1i_nl    &
       +  (dai22(ipoi)+dni22(ipoi)+dqli22(ipoi))*A_noE_2i_nl)&
       * params_b(1,ipoi)/Z_i*params_b(4,ipoi)
!
! total heat fluxes:
!
    fluxes_con_nl(3,ipoi)=(Sb(ipoi)*Q_e_nl - &
    (-Sb(ipoi)*(dae22(ipoi)+dqle22(ipoi))*params_b(1,ipoi)*ddr_params_nl(3,ipoi))) &
    /params_b(3,ipoi)
!
    fluxes_con_nl(4,ipoi)=(Sb(ipoi)*Q_i_nl - &
!colli    (-Sb(ipoi)*(dai22(ipoi)+dni22(ipoi)+dqli22(ipoi)-2.5d0*dqli12(ipoi)) &
     (-Sb(ipoi)*(dai22(ipoi)+dni22(ipoi)+dqli22(ipoi)-2.5d0*dqli21(ipoi)) &
      *params_b(1,ipoi)/Z_i*ddr_params_nl(4,ipoi))) &
     /params_b(4,ipoi)
  enddo
!
  nshift=4
  k=0
  dy=0.d0
!
!$$$$
do iprobe=1,neqset
!
  ibeg=iprobe/nbaleqs-nshift
  iend=iprobe/nbaleqs+nshift
  ibegb=ibeg-1
  iendb=iend+1
  ibegtot=iprobe-nshift*nbaleqs
  iendtot=iprobe+nshift*nbaleqs
  ibeg=max(1,ibeg)
  iend=min(npoi,iend)
  ibegb=max(1,ibegb)
  iendb=min(npoib,iendb)
  ibegtot=max(1,ibegtot)
  iendtot=min(neqset,iendtot)
!
  y_lin(iprobe)=1.d0
!
  do ipoi=ibeg,iend
    do ieq=1,nbaleqs
      i=nbaleqs*(ipoi-1)+ieq
      params_lin(ieq,ipoi) = y_lin(i)
    end do
  end do
!
! Compute fluxes and internal sources:
!
  do ipoi=ibegb,iendb
    do ieq=1,nbaleqs
! radial derivatives of equilibrium parameters at cell boundaries:
      ddr_params(ieq,ipoi) &
         =sum(params_lin(ieq,ipbeg(ipoi):ipend(ipoi))*deriv_coef(:,ipoi))
      params_b_lin(ieq,ipoi)   &
         = sum(params_lin(ieq,ipbeg(ipoi):ipend(ipoi))*reint_coef(:,ipoi))
    enddo
  enddo
!
  Ercov_lin(ibegb:iendb)                                              &
           =sqg_bthet_overc(ibegb:iendb)*params_b_lin(2,ibegb:iendb)  & 
           +(params_b(4,ibegb:iendb)*ddr_params(1,ibegb:iendb)        &
           /params_b(1,ibegb:iendb)+ddr_params(4,ibegb:iendb))        &
           /(Z_i*e_charge)

!
  do ipoi=ibegb,iendb
!
! Thermodynamic forces for zero radial electric field:
    A_noE_1e=ddr_params(1,ipoi)/params_b(1,ipoi)  &
            -1.5d0*ddr_params(3,ipoi)/params_b(3,ipoi)
    A_noE_2e=ddr_params(3,ipoi)/params_b(3,ipoi)
    A_noE_1i=ddr_params(1,ipoi)/params_b(1,ipoi)  &
            -1.5d0*ddr_params(4,ipoi)/params_b(4,ipoi)
    A_noE_2i=ddr_params(4,ipoi)/params_b(4,ipoi)
!
!
! Thermodynamic forces for finite radial electric field:
    A_1e=A_noE_1e+Ercov_lin(ipoi)*e_charge/params_b(3,ipoi)
!ERROR    A_1i=A_noE_1e-Ercov_lin(ipoi)*e_charge*Z_i/params_b(4,ipoi)
    A_1i=A_noE_1i-Ercov_lin(ipoi)*e_charge*Z_i/params_b(4,ipoi) !<-FIXED
!
! particle flux densities:
    gamma_e=-(dae11(ipoi)*A_noE_1e+dae12(ipoi)*A_noE_2e)*params_b(1,ipoi)
    gamma_ql_e=-(dqle11(ipoi)*A_1e+dqle12(ipoi)*A_noE_2e)*params_b(1,ipoi)
    gamma_e=gamma_e+gamma_ql_e
    gamma_i=-(dai11(ipoi)*A_noE_1i+dai12(ipoi)*A_noE_2i)*params_b(1,ipoi)/Z_i
    gamma_ql_i=-(dqli11(ipoi)*A_1i+dqli12(ipoi)*A_noE_2i)*params_b(1,ipoi)/Z_i
    gamma_i=gamma_i+gamma_ql_i
!
! total particle flux:
    fluxes_dif(1,ipoi)=-Sb(ipoi)*ddr_params(1,ipoi)*(dae11(ipoi)              &
                  +dqle11(ipoi)*(1.d0+params_b(4,ipoi)/params_b(3,ipoi)/Z_i))
    fluxes_con(1,ipoi)=(Sb(ipoi)*gamma_e-fluxes_dif(1,ipoi))/params_b(1,ipoi)
!
! toroidal moment flux density divided by mass:
    dfluxvphi=-visca(ipoi)*ddr_params(2,ipoi)*params_b(1,ipoi)/Z_i*gpp_av(ipoi)
! total toroidal moment flux:
    fluxes_dif(2,ipoi)=Sb(ipoi)*dfluxvphi
    fluxes_con(2,ipoi)=0.d0
!
! electron heat flux density:
!colli    Q_e=-(dae12(ipoi)*A_noE_1e+dqle12(ipoi)*A_1e             &
    Q_e=-(dae12(ipoi)*A_noE_1e+dqle21(ipoi)*A_1e             &
       +  (dae22(ipoi)+dqle22(ipoi))*A_noE_2e)               &
       * params_b(1,ipoi)*params_b(3,ipoi)
!
! ion heat flux density:
!colli    Q_i=-(dai12(ipoi)*A_noE_1i+dqli12(ipoi)*A_1i             &
    Q_i=-(dai12(ipoi)*A_noE_1i+dqli21(ipoi)*A_1i             &
       +  (dai22(ipoi)+dni22(ipoi)+dqli22(ipoi))*A_noE_2i)   &
       * params_b(1,ipoi)/Z_i*params_b(4,ipoi)
!
! total heat fluxes:
!
    fluxes_dif(3,ipoi)=-Sb(ipoi)*(dae22(ipoi)+dqle22(ipoi))                 &
                      * params_b(1,ipoi)*ddr_params(3,ipoi)
    fluxes_con(3,ipoi)=(Sb(ipoi)*Q_e-fluxes_dif(3,ipoi))/params_b(3,ipoi)
!
    fluxes_dif(4,ipoi)=-Sb(ipoi)*(dai22(ipoi)+dni22(ipoi)+dqli22(ipoi)      &
!colli             - 2.5d0*dqli12(ipoi))*params_b(1,ipoi)/Z_i*ddr_params(4,ipoi)
             - 2.5d0*dqli21(ipoi))*params_b(1,ipoi)/Z_i*ddr_params(4,ipoi)
    fluxes_con(4,ipoi)=(Sb(ipoi)*Q_i-fluxes_dif(4,ipoi))/params_b(4,ipoi)
!
! Momentum source due to the polarization current:
    polforce(ipoi)=(gamma_e-Z_i*gamma_i)*e_charge*sqg_bthet_overc(ipoi) &
                  /(am*p_mass)
!
! Heat sources due to the radial QL drift in the equilibrium electric field:
    qlheat_e(ipoi)=-Ercov(ipoi)*gamma_ql_e*e_charge
    qlheat_i(ipoi)=Z_i*Ercov(ipoi)*gamma_ql_i*e_charge
!
  enddo
!
! Condition of zero flux at the inner boundary:
  fluxes_dif(:,1)=0.d0
  fluxes_con(:,1)=0.d0
  fluxes_con_nl(:,1)=0.d0
!
!
! Partial time derivatives of equilibrium parameters:
!
  do ipoi=ibeg,iend
!
! Flux divergence:
    do ieq=1,nbaleqs
      dot_params(ieq,ipoi)=-(fluxes_dif(ieq,ipoi+1)-fluxes_dif(ieq,ipoi))    &
                           /(Sc(ipoi)*(rb(ipoi+1)-rb(ipoi)))                 &
                           -(fluxes_con(ieq,ipoi+1)-fluxes_con(ieq,ipoi))    &
                           /(Sc(ipoi)*(rb(ipoi+1)-rb(ipoi)))*params(ieq,ipoi)
      convel=0.5d0*(fluxes_con_nl(ieq,ipoi+1)+fluxes_con_nl(ieq,ipoi))/Sc(ipoi)
! upstream convection:
      if(convel.gt.0.d0) then
!      if(convel.lt.0.d0) then
        dot_params(ieq,ipoi)=dot_params(ieq,ipoi)                             &
          -convel*(params_lin(ieq,ipoi+1)-params_lin(ieq,ipoi))/(rc(ipoi+1)-rc(ipoi))
      else
        if(ipoi.gt.1) then
          dot_params(ieq,ipoi)=dot_params(ieq,ipoi)                            &
            -convel*(params_lin(ieq,ipoi-1)-params_lin(ieq,ipoi))/(rc(ipoi-1)-rc(ipoi))
        else
          dot_params(ieq,ipoi)=dot_params(ieq,ipoi)                            &
            -convel*(params_b_lin(ieq,1)-params_lin(ieq,1))/(rb(1)-rc(1))
!print *,convel
        endif
      endif
    enddo

   !
   ! Add internal sources:
   ! Momentum:
       dot_params(2,ipoi) = dot_params(2,ipoi)     &
                          + 0.5d0*(polforce(ipoi)+polforce(ipoi+1))

   ! Heat into electrons:
       dot_params(3,ipoi) = dot_params(3, ipoi) &
                              + 0.5d0*(qlheat_e(ipoi) + qlheat_e(ipoi+1))

   ! Heat into ions:
       dot_params(4,ipoi) = dot_params(4,ipoi) &
                          + 0.5d0*(qlheat_i(ipoi) + qlheat_i(ipoi+1))


    ! Covert momentum time derivative to time derivative of the rotation frequency:
       dot_params(2,ipoi)=dot_params(2,ipoi)*Z_i/params(1,ipoi) &
                  *2.d0/(gpp_av(ipoi+1)+gpp_av(ipoi))

    ! Convert dot_params from d(nT_{e,i})/dt to d(T_{e,i})/dt:
       dot_params(3,ipoi) = (-params(3,ipoi)*dot_params(1,ipoi)     &
                          + dot_params(3,ipoi)/1.5d0)/params(1,ipoi)

       dot_params(4,ipoi) = (-params(4,ipoi)*dot_params(1,ipoi)     &
                          + dot_params(4,ipoi)/1.5d0)/params(1,ipoi)

  enddo

!
! RHS vector of ODE system
!
  do ipoi=ibeg,iend
    do ieq=1,nbaleqs
      i=nbaleqs*(ipoi-1)+ieq
      dy(i)=dot_params(ieq,ipoi)
    end do
  end do
!
  do i=ibegtot,iendtot
    if(dy(i).ne.0.d0) then
      k=k+1
      if(isw_rhs.eq.1) then
        irow(k)=i
        icol(k)=iprobe
        amat(k)=dy(i)
      endif
    endif
  enddo
!
! Clean:
!
  y_lin(iprobe)=0.d0
  do ipoi=ibeg,iend
    do ieq=1,nbaleqs
      i=nbaleqs*(ipoi-1)+ieq
      params_lin(ieq,ipoi) = y_lin(i)
    end do
  end do
  do ipoi=ibegb,iendb
    do ieq=1,nbaleqs
      ddr_params(ieq,ipoi)=0.d0
    enddo
  enddo
  Ercov_lin(ibegb:iendb)=0.d0
  fluxes_dif(:,ibegb:iendb)=0.d0
  fluxes_con(:,ibegb:iendb)=0.d0
  polforce(ibegb:iendb)=0.d0
  qlheat_e(ibegb:iendb)=0.d0
  qlheat_i(ibegb:iendb)=0.d0
  dot_params(:,ibeg:iend)=0.d0
  do ipoi=ibeg,iend
    do ieq=1,nbaleqs
      i=nbaleqs*(ipoi-1)+ieq
      dy(i)=dot_params(ieq,ipoi)
    end do
  end do
!$$$$
enddo
!
if(isw_rhs.eq.0) then
!  print *,'Number of non-zero elements = ',k,' out of: ',neqset**2
  nz=k
  nsize=neqset
else
!
  call rhs_balance_source(x,y,dy)
!
  rhsvec=dy
endif
!
  end subroutine rhs_balance
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine initialize_rhs(y,dy)

  use grid_mod, only :   neqset
  use matrix_mod

  implicit none
  double precision :: x
  double precision, dimension(neqset) :: y,dy

  x=0.d0
  isw_rhs=0
!
  call rhs_balance(x,y,dy)  
!
  isw_rhs=1
  if(allocated(amat)) deallocate(irow,icol,amat,rhsvec)
  allocate(irow(nz),icol(nz),amat(nz),rhsvec(nsize))
!
  end subroutine initialize_rhs 
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine rhs_balance_source(x,y,dy)

use grid_mod, only :   nbaleqs,neqset,iboutype,npoic,npoib &
                     , params,ddr_params,Sc,Sb,deriv_coef  &
                     , ipbeg,ipend,rb,params_b,reint_coef  &
                     , fluxes_dif,fluxes_con,rc,dot_params &
                     , dae11,dae12,dae22,dai11,dai12,dai22 &
                     , dni22,visca,gpp_av,dery_equisource  &
                     , dqle11,dqle12,dqle21,dqle22         &
                     , dqli11,dqli12,dqli21,dqli22         &
                     , sqg_bthet_overc,Ercov,polforce,qlheat_e,qlheat_i &
                     , params_lin,Ercov_lin,ddr_params_nl,fluxes_con_nl &
                     , init_params, params_b_lin 
use baseparam_mod, only : Z_i,e_charge,am,p_mass,c,btor
use control_mod, only : iwrite
use wave_code_data, only : q,Vth

implicit none

integer :: ipoi,ieq,i,npoi
double precision :: x,A_noE_1e,A_noE_2e,A_noE_1i,A_noE_2i,convel
double precision :: A_noE_1e_nl,A_noE_2e_nl,A_noE_1i_nl,A_noE_2i_nl
double precision :: gamma_e,gamma_i,dfluxvphi,Q_e,Q_i,A_1e,A_1i
double precision :: gamma_e_nl,Q_e_nl,Q_i_nl,A_1e_nl,A_1i_nl
double precision :: gamma_ql_e,gamma_ql_i
double precision, dimension(neqset) :: y,dy,y_lin
!
  if(iboutype.eq.1) then
    npoi=npoic-1
  else
    npoi=npoic
  endif
!npoi=npoic
!
! equilibrium parameters:
!
y_lin = 0.0d0
params_lin = 0.0d0

  params_lin=params

do ipoi=1,npoi
    do ieq=1,nbaleqs
      i=nbaleqs*(ipoi-1)+ieq
      params(ieq,ipoi)=y(i)
      params_lin(ieq,ipoi) = y_lin(i)
    end do
end do

!
! Interpolation:
!
  do ipoi=1,npoib
    do ieq=1,nbaleqs
! radial derivatives of equilibrium parameters at cell boundaries:
      ddr_params_nl(ieq,ipoi) &
         =sum(params(ieq,ipbeg(ipoi):ipend(ipoi))*deriv_coef(:,ipoi))
      ddr_params(ieq,ipoi) &
         =sum(params_lin(ieq,ipbeg(ipoi):ipend(ipoi))*deriv_coef(:,ipoi))
! equilibrium parameters at cell boundaries:
      params_b(ieq,ipoi)   &
         =sum(params(ieq,ipbeg(ipoi):ipend(ipoi))*reint_coef(:,ipoi))
      params_b_lin(ieq,ipoi)   &
         =sum(params_lin(ieq,ipbeg(ipoi):ipend(ipoi))*reint_coef(:,ipoi))
    enddo
  enddo
!
! Compute radial electric field:
!
!   Ercov=sqg_bthet_overc*params_b(2,:)                                 & !OLD
   Ercov=sqg_bthet_overc*(params_b(2,:)-Vth*q/rb)                           &
       +(params_b(4,:)*ddr_params_nl(1,:)/params_b(1,:)+ddr_params_nl(4,:)) &
       /(Z_i*e_charge)

   Ercov_lin=sqg_bthet_overc*params_b_lin(2,:)                              &
        +(params_b(4,:)*ddr_params(1,:)/params_b(1,:)+ddr_params(4,:))      &
       /(Z_i*e_charge)

! Compute diffusion coefficient matrices:
!
  call calc_dequi
!
! Compute fluxes and internal sources:
!
!open(8765,file='torqueden.dat')
  do ipoi=1,npoib
!
! Thermodynamic forces for zero radial electric field:
    A_noE_1e=ddr_params(1,ipoi)/params_b(1,ipoi)  &
            -1.5d0*ddr_params(3,ipoi)/params_b(3,ipoi)
    A_noE_2e=ddr_params(3,ipoi)/params_b(3,ipoi)
    A_noE_1i=ddr_params(1,ipoi)/params_b(1,ipoi)  &
            -1.5d0*ddr_params(4,ipoi)/params_b(4,ipoi)
    A_noE_2i=ddr_params(4,ipoi)/params_b(4,ipoi)
!
    A_noE_1e_nl=ddr_params_nl(1,ipoi)/params_b(1,ipoi)  &
            -1.5d0*ddr_params_nl(3,ipoi)/params_b(3,ipoi)
    A_noE_2e_nl=ddr_params_nl(3,ipoi)/params_b(3,ipoi)
    A_noE_1i_nl=ddr_params_nl(1,ipoi)/params_b(1,ipoi)  &
            -1.5d0*ddr_params_nl(4,ipoi)/params_b(4,ipoi)
    A_noE_2i_nl=ddr_params_nl(4,ipoi)/params_b(4,ipoi)
!
! Thermodynamic forces for finite radial electric field:
    A_1e=A_noE_1e+Ercov_lin(ipoi)*e_charge/params_b(3,ipoi)
!ERROR    A_1i=A_noE_1e-Ercov_lin(ipoi)*e_charge*Z_i/params_b(4,ipoi)
    A_1i=A_noE_1i-Ercov_lin(ipoi)*e_charge*Z_i/params_b(4,ipoi) !<-FIXED
!
    A_1e_nl=A_noE_1e_nl+Ercov(ipoi)*e_charge/params_b(3,ipoi)
!ERROR    A_1i_nl=A_noE_1e_nl-Ercov(ipoi)*e_charge*Z_i/params_b(4,ipoi)
    A_1i_nl=A_noE_1i_nl-Ercov(ipoi)*e_charge*Z_i/params_b(4,ipoi) !<-FIXED

!
! particle flux densities:
    gamma_e=-(dae11(ipoi)*A_noE_1e+dae12(ipoi)*A_noE_2e)*params_b(1,ipoi)
    gamma_ql_e=-(dqle11(ipoi)*A_1e+dqle12(ipoi)*A_noE_2e)*params_b(1,ipoi)
    gamma_e=gamma_e+gamma_ql_e
    gamma_i=-(dai11(ipoi)*A_noE_1i+dai12(ipoi)*A_noE_2i)*params_b(1,ipoi)/Z_i
    gamma_ql_i=-(dqli11(ipoi)*A_1i+dqli12(ipoi)*A_noE_2i)*params_b(1,ipoi)/Z_i
    gamma_i=gamma_i+gamma_ql_i
!
    gamma_e_nl=-(dae11(ipoi)*A_noE_1e_nl+dae12(ipoi)*A_noE_2e_nl      &
              +  dqle11(ipoi)*A_1e_nl+dqle12(ipoi)*A_noE_2e_nl        &
                 )*params_b(1,ipoi)


! total particle flux:
    fluxes_dif(1,ipoi)=-Sb(ipoi)*ddr_params(1,ipoi)*(dae11(ipoi)              &
                  +dqle11(ipoi)*(1.d0+params_b(4,ipoi)/params_b(3,ipoi)/Z_i))
    fluxes_con(1,ipoi)=(Sb(ipoi)*gamma_e-fluxes_dif(1,ipoi))/params_b(1,ipoi)
    fluxes_con_nl(1,ipoi)=(Sb(ipoi)*gamma_e_nl - &
     (-Sb(ipoi)*ddr_params_nl(1,ipoi)*(dae11(ipoi)              &
      +dqle11(ipoi)*(1.d0+params_b(4,ipoi)/params_b(3,ipoi)/Z_i))))/params_b(1,ipoi)
!
! toroidal moment flux density divided by mass:
    dfluxvphi=-visca(ipoi)*ddr_params(2,ipoi)*params_b(1,ipoi)/Z_i*gpp_av(ipoi)
! total toroidal moment flux:
    fluxes_dif(2,ipoi)=Sb(ipoi)*dfluxvphi
    fluxes_con(2,ipoi)=0.d0
    fluxes_con_nl(2,ipoi)=0.d0
!
! electron heat flux density:
!colli    Q_e=-(dae12(ipoi)*A_noE_1e+dqle12(ipoi)*A_1e             &
    Q_e=-(dae12(ipoi)*A_noE_1e+dqle21(ipoi)*A_1e             &
       +  (dae22(ipoi)+dqle22(ipoi))*A_noE_2e)               &
       * params_b(1,ipoi)*params_b(3,ipoi)
!
!colli    Q_e_nl=-(dae12(ipoi)*A_noE_1e_nl+dqle12(ipoi)*A_1e_nl    &
    Q_e_nl=-(dae12(ipoi)*A_noE_1e_nl+dqle21(ipoi)*A_1e_nl    &
       +  (dae22(ipoi)+dqle22(ipoi))*A_noE_2e_nl)            &
       * params_b(1,ipoi)*params_b(3,ipoi)
!
! ion heat flux density:
!colli    Q_i=-(dai12(ipoi)*A_noE_1i+dqli12(ipoi)*A_1i             &
    Q_i=-(dai12(ipoi)*A_noE_1i+dqli21(ipoi)*A_1i             &
       +  (dai22(ipoi)+dni22(ipoi)+dqli22(ipoi))*A_noE_2i)   &
       * params_b(1,ipoi)/Z_i*params_b(4,ipoi)
!
!colli    Q_i_nl=-(dai12(ipoi)*A_noE_1i_nl+dqli12(ipoi)*A_1i_nl    &
    Q_i_nl=-(dai12(ipoi)*A_noE_1i_nl+dqli21(ipoi)*A_1i_nl    &
       +  (dai22(ipoi)+dni22(ipoi)+dqli22(ipoi))*A_noE_2i_nl)&
       * params_b(1,ipoi)/Z_i*params_b(4,ipoi)
!
! total heat fluxes:
!
    fluxes_dif(3,ipoi)=-Sb(ipoi)*(dae22(ipoi)+dqle22(ipoi))                 &
                      * params_b(1,ipoi)*ddr_params(3,ipoi)
    fluxes_con(3,ipoi)=(Sb(ipoi)*Q_e-fluxes_dif(3,ipoi))/params_b(3,ipoi)
    fluxes_con_nl(3,ipoi)=(Sb(ipoi)*Q_e_nl - &
    (-Sb(ipoi)*(dae22(ipoi)+dqle22(ipoi))*params_b(1,ipoi)*ddr_params_nl(3,ipoi))) &
    /params_b(3,ipoi)
!
    fluxes_dif(4,ipoi)=-Sb(ipoi)*(dai22(ipoi)+dni22(ipoi)+dqli22(ipoi)      &
!colli             - 2.5d0*dqli12(ipoi))*params_b(1,ipoi)/Z_i*ddr_params(4,ipoi)
             - 2.5d0*dqli21(ipoi))*params_b(1,ipoi)/Z_i*ddr_params(4,ipoi)
    fluxes_con(4,ipoi)=(Sb(ipoi)*Q_i-fluxes_dif(4,ipoi))/params_b(4,ipoi)
    fluxes_con_nl(4,ipoi)=(Sb(ipoi)*Q_i_nl - &
!colli    (-Sb(ipoi)*(dai22(ipoi)+dni22(ipoi)+dqli22(ipoi)-2.5d0*dqli12(ipoi)) &
     (-Sb(ipoi)*(dai22(ipoi)+dni22(ipoi)+dqli22(ipoi)-2.5d0*dqli21(ipoi)) &
      *params_b(1,ipoi)/Z_i*ddr_params_nl(4,ipoi))) &
     /params_b(4,ipoi)
!
! Momentum source due to the polarization current:
    polforce(ipoi)=(gamma_e-Z_i*gamma_i)*e_charge*sqg_bthet_overc(ipoi) &
                  /(am*p_mass)
!
! Heat sources due to the radial QL drift in the equilibrium electric field:
    qlheat_e(ipoi)=-Ercov(ipoi)*gamma_ql_e*e_charge
    qlheat_i(ipoi)=Z_i*Ercov(ipoi)*gamma_ql_i*e_charge
!
!write(8765,*) rb(ipoi), -e_charge*sqg_bthet_overc(ipoi)            &
!*(dqle11(ipoi)*A_1e_nl+dqle12(ipoi)*A_noE_2e_nl)*params_b(1,ipoi), &
!Z_i*e_charge*sqg_bthet_overc(ipoi)                                 &
!*(dqli11(ipoi)*A_1i_nl+dqli12(ipoi)*A_noE_2i_nl)*params_b(1,ipoi), &
!A_1e_nl,A_noE_2e_nl,dqle11(ipoi),dqle12(ipoi),sqg_bthet_overc(ipoi),&
!A_1i_nl,A_noE_2i_nl,dqli11(ipoi),dqli12(ipoi)
  enddo
!close(8765)
!
! Condition of zero flux at the inner boundary:
  fluxes_dif(:,1)=0.d0
  fluxes_con(:,1)=0.d0
  fluxes_con_nl(:,1)=0.d0
!
!
! Partial time derivatives of equilibrium parameters:
!
  do ipoi=1,npoi
!
! Flux divergence:
    do ieq=1,nbaleqs
      dot_params(ieq,ipoi)=-(fluxes_dif(ieq,ipoi+1)-fluxes_dif(ieq,ipoi))    &
                           /(Sc(ipoi)*(rb(ipoi+1)-rb(ipoi)))                 &
                           -(fluxes_con(ieq,ipoi+1)-fluxes_con(ieq,ipoi))    &
                           /(Sc(ipoi)*(rb(ipoi+1)-rb(ipoi)))*params(ieq,ipoi)
      convel=0.5d0*(fluxes_con_nl(ieq,ipoi+1)+fluxes_con_nl(ieq,ipoi))/Sc(ipoi)
! upstream convection:
!      if(convel.lt.0.d0) then
      if(convel.gt.0.d0) then
        dot_params(ieq,ipoi)=dot_params(ieq,ipoi)                             &
          -convel*(params_lin(ieq,ipoi+1)-params_lin(ieq,ipoi))/(rc(ipoi+1)-rc(ipoi))
      else
        if(ipoi.gt.1) then
          dot_params(ieq,ipoi)=dot_params(ieq,ipoi)                            &
            -convel*(params_lin(ieq,ipoi-1)-params_lin(ieq,ipoi))/(rc(ipoi-1)-rc(ipoi))
        else
             dot_params(ieq,ipoi)=dot_params(ieq,ipoi)                    !    &
!               -convel*(params_b(ieq,1))/(rb(1)-rc(1))
        endif
      endif
    enddo

   !
   ! Add internal sources:
   ! Momentum:
       dot_params(2,ipoi) = dot_params(2,ipoi)     &
                          + 0.5d0*(polforce(ipoi)+polforce(ipoi+1))
   if(iwrite.eq.1) write (11,*) rc(ipoi),dot_params(2,ipoi)

   ! Heat into electrons:
       dot_params(3,ipoi) = dot_params(3, ipoi) &
                              + 0.5d0*(qlheat_e(ipoi) + qlheat_e(ipoi+1))

   ! Heat into ions:
       dot_params(4,ipoi) = dot_params(4,ipoi) &
                          + 0.5d0*(qlheat_i(ipoi) + qlheat_i(ipoi+1))


! Convert momentum time derivative to time derivative of the rotation frequency:
       dot_params(2,ipoi)=dot_params(2,ipoi)*Z_i/params(1,ipoi) &
                  *2.d0/(gpp_av(ipoi+1)+gpp_av(ipoi))

! Convert dot_params from d(nT_{e,i})/dt to d(T_{e,i})/dt:
       dot_params(3,ipoi) = (-params(3,ipoi)*dot_params(1,ipoi)     &
                          + dot_params(3,ipoi)/1.5d0)/params(1,ipoi)

       dot_params(4,ipoi) = (-params(4,ipoi)*dot_params(1,ipoi)     &
                          + dot_params(4,ipoi)/1.5d0)/params(1,ipoi)

if(iwrite.eq.1) write (10,*) rc(ipoi),dot_params(2,ipoi)             &
                           , dot_params(3,ipoi),dot_params(4,ipoi)   &
                           , 0.5d0*(polforce(ipoi)+polforce(ipoi+1)) &
                           , 0.5d0*(qlheat_e(ipoi)+qlheat_e(ipoi+1)) &
                           , 0.5d0*(qlheat_i(ipoi)+qlheat_i(ipoi+1))
enddo

if(iwrite.eq.1) then
   close(10)
   close(21)
end if

!
! RHS vector of ODE system
!
  do ipoi=1,npoi
    do ieq=1,nbaleqs
      i=nbaleqs*(ipoi-1)+ieq
      dy(i)=dot_params(ieq,ipoi)
    end do
  end do
!
  dy=dy+dery_equisource
!
  end subroutine rhs_balance_source
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine calc_dequi
!
  use grid_mod, only : npoib,params_b,dae11,dae12,dae22,dai11,dai12,dai22 &
                     , dni22,visca,rb,cneo
  use baseparam_mod, only : dperp, rsepar
!
  implicit none
!
  integer :: ipoi
  double precision :: rnorm,weight
!
  do ipoi=1,npoib
!    rnorm=rb(ipoi)/rb(npoib)
    rnorm=rb(ipoi)/rsepar
    dae11(ipoi)=dperp*(1.d0-(1.d0-0.2d0)*rnorm**3)
    call localizer(-1.d0,rsepar,rsepar+0.5d0,rb(ipoi),weight)
    dae11(ipoi)=dae11(ipoi)*(1.d0-weight)+weight*1d6
  enddo
!  dae12=1.5d0*dae11
  dae12=1.499999d0*dae11
  dae22=3.75d0*dae11
  dai11=dae11
  dai12=dae12
  dai22=dae22
  visca=dae11
!  visca=100.d0*dae11
!
  dni22=cneo*params_b(1,:)/sqrt(abs(params_b(4,:)))
!
  end subroutine calc_dequi
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine get_dql

  use grid_mod, only : nbaleqs,neqset,iboutype,npoic,npoib      &
                     , params,ddr_params,deriv_coef             &
                     , ipbeg,ipend,rb,params_b,reint_coef       &
                     , rc,sqg_bthet_overc,Ercov                 &
                     , ddr_params_nl,y,mwind                    &
                     , dqle11,dqle12,dqle21,dqle22              &
                     , dqli11,dqli12,dqli21,dqli22              &
                     , de11,de12,de21,de22,di11,di12,di21,di22
  use baseparam_mod, only : Z_i,e_charge,am,p_mass,c,btor,e_mass,ev,rtor
  use control_mod, only : irf,write_formfactors
  use wave_code_data
  use mpi
!DIAG:
  use diag_mod, only : write_diag,iunit_diag,write_diag_b,iunit_diag_b,i_mn_loop
!END DIAG
!
  implicit none
!
  integer :: ierror,np_num,irank,modpernode,imin,imax;
  integer :: ipoi,ieq,i,npoi,i_mn,ierr,mwind_save
  double precision, dimension(:), allocatable :: dummy
  double complex, dimension(npoib) :: amn_psi, amn_theta, amn_theta_cyl
!
  double precision, dimension(npoib) :: spec_weight
!double precision :: r_res, D,rnorm
  double precision, dimension(npoib) :: vT_e, vT_i, nu_e, nu_i
!
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: dqle11_loc
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: dqle12_loc
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: dqle21_loc
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: dqle22_loc
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: dqli11_loc
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: dqli12_loc
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: dqli21_loc
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: dqli22_loc
  double complex, dimension(:), allocatable :: formfactor
!
  allocate(dqle11_loc(npoib))
  allocate(dqle12_loc(npoib))
  allocate(dqle21_loc(npoib))
  allocate(dqle22_loc(npoib))
  allocate(dqli11_loc(npoib))
  allocate(dqli12_loc(npoib))
  allocate(dqli21_loc(npoib))
  allocate(dqli22_loc(npoib))
  allocate(formfactor(npoib))
!
  dqle11_loc = 0.0d0
  dqle12_loc = 0.0d0
  dqle21_loc = 0.0d0
  dqle22_loc = 0.0d0
  dqli11_loc = 0.0d0
  dqli12_loc = 0.0d0
  dqli21_loc = 0.0d0
  dqli22_loc = 0.0d0
!
  if(irf.eq.0) then
    return
  elseif(irf.eq.2) then
    dqle11 = 0.0d0
    dqle12 = 0.0d0
    dqle21 = 0.0d0
    dqle22 = 0.0d0
    dqli11 = 0.0d0
    dqli12 = 0.0d0
    dqli21 = 0.0d0
    dqli22 = 0.0e0
    return
  end if
!
  if(iboutype.eq.1) then
    npoi=npoic-1
  else
    npoi=npoic
  endif
!
! equilibrium parameters:
!
!  do ipoi=1,npoi
!    do ieq=1,nbaleqs
!      i=nbaleqs*(ipoi-1)+ieq
!      params(ieq,ipoi)=y(i)
!    end do
!  end do
!
! Interpolation:
!
  do ipoi=1,npoib
    do ieq=1,nbaleqs
! radial derivatives of equilibrium parameters at cell boundaries:
      ddr_params_nl(ieq,ipoi) &
         =sum(params(ieq,ipbeg(ipoi):ipend(ipoi))*deriv_coef(:,ipoi))
! equilibrium parameters at cell boundaries:
      params_b(ieq,ipoi)   &
         =sum(params(ieq,ipbeg(ipoi):ipend(ipoi))*reint_coef(:,ipoi))
    enddo
  enddo
!
! Smooth input for KILCA
!
!  if(.true.) then
  if(.false.) then
    allocate(dummy(npoib))
    do ieq=1,nbaleqs
      call smooth_array_gauss(npoib,mwind,ddr_params_nl(ieq,:),dummy)
      ddr_params_nl(ieq,:)=dummy
      call smooth_array_gauss(npoib,mwind,params_b(ieq,:),dummy)
      params_b(ieq,:)=dummy
    enddo
    deallocate(dummy)
  endif
!
! End smooth input for KILCA
!
! Compute radial electric field:
!
!   Ercov=sqg_bthet_overc*params_b(2,:)                                 & !OLD
   Ercov=sqg_bthet_overc*(params_b(2,:)-Vth*q/rb)                           &
       +(params_b(4,:)*ddr_params_nl(1,:)/params_b(1,:)+ddr_params_nl(4,:)) &
       /(Z_i*e_charge)
!
call MPI_Comm_rank(MPI_COMM_WORLD,irank,ierror);
if (irank .eq. 0 ) then
  if(write_diag_b) then
    do ipoi=1,npoib
      write (iunit_diag_b,*) rb(ipoi),params_b(1:4,ipoi)
    end do
  endif
end if

!
! Compute diffusion coefficient matrices:

call MPI_Comm_size(MPI_COMM_WORLD,np_num,ierror);
call MPI_Comm_rank(MPI_COMM_WORLD,irank,ierror);

!sum over modes:
if (np_num .gt. dim_mn) then
    print *,' '
    print *, 'Number of processes', np_num, 'is larger than number of modes', dim_mn
    call MPI_finalize(ierror)
    stop
end if

modpernode = ceiling(float(dim_mn)/float(np_num));

imin = modpernode*irank+1;

imax = min(dim_mn, modpernode*(irank+1));

if (irf.eq.1) call update_background_files(path2profs);

if (irf.eq.1) call get_wave_code_data(imin, imax);

if (irf.eq.1) call get_background_magnetic_fields_from_wave_code(flre_cd_ptr(imin), dim_r, r, B0t, B0z, B0);

if (irf.eq.1) call get_collision_frequences_from_wave_code(flre_cd_ptr(imin), dim_r, r, nui, nue);

!
!  nu_e=15.4d-6*params_b(1,:)/sqrt(params_b(3,:)/ev)**3            &
!      *(23.d0-0.5d0*log(params_b(1,:)/(params_b(3,:)/ev)**3))
!  nu_i=1.d-7*params_b(1,:)/sqrt(params_b(4,:)/ev)**3              &
!      *(23.d0-0.5d0*log(2.d0*params_b(1,:)/(params_b(4,:)/ev)**3))
!
  nu_e=nue
  nu_i=nui
!
!initialization before summing up over modes:
  dqle11 = 0.0d0
  dqle12 = 0.0d0
  dqle21 = 0.0d0
  dqle22 = 0.0d0
  dqli11 = 0.0d0
  dqli12 = 0.0d0
  dqli21 = 0.0d0
  dqli22 = 0.0e0

!sum over modes:
  do i_mn = imin,imax
!
    call get_wave_vectors_from_wave_code(flre_cd_ptr(i_mn),dim_r,r,              &
                                         m_vals(i_mn),n_vals(i_mn),ks,kp)
!
    call get_wave_fields_from_wave_code(flre_cd_ptr(i_mn),dim_r,r,               &
                  m_vals(i_mn),n_vals(i_mn),Er,Es,Ep,Et,Ez,Br,Bs,Bp,Bt,Bz)
!
    om_E = ks*c*dPhi0/B0;
!
    vT_e = sqrt(params_b(3,:)/e_mass)
    vT_i = sqrt(params_b(4,:)/p_mass/am)
!
!    if(.true.) then
    i_mn_loop=i_mn
    if(.false.) then
      call calc_transport_coeffs_collisionless(npoib, vT_e, de11, de12, de22)
      de21=de12
      call calc_transport_coeffs_collisionless(npoib, vT_i, di11, di12, di22)
      di21=di12
    else
      if(.true.) then
!      if(.false.) then
        call calc_transport_coeffs_ornuhl(npoib,vT_e,nu_e,de11,de12,de21,de22)
        call calc_transport_coeffs_ornuhl(npoib,vT_i,nu_i,di11,di12,di21,di22)
      else
        call calc_transport_coeffs_ornuhl_drift(1,npoib,de11,de12,de21,de22)
        call calc_transport_coeffs_ornuhl_drift(2,npoib,di11,di12,di21,di22)
      endif
    end if
!
    call get_wave_fields_from_wave_code(vac_cd_ptr(i_mn),dim_r,r,               &
                  m_vals(i_mn),n_vals(i_mn),Bz,Bz,Bz,Bz,Bz, Br, Bz,Bz,Bz,Bz)
!
    formfactor=(1.d0,0.d0)/Br
!
    do ipoi = 1, npoib
      call amn_of_r(+m_vals(i_mn),n_vals(i_mn),r(ipoi),                     &
                    amn_psi(ipoi),amn_theta(ipoi),ierr)
      if(ierr.ne.0) then
!            print *,'amn_of_r error ',ierr
        amn_psi(ipoi)=amn_psi(ipoi-1)
        amn_theta(ipoi)=amn_theta(ipoi-1)
!            stop
      endif
    end do
    amn_theta_cyl = r*rtor/n_vals(i_mn)*Br
    spec_weight = 2.0d0*(abs(amn_theta)**2/abs(amn_theta_cyl)**2)
    
    !diagnostics
    open(855,file='amn_theta.dat')
        do ipoi=1,npoic
          write (855,*) r(ipoi),abs(amn_theta(ipoi)),abs(amn_theta_cyl(ipoi))
        end do
    close(855)
    
if (irank .eq. 0 ) then
!DIAG:
  if(write_diag) then
    open(7000,file='Brvac.dat')
    do ipoi = 1, npoib
       write(7000,*) r(ipoi),abs(Br(ipoi)) 
    enddo
    close(7000)
  endif
!END DIAG
endif
!
    call get_wave_fields_from_wave_code (flre_cd_ptr(i_mn),dim_r,r,              &
                  m_vals(i_mn),n_vals(i_mn),Bz,Bz,Bz,Bz,Bz, Br, Bz,Bz,Bz,Bz)
!
!
    formfactor=Br*formfactor
    if(write_formfactors) then
      do ipoi = 1, npoib
        write(10000+n_vals(i_mn)*1000+m_vals(i_mn),*) r(ipoi),abs(formfactor(ipoi))
      enddo
      close(10000+n_vals(i_mn)*1000+m_vals(i_mn))
    endif

!  spec_weight=1.0d0  ! for DIII-D
!
    dqle11_loc = dqle11_loc + de11*spec_weight
    dqle12_loc = dqle12_loc + de12*spec_weight
    dqle21_loc = dqle21_loc + de21*spec_weight
    dqle22_loc = dqle22_loc + de22*spec_weight
    dqli11_loc = dqli11_loc + di11*spec_weight
    dqli12_loc = dqli12_loc + di12*spec_weight
    dqli21_loc = dqli21_loc + di21*spec_weight
    dqli22_loc = dqli22_loc + di22*spec_weight
!
    call get_current_densities_from_wave_code (flre_cd_ptr(i_mn),dim_r,r,       &
                        m_vals(i_mn),n_vals(i_mn),Jri,Jsi,Jpi,Jre,Jse,Jpe)
!   
  end do
!

call MPI_Allreduce(dqle11_loc,dqle11,npoib,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror);
call MPI_Allreduce(dqle12_loc,dqle12,npoib,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror);
call MPI_Allreduce(dqle21_loc,dqle21,npoib,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror);
call MPI_Allreduce(dqle22_loc,dqle22,npoib,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror);
call MPI_Allreduce(dqli11_loc,dqli11,npoib,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror);
call MPI_Allreduce(dqli12_loc,dqli12,npoib,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror);
call MPI_Allreduce(dqli21_loc,dqli21,npoib,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror);
call MPI_Allreduce(dqli22_loc,dqli22,npoib,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror);

call MPI_Barrier(MPI_COMM_WORLD,ierror);

deallocate(dqle11_loc);
deallocate(dqle12_loc);
deallocate(dqle21_loc);
deallocate(dqle22_loc);
deallocate(dqli11_loc);
deallocate(dqli12_loc);
deallocate(dqli21_loc);
deallocate(dqli22_loc);
deallocate(formfactor)

  call calc_parallel_current_directly
  call calc_ion_parallel_current_directly
!

!if(.false.) then
if(.true.) then
mwind_save=mwind
mwind=30
allocate(dummy(npoib))
call smooth_array_gauss (npoib,mwind,dqle11,dummy)
dqle11=dummy
call smooth_array_gauss (npoib,mwind,dqle12,dummy)
dqle12=dummy
call smooth_array_gauss (npoib,mwind,dqle21,dummy)
dqle21=dummy
call smooth_array_gauss (npoib,mwind,dqle22,dummy)
dqle22=dummy
!
mwind=30
call smooth_array_gauss (npoib,mwind,dqli11,dummy)
dqli11=dummy
call smooth_array_gauss (npoib,mwind,dqli12,dummy)
dqli12=dummy
call smooth_array_gauss (npoib,mwind,dqli21,dummy)
dqli21=dummy
call smooth_array_gauss (npoib,mwind,dqli22,dummy)
dqli22=dummy
mwind=mwind_save
deallocate(dummy)
endif

if (irank .eq. 0 ) then
!
!DIAG:
  if(write_diag) then
    do ipoi = 1, npoib
!      write(iunit_diag,*) r(ipoi),dqle11(ipoi),dqli11(ipoi) &
!                            ,real(Br(ipoi)),imag(Br(ipoi))  &
!                            ,real(Ep(ipoi)),imag(Ep(ipoi))  &
!                            ,real(Er(ipoi)),imag(Er(ipoi))  &
!                            ,om_E(ipoi) &
!                            ,abs(c*ks(ipoi)*Ep(ipoi) - om_E(ipoi)*Br(ipoi))**2
       write(iunit_diag,*) r(ipoi),dqle11(ipoi),dqle12(ipoi) &
                           ,dqle22(ipoi),dqli11(ipoi)        &
                           ,dqli12(ipoi),dqli22(ipoi)        &
                           ,abs(Br(ipoi))                    &
                           ,abs(Br(ipoi)-c*kp(ipoi)*Es(ipoi)/om_E(ipoi)) &
                           ,abs(Br(ipoi)-c*ks(ipoi)*Ep(ipoi)/om_E(ipoi)) &
                           ,abs(Jpe(ipoi)),abs(Jpi(ipoi)) &
                           ,abs(Jpe(ipoi)+Jpi(ipoi)) 
    enddo
    close(iunit_diag)
  endif
!END DIAG
end if

  end subroutine get_dql

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine calc_transport_coeffs_collisionless (dim, vT, D_11, D_12, D_22)

use baseparam_mod
!use wave_code_data, only: om_E, kp, ks, Ep, Br
use wave_code_data, only: om_E, kp, ks, Ep, Br, Bs, Bp, r, B0, Te

integer, intent(in) :: dim
real(8), dimension(dim), intent(in) :: vT
real(8), dimension(dim), intent(out) :: D_11, D_12, D_22
real(8), dimension(dim) :: Z
!complex(8), dimension(dim) :: field_fac
double precision, dimension(dim) :: field_fac

Z = om_E/sqrt(2.0)/kp/vT
!Z = om_E/sqrt(2.0)/sqrt(kp**2+1.2e-5**2)/vT

!field_fac = c*ks*Ep - om_E*Br
field_fac = abs(c*ks*Ep - om_E*Br)**2

!D_11 = sqrt(pi)*vT**2/btor**2*(abs(Z/om_E))**3*exp(-Z**2)*(abs(field_fac))**2
D_11 = sqrt(pi)*vT**2/btor**2*(abs(Z/om_E))**3*exp(-Z**2)*field_fac
D_12 = (1.0 + Z**2)*D_11
D_22 = (1.0 + (1.0 + Z**2)**2)*D_11

!do i=1,dim
!write(333,*) r(i),sqrt(field_fac(i))/abs(c*ks(i)), &
!abs(kp(i)*Bp(i)*Te(i)*ev/e_charge/B0(i))
!write(333,*) r(i),abs(Br(i)),abs(kp(i)*Bs(i))
!enddo
!stop
!do i=1,dim
!write(335,*) r(i),D_11(i),D_12(i),D_12(i),D_22(i)
!enddo
!stop

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calc_transport_coeffs_ornuhl (dim,vT,nu,D_11,D_12,D_21,D_22)

use baseparam_mod
use wave_code_data, only: om_E, kp, Es, Br, B0, r, ks, Ep
use diag_mod, only : i_mn_loop
use grid_mod, only : r_resonant,gg_width,rb

integer, parameter :: mnmax=3
integer, intent(in) :: dim
real(8), dimension(dim), intent(in) :: vT,nu
real(8), dimension(dim), intent(out) :: D_11, D_12, D_21, D_22
double precision, dimension(:),     allocatable :: x1,x2,comfac,d_12a
double precision, dimension(:),     allocatable :: epm2,brm2,epbr_re,epbr_im
double complex,   dimension(:,:,:), allocatable :: symbI

allocate(comfac(dim),d_12a(dim),epm2(dim),brm2(dim),epbr_re(dim),epbr_im(dim))
allocate(x1(dim),x2(dim),symbI(0:mnmax,0:mnmax,dim))
!
!    if  Br=c*kp*Es/om_E diffusion tensor iz zero

comfac=0.5d0/(nu*B0**2)
epm2=c**2*abs(Es)**2
brm2=vT**2*abs(Br)**2
epbr_re=2.d0*c*vT*real(conjg(Es)*Br)
epbr_im=2.d0*c*vT*dimag(conjg(Es)*Br)
!epm2=0.0d0 !c**2*abs(Es)**2
!brm2=1.0d0 !vT**2*abs(Br)**2
!epbr_re=0.0d0 !2.d0*c*vT*real(conjg(Es)*Br)
!epbr_im=0.0d0 !2.d0*c*vT*dimag(conjg(Es)*Br)

x1=kp*vT/nu
x2=-om_E/nu
symbI=(0.d0,0.d0)
do i=1,dim
  if(rb(i).lt.r_resonant(i_mn_loop)-2.d0*gg_width) cycle
  if(rb(i).gt.r_resonant(i_mn_loop)+2.d0*gg_width) cycle
  call getIfunc(x1(i),x2(i),symbI(:,:,i))
enddo

D_11 = comfac*(epm2   *real(symbI(0,0,:))                                     &
     +         epbr_re*real(symbI(1,0,:))                                     &
     +         brm2   *real(symbI(1,1,:)))
D_12 = comfac*(epm2   *real(symbI(0,0,:)+0.5d0*symbI(2,0,:))                  &
     +         epbr_re*real(symbI(1,0,:)+0.25d0*(symbI(3,0,:)+symbI(2,1,:)))  &
     +         brm2   *real(symbI(1,1,:)+0.5d0*symbI(3,1,:)))
D_21 = D_12
D_22 = comfac*(epm2   *real(2.d0*symbI(0,0,:)+symbI(2,0,:)                    &
     +                      0.25d0*symbI(2,2,:))                              &
     +         epbr_re*real(2.d0*symbI(1,0,:)                                 &
     +                      0.5d0*(symbI(3,0,:)+symbI(2,1,:))                 &
     +                      0.25d0*symbI(3,2,:))                              &
     +         brm2   *real(2.d0*symbI(1,1,:)+symbI(3,1,:)                    &
     +                      0.25d0*symbI(3,3,:)))

D_12a=comfac*epbr_im*0.25d0*dimag(symbI(2,1,:)-symbI(3,0,:))

D_12 = D_12 + D_12a
D_21 = D_21 - D_12a

!D_11=comfac*epm2*real(symbI(0,0,:))
!D_12=comfac*epbr_re*real(symbI(1,0,:))
!D_22=comfac*brm2*real(symbI(1,1,:))
!do i=1,dim
!write(7000,*) r(i),D_11(i),D_12(i),D_22(i),D_11(i)+D_12(i)+D_22(i)
!enddo
!D_11=comfac*real(symbI(0,0,:))
!D_12=comfac*real(symbI(1,0,:))
!D_22=comfac*real(symbI(1,1,:))
!do i=1,dim
!write(7001,*) r(i),D_11(i),D_12(i),D_22(i),x1(i),x2(i)
!enddo
!stop

deallocate(x1,x2,symbI)
deallocate(comfac,d_12a,epm2,brm2,epbr_re,epbr_im)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getIfunc(x1,x2,symbI)
  integer, parameter :: mnmax=3
  integer :: m,n
  double precision :: x1,x2,z
  double complex :: denom
  double complex, dimension(0:mnmax,0:mnmax) :: symbI,Imn
!
!  if(.true.) then
  if(.false.) then
! collisionless case:
    symbI=(0.d0,0.d0)
    z=x2/(sqrt(2.d0)*x1)
!
    symbI(0,0)=sqrt(2.d0)*exp(-z**2)/abs(x1)
    symbI(1,0)=symbI(0,0)*x2/x1
    symbI(1,1)=symbI(1,0)*x2/x1
!
    symbI(2,0)=symbI(0,0)*(x2/x1)**2
    symbI(2,1)=symbI(1,0)*(x2/x1)**2
    symbI(3,0)=symbI(2,1)
    symbI(3,1)=symbI(1,1)*(x2/x1)**2
!
    symbI(2,2)=symbI(2,0)*(x2/x1)**2
    symbI(3,2)=symbI(2,1)*(x2/x1)**2
    symbI(3,3)=symbI(3,1)*(x2/x1)**2
  else
! collisional case:
!
    call W2_arr(x1,x2,Imn)
!
    if(.true.) then
!    if(.false.) then
! energy conservation:
      denom=(1.d0,0.d0)-Imn(0,0)+(2.d0,0.d0)*Imn(2,0)-Imn(2,2)
      do m=0,3
        do n=0,3
          symbI(m,n)=Imn(m,n)+(Imn(m,0)-Imn(m,2))*(Imn(n,0)-Imn(n,2))/denom
        enddo
      enddo
    else
      symbI=Imn
    endif
!
  endif
!
end subroutine getIfunc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!

subroutine smooth_array_gauss(dimx, ngauss, y, ys)

implicit none

integer :: dimx, ngauss,i
real(8) :: sgauss,dummy
real(8), dimension(dimx) :: y, ys
double precision, dimension(:), allocatable :: wgauss
integer :: k, j;
real(8) :: sum

sgauss=dfloat(ngauss)/5.d0
allocate(wgauss(-ngauss:ngauss))
do i=-ngauss,ngauss
  wgauss(i)=exp(-(dfloat(i)/sgauss)**2)
enddo
dummy=sum(wgauss)
wgauss=wgauss/dummy

do i=ngauss+1,dimx-ngauss
  ys(i)=sum(wgauss*y(i-ngauss:i+ngauss))
enddo
do i=1,ngauss
  ys(i)=sum(wgauss(1-i:ngauss)*y(1:i+ngauss))/sum(wgauss(1-i:ngauss))
enddo
do i=dimx-ngauss+1,dimx
  ys(i)=sum(wgauss(-ngauss:dimx-i)*y(i-ngauss:dimx))/sum(wgauss(-ngauss:dimx-i))
enddo

end subroutine
!
subroutine equipotentials
!
  use baseparam_mod,  only : c,rtor
  use wave_code_data, only : dim_r,m_vals,n_vals,r,om_E,ks,kp,Es,Br,B0 &
                            ,vac_cd_ptr
  use diag_mod,       only : iunit_diag

  use mpi
!
  implicit none
!
  integer :: ierror,irank;
  integer :: nr,i,m,n,ipoi,ierr
  double complex :: amnp,amnt,ampl
  double precision, dimension(:), allocatable :: psi0,phi0
  double complex,   dimension(:), allocatable :: psi1,phi1,dum,Brv
  integer :: ind = 1
!
  nr=dim_r
  m=m_vals(ind)
  n=n_vals(ind)
!
  allocate(psi0(nr),phi0(nr),psi1(nr),phi1(nr),dum(nr),Brv(nr))
!
  psi0(1)=0.d0
  phi0(1)=0.d0
  ipoi=1
!
  do i=2,nr
    psi0(i)=psi0(i-1)+0.5d0*(r(i)-r(i-1))*(B0(i)*kp(i)+B0(i-1)*kp(i-1))
    phi0(i)=phi0(i-1)+0.5d0*(r(i)-r(i-1))/c                   &
          *(om_E(i)*B0(i)/ks(i)+om_E(i-1)*B0(i-1)/ks(i-1))
    if(kp(i)*kp(i-1).lt.0) ipoi=i
  enddo
!
  call amn_of_r(-m,n,r(ipoi),amnp,amnt,ierr)
  call get_wave_fields_from_wave_code(vac_cd_ptr(ind),nr,r,m,n,                    & 
                                      dum,dum,dum,dum,dum,Brv,dum,dum,dum,dum)

!
  ampl = 2.d0*amnt*dble(n)/(r(ipoi)*rtor*Brv(ipoi))
!  ampl = 2.d0*dble(n)/(r(ipoi)*rtor*Brv(ipoi))  ! DIII-D
!
  psi1=(0.d0,1.d0)*Br*ampl
  phi1=(0.d0,1.d0)*Es/ks*ampl
!
  call MPI_Comm_rank(MPI_COMM_WORLD,irank,ierror);
  if (irank .eq. 0 ) then
  open(iunit_diag,file='equipotentials.dat')
  do i=1,nr
    write(iunit_diag,*) r(i),psi0(i),phi0(i),real(psi1(i)),dimag(psi1(i)) &
                       ,real(phi1(i)),dimag(phi1(i)),abs(Br(i))           &
                       ,abs(Br(i)-c*kp(i)*Es(i)/om_E(i))
  enddo
  close(iunit_diag)
  end if
  
  deallocate(psi0,phi0,psi1,phi1,dum,Brv)
!
end subroutine equipotentials
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine calc_parallel_current_directly

  use grid_mod, only : npoib,rb,params_b,Ercov,ddr_params_nl
  use baseparam_mod, only : Z_i,e_charge,am,p_mass,c,btor,e_mass,ev,rtor
  use wave_code_data
  use mpi
!
  implicit none
!
  integer :: ierror,irank;
  integer :: ipoi,i,iunit,mnmax
  double precision, dimension(:), allocatable :: dummy,x1,x2,vT
  double complex,   dimension(:),     allocatable :: curr_e_par
  double complex,   dimension(:,:,:), allocatable :: symbI
!
  iunit=731
  mnmax=3
  allocate(x1(npoib),x2(npoib),symbI(0:mnmax,0:mnmax,npoib))
  allocate(curr_e_par(npoib),vT(npoib))
!
  vT=sqrt(params_b(3,:)/e_mass)
!
  x1=kp*vT/nue
  x2=-om_E/nue
!
  do i=1,npoib
    call getIfunc(x1(i),x2(i),symbI(:,:,i))
  enddo
!
! Here x1 and x2 are used for A_1 and A_2:
  x2=ddr_params_nl(3,:)/params_b(3,:)
  x1=ddr_params_nl(1,:)/params_b(1,:)+e_charge*Ercov/params_b(3,:)-1.5d0*x2
!
  curr_e_par=e_charge*params_b(1,:)*vT/(nue*B0)                          &
            *(c*Es*((x1+x2)*symbI(1,0,:)+0.5d0*x2*symbI(2,1,:))          &
            +vT*Br*((x1+x2)*symbI(1,1,:)+0.5d0*x2*symbI(3,1,:)))
!
  call MPI_Comm_rank(MPI_COMM_WORLD,irank,ierror);
  if (irank .eq. 0 ) then
  open(iunit,file='par_current_e.dat')
  open(10000,file='cond_e.dat')
  do ipoi=1,npoib
    write(iunit,*) rb(ipoi),real(curr_e_par(ipoi)),dimag(curr_e_par(ipoi)) &
                           ,real(Jpe(ipoi)),dimag(Jpe(ipoi))
    write(10000,*) rb(ipoi),real(symbI(1,0,ipoi)),dimag(symbI(1,0,ipoi)) &
                           ,real(symbI(1,1,ipoi)),dimag(symbI(1,1,ipoi)) &
                           ,real(symbI(2,1,ipoi)),dimag(symbI(2,1,ipoi)) &
                           ,real(symbI(3,1,ipoi)),dimag(symbI(3,1,ipoi)) 

  enddo
  close(10000)
  close(iunit)
  end if
!
  deallocate(x1,x2,symbI,curr_e_par,vT)
!
  end subroutine calc_parallel_current_directly
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine calc_ion_parallel_current_directly

  use grid_mod, only : npoib,rb,params_b,Ercov,ddr_params_nl
  use baseparam_mod, only : Z_i,e_charge,am,p_mass,c,btor,e_mass,ev,rtor
  use wave_code_data
  use mpi
!
  implicit none
!
  integer :: ierror,irank;
  integer :: ipoi,i,iunit,mnmax
  double precision :: ei_charge
  double precision, dimension(:), allocatable :: dummy,x1,x2,vT
  double complex,   dimension(:),     allocatable :: curr_i_par
  double complex,   dimension(:,:,:), allocatable :: symbI
!
  iunit=731
  mnmax=3
  allocate(x1(npoib),x2(npoib),symbI(0:mnmax,0:mnmax,npoib))
  allocate(curr_i_par(npoib),vT(npoib))
!
  ei_charge=-Z_i*e_charge
!
!
  vT=sqrt(params_b(4,:)/p_mass)
!
  x1=kp*vT/nui
  x2=-om_E/nui
!
  do i=1,npoib
    call getIfunc(x1(i),x2(i),symbI(:,:,i))
  enddo
!
! Here x1 and x2 are used for A_1 and A_2:
  x2=ddr_params_nl(4,:)/params_b(4,:)
  x1=ddr_params_nl(1,:)/params_b(1,:)+ei_charge*Ercov/params_b(4,:)-1.5d0*x2
!
  curr_i_par=ei_charge*params_b(1,:)*vT/(nui*B0)                          &
            *(c*Es*((x1+x2)*symbI(1,0,:)+0.5d0*x2*symbI(2,1,:))          &
            +vT*Br*((x1+x2)*symbI(1,1,:)+0.5d0*x2*symbI(3,1,:)))
!
  call MPI_Comm_rank(MPI_COMM_WORLD,irank,ierror);
  if (irank .eq. 0 ) then
  open(iunit,file='par_current_i.dat')
  do ipoi=1,npoib
    write(iunit,*) rb(ipoi),real(curr_i_par(ipoi)),dimag(curr_i_par(ipoi)) &
                           ,real(Jpi(ipoi)),dimag(Jpi(ipoi))
  enddo
  close(iunit)
  end if
!
  deallocate(x1,x2,symbI,curr_i_par,vT)
!
  end subroutine calc_ion_parallel_current_directly



