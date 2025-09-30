!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine magfie_vmec(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
  ! Computes magnetic field module in units of the magnetic code  - bmod,
  ! square root of determinant of the metric tensor               - sqrtg,
  ! derivatives of the logarythm of the magnetic field module
  ! over coordinates                                              - bder,
  ! covariant componets of the unit vector of the magnetic
  ! field direction                                               - hcovar,
  ! contravariant components of this vector                       - hctrvr,
  ! contravariant component of the curl of this vector            - hcurl
  ! Order of coordinates is the following: x(1)=s (normalized toroidal flux),
  ! x(2)=theta (VMEC poloidal angle), x(3)=varphi (geometrical toroidal angle).
  !
  !  Input parameters:
  !            formal:  x(3)             - array of VMEC coordinates
  !  Output parameters:
  !            formal:  bmod
  !                     sqrtg
  !                     bder(3)          - derivatives of $\log(B)$
  !                     hcovar(3)        - covariant components of unit vector $\bh$ along $\bB$
  !                     hctrvr(3)        - contra-variant components of unit vector $\bh$ along $\bB$
  !                     hcurl(3)         - contra-variant components of curl of $\bh$
  !
  !  Called routines: vmec_field

  use libneo_kinds, only : dp
  use math_constants, only : TWOPI

  implicit none

  real(dp), parameter :: hs=1.d-3, ht=hs*TWOPI, hp=ht/5.d0

  real(dp), intent(out) :: bmod,sqrtg
  real(dp) :: s,theta,varphi,A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota, &
                      sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,     &
                      Bcovar_r,Bcovar_vartheta,Bcovar_varphi
  real(dp) :: cjac,bcov_s_vmec,bcov_t_vmec,bcov_p_vmec
  real(dp) :: dhs_dt,dhs_dp,dht_ds,dht_dp,dhp_ds,dhp_dt
  real(dp), dimension(3), intent(in) :: x
  real(dp), dimension(3), intent(out) :: bder,hcovar,hctrvr,hcurl

  ! Begin derivatives over s

  theta=x(2)
  varphi=x(3)
  s=x(1)+hs

  call vmec_field(s,theta,varphi,A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota,     &
                  sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,     &
                  Bcovar_r,Bcovar_vartheta,Bcovar_varphi)

  bmod=sqrt(Bctrvr_vartheta*Bcovar_vartheta+Bctrvr_varphi*Bcovar_varphi)
  bcov_s_vmec=Bcovar_r+dl_ds*Bcovar_vartheta
  bcov_t_vmec=(1.d0+dl_dt)*Bcovar_vartheta
  bcov_p_vmec=Bcovar_varphi+dl_dp*Bcovar_vartheta
  bder(1)=bmod
  dht_ds=bcov_t_vmec/bmod
  dhp_ds=bcov_p_vmec/bmod

  s=x(1)-hs

  call vmec_field(s,theta,varphi,A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota,     &
                  sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,     &
                  Bcovar_r,Bcovar_vartheta,Bcovar_varphi)

  bmod=sqrt(Bctrvr_vartheta*Bcovar_vartheta+Bctrvr_varphi*Bcovar_varphi)
  bcov_s_vmec=Bcovar_r+dl_ds*Bcovar_vartheta
  bcov_t_vmec=(1.d0+dl_dt)*Bcovar_vartheta
  bcov_p_vmec=Bcovar_varphi+dl_dp*Bcovar_vartheta
  bder(1)=(bder(1)-bmod)/(2.d0*hs)
  dht_ds=(dht_ds-bcov_t_vmec/bmod)/(2.d0*hs)
  dhp_ds=(dhp_ds-bcov_p_vmec/bmod)/(2.d0*hs)

  ! End derivatives over s
  !-------------------------
  ! Begin derivatives over theta

  s=x(1)
  theta=x(2)+ht

  call vmec_field(s,theta,varphi,A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota,     &
                  sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,     &
                  Bcovar_r,Bcovar_vartheta,Bcovar_varphi)

  bmod=sqrt(Bctrvr_vartheta*Bcovar_vartheta+Bctrvr_varphi*Bcovar_varphi)
  bcov_s_vmec=Bcovar_r+dl_ds*Bcovar_vartheta
  bcov_t_vmec=(1.d0+dl_dt)*Bcovar_vartheta
  bcov_p_vmec=Bcovar_varphi+dl_dp*Bcovar_vartheta
  bder(2)=bmod
  dhs_dt=bcov_s_vmec/bmod
  dhp_dt=bcov_p_vmec/bmod

  theta=x(2)-ht

  call vmec_field(s,theta,varphi,A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota,     &
                  sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,     &
                  Bcovar_r,Bcovar_vartheta,Bcovar_varphi)

  bmod=sqrt(Bctrvr_vartheta*Bcovar_vartheta+Bctrvr_varphi*Bcovar_varphi)
  bcov_s_vmec=Bcovar_r+dl_ds*Bcovar_vartheta
  bcov_t_vmec=(1.d0+dl_dt)*Bcovar_vartheta
  bcov_p_vmec=Bcovar_varphi+dl_dp*Bcovar_vartheta
  bder(2)=(bder(2)-bmod)/(2.d0*ht)
  dhs_dt=(dhs_dt-bcov_s_vmec/bmod)/(2.d0*ht)
  dhp_dt=(dhp_dt-bcov_p_vmec/bmod)/(2.d0*ht)

  ! End derivatives over theta
  !-------------------------
  ! Begin derivatives over varphi

  theta=x(2)
  varphi=x(3)+hp

  call vmec_field(s,theta,varphi,A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota,     &
                  sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,     &
                  Bcovar_r,Bcovar_vartheta,Bcovar_varphi)

  bmod=sqrt(Bctrvr_vartheta*Bcovar_vartheta+Bctrvr_varphi*Bcovar_varphi)
  bcov_s_vmec=Bcovar_r+dl_ds*Bcovar_vartheta
  bcov_t_vmec=(1.d0+dl_dt)*Bcovar_vartheta
  bcov_p_vmec=Bcovar_varphi+dl_dp*Bcovar_vartheta
  bder(3)=bmod
  dhs_dp=bcov_s_vmec/bmod
  dht_dp=bcov_t_vmec/bmod

  varphi=x(3)-hp

  call vmec_field(s,theta,varphi,A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota,     &
                  sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,     &
                  Bcovar_r,Bcovar_vartheta,Bcovar_varphi)

  bmod=sqrt(Bctrvr_vartheta*Bcovar_vartheta+Bctrvr_varphi*Bcovar_varphi)
  bcov_s_vmec=Bcovar_r+dl_ds*Bcovar_vartheta
  bcov_t_vmec=(1.d0+dl_dt)*Bcovar_vartheta
  bcov_p_vmec=Bcovar_varphi+dl_dp*Bcovar_vartheta
  bder(3)=(bder(3)-bmod)/(2.d0*hp)
  dhs_dp=(dhs_dp-bcov_s_vmec/bmod)/(2.d0*hp)
  dht_dp=(dht_dp-bcov_t_vmec/bmod)/(2.d0*hp)

  ! End derivatives over varphi
  !-------------------------

  varphi=x(3)

  call vmec_field(s,theta,varphi,A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota,     &
                  sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,     &
                  Bcovar_r,Bcovar_vartheta,Bcovar_varphi)

  bmod=sqrt(Bctrvr_vartheta*Bcovar_vartheta+Bctrvr_varphi*Bcovar_varphi)
  cjac=1.d0+dl_dt
  sqrtg=sqg*cjac
  bder=bder/bmod
  bcov_s_vmec=Bcovar_r+dl_ds*Bcovar_vartheta
  bcov_t_vmec=(1.d0+dl_dt)*Bcovar_vartheta
  bcov_p_vmec=Bcovar_varphi+dl_dp*Bcovar_vartheta
  hcovar(1)=bcov_s_vmec/bmod
  hcovar(2)=bcov_t_vmec/bmod
  hcovar(3)=bcov_p_vmec/bmod
  hctrvr(1)=0.d0
  hctrvr(2)=(Bctrvr_vartheta-dl_dp*Bctrvr_varphi)/(cjac*bmod)
  hctrvr(3)=Bctrvr_varphi/bmod
  hcurl(1)=(dhp_dt-dht_dp)/sqrtg
  hcurl(2)=(dhs_dp-dhp_ds)/sqrtg
  hcurl(3)=(dht_ds-dhs_dt)/sqrtg

end subroutine magfie_vmec


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
