module projection_fixture
    implicit none
    real(8) :: last_theta
    real(8), parameter :: alpha=.18d0, beta=.035d0, epsilon=1d-5
end module

module input_files
    implicit none
    character(32) :: gfile='analytic angular map'
end module

module field_sub
    use projection_fixture
    implicit none
contains
    subroutine field_eq(r,phi,z,br,bp,bz,brr,brp,brz,bpr,bpp,bpz,bzr,bzp,bzz)
        real(8), intent(in) :: r,phi,z
        real(8), intent(out) :: br,bp,bz,brr,brp,brz,bpr,bpp,bpz,bzr,bzp,bzz
        real(8) :: jac
        jac=1+alpha*cos(last_theta)+2*beta*cos(2*last_theta)
        br=0; bz=0
        bp=jac*(1+epsilon*(.4d0+cos(17*last_theta)))/r
        brr=0; brp=0; brz=0; bpr=-bp/r; bpp=0; bpz=0
        bzr=0; bzp=0; bzz=0
    end subroutine
end module

subroutine flint_for_Boozer(nstep,nsurfmax,nlabel,ntheta, &
    rmn,rmx,zmn,zmx,raxis,zaxis,sigma,rbeg,rsmall,qsaf, &
    psi_pol,psi_tor_vac,psi_tor_plas,C_const,R_ts,Z_ts,bmod_ts,sqgnorm_ts,Gfunc_ts)
    implicit none
    integer :: nstep,nsurfmax,nlabel,ntheta
    real(8) :: rmn,rmx,zmn,zmx,raxis,zaxis,sigma
    real(8) :: rbeg(nlabel),rsmall(nlabel),qsaf(nlabel),psi_pol(nlabel)
    real(8) :: psi_tor_vac(nlabel),psi_tor_plas(nlabel),C_const(nlabel)
    real(8) :: R_ts(ntheta,nlabel),Z_ts(ntheta,nlabel),bmod_ts(ntheta,nlabel)
    real(8) :: sqgnorm_ts(ntheta,nlabel),Gfunc_ts(ntheta,nlabel)
    rmn=500; rmx=740; zmn=-150; zmx=150; raxis=620; zaxis=0; sigma=1
    rbeg=620; rsmall=100; qsaf=1.7d0; psi_pol=1; psi_tor_vac=1
    psi_tor_plas=1; C_const=1
    R_ts=0; Z_ts=0; bmod_ts=0; sqgnorm_ts=0; Gfunc_ts=0
end subroutine

subroutine spline_magdata_in_symfluxcoord
    use efit_to_boozer_mod
    implicit none
    psitor_max=1; psipol_max=1
end subroutine

subroutine magdata_in_symfluxcoord_ext(inp_label,s,psi,theta,q,dq_ds, &
    C_norm,dC_norm_ds,sqrtg,bmod,dbmod_dtheta, &
    R,dR_ds,dR_dtheta,Z,dZ_ds,dZ_dtheta,G,dG_ds,dG_dtheta)
    use projection_fixture
    implicit none
    integer :: inp_label
    real(8) :: s,psi,theta,q,dq_ds,C_norm,dC_norm_ds,sqrtg,bmod,dbmod_dtheta
    real(8) :: R,dR_ds,dR_dtheta,Z,dZ_ds,dZ_dtheta,G,dG_ds,dG_dtheta
    real(8) :: angle,jac
    last_theta=theta
    angle=theta+alpha*sin(theta)+beta*sin(2*theta)
    jac=1+alpha*cos(theta)+2*beta*cos(2*theta)
    q=1.7d0; dq_ds=0; C_norm=1; dC_norm_ds=0; sqrtg=1; psi=s
    G=q*(angle-theta); dG_ds=0; dG_dtheta=q*(jac-1)
    R=620+90*sqrt(s)*cos(angle)+7*s*cos(2*angle)+4*s*sin(3*angle)
    Z=140*sqrt(s)*sin(angle)+9*s*cos(2*angle)-2.5d0*s*sin(3*angle)
    dR_ds=45/sqrt(s)*cos(angle)+7*cos(2*angle)+4*sin(3*angle)
    dZ_ds=70/sqrt(s)*sin(angle)+9*cos(2*angle)-2.5d0*sin(3*angle)
    dR_dtheta=(-90*sqrt(s)*sin(angle)-14*s*sin(2*angle)+12*s*cos(3*angle))*jac
    dZ_dtheta=(140*sqrt(s)*cos(angle)-18*s*sin(2*angle)-7.5d0*s*cos(3*angle))*jac
    bmod=1+.15d0*cos(angle)+.05d0*sin(2*angle)
    dbmod_dtheta=(-.15d0*sin(angle)+.1d0*cos(2*angle))*jac
end subroutine
