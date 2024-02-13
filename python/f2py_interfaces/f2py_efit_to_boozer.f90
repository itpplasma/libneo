module efit_to_boozer
!
    use efit_to_boozer_mod
    use input_files, only : gfile
!
    implicit none
!
    integer :: nstep,nsurfmax,i,is,it,nsurf,nt,mpol,inp_label,m,iunit
    double precision :: s,theta,hs,htheta,aiota,aJb,B_theta,B_phi,oneovernt,twoovernt,sqrtg00
    double precision :: psi,q,dq_ds,C_norm,dC_norm_ds,sqrtg,bmod,dbmod_dtheta,sigma,     &
                        R,dR_ds,dR_dtheta,Z,dZ_ds,dZ_dtheta,G,dG_ds,dG_dtheta
    double precision :: phi,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
    double precision :: dB_phi_ds,dB_theta_ds,pprime
    double complex   :: four_ampl
!
    double precision, dimension(:), allocatable :: R_oft,Z_oft,al_oft,B_oft
    double precision, dimension(:), allocatable :: Rmn_c,Rmn_s,Zmn_c,Zmn_s,almn_c,almn_s,Bmn_c,Bmn_s
    double complex,   dimension(:), allocatable :: calE,calEm_times_Jb

    contains

    subroutine init()
!
    open(1,file='efit_to_boozer.inp')
    read (1,*) nstep    !number of integration steps
    read (1,*) nlabel   !grid size over radial variable
    read (1,*) ntheta   !grid size over poloidal angle
    read (1,*) nsurfmax !number of starting points between the
                        !magnetic axis and right box boundary
                        !when searching for the separatrix
    read (1,*) nsurf    !number of flux surfaces in Boozer file
    read (1,*) mpol     !number of poloidal modes in Boozer file
    read (1,*) psimax   !psi at plasma boundary
    close(1)
!
    if (allocated(rbeg)) deallocate(rbeg)
    if (allocated(rsmall)) deallocate(rsmall)
    if (allocated(qsaf)) deallocate(qsaf)
    if (allocated(psi_pol)) deallocate(psi_pol)
    if (allocated(psi_tor_vac)) deallocate(psi_tor_vac)
    if (allocated(psi_tor)) deallocate(psi_tor)
    if (allocated(C_const)) deallocate(C_const)
!
    allocate(rbeg(nlabel),rsmall(nlabel),qsaf(nlabel),psi_pol(0:nlabel))
    allocate(psi_tor_vac(nlabel),psi_tor(0:nlabel),C_const(nlabel))
!
    if (allocated(R_spl)) deallocate(R_spl)
    if (allocated(Z_spl)) deallocate(Z_spl)
    if (allocated(bmod_spl)) deallocate(bmod_spl)
    if (allocated(sqgnorm_spl)) deallocate(sqgnorm_spl)
    if (allocated(Gfunc_spl)) deallocate(Gfunc_spl)
!
    allocate(R_spl(0:nspl,0:ntheta,nlabel),Z_spl(0:nspl,0:ntheta,nlabel),bmod_spl(0:nspl,0:ntheta,nlabel))
    allocate(sqgnorm_spl(0:nspl,0:ntheta,nlabel),Gfunc_spl(0:nspl,0:ntheta,nlabel))

    call flint_for_Boozer(nstep,nsurfmax,nlabel,ntheta,      &
                        rmn,rmx,zmn,zmx,raxis,zaxis,sigma, &
                        rbeg,rsmall,qsaf,                  &
                        psi_pol(1:nlabel),psi_tor_vac,     &
                        psi_tor(1:nlabel),C_const,         &
                        R_spl(0,1:ntheta,:),               &
                        Z_spl(0,1:ntheta,:),               &
                        bmod_spl(0,1:ntheta,:),            &
                        sqgnorm_spl(0,1:ntheta,:),         &
                        Gfunc_spl(0,1:ntheta,:))
!
    call spline_magdata_in_symfluxcoord
!
    print *,'Splining done'
!
    nt=ntheta
!
    inp_label=1
    hs=1.d0/dfloat(nsurf)
    oneovernt=1.d0/dfloat(nt)
    twoovernt=2.d0*oneovernt
    htheta=twopi*oneovernt
!
    if (allocated(calE)) deallocate(calE)
    if (allocated(calEm_times_Jb)) deallocate(calEm_times_Jb)
    if (allocated(R_oft)) deallocate(R_oft)
    if (allocated(Z_oft)) deallocate(Z_oft)
    if (allocated(al_oft)) deallocate(al_oft)
    if (allocated(B_oft)) deallocate(B_oft)
    if (allocated(Rmn_c)) deallocate(Rmn_c)
    if (allocated(Rmn_s)) deallocate(Rmn_s)
    if (allocated(Zmn_c)) deallocate(Zmn_c)
    if (allocated(Zmn_s)) deallocate(Zmn_s)
    if (allocated(almn_c)) deallocate(almn_c)
    if (allocated(almn_s)) deallocate(almn_s)
    if (allocated(Bmn_c)) deallocate(Bmn_c)
    if (allocated(Bmn_s)) deallocate(Bmn_s)
!
    allocate(calE(nt),calEm_times_Jb(nt))
    allocate(R_oft(nt),Z_oft(nt),al_oft(nt),B_oft(nt))
    allocate(Rmn_c(0:mpol),Rmn_s(0:mpol),Zmn_c(0:mpol),Zmn_s(0:mpol))
    allocate(almn_c(0:mpol),almn_s(0:mpol),Bmn_c(0:mpol),Bmn_s(0:mpol))

    end subroutine init

    subroutine magdata(inp_label,s,psi,theta,q,dq_ds,C_norm,dC_norm_ds, &
        sqrtg,bmod,dbmod_dtheta,R,dR_ds,dR_dtheta,       &
        Z,dZ_ds,dZ_dtheta,G,dG_ds,dG_dtheta)

        integer, intent(in) :: inp_label
        double precision, intent(inout) :: s,psi,theta
        double precision, intent(out) :: q,dq_ds,C_norm,dC_norm_ds,sqrtg,bmod, &
        dbmod_dtheta,R,dR_ds,dR_dtheta,Z,dZ_ds,dZ_dtheta,G,dG_ds,dG_dtheta

        call magdata_in_symfluxcoord_ext(inp_label,s,psi,theta,q,dq_ds,C_norm,dC_norm_ds, &
        sqrtg,bmod,dbmod_dtheta,R,dR_ds,dR_dtheta,       &
        Z,dZ_ds,dZ_dtheta,G,dG_ds,dG_dtheta)
    end subroutine magdata
!
end module efit_to_boozer
