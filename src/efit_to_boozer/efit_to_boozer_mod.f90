!
  module efit_to_boozer_mod
    integer, parameter :: nspl  = 3 !spline order in poloidal angle interpolation
    integer, parameter :: nplag = 4 !stencil for Largange polynomial interpolation (polynomial order + 1)
    integer, parameter :: nder  = 1 !number of derivatives from Largange polynomial interpolation
    logical :: load=.true.
    integer :: nlabel,ntheta
    double precision :: twopi = atan(1.d0)*8.d0
    double precision :: rmn,rmx,zmn,zmx,raxis,zaxis,h_theta,psipol_max,psitor_max
    double precision, dimension(nplag)        :: R_lag,Z_lag,sqrtg_lag,bmod_lag,G_lag
    double precision, dimension(nplag)        :: dbmod_dt_lag,dR_dt_lag,dZ_dt_lag,dG_dt_lag
    double precision, dimension(0:nder,nplag) :: coef
    double precision, dimension(:),     allocatable :: rbeg,rsmall,qsaf,psi_pol,psi_tor_vac,psi_tor,C_const
    double precision, dimension(:,:,:), allocatable :: R_spl,Z_spl,bmod_spl,sqgnorm_spl,Gfunc_spl
    double precision :: psimax = 1.d300
  end module efit_to_boozer_mod
!
