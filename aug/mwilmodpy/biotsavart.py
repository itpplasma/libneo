import numpy as np
from IPython import embed

rtwopi = 1.0 / (2.*np.pi) 
twopi = (2.*np.pi)
mu0 = 1.25663706144e-6
eps0 = 1.e-8

def K_EllInteg1( x , lx ):                           
#!*****************************************************
#! approximation of elliptic integral of the 1. kind
#! absolute accuracy appr. 2.E-8 !
#! see :HANDBOOK OF MATHEMATICAL FUNCTIONS, M. ABRAMOWITZ AND I.A. STEGUN, eqn. 17.3.34
#! 2 input values to accelerate computation, although only one would be needed (k^2)
#! 1. x  == m1    = 1-m    = 1-k^2   (eqn. 17.2.18 / 17.2.19f)
#! 2. lx == ln(x) = ln(m1) = ln(1-k^2)
    K_EllInteg1 = (((.01451196212*x+.03742563713)*x +.03590092383)*x+.09666344259)*x+1.38629436112  -((((.00441787012*x+.03328355346)*x+.06880248576)*x+.12498593597)*x+.5)*lx       
    return K_EllInteg1

def E_EllInteg2( x , lx ):                            
#!*****************************************************
#! approximation of elliptic integral of the 2. kind
#! absolute accuracy appr. 2.E-8 !
#! see :HANDBOOK OF MATHEMATICAL FUNCTIONS, M. ABRAMOWITZ AND I.A. STEGUN, eqn. 17.3.36
#! 2 input values to accelerate computation, although only one would be needed (k^2)
#! 1. x  == m1    = 1-m    = 1-k^2   (eqn. 17.2.18 / 17.2.19f)
#! 2. lx == ln(x) = ln(m1) = ln(1-k^2)

    E_EllInteg2 = (((.01736506451*x+.04757383546)*x+.06260601220)*x +.44325141463)*x+1.0 -(((.00526449639*x+.04069697526)*x+.09200180037)*x +.24998368310)*x*lx       
    return E_EllInteg2

def BrBzflux_Unit_current_s(ra, za, r, z):
#!*****************************************************************
#! computes flux Br and Bz at points (r,z) for unit current at (ra,za)
#! Br, Bz is taken from Feneberg, Lackner, Martin, Computer Physics Communications 31 (1984) 143-148
#! NOTE: Vertausche in Gln (4) Dz mit Dr und in Gln (5) Dr mit Dz. Sonst ist Br und Bz vertauscht.
    zmza2 = (z-za)**2.
    r2    = r**2.
    ra2   = ra**2.
    rzsum = (r+ra)**2. + zmza2
    s1    = (r-ra)**2. + zmza2
    rsum  = ra*r/rzsum
    acl = np.array(1.0-4.0*rsum)
    acl[acl<eps0] = eps0
    alg   = np.log(acl)
    
    K_EllInteg1_val = K_EllInteg1(acl,alg)
    E_EllInteg2_val = E_EllInteg2(acl,alg)
    
    prefac = mu0 * rtwopi / np.sqrt(rzsum)
    
    dpsi_dz_tmp = prefac * (K_EllInteg1_val - E_EllInteg2_val * (ra2 + r2 + zmza2) / s1) * (z-za)
    dpsi_dz = dpsi_dz_tmp * twopi
    Br = - dpsi_dz_tmp / r

    dpsi_dr_tmp = prefac * ( K_EllInteg1_val + E_EllInteg2_val * (ra2 - r2 - zmza2) / s1) * r
    dpsi_dr = dpsi_dr_tmp * twopi
    Bz = + dpsi_dr_tmp / r
    
    flx  = np.sqrt(rzsum)*((1.0-2.0*rsum)*K_EllInteg1_val-E_EllInteg2_val)*rtwopi
    #embed()
    return Br,Bz,flx


#! ra, za in  [m] coordinates of current filament
#! r, z in  [m] coordinates of grid
#! [A] current in filament at (x_curr,y_curr)
#! flx is equal to G(r,r*) in Lackner Computer Physics Communications 12 (1976) 33-44
def psiFlux_Amp_Current(ra, za,currentIn, r, z):
    sum  = (ra+r)**2 + (z-za)**2
    rsum = ra*r/sum
    acl = np.array(1.0-4.0*rsum)
    acl[acl<eps0] = eps0
    alg  = np.log(acl)
    flx  = np.sqrt(sum)*((1.0-2.0*rsum)*K_EllInteg1(acl,alg)-E_EllInteg2(acl,alg))*rtwopi
    
    flx_const = twopi * mu0* currentIn
    return flx*flx_const

#subroutine Brz_unit_current_s(ra, za, r, z, dpsi_dr, dpsi_dz, Br, Bz)
#!*****************************************************************
#! computes flux Br and Bz at points (r,z) for unit current at (ra,za)
#! Br, Bz is taken from Feneberg, Lackner, Martin, Computer Physics Communicatio#ns 31 (1984) 143-148
#! NOTE: Vertausche in Gln (4) Dz mit Dr und in Gln (5) Dr mit Dz. Sonst ist Br und Bz vertauscht.
#use constants,        only: mu0, twopi, rtwopi                         ! 1/(2 pi)
#use distributions,    only: K_EllInteg1, E_EllInteg2
#real(rkind), intent(in)            :: ra, za             ! coordinates of unit current
#real(rkind), intent(in)            :: r, z               ! coordinates of dpsi_dr, dpsi_dz, Br, Bz
#real(rkind), intent(out), optional :: dpsi_dr, dpsi_dz   ! flux gradients
#real(rkind), intent(out), optional :: Br, Bz             ! poloidal magnetic field components
#real(rkind), parameter    :: eps0 = 1.d-8
#real(rkind)    :: rzsum, rsum, s1
#real(rkind)    :: alg,acl
#real(rkind)    :: K_EllInteg1_val, E_EllInteg2_val, prefac
#real(rkind)    :: zmza2, r2, ra2
#real(rkind)    :: dpsi_dr_tmp, dpsi_dz_tmp
#zmza2 = (z-za)**2
#r2    = r**2
#ra2   = ra**2
#rzsum = (r+ra)**2 + zmza2
#s1    = (r-ra)**2 + zmza2
#rsum  = ra*r/rzsum
#acl   = max(eps0, 1.0d0-4.0d0*rsum)
#alg   = log(acl)
#K_EllInteg1_val = K_EllInteg1(acl,alg)
#E_EllInteg2_val = E_EllInteg2(acl,alg)
#prefac = mu0 * rtwopi / sqrt(rzsum)
#if (present(dpsi_dz) .or. present(Br)) then
#  dpsi_dz_tmp = prefac * (K_EllInteg1_val - E_EllInteg2_val * (ra2 + r2 + zmza2#) / s1) * (z-za)
#  if (present(dpsi_dz)) dpsi_dz = dpsi_dz_tmp * twopi
#  if (present(Br)) Br = - dpsi_dz_tmp / r
#endif
#if (present(dpsi_dr) .or. present(Bz)) then
#  dpsi_dr_tmp = prefac * ( K_EllInteg1_val + E_EllInteg2_val * (ra2 - r2 - zmza2) / s1) * r
#  if (present(dpsi_dr)) dpsi_dr = dpsi_dr_tmp * twopi
#  if (present(Bz)) Bz = + dpsi_dr_tmp / r
#endif
#!write(6,'(a,2e12.4)')"Br/z new=", Br, Bz
#!if (present(Br)) &
#!  Br = prefac * (-K_EllInteg1_val + E_EllInteg2_val * (ra2 + r2 + zmza2) / s1) * (z-za) / r
#!if (present(Bz)) &
#!  Bz = prefac * ( K_EllInteg1_val + E_EllInteg2_val * (ra2 - r2 - zmza2) / s1)
#!write(6,'(a,2e12.4)')"Br/z old=", Br, Bz
#!stop "sub Brz_unit_current_s"
#end subroutine Brz_unit_current_s
#subroutine psv_s(ra,za,r,z,flx)
#!*****************************************************************
#! computes flux "flx" at points (r,z) for unit current at (ra,za)
#! flx is equal to G(r,r*) in Lackner Computer Physics Communications 12 (1976) 33-44
#! n  - number of points p(r,z)
#! /afs/ipp/u/rrf/F90/AUG_equil/Org/equil/solver/psv.f
#use constants,        only: rtwopi                         ! 1/(2 pi)
#use distributions,    only: K_EllInteg1, E_EllInteg2
#real(rkind),   intent(in) :: ra, za
#real(rkind),   intent(in) :: r, z
#real(rkind),   intent(out):: flx
#real(rkind), parameter    :: eps0 = 1.d-8
#real(rkind)    :: sum, rsum
#real(rkind)    :: alg,acl

#sum  = (ra+r)**2 + (z-za)**2
#rsum = ra*r/sum
#acl  = max(eps0, 1.0d0-4.0d0*rsum)
#alg  = log(acl)
#flx  = sqrt(sum)*((1.0d0-2.0d0*rsum)*K_EllInteg1(acl,alg)-E_EllInteg2(acl,alg))*rtwopi
#end subroutine psv_s
#function K_EllInteg1(x,lx)                             
#!*****************************************************
#! approximation of elliptic integral of the 1. kind
#! absolute accuracy appr. 2.E-8 !
#! see :HANDBOOK OF MATHEMATICAL FUNCTIONS, M. ABRAMOWITZ AND I.A. STEGUN, eqn. 17.3.34
#! 2 input values to accelerate computation, although only one would be needed (k^2)
#! 1. x  == m1    = 1-m    = 1-k^2   (eqn. 17.2.18 / 17.2.19f)
#! 2. lx == ln(x) = ln(m1) = ln(1-k^2)
#use f90_kind
#implicit none
#real(rkind), intent(in) :: x, lx
#real(rkind)             :: K_EllInteg1
#K_EllInteg1 = (((.01451196212D0*x+.03742563713D0)*x                     &
#               +.03590092383D0)*x+.09666344259D0)*x+1.38629436112D0     &
#            -((((.00441787012D0*x+.03328355346D0)*x+.06880248576D0)*x   &
#               +.12498593597D0)*x+.5D0)*lx                      
#end function K_EllInteg1
#function E_EllInteg2(x,lx)                             
#!*****************************************************
#! approximation of elliptic integral of the 2. kind
#! absolute accuracy appr. 2.E-8 !
#! see :HANDBOOK OF MATHEMATICAL FUNCTIONS, M. ABRAMOWITZ AND I.A. STEGUN, eqn. 17.3.36
#! 2 input values to accelerate computation, although only one would be needed (k^2)
#! 1. x  == m1    = 1-m    = 1-k^2   (eqn. 17.2.18 / 17.2.19f)
#! 2. lx == ln(x) = ln(m1) = ln(1-k^2)
#use f90_kind
#implicit none
#real(rkind), intent(in) :: x, lx
#real(rkind)             :: E_EllInteg2
#E_EllInteg2 = (((.01736506451D0*x+.04757383546D0)*x+.06260601220D0)*x   &
#               +.44325141463D0)*x+1.0D0                                 &
#             -(((.00526449639D0*x+.04069697526D0)*x+.09200180037D0)*x   &
#               +.24998368310D0)*x*lx       
#end function E_EllInteg2
