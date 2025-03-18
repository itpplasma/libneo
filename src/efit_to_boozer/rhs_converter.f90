!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  module rhs_converter_mod
    double precision :: dz_dphi
  end module rhs_converter_mod
!
!-------------------------------------------------------------------------------
!
  subroutine rhs_converter(phi,y,dy)
!
  use rhs_converter_mod, only : dz_dphi
!
  implicit none
!
  integer, parameter :: ndim = 4
!
  double precision, dimension(ndim) :: y,dy
  double precision :: R,phi,Z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ,   &
                      dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
!
  R=y(1)
  Z=y(2)
!
  call field_eq(R,phi,Z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ   &
               ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
  dy(1)=Br*R/Bp
  dy(2)=Bz*R/Bp
  dy(3)=R/Bp                       !d F_H / d phi
  dy(4)=dy(3)*(Br**2+Bp**2+Bz**2)  !d F_B / d phi
!
  dz_dphi=dy(2)
!
  return
  end subroutine rhs_converter
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
