!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  module rhs_surf_mod
    double precision :: dz_dphi
  end module rhs_surf_mod
!
!-------------------------------------------------------------------------------
module rhs_boozer_sub
    implicit none

contains
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine rhs_axis(phi,y,dy)
!
  use field_sub 
!
  implicit none
!
  double precision, intent(in) :: phi
  double precision, dimension(:), intent(in) :: y
  double precision, dimension(:), intent(out) :: dy
  double precision :: R,Z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ,   &
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
  dy(3)=y(1)
  dy(4)=y(2)
!
  return
  end subroutine rhs_axis
!
!
  subroutine rhs_surf(phi,y,dy)
!
  use field_sub 
  use rhs_surf_mod
!
  implicit none
!
  double precision, intent(in) :: phi
  double precision, dimension(:), intent(in) :: y
  double precision, dimension(:), intent(out) :: dy
  double precision :: R,Z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ,   &
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
  dy(3)=R*dy(2)
  dy(4)=R*Z*Br
  dy(5)=R*(Br**2+Bp**2+Bz**2)/Bp
!
  dz_dphi=dy(2)
!
  return
  end subroutine rhs_surf
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
end module rhs_boozer_sub
