!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine field_mesh3d(rrr,pp,zzz,Brad,Bphi,Bzet,dBrdR,dBrdp,dBrdZ  &
     ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)

  use magfield_mod, only : ierrfield,input_format,npmid,nr,np,nz,npoint,ipoint,&
    rad,phi,zet,Brs,Bzs,Bps,Bx,By,Bz,Br,Bp
  use polylag_3, only : mp,indef, plag3d

  implicit double precision (a-h,o-z)

  dimension xp(mp),yp(mp),zp(mp),fp(mp,mp,mp)
  integer indx(mp), indy(mp), indz(mp)
  character(len=30) :: namedim(3), namevar(3)
  character(len=80) :: bez
  data icall/0/
  save

  !-------first call: read data from disk-------------------------------
  if(icall .eq. 0) then
    icall = 1
    iunit1=77

    open(iunit1,file='MESH3D/field.format')
    read(iunit1,*) input_format
    read(iunit1,*) input_units
    if(input_format.eq.0) then
      read(iunit1,*) nr
      read(iunit1,*) np
      read(iunit1,*) nz
    endif
    close(iunit1)

    if(input_format.eq.0) then

      open(iunit1,file='MESH3D/field.dat')
      read(iunit1,*)
      read(iunit1,*)
      read(iunit1,*) bez
      print *,'description:',bez

      allocate(Bx(nr,np,nz),By(nr,np,nz))
      allocate(Br(nr,np,nz),Bp(nr,np,nz),Bz(nr,np,nz))
      allocate(rad(nr),phi(np),zet(nz))

      read(iunit1,*)
      read(iunit1,*)
      read(iunit1,*) ndim

      read(iunit1,*)
      read(iunit1,*)
      read(iunit1,*) nvar

      read(iunit1,*)
      read(iunit1,*)
      do i=1,ndim
        read(iunit1,*) namedim(i)
      end do
      print *,'dimension:',ndim,namedim(1),namedim(2),namedim(3)
      do i=1,nvar
        read(iunit1,*) namevar(i)
      end do
      print *,'function:',nvar,namevar(1),namevar(2),namevar(3)

      read(iunit1,*)
      read(iunit1,*)
      read(iunit1,*) nwert
      print *,'number of values:',nwert
      read(iunit1,*)
      read(iunit1,*)

      !---Input B      -->T = V*s/m/m
      do i=1,nr
        do j=1,np
          do k=1,nz
            read(iunit1,*) aaa,bbb,ccc,Bx(i,j,k),By(i,j,k),Bz(i,j,k)
            if(i.eq.1 .and. j.eq.1 .and. k.eq.1) then
              rmin = aaa * 100.d0 !cm
              pmin = bbb
              zmin = ccc * 100.d0 !cm
            elseif(i.eq.1 .and. j.eq.2 .and. k.eq.1) then
              delp = bbb
            elseif(i.eq.nr .and. j.eq.np .and. k.eq.nz) then
              rmax = aaa * 100.d0 !cm
              pmax = bbb
              zmax = ccc * 100.d0 !cm
            end if
          end do
        end do
      end do

      close(iunit1)

      do j=1,np
        ptmp = (j-1)*hphi
        sinp = sin(ptmp)
        cosp = cos(ptmp)
        do i=1,nr
          do k=1,nz
            Br(i,j,k) = Bx(i,j,k)*cosp + By(i,j,k)*sinp
            Bp(i,j,k) = -Bx(i,j,k)*sinp + By(i,j,k)*cosp
          end do
        end do
      end do

      deallocate(Bx,By)

      hrad = (rmax - rmin)/(nr-1)
      hphi = (pmax - pmin)/(np-1)
      hzet = (zmax - zmin)/(nz-1)
      np = np + 1
      pmax = pmax + delp
      do i=1,nr
        rad(i) = rmin + hrad*(i-1)
      end do
      do i=1,np
        phi(i) = pmin + hphi*(i-1)
      end do
      do i=1,nz
        zet(i) = zmin + hzet*(i-1)
      end do

      do i=1,nr
        do k=1,nz
          Br(i,np,k) = Br(i,1,k)
          Bp(i,np,k) = Bp(i,1,k)
          Bz(i,np,k) = Bz(i,1,k)
        end do
      end do

      pmin_0=pmin

    elseif(input_format.eq.1) then

      open(iunit1,file='MESH3D/field.dat')
      read(iunit1,*) nr,np,nz
      read(iunit1,*) rmin,rmax
      read(iunit1,*) pmin,pmax
      read(iunit1,*) zmin,zmax
      if(input_units.eq.0) then
        rmin = rmin * 100.d0 !cm
        rmax = rmax * 100.d0 !cm
        zmin = zmin * 100.d0 !cm
        zmax = zmax * 100.d0 !cm
      end if
      np=np+4
      allocate(Br(nr,np,nz),Bp(nr,np,nz),Bz(nr,np,nz))
      allocate(rad(nr),phi(np),zet(nz))
      do i=1,nr
        do j=3,np-2
          do k=1,nz
            read(iunit1,*) Br(i,j,k),Bp(i,j,k),Bz(i,j,k)
          end do
        end do
      end do
      close(iunit1)

      hrad = (rmax - rmin)/(nr-1)
      hphi = (pmax - pmin)/(np-5)
      hzet = (zmax - zmin)/(nz-1)

      do i=1,nr
        rad(i) = rmin + hrad*(i-1)
      end do
      do i=1,np
        phi(i) = pmin + hphi*(i-3)
      end do
      do i=1,nz
        zet(i) = zmin + hzet*(i-1)
      end do

      Br(:,1,:)=Br(:,np-4,:)
      Br(:,2,:)=Br(:,np-3,:)
      Br(:,np-1,:)=Br(:,4,:)
      Br(:,np,:)=Br(:,5,:)
      Bp(:,1,:)=Bp(:,np-4,:)
      Bp(:,2,:)=Bp(:,np-3,:)
      Bp(:,np-1,:)=Bp(:,4,:)
      Bp(:,np,:)=Bp(:,5,:)
      Bz(:,1,:)=Bz(:,np-4,:)
      Bz(:,2,:)=Bz(:,np-3,:)
      Bz(:,np-1,:)=Bz(:,4,:)
      Bz(:,np,:)=Bz(:,5,:)

      pmin_0=pmin-2.d0*hphi

    elseif(input_format.eq.3) then

      print *,'input_format = 3 : sparse field with stellarator symmetry'
      open(iunit1,file='MESH3D/field_sparse.dat')
      read(iunit1,*)
      read(iunit1,*)
      read(iunit1,*)
      read(iunit1,*)
      npoint=0
      do
        read(iunit1,*,end=11)
        npoint=npoint+1
      end do
11    close(iunit1)
      open(iunit1,file='MESH3D/field_sparse.dat')
      read(iunit1,*) nr,np,nz
      read(iunit1,*) rmin,rmax
      read(iunit1,*) pmin,pmax
      read(iunit1,*) zmin,zmax

      pmax=2*pmax-pmin

      npmid=np+2
      np=2*np+3

      allocate(ipoint(nr,np,nz))
      ipoint=0

      allocate(Brs(0:npoint),Bps(0:npoint),Bzs(0:npoint))
      Brs(0)=0.d0
      Bps(0)=0.d0
      Bzs(0)=0.d0

      do i=1,npoint
        read(iunit1,*) ir,ip,iz,Brs(i),Bps(i),Bzs(i)
        ipoint(ir,ip+2,iz)=i
      end do

      close(iunit1)

      do ir=1,nr
        do iz=1,nz
          ipoint(ir,1,iz)=ipoint(ir,5,nz+1-iz)
          ipoint(ir,2,iz)=ipoint(ir,4,nz+1-iz)
        end do
      end do

      do ir=1,nr
        do ip=npmid+1,np
          do iz=1,nz
            ipoint(ir,ip,iz)=ipoint(ir,np+1-ip,nz+1-iz)
          end do
        end do
      end do

      allocate(rad(nr),phi(np),zet(nz))

      hrad = (rmax - rmin)/(nr-1)
      hphi = (pmax - pmin)/(np-5)
      hzet = (zmax - zmin)/(nz-1)

      do i=1,nr
        rad(i) = rmin + hrad*(i-1)
      end do
      do i=1,np
        phi(i) = pmin + hphi*(i-3)
      end do
      do i=1,nz
        zet(i) = zmin + hzet*(i-1)
      end do

      pmin_0=pmin-2.d0*hphi

    else

      print *,'unknown field format'
      stop

    end if

    hrm1 = 1.d0/hrad
    hpm1 = 1.d0/hphi
    hzm1 = 1.d0/hzet
    phi_period=pmax - pmin

  end if
  !------- end first call ----------------------------------------------

  ierrfield=0

  if(pp.gt.pmax) then
    ppp=pp-phi_period*(int((pp-pmax)/phi_period)+1)
  elseif(pp.lt.pmin) then
    ppp=pp+phi_period*(int((pmin-pp)/phi_period)+1)
  else
    ppp=pp
  end if

  call indef(rrr,rmin,hrm1,nr,indx)

  call indef(ppp,pmin_0,hpm1,np,indy)

  call indef(zzz,zmin,hzm1,nz,indz)

  do i=1,mp
    xp(i) = rad(indx(i))
    yp(i) = phi(indy(i))
    zp(i) = zet(indz(i))
  end do

  do k=1,mp
    do j=1,mp
      if(input_format.eq.3) then
        if(indy(j).le.2) then
          sigbr=-1.d0
        elseif(indy(j).le.npmid) then
          sigbr=1.d0
        elseif(indy(j).le.np-2) then
          sigbr=-1.d0
        else
          sigbr=1.d0
        end if
      end if
      do i=1,mp
        if(input_format.eq.3) then
          fp(i,j,k) = Brs(ipoint(indx(i),indy(j),indz(k)))*sigbr
        else
          fp(i,j,k) = Br(indx(i),indy(j),indz(k))
        end if
      end do
    end do
  end do
  call plag3d(rrr,ppp,zzz,fp,hrm1,hpm1,hzm1,xp,yp,zp, &
       polylag,poly1x,poly1y,poly1z)
  Brad = polylag
  dBrdR = poly1x
  dBrdp = poly1y
  dBrdZ = poly1z

  do k=1,mp
    do j=1,mp
      do i=1,mp
        if(input_format.eq.3) then
          fp(i,j,k) = Bps(ipoint(indx(i),indy(j),indz(k)))
        else
          fp(i,j,k) = Bp(indx(i),indy(j),indz(k))
        endif
        if(abs(fp(i,j,k)).le.0.d0) then
          ierrfield=1
          print *,'boundary touched'

          return
        end if
      end do
    end do
  end do
  call plag3d(rrr,ppp,zzz,fp,hrm1,hpm1,hzm1,xp,yp,zp &
       ,polylag,poly1x,poly1y,poly1z)
  Bphi = polylag
  dBpdR = poly1x
  dBpdp = poly1y
  dBpdZ = poly1z

  do k=1,mp
    do j=1,mp
      do i=1,mp
        if(input_format.eq.3) then
          fp(i,j,k) = Bzs(ipoint(indx(i),indy(j),indz(k)))
        else
          fp(i,j,k) = Bz(indx(i),indy(j),indz(k))
        end if
      end do
    end do
  end do
  call plag3d(rrr,ppp,zzz,fp,hrm1,hpm1,hzm1,xp,yp,zp &
       ,polylag,poly1x,poly1y,poly1z)
  Bzet = polylag
  dBzdR = poly1x
  dBzdp = poly1y
  dBzdZ = poly1z

end subroutine field_mesh3d
