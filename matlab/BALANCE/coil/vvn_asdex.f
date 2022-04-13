cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE gbcoil(rI,fI,ZI,BrI,BfI,BZI,
     *BrrI,BrfI,BrZI,BfrI,BffI,BfZI,BZrI,BZfI,BZZI)
      implicit real*8(a-h,o-z),integer(i-n)
 
      common/cfacb/facbc
      save facb,ier
      common/datwc/nnodc
      data ier/0/

      if(ier.eq.1) go to 10
      facb=facbc
      open(36,form='FORMATTED',file='fort.36')
   10 continue       

      x=ri
      y=zi

      CALL GBmodc(X,FI,Y,BX,BF,BY,
     *BRR,BRF,BRZ,BFR,BFF,BFZ,BZR,BZF,BZZ)

      if(ier.eq.0) then
      write(36,*)facb,'=facb',nnodc,'=nnodc'
      write(36,*)x,'=r_st',fi,'=fi_st',y,'=z_st'
      end if

      bri=BX*facb
      Bfi=Bf*facb
      bzi=BY*facb
      BRRi=BRR*facb
      BRfi=BRf*facb
      BRZi=BRZ*facb
      BfRi=BfR*facb
      Bffi=Bff*facb
      Bfzi=Bfz*facb
      BZRi=BZR*facb
      BZfi=BZf*facb
      BZZi=BZZ*facb

      if(ier.eq.0) then
      write(36,*)bri,'=br',bfi,'=bf',bzi,'=bz'
      close(36)
      ier=1 
      end if

      RETURN
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE GBmodc(rI,fI,ZI,BrI,BfI,BZI,
     *BrrI,BrfI,BrZI,BfrI,BffI,BfZI,BZrI,BZfI,BZZI)
cccc w7x version
      implicit real*8(a-h,o-z),integer(i-n)

c 25.07.2012      parameter (npar=10850,ncoil=27)
      parameter (npar=10850,ncoil=16)
      character*24 file2,file3
c 25.07.2012      data file2/'cur_it.dd'/
      data file2/'cur_asd.dd'/
c 25.07.2012      data file3/'co_itm.dd'/
      data file3/'co_asd.dd'/
      
      DIMENSION XO(npar),YO(npar),ZO(npar),cur(npar)
      dimension curco(ncoil),nco(npar)
c     GBX1

      logical prop
      common/datwc/nnodc
c 19.05.2011
      common/datcur/curco
      save
      DATA PROP/.FALSE./
      IF(PROP) GO TO 80
      PROP=.TRUE.

      bfrt=0d0      

      nparx=npar
      CALL datw7xm(XO,YO,ZO,cur,nco,nparx)
      
      npr=1
      K2=nnodc 


      print 10,bfrt
   10 format(' gbrfz from GBXOT   (r*8, GBX1) bfrt=',1pd17.10/)
      print*,'  w7x_g version'
      print*,'curco'
      print*,curco
   80 fd=fi
      rd=ri
      cosf=dcos(fd)
      sinf=dsin(fd)
      Y=rd*sinf
      X=rd*cosf
      Z=ZI
      BX=0d0
      BY=0d0
      BZ=0d0
      BXX=0d0
      BXY=0d0

      BXZ=0d0
      BYY=0d0
      BYZ=0d0
      BZZ=0d0
c 25.12.2008
      Byx=0d0
      Bzx=0d0
      Bzy=0d0
c 25.12.2008 end
      K0=1
      K1=2
      DO 120 N=1,NPR
      BXB=0d0
      BYB=0d0

      BZB=0d0
      BXXB=0d0
      BXYB=0d0
      BXZB=0d0
      BYYB=0d0
      BYZB=0d0
      BZZB=0d0
c 25.12.2008
      ByxB=0d0
      BzxB=0d0
      BzyB=0d0
c 25.12.2008 end
      XMXC=X-XO(K0)
      YMYC=Y-YO(K0)
      ZMZC=Z-ZO(K0)

      R1=dSQRT(XMXC*XMXC+YMYC*YMYC+ZMZC*ZMZC)
      OR1=1d0/R1
      R1X=XMXC*OR1
      R1Y=YMYC*OR1
      R1Z=ZMZC*OR1
      DO 100 K=K1,K2
      AX=XO(K)-XO(K-1)
      XMXC=X-XO(K)
      AY=YO(K)-YO(K-1)
      YMYC=Y-YO(K)

      AZ=ZO(K)-ZO(K-1)
      ZMZC=Z-ZO(K)
      ZPRA=AX*XMXC+AY*YMYC+AZ*ZMZC
      R2=dSQRT(XMXC*XMXC+YMYC*YMYC+ZMZC*ZMZC)
      OR2=1d0/R2
      R2X=XMXC*OR2
      R2Y=YMYC*OR2
      R2Z=ZMZC*OR2
      
      R1PR2=R1+R2
      OR1PR2=1d0/R1PR2
      
      if(cur(k-1).eq.0d0) go to 101
      
      OBCP=1d0/(R2*R1PR2+ZPRA)
      FAZRDA=-R1PR2*OBCP*OR1*OR2
      FLZA=-OBCP
      FLR1=-R2*OBCP+OR1PR2-OR1
      FLR2=-(2d0*R2+R1)*OBCP+OR1PR2-OR2
      FLX=FLR1*R1X+FLR2*R2X+FLZA*AX
      FLY=FLR1*R1Y+FLR2*R2Y+FLZA*AY
      FLZ=FLR1*R1Z+FLR2*R2Z+FLZA*AZ
      BXB1=YMYC*AZ-ZMZC*AY
      BYB1=ZMZC*AX-XMXC*AZ      
      BZB1=XMXC*AY-YMYC*AX
            
c 19.05.2011      fazrda=fazrda*cur(k-1)
      ncoi=nco(k)                
      fazrda=fazrda*cur(k-1)*curco(ncoi)
c 19.05.2011 end
      
      BXB=BXB1*FAZRDA+BXB
      BYB=BYB1*FAZRDA+BYB
      BZB=BZB1*FAZRDA+BZB
      BXXB=FLX*BXB1*FAZRDA+BXXB
      BXYB=(FLY*BXB1+AZ)*FAZRDA+BXYB
      BXZB=(FLZ*BXB1-AY)*FAZRDA+BXZB
      BYYB=FLY*BYB1*FAZRDA+BYYB
      BYZB=(FLZ*BYB1+AX)*FAZRDA+BYZB
      BZZB=FLZ*BZB1*FAZRDA+BZZB
c 25.12.2008
      ByxB=(FLx*ByB1-AZ)*FAZRDA+ByxB
      BzxB=(FLx*BzB1+Ay)*FAZRDA+BzxB
      BzyB=(FLy*BzB1-Ax)*FAZRDA+BzyB
c 25.12.2008 end
      
  101 continue

      R1=R2
      R1X=R2X
      R1Y=R2Y
      R1Z=R2Z
      OR1=OR2
  100 CONTINUE


      BX=BXB
      BY=BYB
      BZ=BZB
      BXX=BXXB
      BXY=BXYB
      BXZ=BXZB
      BYY=BYYB
      BYZ=BYZB
      BZZ=BZZB
c 25.12.2008
      byx=ByxB
      bzx=BzxB
      bzy=BzyB
c 25.12.2008 end
  120 CONTINUE
      cosf2=cosf*cosf
      sinf2=sinf*sinf
      sicof=sinf*cosf
      br=bx*cosf+by*sinf
      bri=br
      bfb=bfrt/rd
      bf=by*cosf-bx*sinf
      bfi=bf+bfb
c 25.12.2008      bxy2sc=bxy*sicof*2d0
      bxybyx=(bxy+byx)*sicof
      byybxx=(byy-bxx)*sicof
      brri=bxx*cosf2+byy*sinf2+bxybyx
      bfr=byx*cosf2-bxy*sinf2+byybxx
      bfri=bfr-bfb/rd
      bzri=bzx*cosf+bzy*sinf
      brfi=(bxy*cosf2-byx*sinf2+byybxx)*rd+bf
      bffi=(byy*cosf2+bxx*sinf2-bxybyx)*rd-br
      bzfi=(bzy*cosf-bzx*sinf)*rd
      brzi=bxz*cosf+byz*sinf
      bfzi=byz*cosf-bxz*sinf
c 25.12.2008 end
      BZI=BZ
      BZZI=BZZ

      RETURN
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      CALL datw7xm(XO,YO,ZO,cur,nco,nparx)
      SUBROUTINE datw7xm(XO,YO,ZO,cur,nco,neli)
      implicit real*8(a-h,o-z),integer(i-n)

      dimension XO(neli),YO(neli),ZO(neli),cur(neli)
c 25.07.2012      parameter (npar=10850,ncoil=27)
      parameter (npar=10850,ncoil=16)
      character*24 file2,file3
c 25.07.2012      data file2/'cur_it.dd'/
      data file2/'cur_asd.dd'/
c 25.07.2012      data file3/'co_itm.dd'/
      data file3/'co_asd.dd'/
      
      dimension curco(ncoil),nco(neli)
      common/datwc/nnodc           
c 19.05.2011
      common/datcur/curco

      open(1, file=file3, status='old')   
c 19.05.2011
      open(2, file=file2, status='old')           
      read(2,*) curco
c 19.05.2011 end
      read(1,*) nnod
      nnodc=nnod
      
      if(nnod.gt.neli) then      
      print*,nnod,'=nnod',neli,'=neli'
      print*,'nnod.gt.neli, stop'     
      stop
      end if
      
      do knod=1,nnod
c 19.05.2011      read(1,*) xo(knod),yo(knod),zo(knod),cur(knod)
      read(1,*) xo(knod),yo(knod),zo(knod),cur(knod),nco(knod)
      end do
      cbc=cur(nnod)
      

      if(cbc.ne.0d0) then
      print*,cbc,'=cbc.ne.0, stop'
      stop
      end if
      
      print*,nnod,'=nnod,  datw7x terminated '
      close(1)
c 25.07.2012
      close(2)
c 25.07.2012 end
      RETURN
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine magfie(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
c
c Computes magnetic field module in units of the magnetic code  - bmod,
c square root of determinant of the metric tensor               - sqrtg,
c derivatives of the logarythm of the magnetic field module
c over coordinates                                              - bder,
c covariant componets of the unit vector of the magnetic
c field direction                                               - hcovar,
c contravariant components of this vector                       - hctrvr,
c contravariant component of the curl of this vector            - hcurl
c Order of coordinates is the following: x(1)=R (big radius), 
c x(2)=phi (toroidal angle), x(3)=Z (altitude).
c
c  Input parameters:
c            formal:  x                -    array of coordinates
c  Output parameters:
c            formal:  bmod
c                     sqrtg
c                     bder
c                     hcovar
c                     hctrvr
c                     hcurl
c
c  Called routines:  GBhs,GBRZd 
c
      double precision x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl
      double precision hr,hf,hz
c
      double precision ri,fii,zi,br,bf,bz,
     *BRR,BRF,BRZ,BFR,BFF,BFZ,BZR,BZF,BZZ,
     *BRK,BZK,BRRK,BRZK,BZRK,BZZK
c
      dimension x(3),bder(3),hcovar(3),hctrvr(3),hcurl(3)
c
      rbig=max(x(1),1d-12)
c
ccccccc computation of gb in cylindrical co-ordinates cccccccc
      ri=rbig
      fii=x(2)
      zi=x(3)

      CALL gbcoil(ri,fii,zi,br,bf,bz,
     *BRR,BRF,BRZ,BFR,BFF,BFZ,BZR,BZF,BZZ)

cccccc end of gb computation cccccccccc
      bmod=dsqrt(br**2+bf**2+bz**2)
      sqrtg=rbig
      hr=br/bmod
      hf=bf/bmod
      hz=bz/bmod
c
      bder(1)=(brr*hr+bfr*hf+bzr*hz)/bmod
      bder(2)=(brf*hr+bff*hf+bzf*hz)/bmod
      bder(3)=(brz*hr+bfz*hf+bzz*hz)/bmod
c
      hcovar(1)=hr
      hcovar(2)=hf*rbig
      hcovar(3)=hz
c
      hctrvr(1)=hr
      hctrvr(2)=hf/rbig
      hctrvr(3)=hz
c
      hcurl(1)=((bzf-rbig*bfz)/bmod+
     +          hcovar(2)*bder(3)-hcovar(3)*bder(2))/sqrtg
      hcurl(2)=((brz-bzr)/bmod+
     +          hcovar(3)*bder(1)-hcovar(1)*bder(3))/sqrtg
      hcurl(3)=((bf+rbig*bfr-brf)/bmod+
     +          hcovar(1)*bder(2)-hcovar(2)*bder(1))/sqrtg
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
