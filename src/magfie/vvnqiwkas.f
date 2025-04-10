!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE tj2vvo(RT0,R0i,L1i,cbfi,BY0i,bf0)
      double precision rt0,r0i,cbfi,by0i,bf0
! 22.09.2010
!      real*8 facr
!      facr=0.01d0
!      facr=1.0d0
! 22.09.2010 end
!c
      INTEGER BMMA,BNC
      LOGICAL NOPRI1,NOPRI2
      double precision RC(6),KBZ,rt,r0,by0
! 08.01.2011
      real*8 facr
      facr=0.01d0
!      facr=1.0d0
! 08.01.2011 end
      BNC=1
      NOPRI1=.false.
      NOPRI2=.FALSE.

      open(17,form='FORMATTED',file='pchsm.d')
      READ(17,340)nopri1,nopri2
  340 FORMAT(7X,L5,6X,7X,L5)

      READ(17,120)NMAX,MMA,NMA,RT,R0,L1,NC
!  120 FORMAT(5X,I2,2X,4X,I3,2X,5X,I2,2X,2(3X,E10.4,2X),2(3X,I3,2X))
  120 FORMAT(4X,I3,2X,4X,I3,2X,5X,I2,2X,2(3X,E10.4,2X),2(3X,I3,2X))
      IF(NC.GT.BNC) NC=BNC
!      PRINT 130,NMAX,BMMA,MMA,NMA,L1,RT,R0,NC
  130 FORMAT(//5HNMAX=,I3,4X,4X,5HBMMA=,I3,4X// &
      4HMMA=,I3,4X,4HNMA=,I3,4X,3HL1=,I3,4X, &
      3HRT=,E14.7,4X,3HR0=,E14.7,4X,3HNC=,I3/)
      READ(17,420)KBZ,RC
  420 FORMAT(5X,E13.7/3(3X,E13.7)/3(3X,E13.7))
      READ(17,180)BY0
  180 FORMAT(4X,E14.7)
      PRINT 250,BY0,KBZ
  250 FORMAT(/4HBY0=,E14.7,3X,4HKBZ=,E14.7/)
      PRINT 480,RC
  480 FORMAT(/3HRC=,6(E14.7,2X))
! 20.09.2010      CALL COINhs(RT,R0,RC,NC,BMMA,MMA,NMA,NMAX,
! 20.09.2010     *KBZ,L1,NOPRI1,NOPRI2)

      close(17)

      RT0=rt/facr
      R0i=r0/facr
      L1i=l1
      cbf=kbz/facr
      cbfi=cbf
      BY0i=by0
      bf0=-cbf*L1/rt0

      return
      END
!ccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE gbtj2(rI,fI,ZI,BrI,BfI,BZI, &
      BrrI,BrfI,BrZI,BfrI,BffI,BfZI,BZrI,BZfI,BZZI)
!
      use magfield_mod, only : ierrfield
!
      implicit real*8(a-h,o-z),integer(i-n)

      COMMON/CRKI1e/ PI,XGAC,RTC,BY0,R0C
! 20.09.2010      save by0c,ier
      save ier,facr
      data ier/0/

      if(ier.eq.1) go to 10
! 20.09.2010      by0c=by0
      facr=0.01d0
!      facr=1.0d0
!      ofacr=1d0/facr
      open(36,form='FORMATTED',file='fort.36')

   10 continue
! 31.08.2008
      figb=fi
! 31.08.2008 end
      x=ri*facr
      y=zi*facr

      CALL field(X,figb,Y,BX,BF,BY, &
      BRR,BRF,BRZ,BFR,BFF,BFZ,BZR,BZF,BZZ)
      if(ierrfield.eq.1) return
      if(ier.eq.0) then
      write(36,*)facr,'=facr'
      write(36,*)x,'=r_st',fi,'=fi_st',y,'=z_st'
      write(36,*)bx,'=br',bf,'=bf',by,'=bz, b_symm'
      end if
! 18.11.2008      CALL GBas(X,figb,Y,BXa,BFa,BYa,
! 18.11.2008     *BRRa,BRFa,BRZa,BFRa,BFFa,BFZa,BZRa,BZFa,BZZa)
! 18.11.2008        bx=bx+bxa
        bri=bx
! 18.11.2008        bf=bf+bfa
        bfi=bf
! 18.11.2008        by=by+bya
        bzi=by
        brri=brr*facr
        brfi=brf
        brzi=brz*facr
        bfri=bfr*facr
        bffi=bff
        bfzi=bfz*facr
        bzri=bzr*facr
        bzfi=bzf
        bzzi=bzz*facr
      if(ier.eq.0) then
! 25.03.2010      write(36,*)bxa,'=br',bfa,'=bf',bya,'=bz, b_asym'
      write(36,*)bx,'=br',bf,'=bf',by,'=bz'
      write(36,*)bzi,'=bzi'

      close(36)
      ier=1
      end if
! 30.08.2008 end

      RETURN
      END
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! .rkin1

