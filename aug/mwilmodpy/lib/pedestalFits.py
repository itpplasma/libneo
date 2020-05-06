#!/usr/bin/env python

import numpy as np
from scipy.optimize import leastsq
from IPython import embed
from scipy import odr


def twoLineFunc( params, x ):
    a0 = params[0]
    a1 = params[1]
    a2 = params[2]
    a3 = params[3]
    f=np.zeros_like(x)
    f[x<= a0] = a2*(a0-x[x<= a0])+a1
    f[x > a0] = a3*(x[x > a0]-a0)+a1
    return f

def twoLineDerivFunc( params, x ):
    a0 = params[0]
    a1 = params[1]
    a2 = params[2]
    a3 = params[3]
    df=np.zeros_like(x)
    df[x<= a0] = -a2
    df[x > a0] = a3
    return df

def twoLineFuncResi( params,x, f ):
    return  f - twoLineFunc( params, x )  


def mtanhFunc( paramIn, x ):

    params = np.zeros((10),dtype='float') 
    params[:paramIn.size] = paramIn

    a0 = params[0]
    a1 = params[1]
    a2 = params[2]
    a3 = params[3]
    a4 = params[4]
    a5 = params[5]
    a6 = params[6]
    a7 = params[7]
    a8 = params[8]
    a9 = params[9]    
 

    z = (a2-x)/a3
    p1 = (1+a4*z+a5*z**2.+a6*z**3.)
    p2 = (1+a7*z+a8*z**2.+a9*z**3.)
    f = (a0+a1)/2.+ (a1-a0)/2.*(p1*np.exp(z)-p2*np.exp(-z))/(np.exp(z)+np.exp(-z))

    return f

def mtanhDerivFunc( paramIn, x ):

    params = np.zeros((10),dtype='float') 
    params[:paramIn.size] = paramIn

    a0 = params[0]
    a1 = params[1]
    a2 = params[2]
    a3 = params[3]
    a4 = params[4]
    a5 = params[5]
    a6 = params[6]
    a7 = params[7]
    a8 = params[8]
    a9 = params[9]    
   
    z = (a2-x)/a3
    p1 = (1+a4*z+a5*z**2.+a6*z**3.)
    p2 = (1+a7*z+a8*z**2.+a9*z**3.)

    zp = -1./a3
    p1p = (a4 + 2.*a5*z+ 3.*a6*z**2.)*zp
    p2p = (a7 + 2.*a8*z+ 3.*a9*z**2.)*zp

    df =  (a1-a0)/2.*(  (p1p*np.exp(z)-p2p*np.exp(-z) )*(np.exp(z)+np.exp(-z))+ 2.*zp*(p1+p2)  )/(np.exp(z)+np.exp(-z))**2.

    return df





def mtanhFuncResi( params, x, f, boundCond=0,  weight=1., xOut=1.15, fOut = 0.0, xIn = 0.85, fIn= 10000. ):
    bound = np.array(0.0)
    #print np.shape(x),f - mtanhFunc( params, x )
    if boundCond==0:
        return f - mtanhFunc( params, x )
    elif boundCond==1:
        bound =( fOut - mtanhFunc( params, xOut ) )
    elif boundCond==2:
        bound = np.abs( fOut - mtanhFunc( params, xOut ) ) + np.abs( fIn - mtanhFunc( params, xIn ) )
    elif boundCond==3:
        bound = ( fOut - mtanhDerivFunc( params, xOut ) )
                 ##f(xout)=0.0,df/dx(xout)=0.0
    elif boundCond==4:
        bound =  np.abs(mtanhFunc( params, xOut )) + np.abs(mtanhDerivFunc( params, xOut ))
    elif boundCond==5:
        if (mtanhDerivFunc( params, xIn ) <= 0):
            add = 0.0
        else:
            add = 1.e6
        bound =  np.abs(mtanhFunc( params, xOut )) + np.abs(mtanhDerivFunc( params, xOut )) + add
        ### f(out)=0.0,df/dx(out) = 0, df/dx(in) = 0
    elif boundCond==6:
        bound =  np.abs(mtanhFunc( params, xOut )) + np.abs(mtanhDerivFunc( params, xOut )) + np.abs(mtanhDerivFunc( params, xIn ))

    return np.abs(f - mtanhFunc( params, x )) + np.abs(bound)
    #print mtanhFunc( params, xOut )
    #if np.abs(bound) < 10.:
    #    return  f - mtanhFunc( params, x ) 
    #else:
    #    Out = np.zeros_like(x)
    #    Out[:]=1.e8
    #    return Out

class Fit:

    def __init__( self, xIn, yIn, initCondIn=None, twoLine=False, mtanh=False):
        self.Status = False  
        
        if ((twoLine == False) & (mtanh == False)):
            twoLine=True

        if twoLine == False:
            mtanh=True
        else:
            mtanh=False
        
        if twoLine:
            nCond=4
        else:
            nCond=10
        
        self.twoLine = twoLine
        self.mtanh = mtanh

        if np.all(initCondIn)==None:
            self.initCond = np.zeros((nCond))
        elif(np.size(initCond)!=4):
            self.initCond = np.zeros((nCond))
        else:
            self.initCond = initCondIn


        if np.size(xIn)!=np.size(yIn):
            self.Status = False
            print 'not equal sized input'

        xMean = xIn.mean()
            #do initial conditions
        if twoLine:
            self.initCond[0] = xMean
            self.initCond[1] = np.interp(xMean,xIn,yIn)
            self.initCond[2] = np.abs(np.mean(np.diff(yIn[xIn <= xMean])/np.diff(xIn[xIn <= xMean])))
            self.initCond[3] = np.abs( np.mean( np.diff(yIn[xIn > xMean])/np.diff(xIn[xIn > xMean])))
        else:
            self.initCond[0] = yIn.min()
            self.initCond[1] = yIn.mean()+yIn.std()
            self.initCond[2] = xMean
            self.initCond[3] = xIn.std()

        self.xData = xIn
        self.yData = yIn
        self.doFit()
        self.Status = True
 
    def __del__( self ):

        del self.initCond

        if self.Status:
            del self.xData
            del self.yData
            del self.params
            del self.error

        del self.Status

    def doFit(self):
        
        if self.twoLine:
            fitOut=leastsq(twoLineFuncResi,self.initCond ,args=(self.xData, self.yData),full_output=1)
        else:
            if ((self.xData.mean()+self.xData.std() > 1.0) | np.size(np.where(self.xData>1.0)[0])):
                fitOut=leastsq(mtanhFuncResi,self.initCond ,args=(self.xData, self.yData),full_output=1)
            else:
                #embed()
                #myFunc = odr.Model(mtanhFunc)
                #myData = odr.Data(self.xData, self.yData)
                #myRun = odr.ODR(data, func, beta0=self.initCond)
                #myOut = myRun.run()
                #self.params = myOut.beta
                #myRun = odr.ODR(data, func, beta0=self.params)
                #myOut = myRun.run()
                ##no SOL data ot at least 5 data points must lie in the SOL
                fitOut=leastsq(mtanhFuncResi,self.initCond[:5] ,args=(self.xData, self.yData,4),full_output=1)
                

        self.params = fitOut[0]
       
        if self.twoLine:
            s_sq = (twoLineFuncResi(self.params,self.xData,self.yData)**2).sum()/(self.yData.size-self.params.size)
        else:
            s_sq = (mtanhFuncResi(self.params,self.xData,self.yData)**2).sum()/(self.yData.size-self.params.size)

        
        if np.all(fitOut[1]) != None:
            pcov = fitOut[1] * s_sq

            error = np.zeros_like(self.params) 
            for i in np.arange(self.params.size):
                error[i] = (np.absolute(pcov[i][i])**0.5)

            self.error = error
        else:
            self.error = np.zeros_like(self.params)

    def getFunc(self,xIn):
        if self.Status == True:
            if self.twoLine:
                return twoLineFunc(self.params,xIn)
            else:
                return mtanhFunc(self.params,xIn)

    def deriv(self,xIn):
        if self.Status == True:
            if self.twoLine:
                return twoLineDerivFunc(self.params,xIn)
            else:
                return mtanhDerivFunc(self.params,xIn)

    def __call__( self,xIn):
        return self.getFunc(xIn)









"""
;+
; NAME:
;       PROF_MTANH
; PURPOSE:
;       EVALUATE A MODIFIED HYPERBOLIC TANGENT FUNCTION AND
;       OPTIONALLY RETURN THE VALUE OF ITS PARTIAL DERIVATIVES.
;       NORMALLY, THIS FUNCTION IS USED BY CURVEFIT TO FIT THE
;       EDGE PEDESTAL PROFILES IN H-MODES.  BASED ON THE IDL
;       ROUTINE 'FUNCT'.
;
; CATEGORY:
;       E2 - CURVE AND SURFACE FITTING.
; CALLING SEQUENCE:
;       PROF_MTANH,X,A,F,PDER
; INPUTS:
;       X = VALUES OF INDEPENDENT VARIABLE.
;       A = PARAMETERS OF EQUATION DESCRIBED BELOW.
; OUTPUTS:
;       F = VALUE OF FUNCTION AT EACH X(I).
;
; OPTIONAL OUTPUT PARAMETERS:
;       PDER = (N_ELEMENTS(X),10) ARRAY CONTAINING THE
;               PARTIAL DERIVATIVES.  P(I,J) = DERIVATIVE
;               AT ITH POINT W/RESPECT TO JTH PARAMETER.
; COMMON BLOCKS:
;       MTANH: for passing parameters
;           NPARA: number of free parameters - must be consistent
;                  with boundary conditions (not checked (yet?))
;          WHPARA: indices of free parameters - must be consistent
;                  with boundary conditions (not checked (yet?))
;           INORD: order of polynomial for fitting inside (lower x)
;                  of pedestal region 0<INORD<3
;            INBC: 0: free
;                  1: df/dx(x=0) = 0 (assumes a(2)/a(3)>>1)
;          OUTORD: order of polynomial for fitting outside (higher x)
;                  of pedestal region 0<OUTORD<3
;           OUTBC: 0: free
;                  1: f(x=xmax) = 0     (assumes zmax=(a(2)-xmax)/a(3)<<-1)
;                  2: df/dx(x=xmax) = 0 (")
;                  3: f(x=xoutbc) = foutbc by fixing a(2), the symmetry
;                     point
;          XOUTBC: see OUTBC=3 above
;          FOUTBC: set OUTBC=3 above
;               A: vector for returning entire set of fit parameters,
;                  including fixed & derived values
;       OUTBC: for passing parameters into the root function
;          XOUTBCFUNC : the x-value at which the boundary condition is
;                       applied (copy of MTANH.XOUTBC)
;          FOUTBCFUNC : the function value at the boundary position
;                       (copy of MTANH.FOUTBC)
;          AOUTBCFUNC : full set of parameters for the function
;                       evaluation (AOUTBCFUNC[2] replaced as required)
; SIDE EFFECTS:
;       NONE.
; RESTRICTIONS:
;       NONE.
; PROCEDURE:
;       F = (A(0)+A(1))/2 + (A(1)-A(0))/2*MTANH(A(4-9),Z) (HERE A(0),A(1) +VE)
;       Z = (A(2)-X)/A(3)
;       MTANH(A(4-9),Z) = [(1+A(4)*Z+A(5)*Z^2+A(6)*Z^3)*EXP( Z) $
;                         -(1+A(7)*Z+A(8)*Z^2+A(9)*Z^3)*EXP(-Z)] $
;                         /[EXP(Z)+EXP(-Z)]
; MODIFICATION HISTORY:
;       Created as a modification of MTANHFIT, 8/4/01 LDH
;       Add constraint at position in pedestal region, 23/7/02 LDH
;       Modified to pass parameters as a structure, 17/9/02 LDH
;-

FUNCTION outbcfunc, a2

COMMON OUTBC, xoutbcfunc, foutbcfunc, aoutbcfunc

    big = 1.0e1
    nx  = 1L

    a = aoutbcfunc

; Loop over a vector of values for a2
    na2 = N_ELEMENTS(a2)
    f = DBLARR(na2)
    FOR i = 0,na2-1 DO BEGIN

; Replace a[2]
      a[2] = a2[i]

      z = (a(2)-xoutbcfunc)/a(3)
      whok = WHERE(z LE +big AND z GE -big)
      whpl = WHERE(z GT +big)
      whmn = WHERE(z LT -big)
      IF whok(0) NE -1 THEN zok = z(whok)
      IF whpl(0) NE -1 THEN zpl = z(whpl)
      IF whmn(0) NE -1 THEN zmn = z(whmn)

      p1    = 1.0 + a(4)*z + a(5)*z^2 + a(6)*z^3
      p2    = 1.0 + a(7)*z + a(8)*z^2 + a(9)*z^3
      numer = DBLARR(nx) & denom = DBLARR(nx)
      IF whok(0) NE -1 THEN BEGIN
        numer(whok) = p1(whok)*exp(zok)-p2(whok)*exp(-zok)
        denom(whok) = exp(zok)+exp(-zok)
      ENDIF
      IF whpl(0) NE -1 THEN BEGIN
        numer(whpl) = p1(whpl)
        denom(whpl) = 1.0
      ENDIF
      IF whmn(0) NE -1 THEN BEGIN
        numer(whmn) = -p2(whmn)
        denom(whmn) = 1.0
      ENDIF
      mtanh = numer / denom
      f[i]=0.5*( (a(0)+a(1)) + (a(1)-a(0))*mtanh )
    ENDFOR

    f = f-foutbcfunc
    IF na2 EQ 1 THEN f = f[0]
    
    RETURN, f
END


PRO prof_mtanh,X,APARA,F,pderp,dY=dY,PRINT=PRINT

COMMON MTANH, mtanh_params
COMMON OUTBC, xoutbcfunc, foutbcfunc, aoutbcfunc

    ON_ERROR,0
    IF TOTAL(FINITE(apara)) LT N_ELEMENTS(apara) THEN STOP

; local copies of parameters
    npara  = mtanh_params.npara
    whpara = mtanh_params.whpara[0:npara-1]
    inord  = mtanh_params.inord
    inbc   = mtanh_params.inbc
    outord = mtanh_params.outord
    outbc  = mtanh_params.outbc
    xoutbc = mtanh_params.xoutbc
    foutbc = mtanh_params.foutbc
    width  = mtanh_params.width
    a      = mtanh_params.a

    type = size(apara)
    type = type[type[0]+1]
    double = type EQ 5
    IF double THEN big = 1.0D1 ELSE big = 1.0E1
    res = machar(DOUBLE=double)
    eps = sqrt(res.eps)
    nx = N_ELEMENTS(x)

    FOR i=0,npara-1 DO BEGIN
      a(whpara(i)) = APARA(i)
    ENDFOR

; for fixed boundary conditions at xmax, set highest order coefficient
; - here we're assuming zmax = (a(2)-xmax)/a(3) << -1
    IF outbc EQ 1 THEN BEGIN ; f(xoutbc) = foutbc
      xmax = xoutbc
      zmax = (a(2)-xmax)/a(3)
      IF outord EQ 1 THEN $
        a(7) = 2.*(a(0)-foutbc)/(zmax  *(a(1)-a(0))) $
      ELSE IF outord EQ 2 THEN $
        a(8) = 2.*(a(0)-foutbc)/(zmax^2*(a(1)-a(0))) - a(7)/zmax $
      ELSE IF outord EQ 3 THEN $
        a(9) = 2.*(a(0)-foutbc)/(zmax^3*(a(1)-a(0))) - a(7)/zmax^2 - a(8)/zmax
    ENDIF ELSE IF outbc EQ 2 THEN BEGIN ; df/dx(xoutbc) = 0
      xmax = xoutbc
      zmax = (a(2)-xmax)/a(3)
      IF outord EQ 1 THEN $
        a(7) = 0.0 $
      ELSE IF outord EQ 2 THEN $
        a(8) = -a(7)/(2.*zmax  ) $
      ELSE IF outord EQ 3 THEN $
        a(9) = -a(7)/(3.*zmax^2) - 2.*a(8)/(3.*zmax)
    ENDIF ELSE IF outbc EQ 3 THEN BEGIN ; f(x) = f
      a[0] = 0.0    ; set these here so they are available for inbc=1
      IF SIZE(aoutbcfunc,/TYPE) EQ 0 THEN $
        a[2] = 1.0 $
      ELSE IF aoutbcfunc[2] NE 0.0 THEN $
        a[2] = aoutbcfunc[2] $
      ELSE $
        a[2] = 1.0
      a[3] = width
    ENDIF

; for zero slope on axis, set highest order coefficient
; If outbc=3, we don't know a[2] yet, use value from last
; iteration (and start with 1.0 for first iteration)

    IF inbc THEN $
      IF inord EQ 3 THEN $
        a(6) = -a(3)/(3.*a(2))*(a(3)*a(4)/a(2)+2.*a(5)) $
      ELSE IF inord EQ 2 THEN $
        a(5) = -a(3)*a(4)/(2.*a(2)) $
      ELSE IF inord EQ 1 THEN $
        a[4] = 0.0

; for fixed boundary conditions in pedestal region:
; set the pedestal symmetry point a[2] using the b.c.
; and fix the SOL value a[0] to 0.0
; In addition: a special which also sets the pedestal
; width a[3] to mtanh_params.width
; a[2] must be set last because we need all the other a's
    IF outbc EQ 3 THEN BEGIN  ; f(xoutbc) = foutbc and a[3] = width
      a[0] = 0.0
      a[3] = width
      xoutbcfunc = xoutbc
      foutbcfunc = foutbc
      aoutbcfunc = a
      XROOT = [xoutbc,xoutbc-3.0*a[3],xoutbc+3.0*a[3]]
;      XROOT = [xoutbc,0.0,xoutbc+3.0*a[3]]
;      IF outbcfunc(xroot[1])*outbcfunc(xroot[2]) GT 0.0 THEN BEGIN
;        print,xroot,outbcfunc(xroot)
;        stop
      IF a[1] LT 2.0*foutbc THEN $  ; relax boundary condition
        a[2] = 1.0 $
      ELSE $
        a[2] = FX_ROOT(XROOT,'OUTBCFUNC',/DOUBLE)
; Update aoutbcfunc (for partial derivatives)
      aoutbcfunc[2] = a[2]
    ENDIF

; evaluate function: negative x --> f=0
    IF double THEN f = DBLARR(nx) ELSE f = FLTARR(nx)
    whpos = WHERE(x GE 0.0,npos)
    IF npos EQ 0 THEN STOP,'No positive x values in prof_mtanh!'

    z=(a(2)-x[whpos])/a(3)
    whok = WHERE(z LE +big AND z GE -big)
    whpl = WHERE(z GT +big)
    whmn = WHERE(z LT -big)
    IF whok(0) NE -1 THEN zok = z(whok)
    IF whpl(0) NE -1 THEN zpl = z(whpl)
    IF whmn(0) NE -1 THEN zmn = z(whmn)

    p1    = 1.0 + a(4)*z + a(5)*z^2 + a(6)*z^3
    p2    = 1.0 + a(7)*z + a(8)*z^2 + a(9)*z^3
    IF double THEN numer = DBLARR(npos) ELSE numer = FLTARR(npos)
    IF double THEN denom = DBLARR(npos) ELSE denom = FLTARR(npos)
    IF whok(0) NE -1 THEN BEGIN
      numer(whok) = p1(whok)*exp(zok)-p2(whok)*exp(-zok)
      denom(whok) = exp(zok)+exp(-zok)
    ENDIF
    IF whpl(0) NE -1 THEN BEGIN
      numer(whpl) = p1(whpl)
      denom(whpl) = 1.0
    ENDIF
    IF whmn(0) NE -1 THEN BEGIN
      numer(whmn) = -p2(whmn)
      denom(whmn) = 1.0
    ENDIF
    mtanh = numer / denom
    f[whpos] = 0.5*( (a(0)+a(1)) + (a(1)-a(0))*mtanh )

    IF KEYWORD_SET(dY) THEN BEGIN
      IF double THEN dY = DBLARR(nx) ELSE dY = FLTARR(nx)
      dp1 = a(4) + 2.0*a(5)*z + 3.0*a(6)*z^2
      dp2 = a(7) + 2.0*a(8)*z + 3.0*a(9)*z^2
      IF double THEN dmdz = DBLARR(npos) ELSE dmdz = FLTARR(npos)
      IF whok(0) NE -1 THEN $
        dmdz(whok) = (dp1(whok)*exp(2.0*zok) $
         + 2.0*(p1(whok)+p2(whok))+dp1(whok)-dp2(whok) $
         - dp2(whok)*exp(-2.0*zok)) $
        / denom(whok)^2
      IF whpl(0) NE -1 THEN $
        dmdz(whpl) =  dp1(whpl)
      IF whmn(0) NE -1 THEN $
        dmdz(whmn) = -dp2(whmn)
      dY[whpos] = -0.5*(a(1)-a(0))/a(3)*dmdz
    ENDIF

    IF N_PARAMS(0) GT 3 THEN BEGIN   ;NEED PARTIAL?
      IF double THEN pderp = DBLARR(nx,npara) $
                ELSE pderp = FLTARR(nx,npara)
      dp1 = a(4) + 2.0*a(5)*z + 3.0*a(6)*z^2
      dp2 = a(7) + 2.0*a(8)*z + 3.0*a(9)*z^2
      IF double THEN dmdz = DBLARR(npos) ELSE dmdz = FLTARR(npos)
      IF whok(0) NE -1 THEN $
        dmdz(whok) = (dp1(whok)*exp(2.0*zok) $
         + 2.0*(p1(whok)+p2(whok))+dp1(whok)-dp2(whok) $
         - dp2(whok)*exp(-2.0*zok)) $
        / denom(whok)^2
      IF whpl(0) NE -1 THEN $
        dmdz(whpl) =  dp1(whpl)
      IF whmn(0) NE -1 THEN $
        dmdz(whmn) = -dp2(whmn)
      IF double THEN pder = DBLARR(npos,10) $
                ELSE pder = FLTARR(npos,10)
      pder(*,0:3) = [ [ 0.5*(1.-mtanh) ], $
    	              [ 0.5*(1.+mtanh) ], $
                      [ 0.5*(a(1)-a(0))/a(3)*dmdz ], $
    	              [-0.5*(a(1)-a(0))/a(3)*z*dmdz ] ]
      IF whok(0) NE -1 THEN $
        pder(whok,4:9) = 0.5*(a(1)-a(0)) * $
                         [ [ zok  *exp(+zok)/denom(whok) ], $
                           [ zok^2*exp(+zok)/denom(whok) ], $
                           [ zok^3*exp(+zok)/denom(whok) ], $
                           [-zok  *exp(-zok)/denom(whok) ], $
                           [-zok^2*exp(-zok)/denom(whok) ], $
                           [-zok^3*exp(-zok)/denom(whok) ] ]
      IF whpl(0) NE -1 THEN $
        pder(whpl,4:6) = 0.5*(a(1)-a(0)) * $
                         [ [ zpl   ], $
                           [ zpl^2 ], $
                           [ zpl^3 ] ]
      IF whmn(0) NE -1 THEN $
        pder(whmn,7:9) = 0.5*(a(1)-a(0)) * $
                         [ [-zmn   ], $
                           [-zmn^2 ], $
                           [-zmn^3 ] ]

      IF outbc EQ 1 THEN BEGIN              ; f(xoutbc)=foutbc
        IF outord EQ 1 THEN BEGIN
          pder(*,0) = pder(*,0) $
            + pder(*,7)*  2.*(a(1)-foutbc)/(zmax*(a(1)-a(0))^2)
          pder(*,1) = pder(*,1) $
            + pder(*,7)*(-2.*(a(0)-foutbc))/(zmax*(a(1)-a(0))^2)
          pder(*,2) = pder(*,2) $
            + pder(*,7)*(2.*(a(0)-foutbc))/(zmax^2*(a(1)-a(0))*a(3))
          pder(*,3) = pder(*,3) $
            + pder(*,7)*(-2.*(a(0)-foutbc)*(a(2)-xmax)) $
                       /(zmax^2*(a(1)-a(0))*a(3)^2)
        ENDIF ELSE IF outord EQ 2 THEN BEGIN
          pder(*,0) = pder(*,0) $
            + pder(*,8)*  2.*(a(1)-foutbc)/(zmax^2*(a(1)-a(0))^2)
          pder(*,1) = pder(*,1) $
            + pder(*,8)*(-2.*(a(0)-foutbc))/(zmax^2*(a(1)-a(0))^2)
          pder(*,2) = pder(*,2) $
            + pder(*,9)*(-4.*(a(0)-foutbc)/(zmax^3*(a(1)-a(0))) $
                         +a(7)/zmax^2)/a(3)
          pder(*,3) = pder(*,3) $
            + pder(*,9)*(-4.*(a(0)-foutbc)/(zmax^3*(a(1)-a(0))) $
                         +a(7)/zmax^2)*(-1.*(a(2)-xmax))/a(3)^2
          pder(*,7) = pder(*,7) $
            + pder(*,8)*(-1./zmax)
        ENDIF ELSE IF outord EQ 3 THEN BEGIN
          pder(*,0) = pder(*,0) $
            + pder(*,9)*  2.*(a(1)-foutbc)/(zmax^3*(a(1)-a(0))^2)
          pder(*,1) = pder(*,1) $
            + pder(*,9)*(-2.*(a(0)-foutbc))/(zmax^3*(a(1)-a(0))^2)
          pder(*,2) = pder(*,2) $
            + pder(*,9)*(-6.*(a(0)-foutbc)/(zmax^4*(a(1)-a(0))) $
                         +2.*a(7)/zmax^3+a(8)/zmax^2)/a(3)
          pder(*,3) = pder(*,3) $
            + pder(*,9)*(-6.*(a(0)-foutbc)/(zmax^4*(a(1)-a(0))) $
                         +2.*a(7)/zmax^3+a(8)/zmax^2) $
                       *(-1.*(a(2)-xmax))/a(3)^2
          pder(*,7) = pder(*,7) $
            + pder(*,9)*(-1./zmax^2)
          pder(*,8) = pder(*,8) $
            + pder(*,9)*(-1./zmax)
        ENDIF
      ENDIF ELSE IF outbc EQ 2 THEN BEGIN   ; df/dx(xmax)=0
        IF outord EQ 2 THEN BEGIN
          pder(*,2) = pder(*,2) $
            + pder(*,8)*(a(7)/(2.*a(3)*zmax^2))
          pder(*,3) = pder(*,3) $
            + pder(*,8)*(-1.*a(7)*(a(2)-xmax)/(2.*a(3)^2*zmax^2))
          pder(*,7) = pder(*,7) $
            + pder(*,8)*(-1./(2.*zmax))
        ENDIF ELSE IF outord EQ 3 THEN BEGIN
          pder(*,2) = pder(*,2) $
            + pder(*,9)*(2./(3.*a(3)*zmax^2)*(a(7)+a(8)/zmax))
          pder(*,3) = pder(*,3) $
            + pder(*,9)*(-2.*(a(2)-xmax)) $
                       /(3.*a(3)^2*zmax^2)*(a(7)+a(8)/zmax)
          pder(*,7) = pder(*,7) $
            + pder(*,9)*(-1./(3.*zmax^2))
          pder(*,8) = pder(*,8) $
            + pder(*,9)*(-2./(3.*zmax))
        ENDIF
      ENDIF

      IF inbc THEN BEGIN                   ; df/dx(0)=0
        IF inord EQ 3 THEN BEGIN
          pder(*,2) = pder(*,2) $
            + pder(*,6)*(2.*a(3)/(3.*a(2)^2)*(a(3)*a(4)/a(2)+a(5)))
          pder(*,3) = pder(*,3) $
            + pder(*,6)*(-2./(3.*a(2))*(a(3)*a(4)/a(2)+a(5)))
          pder(*,4) = pder(*,4) $
            + pder(*,6)*(-1.*a(3)^2/(3.*a(2)^2))
          pder(*,5) = pder(*,5) $
            + pder(*,6)*(-2.*a(3)/(3.*a(2)))
        ENDIF ELSE IF inord EQ 2 THEN BEGIN
          pder(*,2) = pder(*,2) $
            + pder(*,5)*(a(3)*a(4)/(2.*a(2)^2))
          pder(*,3) = pder(*,3) $
            + pder(*,5)*(-1.*a(4)/(2.*a(2)))
          pder(*,4) = pder(*,4) $
            + pder(*,5)*(-1.*a(3)/(2.*a(2)))
        ENDIF
      ENDIF

      IF outbc EQ 3 THEN BEGIN  ; f(xoutbc) = foutbc
        IF a[1] GE 2.0*foutbc THEN BEGIN
          FOR i=0,npara-1 DO BEGIN
            j = whpara[i]
            inc = eps*aoutbcfunc[j]
            IF inc EQ 0.0 THEN inc = eps
            aoutbcfunc[j] = aoutbcfunc[j]+inc
            XROOT = [xoutbc,xoutbc-3.0*a[3],xoutbc+3.0*a[3]]
            a2eps = FX_ROOT(XROOT,'OUTBCFUNC',/DOUBLE)
            aoutbcfunc[j] = aoutbcfunc[j]-inc
            pder[*,j] = pder[*,j] $
                       + pder[*,2]*(a2eps-aoutbcfunc[2])/inc
          ENDFOR
        ENDIF
      ENDIF

      pderp[whpos,*] = pder[*,whpara[0:npara-1]]

    ENDIF

    IF KEYWORD_SET(PRINT) THEN print, a

; Pass full parameter set back through common
    mtanh_params.a = a

    RETURN
END  ; prof_mtanh
"""
