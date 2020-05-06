"""

    Routine to calculate theta star


"""

__author__='Matthias Willensdorfer (mwillens@ipp.mpg.de)'
__date__='2 March 2015'
__version__='1.0'


#import matplotlib.pylab as plt
from numpy import *
from IPython import embed
import kk_mwillens as kk
import matplotlib.pylab as plt
import dd
from scipy import interpolate
from scipy.integrate import cumtrapz
import EQU

class thetaStar:

    def __init__( self ,  shot = None, time = None, Experiment = 'AUGD', Diagnostic = 'EQI', Edition=0 ):
        self.Status = False
        self.Status_Func = False
        self.Status_wls = False
        if ((shot != None) & (time != None)):
            self.exp = Experiment
            self.diag = Diagnostic
            self.ed = Edition
            self.shot = shot
            self.time = time
            self.Status = True


    def define(self, nrhop=1000,nPol=3601, origin='HFS'):
        if self.Status == True:
            rhopIn = arange(0.0,1.0,1.0/nrhop)
            self.origin = origin
            if self.origin=='HFS':
                tomas=True
            elif self.origin=='LFS':
                tomas=False
            magr,magz = self.get_surface_rz(rhop=rhopIn,nPolAngle=nPol,tomas=tomas)
            theta,theta_star,self.R0,self.z0 = self.tomas(magr,magz)
    # theta_star_out = ((2.*pi*theta_star.T / theta_star[:,-1])).T
        #    theta_out = ((2.*pi*theta.T / theta[:,-1])).T  
                        
            self.theta_star = (2.*pi*theta_star.T / theta_star[:,-1]).T-pi
            self.theta = (2.*pi*theta.T / theta[:,-1]).T-pi

            #embed()
            
            rhop=reshape(repeat(rhopIn[1:],nPol),(nrhop-1,nPol))
            #self.funcRBi=interpolate.RectBivariateSpline(rhop.mean(axis=1), theta.mean(axis=0), theta_star)
            
            self.func2D = []
            self.rhopFunc = rhopIn[1:]
            for i in arange(0,size(self.rhopFunc)):
                self.func2D.append( interpolate.interp1d(self.theta[i],self.theta_star[i],bounds_error=True) ) 
            self.Status_Func = True


# f=interpolate.interp2d(rhop[::10],theta[::10],theta_star[::10])          
#f=interpolate.RectBivariateSpline(equR, equz, equPSI)
    def getf_theta2thetaStar( self, rhopIn=None, thetaIn=None ):

        if self.Status_Func:
            if (rhopIn == None) | (thetaIn == None) :
                print 'no input or None input'
                return None
            if size(rhopIn) == 1:
                idx = argmin( abs( self.rhopFunc - rhopIn ) )
                return self.func2D[idx](thetaIn)
            elif size(rhopIn) > 1 :
                nPoints = size(rhopIn)
                if nPoints == size(thetaIn):                
                    thetaStar=zeros((nPoints))
                    for i in arange(nPoints):
                        idx = argmin( abs( self.rhopFunc - rhopIn[i] ) )
                        #try:
                        thetaStar[i] = self.func2D[idx](thetaIn[i])
                        #except:
                        #    embed()
                    return thetaStar
                else:
                    print 'theta and rhop not the same size'
        else:
            print 'no function defined'
            return None

    def get_theta(self, R=None, z=None):
        
        if self.Status_Func:
            if (R != None) & (z != None):
                if size(R)==size(z):
                    if self.origin=='HFS':
                        return arctan2(z-self.z0,R-self.R0) 
                    elif self.origin=='LFS':
                        return arctan2(z-self.z0,self.R0-R) 
    
                else:
                    print 'R and z input must have the same size'
                    return None
            else:
                print 'no input or None input'
                return None
        else:
            print 'no function defined'
            return None
 


    def get_surface_rz( self, rz=None, q=None, rhop=None, nPolAngle=3601,tomas=False):
   
        if self.Status == True:
            try: 
                import kk_mwillens as kk
            except:
                print 'module import not possible'
                return None
    
            rhopIn=rhop

            if (rhopIn == None) & (q == None) & (rz == None):
                print 'no input chosen'
                return None

   
            if (q != None):
        ## q is given and thop will be calculated
                print 'q it is...'
                if (size(q) < 1):
                    print 'rz must have one value at least'
                    return None
    
            #first read magnetic axis and separatrix
                output = kk.KK().kkrhorz( self.shot, self.time, [0.0,1.0],exp=self.exp, diag=self.diag )
                RmagAxis_r = output.r[0]
                RmagAxis_z = output.z[0]
                Sep_r = output.r[1]
                Sep_z = output.z[1]
                 
            #generate
                runningR =  arange( RmagAxis_r, Sep_r,0.0005)
                runningz =  repeat( RmagAxis_z, size(runningR))
            #embed()
                output = kk.KK().kkrzq( self.shot, self.time, runningR, runningz,exp=self.exp, diag=self.diag )
                f = scipy.interpolate.interp1d(-output.q,output.rhop, bounds_error=False, fill_value = float('NaN'))	
                print 'redefinition of rhop'
                rhop = f(q).tolist()

    
            elif (rz != None):
                print 'rz it is...'
                if (size(rz) != 2):
                     print 'rz must have two values'
                     return None

                output = kk.KK().kkrzptfn(self.shot, self.time, [rz[0]], [rz[1]],exp=self.exp, diag=self.diag )
                rhop =  squeeze(output.rho_p)
    
            elif (rhopIn != None):
                print 'rhop it is...'
                
                #rhopIn = array(rhopIn)
                if ( size(rhopIn) < 1) :
                    print rhopIn
                    print 'at least one rhop value'
                    return None
                elif (size(rhopIn) == 1):
                    rhop=array(rhopIn)
                else: 
                    idx=squeeze(where((rhopIn > 0.0 ) & ( rhopIn < 1.0 )))
                    if size(idx)<1:
                        print 'rhop must be between 0 and 1'
                        return None
                    else:
                        rhop=rhopIn[idx]

            if tomas:
                polAngle =  linspace(0.,360.,nPolAngle)-180.
            else:
                polAngle =  linspace(0.,360.,nPolAngle)
#arange(nPolAngle)*360./float(nPolAngle)-180.

            #get rz of one surface
            #print 'get rz of surface of ', rhop
            if size(rhop) == 1: 
                rhop = [rhop]
        
            output = kk.KK().kkrhorz( self.shot, self.time, rhop, angle = polAngle,exp=self.exp, diag=self.diag )
            surf_r =  squeeze(output.r)
            surf_z =  squeeze(output.z)
    #  out = kk.KK().kkrhopto( shot, time, [rhop],,exp=exp_eq, diag=dia_eq)
            return surf_r,surf_z
           


#routine from you should use KKRHORZ to find the flux surfaces
    def tomas(self,magr,magz):

        if self.Status == True:

            if size(shape(magr))<1:
                print 'for this method 3 surfaces must be given'
        
            n_theta, n_rho = magr.shape
 
            rho = linspace(1./n_rho,1,n_rho)
    #        theta = linspace(0,2*pi,n_theta,endpoint=False)
    
    
            magr, magz = copy(magr.T), copy(magz.T)
    
            r0,z0 = magr[0].mean(), magz[0].mean()
        #calculate gradient of Phi with resoect to R and z
            magr[0] +=  (magr[1]-magr[0])/100
            magz[0] +=  (magz[1]-magz[0])/100


            drdrho,drtheta = gradient(magr)
            dzdrho,dztheta = gradient(magz)
            dpsidrho,dpsitheta = gradient(tile(rho**2, (n_theta,1)).T )

            grad_rho = dstack((drdrho,dzdrho,dpsidrho ))
            grad_theta = dstack((drtheta,dztheta,dpsitheta))
            normal = cross(grad_rho,grad_theta,axis=-1)



            dpsi_dr = -normal[:,:,0]/(normal[:,:,2]+1e-8) #Bz
            dpsi_dz = -normal[:,:,1]/(normal[:,:,2]+1e-8) #Br

    #WARNING not defined on the magnetics axis

            dtheta_star =((magr-r0)**2+(magz-z0)**2)/(dpsi_dz*(magz-z0)+dpsi_dr*(magr-r0))/magr
            if self.origin=='HFS':
                theta = arctan2(magz-z0,magr-r0)
            elif self.origin=='LFS':
                theta = arctan2(magz-z0,r0-magr)

            theta = unwrap(theta-theta[:,(0,)],axis=1)

           
        #definition of the thetat star by integral
            theta_star = cumtrapz(dtheta_star,theta,axis=1,initial=0)
            #correction = (n_theta-1.)/n_theta
          #  if all(magr[:,0]==magr[:,-1]) and all(magz[:,0]==magz[:,-1]):
          #      correction = 1
          #  theta_star/= theta_star[:,(-1,)]/(2*pi)/correction     #normalizeto 2pi

            #theta_star 
            
            theta_star_out = ((2.*pi*theta_star.T / theta_star[:,-1])).T
            theta_out = ((2.*pi*theta.T / theta[:,-1])).T
           
          #  embed()

            return theta, theta_star, r0, z0



    def define_wls(self, rhopIn=0.8,nPol=3601):    
        if self.Status == True:
            magr,magz = self.get_surface_rz(rhop=[rhopIn],nPolAngle=nPol)
            theta_star,theta,self.R0wls,self.z0wls = self.wls(magr,magz)
            self.func_wls = interpolate.interp1d(theta,theta_star)
            self.Status_wls = True

    def getwls_theta2thetaStar( self, theta=None ):

        if self.Status_wls:
            if (theta != None) :
                return self.func_wls(theta)
            else:
                print 'no input or None input'
                return None
        else:
            print 'no function defined'
            return None


    def wls( self,  Rq_in, zq_in):

        if self.Status == True:
       #Define gemeometric axis
            Rgeo = 0.5*(max(Rq_in) + min(Rq_in))
            zgeo = 0.5*(max(zq_in) + min(zq_in))
       
        # poloidal angle (geometrical)
            theta_in = arctan2(zq_in-zgeo, Rq_in-Rgeo)
            
# rearrange to start at inner midplane
    
            theta,indices = unique(theta_in, return_index=True)
            theta_out=theta

            zq = zq_in[indices]
            Rq = Rq_in[indices]
    #take first and last R to get R of inner midplane
            Rimp = 0.5*(Rq[0]+Rq[-1])

#add (interpolated) point at theta = -\pi
            if theta_in[0] > -pi:
                theta = append(-pi,theta)
                zq = append(zgeo,zq)
                Rq = append(Rimp,Rq)

#same for pi
            if theta_in[-1] < pi:
                theta = append(theta,pi)
                theta_out=theta[:-1]
                zq = append(zq,zgeo)
                Rq = append(Rq,Rimp)  

        #make new coordinate system
            zz0 = zq - zgeo
            RR0 = Rq - Rgeo
            r = sqrt(zz0*zz0 + RR0*RR0)

            EQ = EQU.EQU()
            EQ.Load(self.shot,self.exp, self.diag)         
            equR,equz,equPSI=EQ(self.time) 
            
            #equRR, equzz = meshgrid(equR, equz)
            #f = interpolate.interp2d(equRR, equzz, equPSI, kind='cubic')
            f=interpolate.RectBivariateSpline(equR, equz, equPSI)#,bbox=[Rq.min(), Rq.max(), zq.min(), zq.max()])
            [psiq,dpsidR,dpsidz] = [f.ev(Rq,zq),f.ev(Rq,zq,dx=1),f.ev(Rq,zq,dy=1)]
##// g^{-1/2} R
#// = (\nabla \theta) \cdot (\nabla \Psi \times \nabla \Phi) R
#// = (1/r * e_\theta) \cdot (\vec{B}_p) R
#// = (1/r^2) (z_0-z, R-R_0) \cdot (-dPsi/dz, dPsi/dR) 
            gm12Rq = (zz0*dpsidz + RR0*dpsidR)/ r / r
            y = 1.0/ gm12Rq / Rq

            dx = theta[1:]-theta[0:-1]
            theta_star = cumtrapz(0.5*dx*(y[0:-1]+y[1:]),initial=0)
 
            
           # theta_star = (2*pi*theta_star / theta_star[-1]) - pi
            theta_star = (2.*pi*theta_star / theta_star[-1]) - pi
            theta_out = theta_out+pi
            theta_out = (2.*pi*theta_out / theta_out[-1]) - pi
            return theta_star,theta_out,Rgeo,zgeo

    
"""
// field_line_angle.sci

// For equilibrium 'equ' and given surface with poloidal flux 'psi'
// calculate poloidal angle 'theta_star' in straight-field-line system
// as a function of geometrical poloidal angle 'theta'.

// based on PhD thesis by M Schittenhelm, ch. 3.3
// see also D'haeseleer et al., Flux coordinates and magnetic field structure

// 04-Sep-2005 wls return (Rgeo,zgeo) of selected surface
// 06-Sep-2005 wls use flux_contour subroutine
// 21-Jul-2009 wls formulation modified for numerical stability
// 22-Jul-2009 wls fixed theta, theta star range -%pi ... %pi


function [theta_star, theta, r, Rgeo, zgeo] = field_line_angle(equ, psi)

// (R,z) coordinates of flux surface
  if exists("flux_contour")==0 then getf("flux_contour.sci"); end
  [Rq,zq]=flux_contour(equ,psi);

// minor axis 
  Rgeo = 0.5*(max(Rq) + min(Rq));
  zgeo = 0.5*(max(zq) + min(zq));

// poloidal angle (geometrical)
  theta_ = atan(zq-zgeo, Rq-Rgeo);

// rearrange to start at inner midplane
  [theta,si] = unique(theta_);
  zq = zq(si); Rq=Rq(si);
  Rimp = 0.5*(Rq(1)+Rq($));
  
// add (interpolated) point at theta = -\pi
  if theta(1) > -%pi then
    theta = [-%pi theta];
    zq = [zgeo zq];
    Rq = [Rimp Rq];
  end

// ditto. at theta = \pi
  if theta($) < %pi then
    theta = [theta %pi];
    zq = [zq zgeo];
    Rq = [Rq Rimp];
  end

  zz0 = zq - zgeo;
  RR0 = Rq - Rgeo;
  r = sqrt(zz0.^2 + RR0.^2);

// flux (test) and flux R,z-derivatives on surface
//  global C;
  C = splin2d(equ.R, equ.z, equ.PsiRz);
  [psiq,dpsidR,dpsidz] = interp2d(Rq,zq,equ.R,equ.z,C);

// g^{-1/2} R
// = (\nabla \theta) \cdot (\nabla \Psi \times \nabla \Phi) R
// = (1/r * e_\theta) \cdot (\vec{B}_p) R
// = (1/r^2) (z_0-z, R-R_0) \cdot (-dPsi/dz, dPsi/dR) 

  gm12Rq = (zz0.*dpsidz + RR0.*dpsidR) ./ r ./ r;
  y = 1.0 ./ gm12Rq ./ Rq;
// trapezoidal integration
  dx = theta(2:$)-theta(1:$-1);
  theta_star = cumsum([0 0.5*dx.*(y(1:$-1)+y(2:$))]);
  theta_star = (2*%pi*theta_star / theta_star($)) - %pi;

endfunction

"""
