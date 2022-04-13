


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
        ## if magnetic axis is used only used for wls
        self.Status_mag0 = False

        if ((shot != None) & (time != None)):
            self.exp = Experiment
            self.diag = Diagnostic
            self.ed = Edition
            self.shot = shot
            self.time = time
            self.Status = True


    def define(self, rhoIn=None,nrhop=1000,nPol=361, origin='HFS'):

        if self.Status == True:
            if all(rhoIn)==None:
                rhopIn = linspace(0.0,1.0,nrhop)
            else:
                rhopIn = rhoIn
            
            if origin=='HFS':
                tomas=True
            elif origin=='LFS':
                tomas=False

            magr,magz = self.get_surface_rz(rhop=rhopIn,nPolAngle=nPol,tomas=tomas)
            theta,theta_star,self.R0,self.z0 = self.tomas(magr,magz,rhopIn = rhopIn)
                   
            #self.theta_star = (2.*pi*theta_star.T / theta_star[:,-1]).T-pi
            #self.theta = (2.*pi*theta.T / theta[:,-1]).T-pi
            
            #embed()
            if origin == 'LFS':
                ind1 = int(mean(argmin(abs(theta-pi), axis=1)))
                self.theta =concatenate((theta[:,ind1:]-2.*pi,theta[:,:ind1]),axis=1)
                self.theta_star = concatenate((theta_star[:,ind1:]-2.*pi,theta_star[:,:ind1]),axis=1)
                self.R = concatenate(( (magr.T)[:,ind1:],(magr.T)[:,:ind1]),axis=1)
                self.z = concatenate(( (magz.T)[:,ind1:],(magz.T)[:,:ind1]),axis=1)
            else:
                self.theta = theta - pi
                self.theta_star = theta_star - pi
                self.R,self.z = magr.T, magz.T
            
            #rhop=reshape(repeat(rhopIn[1:],nPol),(nrhop-1,nPol))
               
            self.func2D = []
            
            self.rhop = rhopIn
            self.origin = origin
            
            for i in arange(0,size(self.rhop)):
                self.func2D.append( interpolate.interp1d(self.theta[i],self.theta_star[i],bounds_error=True) ) 
            self.Status_Func = True

### Get Theta star for everything
    def getf_theta2thetaStar( self, rhopIn=None, thetaIn=None, wls = False ):

        if ( (self.Status_Func & (wls==False)) | (self.Status_wls & wls)) :
            if (all(rhopIn) == None) | (all(thetaIn) == None) :
                print 'no input or None input'
                return None

            if not(wls):
                rhopFunc = self.rhop
                func2D = self.func2D
            else:
                rhopFunc = self.rhopwls
                func2D = self.func2Dwls

            if size(rhopIn) == 1:
                idx = argmin( abs( rhopFunc - rhopIn ) )
                return func2D[idx](thetaIn)
            elif size(rhopIn) > 1 :
                nPoints = size(rhopIn)
                if nPoints == size(thetaIn):                
                    thetaStar=zeros((nPoints))
                    for i in arange(nPoints):
                        idx = argmin( abs( rhopFunc - rhopIn[i] ) )
                        thetaStar[i] = func2D[idx](thetaIn[i])
                    return thetaStar
                else:
                    print 'theta and rhop not the same size'
        else:
            print 'no function defined'
            return None

    def get_theta(self, R=None, z=None, wls=False):
   
####check input
        if ( (all(R) == None) | (all(z) == None)):
            print 'no input or None input'
            return None    
####check size
        if ( size(R)!=size(z) ) :
            print 'R and z input must have the same size'
            return None

        if (self.Status_Func & (wls==False) ):
            if self.origin=='HFS':
                return arctan2(z-self.z0,-R+self.R0) 
            elif self.origin=='LFS':
                return arctan2(z-self.z0,R-self.R0)     
               
            ###use theta
        elif (self.Status_wls & wls ):
            ## if only one z0 and R0 are given e.g.> magnetic axis was used
            if ( ( size(self.z0wls) == 1) & (size(self.R0wls) == 1) ) :
                return arctan2(z-self.z0wls,R-self.R0wls)    
###super langsam!!!!
            elif ( ( size(self.z0wls)>1) & (size(self.R0wls)>1 ) ):
                ##get rhop valuse
                output = kk.KK().kkrzptfn( self.shot, self.time, R, z,exp=self.exp, diag=self.diag )
                rhopIn = squeeze( output.rho_p )
                rhopFunc = self.rhopwls
                
                nPoints = size(rhopIn)              
                theta=zeros((nPoints))
                for i in arange(nPoints):
                    idx = argmin( abs( rhopFunc - rhopIn[i] ) )
                    theta[i] = arctan2(z[i]-self.z0wls[idx],R[i] - self.R0wls[idx])   
                return theta
            else:
                print 'something wrent wrong with R0 and z0 using wls'
                return None
        else:
            print 'no appropriate function defined'
            return None


 
    def get_rz_theta_thetaS( self, rz=None, q=None, rhop=None):

        if self.Status_Func:
            rhopIn=rhop
            
            if (all(rhopIn) == None) & (all(q) == None) & (all(rz) == None):
                print 'no input chosen'
                return None                

            if (all(q) != None):
        ## q is given and thop will be calculated
                print 'q it is...'
                if (size(q) < 1):
                    print 'rz must have one value at least'
                    return None
    
            #first read magnetic axis and separatrix
                output = kk.KK().kkrhorz( self.shot, [self.time], [0.0,1.0],exp=self.exp, diag=self.diag, ed=self.ed )
                
                RmagAxis_r = squeeze(output.r)[0]
                RmagAxis_z = squeeze(output.z)[0]
                Sep_r = squeeze(output.r)[1]
                Sep_z = squeeze(output.z)[1]
                 
            #generate
                runningR =  arange( RmagAxis_r, Sep_r,0.0005)
                runningz =  repeat( RmagAxis_z, size(runningR))
            #embed()
                output = kk.KK().kkrzq( self.shot, self.time, runningR, runningz,exp=self.exp, diag=self.diag, ed=self.ed )
                f = interpolate.interp1d(-output.q,output.rho_p, bounds_error=False, fill_value = float('NaN'))	
                print 'redefinition of rhop'
                rhop = f(q).tolist()
                print rhop

            elif (all(rz) != None):
                print 'rz it is...'
                if (size(rz) != 2):
                     print 'rz must have two values'
                     return None

                output = kk.KK().kkrzptfn(self.shot, self.time, [rz[0]], [rz[1]],exp=self.exp, diag=self.diag, ed=self.ed )
                rhop =  squeeze(output.rho_p)
    

            elif (all(rhopIn) != None):
                print 'rhop it is...'
                
                #rhopIn = array(rhopIn)
                if ( size(rhopIn) < 1) :
                    print rhopIn
                    print 'at least one rhop value'
                    return None
                elif (size(rhopIn) == 1):
                    rhop=array(rhopIn)
                else: 
                    idx=squeeze(where((rhopIn >= 0.0 ) & ( rhopIn <= 1.0 )))
                    if size(idx)<1:
                        print 'rhop must be between 0 and 1'
                        return None
                    else:
                        rhop=rhopIn[idx]

            if size(rhop) == 1:
                idx = argmin( abs( self.rhop - rhop ) )
                return self.R[idx],self.z[idx],self.theta[idx],self.theta_star[idx],self.rhop[idx]
            elif size(rhop) > 1 :
                nPoints = size(rhopIn)    
                R_out = zeros((nPoints))
                z_out = zeros((nPoints))
                rhop_out = zeros((nPoints))
                theta_out = zeros((nPoints))
                thetaStar_out=zeros((nPoints))
                for i in arange(nPoints):
                    idx = argmin( abs( self.rhop - rhop[i] ) )
                    R_out[i] = self.R[idx]
                    z_out[i] = self.z[idx]
                    rhop_out[i] = self.rhop[idx]
                    theta_out[i] = self.theta[idx]   
                    thetaStar_out[i] = self.theta_star[idx] 
                return R_out, z_out, theta_out, thetaStar_out, rhop_out



    def get_surface_rz( self, rz=None, q=None, rhop=None, nPolAngle=3601,thetaIn=None,tomas=False):
   
        if self.Status == True:
            try: 
                import kk_mwillens as kk
            except:
                print 'module import not possible'
                return None
    
            rhopIn=rhop

            if (all(rhopIn) == None) & (all(q) == None) & (all(rz) == None):
                print 'no input chosen'
                return None

   
            if (all(q) != None):
        ## q is given and thop will be calculated
                print 'q it is...'
                if (size(q) < 1):
                    print 'rz must have one value at least'
                    return None
    
            #first read magnetic axis and separatrix
                output = kk.KK().kkrhorz( self.shot, self.time, [0.0,1.0],exp=self.exp, diag=self.diag, ed=self.ed )
                RmagAxis_r = output.r[0]
                RmagAxis_z = output.z[0]
                Sep_r = output.r[1]
                Sep_z = output.z[1]
                 
            #generate
                runningR =  arange( RmagAxis_r, Sep_r,0.0005)
                runningz =  repeat( RmagAxis_z, size(runningR))
            #embed()
                output = kk.KK().kkrzq( self.shot, self.time, runningR, runningz,exp=self.exp, diag=self.diag, ed=self.ed )
                f = interpolate.interp1d(-output.q,output.rho_p, bounds_error=False, fill_value = float('NaN'))	
                print 'redefinition of rhop'
                rhop = f(q).tolist()
                print rhop
    
            elif (all(rz) != None):
                print 'rz it is...'
                if (size(rz) != 2):
                     print 'rz must have two values'
                     return None

                output = kk.KK().kkrzptfn(self.shot, self.time, [rz[0]], [rz[1]],exp=self.exp, diag=self.diag, ed=self.ed )
                rhop =  squeeze(output.rho_p)
    
            elif (all(rhopIn) != None):
                print 'rhop it is...'
                
                #rhopIn = array(rhopIn)
                if ( size(rhopIn) < 1) :
                    print rhopIn
                    print 'at least one rhop value'
                    return None
                elif (size(rhopIn) == 1):
                    rhop=array(rhopIn)
                else: 
                    idx=squeeze(where((rhopIn >= 0.0 ) & ( rhopIn <= 1.0 )))
                    if size(idx)<1:
                        print 'rhop must be between 0 and 1'
                        return None
                    else:
                        rhop=rhopIn[idx]

            if tomas:
                if (all(thetaIn) == None):
                    polAngle =  linspace(0.,360.,nPolAngle)-180.
                else:
                    polAngle = thetaIn
            else:
                if (all(thetaIn) == None):
                    polAngle =  linspace(0.,360.,nPolAngle)
                else:
                    polAngle = thetaIn

            

            #get rz of one surface
            #print 'get rz of surface of ', rhop
            if size(rhop) == 1: 
                rhop = [rhop]

            #embed()

            if size(shape(polAngle)) == 1:
                output = kk.KK().kkrhorz( self.shot, [self.time], rhop, angle = polAngle,exp=self.exp, diag=self.diag, ed=self.ed )
                surf_r =  squeeze(output.r)
                surf_z =  squeeze(output.z)

            elif size(shape(polAngle)) > 1 :
                idxRhop = where(array(shape(polAngle)) == rhop.size)[0]
                if idxRhop.size == 1:
                    polAngle = swapaxes(polAngle,idxRhop,0)
                    output = kk.KK().kkrhorz( self.shot, [self.time], rhop, angle = polAngle,exp=self.exp, diag=self.diag, ed=self.ed)
                    surf_r =  squeeze(output.r)
                    surf_z =  squeeze(output.z)
                elif idxRhop.size == 2:
                    output = kk.KK().kkrhorz( self.shot, [self.time], rhop, angle = polAngle,exp=self.exp, diag=self.diag, ed=self.ed )
                    surf_r =  squeeze(output.r)
                    surf_z =  squeeze(output.z)
                else:
                    print 'input theta did not have the the same index size as Input rhop'
            else:
                print 'something went wront'
     
            
   #  out = kk.KK().kkrhopto( shot, time, [rhop],,exp=exp_eq, diag=dia_eq)
            return surf_r,surf_z
           


#routine from tomas you should use KKRHORZ to find the flux surfaces, uses magnetic axis as center
    def tomas(self,magr,magz, rhopIn=None, origin = None):

        if self.Status == True:

            if size(shape(magr))<1:
                print 'for this method 3 surfaces must be given'
        
            if origin == None:
                origin = 'HFS'
            

            n_theta, n_rho = magr.shape
            
            if all(rhopIn) == None:
                rho = linspace(1./n_rho,1,n_rho)
            else:
                print 'use rhop from input'
                #embed()
                rho = rhopIn
      
            magr, magz = copy(magr.T), copy(magz.T)
    
            #to avoid singularity not exactly the center
            r0,z0 = magr[0].mean(), magz[0].mean()
        #calculate gradient of Phi with resoect to R and z
            magr[0] +=  (magr[1]-magr[0])/100
            magz[0] +=  (magz[1]-magz[0])/100

#gradient in rho and theta
            drdrho,drtheta = gradient(magr)
            dzdrho,dztheta = gradient(magz)
            dpsidrho,dpsitheta = gradient(tile(rho**2, (n_theta,1)).T )
      #      dpsidrho2,dpsitheta2 = gradient(tile(rho2**2, (n_theta,1)).T )

            grad_rho = dstack((drdrho,dzdrho,dpsidrho ))
       #     grad_rho2 = dstack((drdrho,dzdrho,dpsidrho2 ))

            grad_theta = dstack((drtheta,dztheta,dpsitheta))
       #     grad_theta2 = dstack((drtheta,dztheta,dpsitheta2))
            normal = cross(grad_rho,grad_theta,axis=-1)
       #     normal2 = cross(grad_rho2,grad_theta2,axis=-1)


            dpsi_dr = -normal[:,:,0]/(normal[:,:,2]+1e-8) #Bz
            dpsi_dz = -normal[:,:,1]/(normal[:,:,2]+1e-8) #Br

       #     dpsi_dr2 = -normal2[:,:,0]/(normal2[:,:,2]+1e-8) #Bz
        #    dpsi_dz2 = -normal2[:,:,1]/(normal2[:,:,2]+1e-8) #Br

    #WARNING not defined on the magnetics axis

            dtheta_star =((magr-r0)**2+(magz-z0)**2)/(dpsi_dz*(magz-z0)+dpsi_dr*(magr-r0))/magr
        #    dtheta_star2 =((magr-r0)**2+(magz-z0)**2)/(dpsi_dz2*(magz-z0)+dpsi_dr2*(magr-r0))/magr

            theta = arctan2(magz-z0,r0-magr)

            theta = unwrap(theta-theta[:,(0,)],axis=1)

           
        #definition of the thetat star by integral
            theta_star = cumtrapz(dtheta_star,theta,axis=1,initial=0)
            
            theta_star_out = (2.*pi*(theta_star.T - theta_star[:,0])/ (theta_star[:,-1]- theta_star[:,0])).T
            theta_out = (2.*pi*(theta.T - theta[:,0])/ (theta[:,-1]- theta[:,0])).T
     

            return theta_out, theta_star_out, r0, z0


## input rhopIn = arange(0.0,1.0,1.0/nrhop)
    def define_wls(self, rhopIn=0.8,nPol=3601, magnAxis = False):    

        if self.Status == True:
            rhop = array(rhopIn)
            ##embed()
            ## get contours
            magr,magz = self.get_surface_rz(rhop=rhop, nPolAngle=nPol)

            ##use magnetic axis to get R0 and z0
            if magnAxis:
                            #first read magnetic axis and separatrix
                print 'read magnetic axis'
                output = kk.KK().kkrhorz( self.shot, self.time, [0.0,1.0],exp=self.exp, diag=self.diag )
                Rmag = output.r[0]
                zmag = output.z[0]
            else:
                Rmag = None
                zmag = None

### use EQulirium
            EQ = EQU.EQU()
            EQ.Load(self.shot,self.exp, self.diag)         
            equR,equz,equPSI=EQ(self.time) 
            EQ.Unload()

            #equRR, equzz = meshgrid(equR, equz)
            #f = interpolate.interp2d(equRR, equzz, equPSI, kind='cubic')
            self.fPsi=interpolate.RectBivariateSpline(equR, equz, equPSI)#,bbox=[Rq.min(), Rq.max(), zq.min(), zq.max()])

            self.rhopwls = rhop
            
            if size(self.rhopwls) == 1:
               
                theta_star,theta,self.R0wls,self.z0wls = self.wls(magr,magz,R0=Rmag,z0=zmag)

                self.func2Dwls = interpolate.interp1d(theta,theta_star) 
                self.Status_wls = True
            elif size(self.rhopwls) > 1:
                self.func2Dwls = []
                self.R0wls = zeros_like(self.rhopwls)
                self.z0wls = zeros_like(self.rhopwls)

                for i in arange(0,size(self.rhopwls)):
                    theta_star,theta,self.R0wls[i],self.z0wls[i] = self.wls(magr[:,i],magz[:,i],R0=Rmag,z0=zmag)
                    self.func2Dwls.append( interpolate.interp1d(theta,theta_star,bounds_error=True) )     
 
                self.Status_wls = True

            if magnAxis:
                self.R0wls = Rmag
                self.z0wls = zmag
    ### 
    def wls( self,  Rq_in, zq_in, R0=None,z0=None):

        if self.Status == True:
       #Define gemeometric axis
            if ( all(R0) == None) | (all(z0) == None):
                Rgeo = 0.5*(max(Rq_in) + min(Rq_in))
                zgeo = 0.5*(max(zq_in) + min(zq_in))
            else:
                #Rgeo and zgeo are not the geom. axis anymore
                Rgeo = R0
                zgeo = z0
       
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

 
            [psiq,dpsidR,dpsidz] = [self.fPsi.ev(Rq,zq),self.fPsi.ev(Rq,zq,dx=1),self.fPsi.ev(Rq,zq,dy=1)]
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
