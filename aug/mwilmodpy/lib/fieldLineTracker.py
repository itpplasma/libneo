#!/usr/bin/env python

import dd
import numpy
import kk
from IPython import embed
import EQU
import matplotlib.pylab as plt


class fieldLineTracker(object):




    

    def trace_field_line(shot,tshot,r_start,z_start,phi_start,npoints,lenght,phi_plane=None, edition=0,exper='AUGD',diag='EQI',plotit=False,direction='both'):
                    
#   -> shot         : shotnumber
#   -> tshot         : time in discharge
#   -> r_start      : r-coordinate in m at middle point of magnetic field line
#   -> z_start      : z-coordinate in m at middle point of magnetic field line 
#   -> phi_start    : phi-coordinate in degrees at middle point of magnetic fieldline [radians]
#   -> npnts        : number of points along magnetic field line (don't worry aboutthe discretization error)
#   -> lenght       : lenght of the field line [m]
    #phi_plane      : compute R,Z crossection for some plane at ohi = phi_plane 
    #direction     : forward, backward, both - how to trace the fieldlines
#   keywords for FPP/EQU diagnostic
#

        if direction == 'both':
            npoints/= 2
            from scipy.integrate import odeint
            from scipy.linalg import norm
    
        Binterp = kkrzbrzt(shot,tshot, exp=exper,diag=diag,ed=edition)
    #out = kk.kkrzBrzt(shot, tshot, 1.6,0,exp=exper, diag=diag, ed=edition )
    
    #NOTE  kkrzBrzt is using stupid linear interpolation and it is slow for EQH. 
    
        def fun(y,t0, Binterp ):

            r,z,theta = y
        
            B_r,B_z,B_t = Binterp(r,z)
            
            B = norm([B_r,B_z,B_t])*sign(t0)

            dr = -B_r/B
            dz = -B_z/B
            dphi = -B_t/B/r
        
            return squeeze(array((dr,dz,dphi)))
    
        Rlines, Zlines, phi_lines = [],[],[]
        if direction in ['backward','both']:
            T = linspace(-lenght/2.,0, npoints, endpoint=False)
        
            Rline,Zline,phi_line =  odeint(fun, (r_start,z_start,phi_start),T,args =
(Binterp,) ).T

            Rlines.append(Rline[::-1])
            Zlines.append(Zline[::-1])
            phi_lines.append(phi_line[::-1])
    
    
        if direction in ['forward','both']:
            T = linspace(lenght/2.,0, npoints, endpoint=False)[::-1]
            
            Rline,Zline,phi_line =  odeint(fun, (r_start,z_start,phi_start),T,args =
(Binterp,)).T

            Rlines.append(Rline)
            Zlines.append(Zline)
            phi_lines.append(phi_line)
    
            Rline = hstack(Rlines)
            Zline = hstack(Zlines)
            phi_line = hstack(phi_lines)

    
    
            
    if phi_plane != None:
        #n_rounds = int(ceil((phi_line[-1]-phi_line[0])/(2*pi)))+1
        #phi_plane2 = arange(n_rounds)*2*pi+phi_plane
        
        #print phi_plane, n_rounds
        #print phi_plane2
        #arange(ceil((phi_line[0]-phi_plane)/(2*pi)), floor((phi_line[-1]-phi_plane)/(2*pi)))
        
        
        phi_plane_ = 2*pi*arange(floor((phi_line[0]-phi_plane)/(2*pi)),
ceil((phi_line[-1]-phi_plane)/(2*pi)))+phi_plane


        Rplane = interp(phi_plane_,phi_line,Rline )
        Zplane = interp(phi_plane_,phi_line,Zline )
        


    if plotit:
        fig = figure(figsize=(10,5))
        ax = fig.add_subplot(121)

        kgc = kk.kkGCd0(shot)
        comp = kgc.gc_x.iterkeys
        for key in comp():
            ax.plot(kgc.gc_x[key], kgc.gc_y[key], 'k-', lw=0.5)
        #plot(Rplane,Zplane,',')
        #ax.set_xlim(0.9, 3)
        #ax.set_ylim(-1.5, 1.5)
        
        #show()
        if phi_plane != None:
            ax.plot(Rplane, Zplane,'bo')
            ind = phi_plane_>phi_start
            for i,(r,z) in enumerate(zip(Rplane[ind], Zplane[ind])):
                text(r+0.01,z+0.01,str(i+1))
            for i,(r,z) in enumerate(zip(Rplane[~ind][::-1], Zplane[~ind][::-1])):
                text(r+0.01,z+0.01,str(-i-1)) 
            ax.plot(r_start,z_start,'rx')
            ax.text(r_start+0.01,z_start+0.01,'start')
            #ax.plot(Rline*cos(phi_line-phi_plane),Zline,'b',lw=.2)

            
        else:
            ax.plot(Rline,Zline)
            ax.plot(Rline[[0,-1]],Zline[[0,-1]],'ro')
            ax.plot(r_start,z_start,'bo')

        #ax.plot(Rline[-1],Zline[-1],'ro')
        ax.set_xlabel('R [m]')
        ax.set_ylabel('z [m]')

        ax.axis('equal')
        #ax.set_xlim(0.9, 3)
        #ax.set_ylim(-1.5, 1.5)

     


        ax = fig.add_subplot(122)
        phi = linspace(0,2*pi,200)
        
        ax.plot(2.24*cos(phi),2.24*sin(phi),'k')
        ax.plot(0.96*cos(phi),0.96*sin(phi),'k')
        if  phi_plane != None:
            ax.plot([0.96*cos(phi_plane),2.24*cos(phi_plane)],[0.96*sin(phi_plane),2.24*sin(phi_plane)],'g--')
            text(cos(phi_start)*r_start+0.05,sin(phi_start)*r_start+0.05,'start')

        else:
            ax.plot(cos(phi_line[[0,-1]])*Rline[[0,-1]],sin(phi_line[[0,-1]])*Rline[[0,-1]],'ro'
)

        ax.plot(cos(phi_line)*Rline,sin(phi_line)*Rline )
        #plot(cos(phi_line)*Rline,sin(phi_line)*Rline )
        ax.plot(cos(phi_start)*r_start,sin(phi_start)*r_start,'ob' )

        
        #ax.plot(Rline[[0,-1]],Zline[[0,-1]],'ro')
        #ax.plot(Rline[-1],Zline[-1],'ro')
        ax.set_xlabel('X [m]')
        ax.set_ylabel('Y [m]')
        ax.axis('equal')
        show()
        

        
    if phi_plane != None:

        return Rplane,Zplane,phi_plane


        
        
        
    return Rline,Zline,phi_line



class B_interp:
    def __init__(self,PFM,R,Z,Bt):
        self.Phi = PFM
        self.R = R
        self.Z = Z
        
        dr = (R[-1]-R[0])/(len(R)-1)
        dz = (Z[-1]-Z[0])/(len(Z)-1)

        Br = -np.diff(self.Phi,axis=1)/dz/(2*np.pi*self.R[:,None])
        self.Br = (Br[1:,:]+Br[:-1,:])/2
        Bz = 
        np.diff(self.Phi,axis=0)/dr/(2*np.pi*.5*(self.R[1:]+self.R[:-1])[:,None])
        self.Bz = (Bz[:,1:]+Bz[:,:-1])/2
        self.Bt = Bt
        
        self.scaling = np.array([dr,dz])
        self.offset  = np.array([R[0]+dr/2,Z[0]+dz/2])
            
            
    def __call__(self,r,z):
            
        idx = (array([r,z])-self.offset)/self.scaling

        Br = map_coordinates(self.Br,idx[:,None],mode='nearest',order=2,prefilter=True)
        Bz =map_coordinates(self.Bz,idx[:,None],mode='nearest',order=2,prefilter=True)
        Bt = self.Bt*1.65/r  #just aproximation
            
        return Br,Bz,Bt
        
    #return B_interp(PFM,R,Z,Bt)
