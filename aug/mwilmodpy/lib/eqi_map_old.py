#!/usr/bin/env python 
# -*- coding: utf-8 -*-

#Fast routines for equlibrium mapping (beta version)
__authors__ = 'Tomas Odstrcil'
__email__ = "tomas.odstrcil@ipp.mpg.de"
__version__ = '0.1'
__date__ = '25.02.2014'

""" ### Disclaimer: I'm not responsible for any harm done by this script  """

import numpy as np
import os,sys

#Load specific Libraries for kk and dd
import imp

import dd_alex as dd
from scipy.ndimage.interpolation import map_coordinates
import time
from scipy.interpolate import UnivariateSpline, interp1d
from scipy import signal
from IPython import embed

dd = dd.shotfile()


def MovingAveradge(sig, n,axis=-1):
    #Fast algorithm for calculation of the moving averadge
    sig.swapaxes(axis, -1)
    n = int(n)


    sig.cumsum(axis=-1,out=sig)
    right_pad  = np.ones(n/2)*sig[...,-1][...,None]
    left_pad = np.zeros(np.shape(sig[...,0])+((n+1)/2,))
    

    cs = np.concatenate((sig[...,n/2:],right_pad), axis=-1)
    cs-= np.concatenate((left_pad,sig[...,:-n/2]), axis=-1)
    cs *=1./n
    edge = 1-np.arange(n/2+1.)/n
    cs[...,:(n+1)/2] /= edge[-2+n%2::-1]
    cs[...,-(n+1)/2:]/= edge
    cs.swapaxes( -1,axis)
    

    return cs

def MovingAveradgeFast(sig, n,axis=-1):
    #Fast inplace algorithm for calculation of the moving averadge
    #return only data from n/2 to end-n/2!!!
    sig.swapaxes(axis, -1)

    sig.cumsum(axis=-1,out=sig)

    sig[...,:-n]-= sig[...,n:]
    sig = sig[...,:-n]

    sig *=-1./n

    sig.swapaxes( -1,axis)
    
    return sig

def interp_core(input):    return interp1d(*input[:-1])(input[-1])

class eqi_map:
    def __init__(self):
        self.shot = None
        self.diag = None
        self.ed = None
        self.exp = None
        self.dd_open=False
        self.ed_open=0
     
    def Open(self,diag,shot,experiment='AUGD',edition=0):
        if diag==self.diag and edition==self.ed and experiment==self.exp and shot==self.shot:
            if self.dd_open:
                return True
            else:
                self.diag,self.ed,self.exp,self.shot = diag,edition,experiment,shot
                self.dd_open = dd.Open(diag,int(shot),experiment=experiment,edition=edition)
                return self.dd_open
        else:
            if self.dd_open:    dd.close()
            self.diag,self.ed,self.exp,self.shot = diag,edition,experiment,shot
            self.dd_open=dd.Open(diag,int(shot),experiment=experiment,edition=edition)
            return self.dd_open
    def Close(self):
        self.dd_open=False
        dd.Close()
        



    def kkrhopto(self, rho,tin,shot,diag='EQH', exp='AUGD', ed=0,integ_time=0,rho_lbl='rho_pol',extrapolate=False):
        """Fast and smooth!  mapping from rho_pol to rho_tor
        
        Parameters
        ----------
        rho : ndarray
            rho_pol/roh_tor coordinates, 1D (time constant) or 2D (time variable) of size (nt,nx)
        tin : 1darray
            T coordinates, float, 1D (time constant) 
            
        shot: int
            number of the mapped shot
        diag: str
            diagnsotics used for mapping (EQI,EQH,...)
        exp: str
            experiment (AUGD)
        ed:  int
            diag edition
        integ_time:float
            smooth the position by boxcar method (averadge position over integration window of the diagnostic)
        rho_lbl:  str ''rho_pol, 'rho_tor' 
            input profile
        extrapolate:bool:
            extrapolate rho_tor out of 1
        
        Returns
        -------
        rho_tor : 2d array
        Magnetics flux rho_tor coordinates of the points
        
        """
    
        rho = np.atleast_2d(rho)
        self.shot = int(shot)
        tin = np.atleast_1d(tin)
        
        if tin.ndim==2:
            tin = tin.mean(1) #BUG 
        ntin = np.size(tin)
        
        if np.size(rho,0) == 1:
            rho = np.tile(rho, (ntin,1))
        
        stat = False
        if diag==self.diag and ed==self.ed and exp==self.exp and shot==self.shot:
            stat = True
        
        if not self.Open(diag,self.shot,experiment=exp,edition=ed):
            print 'opening failured, no diag: '+diag
            return rho
        

        if not hasattr(self,'tvec') or not stat:  
            self.tvec = dd.GetSignal('time')
        tvec = self.tvec
        #value of the toroidal flux at the different poloidal flux surfaces
        if not hasattr(self,'PFL') or not stat:  
            self.PFL = dd.GetSignalGroup('PFL')
        PFL = np.copy(self.PFL)
        
         #toroidal flux label vs PFL 
        if not hasattr(self,'TFLx') or not stat:  
            self.TFLx = dd.GetSignalGroup('TFLx')
        TFLx = np.copy(self.TFLx)
        #Poloidal Flux values (0:LPFx) 
        if not hasattr(self,'PFxx') or not stat:  
            self.PFxx = dd.GetSignalGroup('PFxx')
        PFxx = np.copy(self.PFxx)
        #0,: is mag axis (tvec)
        #1.: is separatrix (tvec)

        #dd.Close()
        

        
        nl = np.size(PFL,0)
        dt = (self.tvec[-1]-self.tvec[0])/len(self.tvec)
    
        ind = (self.tvec >= np.min(tin)-dt-integ_time)&(self.tvec <= np.max(tin)+dt+integ_time)
        


        mag = PFxx[0,:]
        
        orientation = np.sign(PFxx[0,:].mean())
        ikCAT = np.argmax(PFxx[1:,:]*orientation,axis=0)+1
        sep = PFxx[ikCAT,np.arange(len(mag))]

    
        
        rho_out = np.zeros_like(rho)
    
        if integ_time!= 0:
            NF = round(integ_time/dt)
            if NF!= 0:
                TFLx= MovingAveradge(TFLx,NF)
                PFL = MovingAveradge(PFL,NF)
                mag = MovingAveradge(mag,NF)
                sep = MovingAveradge(sep,NF)
    

        i0 = -1
        for jt,t in enumerate(tin):
            i = np.argmin(abs(self.tvec-t))  
            if i0 != i:  #already calculated in the previous step

                sort_wh=np.argsort(PFL[:,i])


                mg = np.interp(mag[i], PFL[sort_wh,i],TFLx[sort_wh,i])
                #sp = np.interp(sep[i], PFL[sort_wh,i],TFLx[sort_wh,i])
                sp = np.min(TFLx[:,i]) #best estimate of sep in rho_tor
                if abs(sp-mg)<1e-4:
                    continue
                tfl = (TFLx[sort_wh,i]-mg)/(sp-mg)
                pfl = (PFL[sort_wh,i]-mag[i])/(sep[i]-mag[i])

                
                pfl[pfl < 0] = 0#round errors
                tfl[tfl < 0] = 0#round errors
                ind = (pfl < 1-1e-3) & (tfl > 1e-3) 
                ind[-1] = True
                pfl,tfl =  np.sqrt(pfl[ind]), np.sqrt(tfl[ind])
                try:
                    if rho_lbl == 'rho_pol':
                        s = UnivariateSpline(pfl[::-1],tfl[::-1],k=3, s=1e-4,bbox=[0, 1],) 
                    elif rho_lbl == 'rho_tor':
                        s = UnivariateSpline(tfl[::-1],pfl[::-1],k=3, s=1e-4,bbox=[0, 1],) 
                except:
                    continue
            
            i0 = i
            #TODO zkontrolovat jestli extrapolace funguje, není tam záporná derivace? 
            
            edge = rho[jt,:] >= 1
            rho_out[jt,:] = s(rho[jt,:])
            if any(np.isnan(rho_out[jt,:])):  #UnivariateSpline has failured
                rho_out[jt,:] = np.interp(rho[jt,:], pfl[::-1],tfl[::-1])
            if not extrapolate:
                rho_out[jt,edge] = 1
            
        if not extrapolate:
            rho_out[rho_out>1] = 1  #round errors
            
        rho_out[rho_out<0] = 0  #round errors
        
        
        return rho_out
    

        import matplotlib.pylab as plt

        rhotor2 = np.zeros_like(rhotor)
        for jt,t in enumerate(tin):
            pto=kk.KK().kkrhopto(shot,t,rho[jt,:],ed=ed)
            rhotor2[jt,:] = pto.rho_t
            ed = pto.ed
        plt.title('kkrhopto')
        plt.plot(rhotor2,'r',linewidth = 0.2)
        plt.plot(rhotor,'b--',linewidth = 0.2)
        plt.show()
            


    def kkrhoPVo(self,shot,tin,rho,diag='EQH', exp='AUGD', ed=0,rho_lbl='rho_pol',strahl_def=True): 
        """Fast mapping from rho_tor to r_V, Normalized volume flux radius 
            https://www.aug.ipp.mpg.de/aug/local/aug_only/flcoord/flcoord_2.html
        
        Parameters
        ----------
        rho : ndarray
            rho coordinates (rho_pol or rho_tor), 
            1D array
        tin : 1darray
            T coordinates, float, 1D (time constant) 
            
        shot: int
            number of the mapped shot
        diag: str
            diagnsotics used for mapping (EQI,EQH,...)
        exp: str
            experiment (AUGD)
        ed:  int
            diag edition
            
        rho_lbl:str  choose 'rho_tor' or 'rho_pol'

        
        Returns
        -------
        r_V : 2d array
        Normalized volume flux radius 
        
        """
        #tin = atleast_1d(tin)
        
        if rho_lbl in ('r_V','r_V_strahl'):
            raise Exception('can not map r_V -> r_V')
        
        
        
        tin = np.squeeze(tin)
        #rho = np.asarray(rho)
        rho = np.atleast_2d(rho)
        #if rho.ndmin== 1:
            #rho = np.r_[0,rho,1]
        #elif rho.ndmin == 2:
            #rho = np.c_[tin*0,rho,tin*0+1]
        rho = np.squeeze(np.c_[[0,]*rho.shape[0],rho,[1,]*rho.shape[0]])
            
        _,V,_,_ = self.PlasmaVolume(shot,rho,tin,diag=diag,exp=exp,ed=ed,rho_lbl=rho_lbl)
        V0,Vs,V = V[:,0],V[:,-1],V[:,1:-1]


        if strahl_def:

            if not self.Open('GQ'+diag[-1],shot, experiment=exp, edition=ed):        
                return 

            R = dd.GetSignal('Rmag')
            tvec = dd.GetTimebase('Rmag')
            R = np.atleast_1d(np.interp(tin,tvec,R))
            r_V = np.sqrt(abs(V-V0[:,None])/(2*np.pi**2*R[:,None]))  #abs is there becaause of very small numerical errors
        else:
            r_V = np.sqrt(abs(V-V0[:,None])/(Vs-V0)[:,None],out=V)#abs is there becaause of very small numerical errors
        
     
        
        
        
        return r_V
                   
        
    def kkeqrinv(self,shot,tin,rho,diag='EQH', exp='AUGD', ed=0,rho_lbl='rho_pol'):
        iR = self.remapFromPFS(shot,tin,rho,'Rinv',diag=diag, exp=exp, ed=ed,rho_lbl=rho_lbl)
        iR_2 = self.remapFromPFS(shot,tin,rho,'R2inv',diag=diag, exp=exp, ed=ed,rho_lbl=rho_lbl)

        return iR,iR_2
    
    
        z_grid = np.linspace(-.90,.90,400)
        r_grid = np.linspace(.1,2.3,400)
        R,Z = np.meshgrid(r_grid, z_grid)
        M =  self.map_RZ2rho_matrix(shot,tin,rho,r_grid,z_grid,rho_lbl=rho_lbl)
        iR2 = np.vstack([m*(1/R.ravel()) for m in M])
        iR_22 = np.vstack([m*(1/R.ravel()**2) for m in M])
        
        #from matplotlib.pylab import *
        plot(iR.T,'--');plot(iR_2.T,'--');plot(iR_22.T);plot(iR2.T);savefig('iR.png');clf()
        exit()
    

    def remapFromPFS(self,shot,tin,rho,quantity,diag='EQH', exp='AUGD', ed=0,rho_lbl='rho_pol'):
        
        """Fast calculation of any PFL quantity from EQH
        
        Parameters
        ----------
        rho : ndarray
            rho coordinates (rho_pol or rho_tor), 
            1D (time constant) or 2D (time variable) of size (nt,nx)
        tin : 1darray
            T coordinates, float, 1D (time constant) 
            
        shot: int
            number of the mapped shot
        quantity: str
            name of the quantity, like  Qpembedsi, B2ave,...
        diag: str
            diagnsotics used for mapping (EQI,EQH,...)
        exp: str
            experiment (AUGD)
        ed:  int
            diag edition
            
        rho_lbl:str  choose 'rho_tor' or 'rho_pol'
    
        
        Returns
        -------
        Quantity : 2d array
         profile at times tin and positions rho_tor
        """
        
        
        
        stat = False
        if diag==self.diag and ed==self.ed and exp==self.exp and shot==self.shot:
            stat = True
        if not self.Open(diag,shot,experiment=exp,edition=ed):
            print 'no '+diag
            return 
  

        if not hasattr(self,'tvec') or not stat:  
            self.tvec = dd.GetSignal('time')
        tvec = np.copy(self.tvec)   
            
        #index vector for monotonously increasing time 
        if not hasattr(self,'ixti') or not stat:  
            self.ixti = np.int_(dd.GetSignal('ixti')-1)
        ixti = np.copy(self.ixti)
        
        #value of the toroidal flux at the different poloidal flux surfaces
        if not hasattr(self,'PFL') or not stat:  
            self.PFL = dd.GetSignalGroup('PFL')
        PFL = np.copy(self.PFL)
        
        if not hasattr(self,'TFLx') or not stat:  
            self.TFLx = dd.GetSignalGroup('TFLx')
        TFLx = np.copy(self.TFLx)
        #Poloidal Flux values (0:LPFx) 
        if not hasattr(self,'PFxx') or not stat:  
            self.PFxx = dd.GetSignalGroup('PFxx')
        PFxx = np.copy(self.PFxx)
        
        #toroidal flux label vs PFL 
        if not hasattr(self,'PFxx') or not stat:  
            self.PFxx = dd.GetSignalGroup('PFxx')
        PFxx = np.copy(self.PFxx)
        #0,: is mag axis (tvec)
        #1.: is separatrix (tvec)

        #q_value vs PFL 
        if not hasattr(self,quantity) or not stat:  
            setattr(self,quantity, dd.GetSignalGroup(quantity))
            #self.Qpsi = dd.GetSignalGroup(quantity)
        Qpsi = np.copy(getattr(self, quantity))

        #dd.Close()
        
        
        
        rho = np.atleast_2d(rho)
        tin = np.atleast_1d(tin)
        
        nl = np.size(PFL,0)
        dt = (tvec[-1]-tvec[0])/len(tvec)
    
        ind = (tvec >= np.min(tin)-dt) &  (tvec <= np.max(tin)+dt)

        ixti = ixti[ind]
        tvec = tvec[ixti]
        TFLx = TFLx[:,ixti]
        PFL = PFL[:,ixti]
        mag = PFxx[0,ixti]
        sep = PFxx[1,ixti]
        Qpsi = Qpsi[:,ixti]
        
        
        ntin = np.size(tin)
        
        if np.size(rho,0) == 1:
            rho = np.tile(rho, (ntin,1))
        
        Q = np.zeros_like(rho)
        i0 = -1
        for jt,t in enumerate(tin):
            
            i = np.argmin(abs(tvec-t))  
            if i0 == i:  #already calculated in the previous step
                Q[jt,:] = Q[jt-1,:]
                continue
                
            i0 = i
            
            sort_wh=np.argsort(PFL[:,i])


            mg = np.interp(mag[i], PFL[sort_wh,i],TFLx[sort_wh,i])
            sp = np.min(TFLx[:,i]) #best estimate of sep in rho_pol
                
            tfl = (TFLx[sort_wh,i]-mg)/(sp-mg)
            pfl = (PFL[sort_wh,i]-mag[i])/(sep[i]-mag[i])
            q = Qpsi[sort_wh,i]

            
            pfl[pfl < 0] = 0#round errors
            tfl[tfl < 0] = 0#round errors
            ind = (pfl < 1-1e-3) & (tfl > 1e-3)  #for pfl>1 is not tfl defined and sometimes tfl is zero even if it should not be
            ind[-1] = True

            pfl,tfl,q =  np.sqrt(pfl[ind]), np.sqrt(tfl[ind]),q[ind]
                
        
            
            #smooth fit of the noisy function  
            try:
                if rho_lbl == 'rho_tor':
                    s = UnivariateSpline(tfl[::-1],q[::-1],k=3, s=1e-4,bbox=[0, 1])
                elif rho_lbl == 'rho_pol':
                    s = UnivariateSpline(pfl[::-1],q[::-1],k=3, s=1e-4,bbox=[0, 1])  
                    
            except:
                print 'equilibrium mapping failured at %.4f'%t
                continue
            

            Q[jt,:] = s(rho[jt,:])
            if any(np.isnan(rho[jt,:])):  #UnivariateSpline has failured
                Q[jt,:] = np.interp(rho[jt,:], pfl[::-1],tfl[::-1])
        

        return Q
        
        
        


    def kkrhotpq(self, shot,tin,rho,diag='EQH', exp='AUGD', ed=0,rho_lbl='rho_pol'):
        """Fast calculation of safety factor profile 
        
        Parameters
        ----------
        rho : ndarray
            rho coordinates (rho_pol or rho_tor), 
            1D (time constant) or 2D (time variable) of size (nt,nx)
        tin : 1darray
            T coordinates, float, 1D (time constant) 
            
        shot: int
            number of the mapped shot
        diag: str
            diagnsotics used for mapping (EQI,EQH,...)
        exp: str
            experiment (AUGD)
        ed:  int
            diag edition
            
        rho_lbl:str  choose 'rho_tor' or 'rho_pol'
    
        
        Returns
        -------
        Q : 2d array
        safety profile at times tin and positions rho_tor
        
        """
        #stat = False
        #if diag==self.diag and ed==self.ed and exp==self.exp and shot==self.shot:
            #stat = True
        #if not self.Open(diag,shot,experiment=exp,edition=ed):
            #print 'no '+diag
            #return 
  

        #if not hasattr(self,'tvec') or not stat:  
            #self.tvec = dd.GetSignal('time')
        #tvec = np.copy(self.tvec)   
            
        ##index vector for monotonously increasing time 
        #if not hasattr(self,'ixti') or not stat:  
            #self.ixti = np.int_(dd.GetSignal('ixti')-1)
        #ixti = np.copy(self.ixti)
        
        ##value of the toroidal flux at the different poloidal flux surfaces
        #if not hasattr(self,'PFL') or not stat:  
            #self.PFL = dd.GetSignalGroup('PFL')
        #PFL = np.copy(self.PFL)
        
        #if not hasattr(self,'TFLx') or not stat:  
            #self.TFLx = dd.GetSignalGroup('TFLx')
        #TFLx = np.copy(self.TFLx)
        ##Poloidal Flux values (0:LPFx) 
        #if not hasattr(self,'PFxx') or not stat:  
            #self.PFxx = dd.GetSignalGroup('PFxx')
        #PFxx = np.copy(self.PFxx)
        
        ##toroidal flux label vs PFL 
        #if not hasattr(self,'PFxx') or not stat:  
            #self.PFxx = dd.GetSignalGroup('PFxx')
        #PFxx = np.copy(self.PFxx)
        ##0,: is mag axis (tvec)
        ##1.: is separatrix (tvec)

        ##q_value vs PFL 
        #if not hasattr(self,'Qpsi') or not stat:  
            #self.Qpsi = dd.GetSignalGroup('Qpsi')
        #Qpsi = np.copy(self.Qpsi)

        ##dd.Close()
        
        
        
        #rho = np.atleast_2d(rho)
        #tin = np.atleast_1d(tin)
        
        #nl = np.size(PFL,0)
        #dt = (tvec[-1]-tvec[0])/len(tvec)
    
        #ind = (tvec >= np.min(tin)-dt) &  (tvec <= np.max(tin)+dt)

        #ixti = ixti[ind]
        #tvec = tvec[ixti]
        #TFLx = TFLx[:,ixti]
        #PFL = PFL[:,ixti]
        #mag = PFxx[0,ixti]
        #sep = PFxx[1,ixti]
        #Qpsi = Qpsi[:,ixti]
        
        
        #ntin = np.size(tin)
        
        #if np.size(rho,0) == 1:
            #rho = np.tile(rho, (ntin,1))
        
        #Q = np.zeros_like(rho)
        #i0 = -1
        #for jt,t in enumerate(tin):
            
            #i = np.argmin(abs(tvec-t))  
            #if i0 == i:  #already calculated in the previous step
                #Q[jt,:] = Q[jt-1,:]
                #continue
                
            #i0 = i
            
            #sort_wh=np.argsort(PFL[:,i])


            #mg = np.interp(mag[i], PFL[sort_wh,i],TFLx[sort_wh,i])
            #sp = np.min(TFLx[:,i]) #best estimate of sep in rho_pol
                
            #tfl = (TFLx[sort_wh,i]-mg)/(sp-mg)
            #pfl = (PFL[sort_wh,i]-mag[i])/(sep[i]-mag[i])
            #q = Qpsi[sort_wh,i]

            
            #pfl[pfl < 0] = 0#round errors
            #tfl[tfl < 0] = 0#round errors
            #ind = (pfl < 1-1e-3) & (tfl > 1e-3)  #for pfl>1 is not tfl defined and sometimes tfl is zero even if it should not be
            #ind[-1] = True

            #pfl,tfl,q =  np.sqrt(pfl[ind]), np.sqrt(tfl[ind]),q[ind]
                
        
            
            ##smooth fit of the noisy function  
            #try:
                #if rho_lbl == 'rho_tor':
                    #s = UnivariateSpline(tfl[::-1],q[::-1],k=3, s=1e-3,bbox=[0, 1])
                #elif rho_lbl == 'rho_pol':
                    #s = UnivariateSpline(pfl[::-1],q[::-1],k=3, s=1e-3,bbox=[0, 1])  
                    
            #except:
                #print 'equilibrium mapping failured at %.4f'%t
                #continue
            

            #Q[jt,:] = s(rho[jt,:])
            #if any(np.isnan(rho[jt,:])):  #UnivariateSpline has failured
                #Q[jt,:] = np.interp(rho[jt,:], pfl[::-1],tfl[::-1])
        

        #return Q
        
        Q = self.remapFromPFS(shot,tin,rho,'Qpsi',diag=diag, exp=exp, ed=ed,rho_lbl=rho_lbl)
        return Q

        Q2 = np.zeros_like(rho)
        for jt,t in enumerate(tin):
            #print t, np.shape(tin)
            tpq=kk.KK().kkrhotpq(int(shot),t,rho[jt,:],ed=ed)
            Q2[jt,:]=-tpq.qp
            ed = tpq.ed
            #plt.plot(tpq.ft,'r--')
            #plt.plot(tpq.fp,'b--')
            #plt.plot(tpq.qp,'y--')

            #plt.plot(rho[jt,:],'k')
            #plt.show()


        #print rhopol
        plt.title('kkrhopto')
        plt.plot(-Q,'r',linewidth = 0.2)
        plt.plot(Q2,'b--',linewidth = 0.2)
        plt.ylim(0,10)
        plt.show()
            



    def kkrzbrzt(self,shot,tin,rin,zin,exp='AUGD',diag='EQH',ed=0,dR=0, dZ=0):
        """calculate a Br and Bz profiles
        
        Parameters
        ----------
        rin : ndarray
            R coordinates 
            1D (time constant) or 2D (time variable) of size (nt,nx)
        zin : ndarray
            Z coordinates 
            1D (time constant) or 2D (time variable) of size (nt,nx)   
            
        tin : 1darray
            T coordinates, float, 1D (time constant) 
            
        shot: int
            number of the mapped shot
        diag: str
            diagnsotics used for mapping (EQI,EQH,...)
        exp: str
            experiment (AUGD)
        ed:  int
            diag edition     
        dR:  float
            shift the  equilibrium to right by dR [m]
        dZ: float 
            shift the  equilibrium up by dZ [m]
        
        
        Returns
        -------
        interpBr : ndarray
            profile of Br on the grid
        interpBr : ndarray
            profile of Bz on the grid


        
        """
        
        stat = False
        if diag==self.diag and ed==self.ed and exp==self.exp and shot==self.shot:
            stat = True
        #import matplotlib.pylab as plt

        if not self.Open(diag,shot,experiment=exp,edition=ed):
            print 'no '+diag
            return 
        
     

        if not hasattr(self,'tvec') or not stat:  
            self.tvec = dd.GetSignal('time')
        tvec = np.copy(self.tvec)  
        
        #radial mesh (0:M) 
        if not hasattr(self,'Rmesh') or not stat:  
            self.Rmesh = dd.GetSignalGroup('Ri')
        Rmesh = np.copy(self.Rmesh)  
        
        #vertical mesh (0:N) 
        if not hasattr(self,'Zmesh') or not stat:  
            self.Zmesh = dd.GetSignalGroup('Zj')
        Zmesh = np.copy(self.Zmesh)  
        
        NR = int(dd.GetParameter('PARMV','M'))
        NZ = int(dd.GetParameter('PARMV','N'))


        #poloidal flux label (0:Lpf) 
        if not hasattr(self,'PFM') or not stat:  
            self.PFM = dd.GetSignalGroup('PFM')
        PFM = np.copy(self.PFM[:NR+1,:NZ+1,:])  
     
    
      
        #dd.Close()
        
        T = time.time()
        tin = np.atleast_1d(np.double(np.copy(tin)))
        rin = np.atleast_2d(rin)
        zin = np.atleast_2d(zin)
        nt = np.size(tin)
        
    
        if np.size(rin,0)!= nt:
            rin = np.tile(rin,(nt,1))
        if np.size(zin,0)!= nt:
            zin = np.tile(zin,(nt,1))   
    
        ntin = np.size(tin,-1)
        nrin = np.size(rin,1)
        
        interpBr = np.empty((nt, nrin))
        interpBz = np.empty((nt, nrin))

        i0 = -1
        
        
        for jt in xrange(ntin):
            t = tin[jt]
            if np.isnan(t): #calculated in the previous step
                continue
            
            
            ind = (tin == t)
    
            
            i = np.argmin(abs(tvec-t))  
            
            if i0 != i:  #already calculated in the previous step
    
                i0 = i
            
                R = Rmesh[:,i]
                Z = Zmesh[:,i]
                dr = (R[-1]-R[0])/(len(R)-1)
                dz = (Z[-1]-Z[0])/(len(Z)-1)
                Phi = PFM[...,i]

                Br = -np.diff(Phi,axis=1)/dz/(2*np.pi*R[:,None])
                Br = (Br[1:,:]+Br[:-1,:])/2
                Bz =  np.diff(Phi,axis=0)/dr/(2*np.pi*.5*(R[1:]+R[:-1])[:,None])
                Bz = (Bz[:,1:]+Bz[:,:-1])/2

            
                
                scaling = np.array([dr,dz])
                offset  = np.array([R[0]+dr/2,Z[0]+dz/2])
            
            
            
            #coords = np.dstack((rin[ind,:],zin[ind,:]))
            #coords = np.rollaxis(coords, -1)
            coords = np.array((rin[ind,:]-dR,zin[ind,:]-dZ))  #BUG shifting the whole grid equally
            
        
        
            idx = (coords-offset[:,None,None]) / scaling[:,None,None]

            interpBr[ind,:] = map_coordinates(Br,idx,mode='nearest',order=1,prefilter=True)
            interpBz[ind,:] = map_coordinates(Bz,idx,mode='nearest',order=1,prefilter=True)
            tin[ind] = np.nan
    
        return interpBr, interpBz

            



    def kkrzptfn(self,shot,tin,rin,zin,exp='AUGD',diag='EQH',  ed=0,rho_lbl='rho_pol',integ_time=0,dR=0,dZ=0):
        """Equilibrium mapping routine, fast for large number of points
        
        Parameters
        ----------
        rin : ndarray
            R coordinates, 1D (time constant) or 2D (time variable) of size (nt,nx)
        zin : ndarray
            Z coordinates, 1D (time constant) or 2D (time variable) of size (nt,nx)
        tin : ndarray
            T coordinates, float, 1D (time constant) or 2D (time variable) of size (nt,nx)  
            
        shot: int
            number of the mapped shot
        diag: str
            diagnsotics used for mapping (EQI,EQH,...)
        exp: str
            experiment (AUGD)
        ed:  int
            diag edition
        rho_lbl: str
            mapped coordinates - rho_pol or rho_tor, r_V, r_V_strahl
        integ_time:float
            smooth the position by boxcar method (averadge position over integration window of the diagnostic)
        dR:  float
            shift the  center of equilibrium to right by dR [m]
        dZ: float 
            shift the center of equilibrium up by dZ [m]
        
        Returns
        -------
        rho : 2d array
        Magnetics flux coordinates of the points
        
        """
      
        #BUG possible issue with non-equal time spacing, (variating of the sampling frequency), was not taken in account!
        

        t = time.time()
        #print 'shot', shot
        
        stat = False
        if diag==self.diag and ed==self.ed and exp==self.exp and shot==self.shot:
            stat = True
        
        if not self.Open(diag,shot, experiment=exp,edition=ed):
            print 'no '+diag
            return 
        self.ed_open=dd.edition

        if not hasattr(self,'tvec') or not stat:  
            self.tvec = dd.GetSignal('time')
        tvec_ = np.copy(self.tvec)  
        
        #index vector for monotonously increasing time 
        if not hasattr(self,'ixti') or not stat:  
            self.ixti = np.int_(dd.GetSignal('ixti')-1)
        ixti = np.copy(self.ixti)  
        
        #radial mesh (0:M) 
        if not hasattr(self,'Rmesh') or not stat:  
            self.Rmesh = dd.GetSignalGroup('Ri')
        Rmesh = np.copy(self.Rmesh)  
        
        #vertical mesh (0:N) 
        if not hasattr(self,'Zmesh') or not stat:  
            self.Zmesh = dd.GetSignalGroup('Zj')
        Zmesh = np.copy(self.Zmesh)  
       
        #poloidal flux label (0:Lpf) 
        if not hasattr(self,'PFM') or not stat:  
            self.PFM = dd.GetSignalGroup('PFM')
            self.PFM = self.PFM.astype(np.float32, copy=False)
        PFM_ = np.copy(self.PFM)  
        
        #Poloidal Flux values (0:LPFx) 
        if not hasattr(self,'PFxx') or not stat:  
            self.PFxx = dd.GetSignalGroup('PFxx')
        PFxx_ = np.copy(self.PFxx)          
        #0,: is mag axis (tvec)
        #1.: is separatrix (tvec)
        
        if not hasattr(self,'ikCAT') or not stat:  
            self.ikCAT = dd.GetSignalGroup('ikCAT')
        ikCAT_ = np.copy(self.ikCAT)   
        
        #dd.Close()
        
        print 'load ',time.time()-t

        t = time.time()
        
        dt = np.mean(np.diff(tvec_))
        #print tvec_.shape, ixti.shape
        

        ind=(tvec_>=np.min(tin)-dt*2.5-integ_time)&(tvec_<=np.max(tin)+dt*2.5+integ_time)
        

        Rmesh = Rmesh[:,1]
        nr = len(Rmesh)
        Zmesh = Zmesh[:,1]    
        nz =len(Zmesh)
        nl = np.size(PFxx_,0)


        if not any(ind):
            print 'Out of time range'
            return 

        #
        #BUG make a tvec equally spaced 
        dt = np.median(np.diff(tvec_[ixti]))

        index = np.cumsum(np.int_(np.round(np.diff(tvec_[ixti])/dt)))
        index = np.r_[0,index]

        
        index = index[ind]
        index-= index[0]
        nt = index[-1]+1
        ixti = ixti[ind]

        tvec = np.zeros(nt)
        tvec[index] = tvec_[ixti]

        PFM = np.zeros((nr,nz,nt),dtype=np.single)
        PFM[:,:,index] = PFM_[:nr,:nz,ixti]
        PFxx = np.zeros((nl,nt))
        PFxx[:,index] = PFxx_[:,ixti]
 #       ikCAT = np.zeros((nt))
 #       ikCAT[:,index] = ikCAT_[:,ixti]
        ind = np.where(tvec == 0)[0]


        orientation = np.sign(PFxx[0,:].mean())
        
        for m in ind:
            for j1 in xrange(m-1,-1,-1):
                if not j1 in ind: break
            for j2 in xrange(m,nt):
                if not j2 in ind: break
            #linear interpolation
            #PFM[...,m] =  ( PFM[...,j1]*(j2-m)+ PFM[...,j2]*(m-j1))/float(j2-j1)
            #PFxx[...,m] = (PFxx[...,j1]*(j2-m)+PFxx[...,j2]*(m-j1))/float(j2-j1)
            #nearest neighbour 
            PFM[...,m]  = PFM[...,j1]  if (m-j1)<(j2-m) else PFM[...,j2]
            PFxx[...,m] = PFxx[...,j1] if (m-j1)<(j2-m) else PFxx[...,j2]
            
            tvec[m] = (tvec[j1]*(j2-m)+tvec[j2]*(m-j1))/float(j2-j1)
    
        
        mag = PFxx[0,:]

        #find the edge value of the poloidal flux 
        #print 'ikCAT', ikCAT
        ikCAT = np.argmax(PFxx[1:,:]*orientation,axis=0)+1
        #print 'ikCAT',ikCAT
        sep = PFxx[ikCAT,np.arange(len(mag))]

        PFM-= mag
        PFM*=1/(sep-mag)
        rp = PFM
            
        if integ_time!= 0:
            NF = max(1.,round(integ_time/dt))

            rp   = MovingAveradgeFast(rp,NF)
            tvec = MovingAveradgeFast(tvec,NF)
    
        
        
        dr = (Rmesh[-1]-Rmesh[0])/nr
        dz = (Zmesh[-1]-Zmesh[0])/nz

        rin = np.atleast_2d(rin)
        zin = np.atleast_2d(zin)
        tin = np.atleast_2d(tin)
        ntin = np.size(tin,-1)
        nrin = np.size(rin,1)
        
        if np.size(rin,0) == 1:
            rin = np.tile(rin, (ntin,1))
        if np.size(zin,0) == 1:
            zin = np.tile(zin, (ntin,1))
        if np.size(tin,0) == 1:
            tin_ = np.tile(tin, (nrin,1)).T
        
        
        if rin.shape!= zin.shape:
            raise Exception('Wrong shape of rin or zin')
        if np.size(rin,0) != ntin:
            raise Exception('Wrong shape of rin %s'%str(rin.shape))
        if np.size(zin,0) != ntin:
            raise Exception('Wrong shape of zin %s'%str(zin.shape))
        if np.size(rin,0) != np.size(zin,0):
            raise Exception('Not equal shape of zin and rin %s,%s'%(str(zin.shape), str(zin.shape)))
        
        #coords = np.dstack((rin,zin,tin))
        #coords = np.rollaxis(coords, 2)
        coords = np.array((rin,zin,tin_))

        scaling = np.array([dr,dz,dt])
        offset = np.array([Rmesh[0]+dr/2,Zmesh[0]+dz/2,tvec[0]+0*dt/2])
        
    
        idx = (coords-offset[:,None,None]) / scaling[:,None,None]
        

        rho = np.empty((ntin,nrin),dtype=np.single)
        map_coordinates(rp, idx, mode='nearest', order=2,prefilter=True,output=rho)
        
        
        
        
        rho[rho<0] = 0
        rho = np.sqrt(rho,out=rho)*0.993  # must by multiplied 
        
        if dR!= 0 or dZ != 0:
            if dR!= 0:  rin-= dR*(1-rho)
            if dZ != 0: zin-= dZ*(1-rho)
            #coords = np.dstack((rin,zin,tin))
            #coords = np.rollaxis(coords, 2)
            coords = np.array((rin,zin,tin_))
            idx = (coords-offset[:,None,None]) / scaling[:,None,None]
            map_coordinates(rp, idx, mode='nearest', order=2,prefilter=True,output=rho)
            rho[rho<0] = 0
            rho = np.sqrt(rho,out=rho)*0.993  # must by multiplied 
        
        del rp
        

       
        
        print 'map ',time.time()-t
        t = time.time()
        if rho_lbl == 'rho_tor': 
            rho = self.kkrhopto(rho,tin_,shot,diag, exp, ed,integ_time=integ_time)
        if rho_lbl == 'r_V':
            rho = self.kkrhoPVo(shot,tin,rho,diag, exp, ed,rho_lbl='rho_pol',strahl_def=False)
        if rho_lbl == 'r_V_strahl':
            rho = self.kkrhoPVo(shot,tin,rho,diag, exp, ed,rho_lbl='rho_pol',strahl_def=True)


        if rho_lbl!= 'rho_pol':
            print 'remapping from rho_pol to '+rho_lbl, time.time()-t
        
        
        return rho




        
        t = time.time()

        #import matplotlib.pylab as plt 
        #old version
        #rhop2 = np.zeros((ntin,nrin))
        #rhot2 = np.zeros((ntin,nrin))
        #for j in range(ntin) :
            #ptfn=kk.KK().kkrzptfn(int(shot),tin[j,0],rin[j,:],zin[j,:],exp=exp,diag=diag,ed=ed)
            #rhop2[j,:] = ptfn.rho_p                       
            #rhot2[j,:] = ptfn.rho_t
        
        ptfn=kk.KK().kkrzptfn(int(shot),tin[:,0],rin[0,:],zin[0,:],exp=exp,diag=diag,ed=ed)  
        rhop2 = ptfn.rho_p                       
        rhot2 = ptfn.rho_t  
            
        print time.time()-t
        #exit()
        import matplotlib.pylab as plt
        
        fig = plt.figure('mag')
        ax = fig.add_subplot(111)
        ax.plot(tin, rho,'r',linewidth = 0.2)
        
        
        if rho_lbl=='rho_pol':
            ax.plot(tin, rhop2,'b--',linewidth = 0.2)
        else:
            ax.plot(tin,  rhot2,'b--',linewidth = 0.2)

        plt.xlim(tin.min(), tin.max())
        plt.show()

        return rhop2








    def map_rho2RZ_matrix(self,shot,tin,rho,r_grid,z_grid,rho_lbl='rho_tor',diag='EQH',exp='AUGD',ed=0,dR=0,dZ=0):
        """Calculate a matrix of the map from some flux averadged quantity to the 
            2d RZ coordinate
        
        Parameters
        ----------
        tin : ndarray
            time coordinate, constant or 1D vector
        r_grid : ndarray
            equally sampled 1D grid vector of R coordinate of the matrix
        z_grid : ndarray
            equally sampled 1D grid vector of Z coordinate of the matrix
        rho : ndarray
            1D  vector of rho coordinates for mapping
        shot: int
            number of the mapped shot
        diag: str
            diagnsotics used for mapping (EQI,EQH,...)
        exp: str
            experiment (AUGD)
        ed:  int
            diag edition
        rho_lbl: str
            mapped coordinates - rho_pol or rho_tor
        dR:  float
            shift the  center of equilibrium to right by dR [m]
        dZ: float 
            shift the center of equilibrium up by dZ [m]
        
    
        Returns
        -------
        maps : list of 2D arrays
        list of sparse matrixes for every time point,  matrix has a shape (nrho, nr*nt)
        
        """

        print 'map_rho2RZ_matrix'
        
        import scipy.sparse as sp

        
        nrho =  np.size(rho)
        nt   = np.size(tin)

        test_points = np.eye(nrho)
        


        
        R,Z,data_rho_grid = self.map_rho2RZ(shot,tin,rho,test_points,r_grid,z_grid,rho_lbl,diag,exp,ed,0,dR,dZ)

    
        M = data_rho_grid.reshape(nt,-1,nrho)
        
        try:
            import multiprocessing
            
            numcores = multiprocessing.cpu_count()
            pool = multiprocessing.Pool(numcores)
            M = pool.map(sp.csc_matrix, [M[i,...] for i in range(nt)])
                    
            pool.close()
            pool.join()
            
        except:
            M = [sp.csc_matrix(M[i,...]) for i in range(nt)]


        #pcolor(r_grid, z_grid, (M[0]*sin(100*theta)).reshape(200,201))
        return M


    def map_theta2RZ_matrix(self,shot,tin,theta,r_grid,z_grid,diag='GQH',exp='AUGD',ed=0,dR=0,dZ=0):
        """Calculate a matrix of the map from some angular function to the 
            2d RZ coordinate
        
        Parameters
        ----------
        tin : ndarray
            time coordinate, constant or 1D vector
        r_grid : ndarray
            equally sampled 1D grid vector of R coordinate of the matrix
        z_grid : ndarray
            equally sampled 1D grid vector of Z coordinate of the matrix
        theta : ndarray
            1D  vector of rho coordinates for mapping
        shot: int
            number of the mapped shot
        diag: str
            diagnsotics used for mapping (EQI,EQH,...)
        exp: str
            experiment (AUGD)
        ed:  int
            diag edition

        dR:  float
            shift the  center of equilibrium to right by dR [m]
        dZ: float 
            shift the center of equilibrium up by dZ [m]
        Returns
        -------
        maps : list of 2D arrays
        list of sparse matrixes for every time point,  matrix has a shape (ntheta, nr*nt)
        
        """

        print 'map_rho2RZ_matrix'
        
        import scipy.sparse as sp

        
        ntheta =  np.size(theta)
        nt   = np.size(tin)

        test_points = np.eye(ntheta)
        

       
        R,Z,data_theta_grid = self.map_theta2RZ(shot,tin,theta,test_points,r_grid,z_grid,diag,exp,ed,fill_value=0,dR=dR,dZ=dZ)

    
        M = data_theta_grid.reshape(nt,-1,ntheta)
        
        try:
            import multiprocessing
            
            numcores = multiprocessing.cpu_count()
            pool = multiprocessing.Pool(numcores)
            M = pool.map(sp.csc_matrix, [M[i,...] for i in range(nt)])
                    
            pool.close()
            pool.join()
            
        except:
            M = [sp.csc_matrix(M[i,...]) for i in range(nt)]
            
            
        
        

        return M










    def map_rho2RZ(self, shot,tin,rho,data_rho,r_grid,z_grid,rho_lbl='rho_tor',diag='EQH',exp='AUGD',ed=0,fill_value=np.nan,dR=0,dZ=0,integ_time=0):
        """Mapping routine from 1D averadged quantity to the 2D R,Z coordinates
        
        Parameters
        ----------
        tin : ndarray
            time coordinate, constant or 1D vector
        r_grid : ndarray
            equally sampled 1D grid vector of R coordinate of the matrix
        z_grid : ndarray
            equally sampled 1D grid vector of Z coordinate of the matrix
        data_rho : ndarray
            2D  vector of rho coordinates for mapping (nt,nr)
        rho : 1darray
            1D  vector of of the mapped function of rho
        shot: int
            number of the mapped shot
        diag: str
            diagnsotics used for mapping (EQI,EQH,...)
        exp: str
            experiment (AUGD)
        ed:  int
            diag edition
        rho_lbl: str
            mapped coordinates - rho_pol or rho_tor
        fill_value:float
            how to fill R,Z poinst not covered by rho
        dR:  float
            shift the  center of equilibrium to right by dr [m]
        dZ: float 
            shift the center of equilibrium up dz [m]
        integ_time: float
            integration time of the used diagnostic - data will be smoothed over this time
        Returns
        -------
        R: 2d array of the grid coordinates
        Z: 2d array of the grid coordinates
        data_rho_grid: mapped quantity
        
        """
        
        print 'map_rho2RZ'


        tin = np.atleast_1d(tin)
        nt = np.size(tin)
        data_rho = np.atleast_3d(data_rho.T).T
        R,Z = np.meshgrid(r_grid,z_grid)
        
        
        
        
        rho_grid = self.kkrzptfn(shot,tin,R.ravel(),Z.ravel(),exp=exp,diag=diag,ed=ed,rho_lbl=rho_lbl,integ_time=0,dR=dR,dZ=dZ)
        #rho_grid = self.kkrzptfn(shot,tin,R.ravel(),Z.ravel(),diag,rho_lbl,0,dR,dZ)
        rho_grid = rho_grid.reshape((nt,len(z_grid),len(r_grid)))
        
        #short version
        #data_rho_grid=interp1d(rho,data_rho,'linear',1,False,False,fill_value)(rho_grid)

        #long by multiprocessing  
        #print 
        #print data_rho.shape
        if nt==data_rho.shape[1]:
            arg = [(np.single(rho),np.single(data_rho[:,i,:]),'linear',1,False,False,fill_value, rho_grid[i,...]) for i in range(nt)]
        else:
            arg = [(np.single(rho),np.single(data_rho),'linear',1,False,False,fill_value, rho_grid[i,...]) for i in range(nt)]
    
        import multiprocessing
        
        numcores = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(numcores)
        out = pool.map(interp_core, arg)
                
        pool.close()
        pool.join()

        data_rho_grid = np.concatenate(out, axis=0)


        return R,Z,data_rho_grid



    def map_theta2RZ(self,shot,tin,theta,data_theta,r_grid,z_grid,diag='GQH',exp='AUGD',ed=0,fill_value=np.nan,dR=0,dZ=0):
        """Mapping routine from 1D angular quantity to the 2D R,Z coordinates
        
        Parameters
        ----------
        tin : ndarray
            time coordinate, constant or 1D vector
        r_grid : ndarray
            1D grid vector of R coordinate of the matrix
        z_grid : ndarray
            1D grid vector of Z coordinate of the matrix
        data_theta : ndarray
            nD  vector of theta coordinates for mapping  (periodic function)
        theta : 1darray
            nD  vector of of the mapped function of theta  rfom 0 to 2*pi
        shot: int
            number of the mapped shot
        diag: str
            diagnsotics used for mapping (EQI,EQH,...)
        exp: str
            experiment (AUGD)
        ed:  int
            diag edition
        dR:  float
            shift the  center of equilibrium to right by dr [m]
        dZ: float 
            shift the center of equilibrium up dz [m]

        fill_value:float
            how to fill R,Z poinst not covered by rho
        
        Returns
        -------
        R: 2d array of the grid coordinates
        Z: 2d array of the grid coordinates
        data_rho_grid: mapped quantity
        
        """
        print 'map_theta2RZ'

        if not self.Open(diag,shot,experiment=exp,edition=ed):
            return 
        
        R = dd.GetSignal('Rmag')
        Z = dd.GetSignal('Zmag')
        tvec = dd.GetTimebase('Rmag')
        #dd.Close()
        
        tin = np.atleast_1d(tin)
        nt = np.size(tin)
        
        R = np.interp(tin,tvec,R)+dR
        Z = np.interp(tin,tvec,Z)+dZ

        theta_grid = -np.arctan2((z_grid[None,:,None]-Z[:,None,None]), -(r_grid[None,None,:]-R[:,None,None]))+np.pi


        #data_theta = np.atleast_2d(data_theta)
        data_theta = np.atleast_3d(data_theta.T).T

        #short version
        #data_theta_grid=interp1d(theta,data_theta,'linear',1,False,False,fill_value)(theta_grid)

    
        print theta_grid.shape
        #long by multiprocessing  
        arg = [(theta,data_theta,'linear',1,False,False,fill_value,theta_grid[i,...]) for i in range(nt)]

        import multiprocessing

        print 'interp'
        numcores = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(numcores)
        out = map(interp_core, arg)
                
        pool.close()
        pool.join()

        data_theta_grid = np.concatenate(out, axis=0)


        return R,Z,data_theta_grid







    def map_RZ2rho_matrix(self,shot,tin,rho_out,r_grid,z_grid,rho_lbl='rho_tor',diag='EQH',exp='AUGD',ed=0,dR=0,dZ=0):
        """Calculate a matrix of the map from the flux averadged quantity to the 
            2d RZ coordinate
        
        Parameters
        ----------
        tin : ndarray
            time coordinate, constant or 1D vector
        r_grid : ndarray
            equally sampled 1D grid vector of R coordinate of the matrix
        z_grid : ndarray
            equally sampled 1D grid vector of Z coordinate of the matrix
        rho_out : ndarray
            1D equally sampled vector of rho coordinates for mapping
        shot: int
            number of the mapped shot
        diag: str
            diagnsotics used for mapping (EQI,EQH,...)
        exp: str
            experiment (AUGD)
        ed:  int
            diag edition
        rho_lbl: str
            mapped coordinates - rho_pol or rho_tor
        dR:  float
            shift the  center of equilibrium to right by dr [m]
        dZ: float 
            shift the center of equilibrium up dz [m]
    
        Returns
        -------
        maps : list of 2D arrays
        list of sparse matrixes for every time point,  matrix has a shape ( nr*nt,nrho)
        
        """
        #BUG sometimes wrong in the zero!!
        print 'map_RZ2rho_matrix'
        
        
        import matplotlib.pylab as plt

        tin = np.atleast_1d(tin)
    
        nt = np.size(tin)
        nrho = np.size(rho_out)
        nr = len(r_grid)
        nz = len(z_grid)

        R,Z = np.meshgrid(r_grid,z_grid)
        R,Z = R.ravel(),Z.ravel()
        rho =  self.kkrzptfn(shot,tin,R,Z,exp,diag,ed,rho_lbl,dR=dR,dZ=dZ)
        ind0 = np.argmin(rho,axis=1)
        R0,Z0 = R.flat[ind0][:,None],Z.flat[ind0][:,None]
        
        self.Close()
        if not dd.Open('GQ'+diag[-1],shot, experiment=exp, edition=ed):  
            return 

        R0 = dd.GetSignal('Rmag')
        Z0 = dd.GetSignal('Zmag')
        RZtvec = dd.GetTimebase('Rmag')
        R0 = np.atleast_1d(np.interp(tin,RZtvec,R0))[:,None]
        Z0 = np.atleast_1d(np.interp(tin,RZtvec,Z0))[:,None]
        dd.Close()
            
        
        import scipy.sparse as sp

        
        Br,Bz = self.kkrzbrzt(shot,tin,R,Z,exp,diag,ed,dR=dR,dZ=dZ)
        Bp = np.hypot(Br,Bz)

        
        weight = np.empty((nt, nr*nz,2)) #lower and upper weight of the rho
        drho = (rho_out.max())/(nrho-1)
        idx_rho = rho/drho
        weight[...,0] = (np.floor(idx_rho)+1-idx_rho)
        weight[...,1] = (idx_rho-np.floor(idx_rho))
        index_p = np.tile(np.arange(nr*nz,dtype=int), (2,1)).T
        index_rho = np.int_(np.dstack((np.floor(idx_rho),np.floor(idx_rho)+1)))
        
        weight[index_rho>=nrho] = 0
        index_rho[index_rho>=nrho] = nrho-1

        matrixes = []
        

        #ugly fix of the zero division in the core - weight for LFS and HFS in the core should by very similar
        #=> it should by OK
       
        r = np.hypot(R-R0,Z-Z0)+0.01
        J = (-Bz*(R-R0)+Br*(Z-Z0))/r  #BUG how is ist possible that J is sometimes negatice???
        J = np.hypot(J,0.1)
        
        #r = np.hypot(R-R0,Z-Z0)
        #J = (-Bz*(R-R0)+Br*(Z-Z0))/r  #BUG how is ist possible that J is sometimes negatice???    
        
        #weight*=(R/J)[:,:,None] #BUG!!!
        
        #BUg is it right? lower weight in the HFS?? 
        
   
        for it in xrange(nt):
            index_p   = index_p.ravel()
            irho = index_rho[it,...].ravel()
            w = weight[it,...].ravel()
            
            #w[w<0] = 1

            M = sp.csc_matrix((w,(index_p,irho)),shape=(nr*nz,nrho)).T
            norm = np.asarray(M.sum(1))[:,0]+1e-10  #borerct normalization of the flux surfaces 
          
            M = sp.diags(1/norm,0)*M #normalization 
            matrixes.append(M)
    
        return matrixes
            



    #def map_RZ2rho(self,shot,tin,r_grid,z_grid,Fgrid,rho_out,rho_lbl='rho_tor',diag='EQH'
                #,exp='AUGD',ed=0):
        #print 'map_RZ2rho'
        ##BUG not working 
        ##jak udělat matici, co to zobrazí?? pro jeden čas.. 
        ##calculate flux averadge from the quantity on the grid
        
        ##Need to get R,Z coordinates for constant flux 
        #tin = np.atleast_1d(tin)
        #Fgrid = np.atleast_3d(Fgrid)


        ##pořešit vstupní dimenze? 
        #nt = np.size(tin)
        #nrho = np.size(rho_out)
        #flux_surfaces = self.kkeqpsp(shot,tin,rho_out,diag,exp,ed,rho_lbl)
        #flux_surfs = []
        #times = []
        #indexes = []
        #i = 0
        #for t,FStime in zip(tin,flux_surfaces):
            #Tindexes = []
            #for FS in FStime:
                ##print len(FS)
                #flux_surfs.append(FS)
                #Tindexes.append(slice(i,i+len(FS)))
                #i+= len(FS)
                #times.append(np.ones(len(FS))*t)
            
            #indexes.append(Tindexes)
            
        ##calculate Br,Bz on these RZ points of the flux surfaces 

        #times = np.hstack(times)
        #flux_surfs = np.vstack(flux_surfs)
        #Rflux = flux_surfs[:,0][:,None]
        #Zflux = flux_surfs[:,1][:,None]

        #Br,Bz = self.kkrzbrzt(shot,times,Rflux,Zflux,exp,diag,ed)
        
        #import matplotlib.pylab as plt

        #Bp = np.squeeze(np.hypot(Br,Bz))
        #plt.plot(Bp)
        #plt.show()
        
        #dr = (r_grid[-1]-r_grid[0])/len(r_grid)
        #dz = (z_grid[-1]-z_grid[0])/len(z_grid)
        #scaling = np.array([dr,dz])
        #offset = np.array([r_grid[0]+dr/2,z_grid[0]+dz/2])
        
        #FSA = np.empty((nt,nrho))
        #print 3

        #for it,tslices in enumerate(indexes):
            #R_inter = []
            #Z_inter = []
            #dlt_b = []
            #for ir,rslice in enumerate(tslices):
                ##print it, ir, rslice
                #if rslice.start==rslice.stop:  #empty - surface was not found
                    #dlt_b.append(None)
                    #continue
                    

                #surfBp = Bp[rslice]
                #R = Rflux[rslice,0]
                ##print R.shape, R[0]
                #R_inter.append(R)
                #R = np.r_[R,R[0]]
                #Z = Zflux[rslice,0]
                #Z_inter.append(Z)
                #Z =  np.r_[Z,Z[0]]
                ##Now need dl/Bp for each

                #dlt_b.append(np.hypot(np.diff(R),np.diff(Z))/surfBp)
            
            ##print np.hstack(R_inter).shape, np.hstack(Z_inter).shape
            
            ##plt.plot(np.hstack(R_inter),np.hstack(Z_inter))
            ##plt.show()
            #coords = np.vstack((np.hstack(R_inter),np.hstack(Z_inter)))
            ##print coords.shape
            ##exit()
            ##plt.imshow(Fgrid[it,:,:])
            ##plt.show()
            #idx = (coords-offset[:,None]) / scaling[:,None]
            
            ##print  Fgrid[it,...].shape
            
            
            #Fsurf = map_coordinates(Fgrid[it,...],idx, mode='nearest', order=1,prefilter=True)
            ##print Fsurf.shape
            ##print sum([len(w) for w in dlt_b if w!= None])
            ##print sum([len(w) for w in R_inter])
            ##finally the flux surface averadge
            #for ir,w in enumerate(dlt_b):
                #if w == None:  #this flux surface was not found
                    #FSA[it,ir] = np.nan
                #else:
                    #FSA[it,ir] = np.average(Fsurf[:len(w)],weights=w)
                    #FSA[it,ir] = np.mean(Fsurf[:len(w)]) #BUG 

                    #Fsurf = Fsurf[len(w):]
                    


            ##print 'Fsurf',len(Fsurf)
            
            ##exit()
        #print 4

        ##solve the lost magnetic surfaces
        #nanind = np.where(np.isnan(FSA))
        #for it,ir in zip(nanind[0],nanind[1]):
            #ind = np.isfinite(FSA[it,:])
            #i = np.argmin(abs(rho_out[ind]-rho_out[ir]))
            #FSA[it,ir] = FSA[it,ind][i]
    
        ##import matplotlib.pylab as plt
        ##plt.plot(tin,FSA)
        ##plt.show()
        ##pořešit nany 
        #return FSA




    def kkeqpsp(self,shot,tin,rhoin,diag='EQH',exp='AUGD',ed=0,rho_lbl='rho_pol',dR=0,dZ=0):
        """Get R,Z coordinates of a flux surfaces
        
        Parameters
        ----------
    
        rhoin : 1darray,float
            rho coordinates of the searched flux surfaces
        tin : 1darray,float
            time 
            
        shot: int
            the shot number
        diag: str
            diagnsotics used for mapping (EQI,EQH,...)
        exp: str
            experiment (AUGD)
        ed:  int
            diag edition
        rho_lbl: str
            mapped coordinates - rho_pol or rho_tor
        dR:  float
            shift the  center of equilibrium to right by dr [m]
        dZ: float 
            shift the center of equilibrium up dz [m]
    
            
    
        Returns
        -------
        rho : list of lists of arrays [npoinst,2]
            list of times containg list of surfaces for different rho 
            and every surface is decribed by 2d array [R,Z]
        
        """
        stat = False
        if diag==self.diag and ed==self.ed and exp==self.exp and shot==self.shot:
            stat = True
        if not self.Open(diag,shot,experiment=exp,edition=ed):
            print 'no '+diag
            return 
        
  

        if not hasattr(self,'tvec') or not stat:  
            self.tvec = dd.GetSignal('time')
        tvec = np.copy(self.tvec)  
     
        #radial mesh (0:M) 
        if not hasattr(self,'Rmesh') or not stat:  
            self.Rmesh = dd.GetSignalGroup('Ri')
        Rmesh = np.copy(self.Rmesh)  
        
        #vertical mesh (0:N) 
        if not hasattr(self,'Zmesh') or not stat:  
            self.Zmesh = dd.GetSignalGroup('Zj')
        Zmesh = np.copy(self.Zmesh)  
       
        #poloidal flux label (0:Lpf) 
        if not hasattr(self,'PFM') or not stat:
            #print 'load PFM'
            self.PFM = dd.GetSignalGroup('PFM')
            self.PFM = self.PFM.astype(np.float32, copy=False)
        PFM = np.copy(self.PFM)  
        
        if not hasattr(self,'PFxx') or not stat:  
            self.PFxx = dd.GetSignalGroup('PFxx')
        PFxx = np.copy(self.PFxx)  
        #0,: is mag axis (tvec)
        #1.: is separatrix (tvec)
    

        rhoin = np.atleast_1d(rhoin)
        
        if rho_lbl=='rho_tor':
            rhoin = self.kkrhopto(rhoin,tin,shot,diag, exp, ed)
        elif rho_lbl == 'r_V':
            rhoin = self.kkrhoPVo(shot,tin,rhoin,diag, exp, ed,rho_lbl='rho_pol',strahl_def=False)
        elif rho_lbl == 'r_V_strahl':
            rhoin = self.kkrhoPVo(shot,tin,rhoin,diag, exp, ed,rho_lbl='rho_pol',strahl_def=True)
        else:
            rhoin = np.tile(rhoin, (np.size(tin),1))
        #print rhoin.shape
            
        orientation = np.sign(PFxx[0,:].mean())

        mag = PFxx[0,:]

        #find the edge value of the poloidal flux 
        ikCAT = np.argmax(PFxx[1:,:]*orientation,axis=0)+1
        
        sep = PFxx[ikCAT,np.arange(len(mag))]
        
        import matplotlib._cntr as cntr

        contours = []
        tin = np.atleast_1d(tin)
        
        Rmesh = Rmesh[:,1]
        nr = len(Rmesh)
        Zmesh = Zmesh[:,1]    
        nz =len(Zmesh)
  

        i0 = -1
        for it, t in enumerate(tin):        
            i = np.argmin(abs(tvec-t))  
            if i0 == i:  #already calculated in the previous step
                contours.append(contours[-1])
                continue
                
            i0 = i
            Flux = rhoin[it,:]**2*(sep[i]-mag[i])+mag[i]
        
            R,Z = np.meshgrid(Rmesh,Zmesh)

            #contour searching routine from matplotlib
            try:
                c = cntr.Cntr(R,Z, PFM[:nr,:nz,i].T)
            except Exception as e:
                print 'cntr.Cntr',R.shape, Z.shape, PFM[:nr,:nz,i].T.shape, PFM.shape
                raise e 
                
            rho_contours = []

            

            for r,f in zip(rhoin[it,:],Flux):
                nlist = c.trace(level0=f,level1=f,nchunk=0)
                lines = nlist[:len(nlist)//2]
                if len(lines) == 0:
                    #lines = [np.array((None,None),ndmin=2),]
                    lines = [np.empty((0,2)),]
                line = []
                #choose the longest line
                for l in lines:
                    if len(l)>=len(line):
                        line = l+np.array((dR,dZ))[None,:]*(1-r)
                rho_contours.append(line)    
                #rho_contours.append(np.vstack(lines))  #BUG use only the first line?? or the largest one? 
                
            contours.append(rho_contours)
            
        
        
        return contours

        #original version 
        import matplotlib.pylab as plt

        
        #print  time.time()-T
        #flux = rho**2

        T = time.time()

        contours2 = []
        KK = kk.KK()
        
        for it,t in enumerate(tin):
            print t
            it = np.argmin(abs(tvec-t))
            Flux = rhoin[it,:]**2*(sep[it]-mag[it])+mag[it]
            rho_contours = []

            for r,f in zip(rhoin[it,:],Flux):
                pto=KK.kkrhopto(shot,t,[r,],exp=exp,diag=diag,ed=ed)
                psp=KK.kkeqpsp(shot,t,pto.pf,exp=exp,diag=diag,ed=ed)
                ed = psp.ed
                
                cont = np.c_[psp.r_surf,psp.z_surf ]
            
                rho_contours.append(cont)
            contours2.append(rho_contours)
    
        
        for c1,c2 in zip(contours, contours2):
            for r1,r2,r in zip(c1,c2,rhoin[it,:]):
            
                plt.plot(r1[:,0],r1[:,1],'k',linewidth=0.3)
                plt.plot(r2[:,0],r2[:,1],'b--',linewidth=0.3)
                
        plt.xlim(Rmesh[:,it].min(),Rmesh[:,it].max())
        plt.ylim(Zmesh[:,it].min(),Zmesh[:,it].max())

                
        plt.show()

    
        
        print  time.time()-T


    def PlasmaVolume(self,shot,rhoin, tin, diag='EQH', exp='AUGD', ed=0,rho_lbl='rho_pol'):
        """calculate the volume unclosed by flux surface rhoin
        
        Parameters
        ----------
    
        rhoin : 1darray,float or 2d array (nt,nrho)
            rho coordinates of the searched flux surfaces
            
        tin : 1darray,float
            time 
            
        shot: int
            the shot number
        diag: str
            diagnsotics used for mapping (EQI,EQH,...)
        exp: str
            experiment (AUGD)
        ed:  int
            diag edition
        rho_lbl: str
            mapped coordinates - rho_pol or rho_tor

    
            
    
        Returns
        -------
        rho: poloidal flux surfaces where V was calculated 
        V : Volumes inside the flux surfaces
        rho_b: poloidal flux surfaces where dV was calculated (rho=rhoin_pol[1:])
        dV : dV/drho
        
        """
        stat = False
        if diag==self.diag and ed==self.ed and exp==self.exp and shot==self.shot:
            stat = True
        if not self.Open(diag,shot,experiment=exp,edition=ed):
            print 'no '+diag
            return 
        
       
        if not hasattr(self,'time') or not stat:  
            self.time = dd.GetSignal('time')
        tvec = np.copy(self.time)
        nt = len(tvec)
        
        #value of the poloidal flux at the different poloidal flux surfaces
        if not hasattr(self,'PFL') or not stat:  
            self.PFL = dd.GetSignalGroup('PFL')
        PFL = np.copy(self.PFL)[:,:nt]

        #Poloidal Flux values (0:LPFx) 
        if not hasattr(self,'PFxx') or not stat:  
            self.PFxx = dd.GetSignalGroup('PFxx')
        PFxx = np.copy(self.PFxx)[:,:nt]
        #0,: is mag axis (tvec)
        #1.: is separatrix (tvec)
        
        if not hasattr(self,'Vol') or not stat:  
            self.Vol = dd.GetSignalGroup('Vol')
        Vol = np.copy(self.Vol).T
     
        #dd.Close()
        
   
        V  = Vol[:nt,  ::2]
        dV = Vol[:nt, 1::2]
        nl = np.size(PFL,0)
        
        
        mag = PFxx[0,:]
        orientation = np.sign(PFxx[0,:].mean())

        #find the edge value of the poloidal flux 
        ikCAT = np.argmax(PFxx[1:,:]*orientation,axis=0)+1
        
        sep = PFxx[ikCAT,np.arange(len(mag))]

    
        rho = (PFL-mag)/(sep-mag)
        rho[rho <0] = 0
        rho = np.sqrt(rho).T

        ##BUG very uneffective!!

        
        
        rhoin = np.atleast_1d(rhoin)
        tin= np.atleast_1d(tin)
        #print rhoin.size
        
        if rhoin.ndim < 2:
            rhoin = np.tile(rhoin, (np.size(tin),1))
            
        #print rhoin.shape
        V_new = np.zeros((np.size(tin),np.size(rhoin,-1)))
        
        
        j0 = -1
        for i,t in enumerate(tin):
            j = np.argmin(abs(tvec-t))  
            if j0 == j:  #already calculated in the previous step
                V_new[i,:] = V_new[i-1,:]
                continue
            j0 = j
            
            ind = np.argsort(rho[j,:])
            ind = ind[rho[j,ind] <0.98]
            rho_tmp =  rho[j,ind]
            
            if rho_lbl == 'rho_tor':
                #print '->rho_tor'
                rho_tmp = self.kkrhopto(rho_tmp,t,shot,diag, exp, ed)
            if rho_lbl == 'r_V':
                #print '->r_V'
                rho_tmp = self.kkrhoPVo(shot,t,rho_tmp,diag, exp, ed,strahl_def=False)
            if rho_lbl == 'r_V_strahl':
                #print '->r_V_strahl'
                rho_tmp = self.kkrhoPVo(shot,t,rho_tmp,diag, exp, ed,strahl_def=True)
            if rho_lbl == 'rho_pol':
                pass

            #print rho[j,::-1]
            
            
            
            r = np.r_[0, np.squeeze(rho_tmp)]
            v = np.r_[0, V[j,ind]]
            V_new[i,:] = UnivariateSpline(r,np.sqrt(v),k=3, s=1e-2,bbox=[-1e-3, 1])(rhoin[i,...])**2

            #V_new[i,:] = interp1d(r,v,kind='cubic',bounds_error=False,fill_value=v.max())(rhoin)     
        #V_new[:,rhoin>1] = 0 #BUG proč to tu bylo??? 

        
        dV = np.diff(V_new,axis=1)#/np.diff(rhoin)[np.newaxis,:]
        #rho = rhoin[1:]
        dV[dV<= 0] = 0 
        
        #plt.plot(rhoin[1:],dV)
        #plt.plot(rhoin,V_new.T)

        #plt.show()
        return rhoin, V_new,rhoin[1:], dV
        
    
    def MapLOStoRho(self,shot,rhoin, tin,LOScoords, diag='EQH', exp='AUGD', ed=0,rho_lbl='rho_pol',dR=0,dZ=0):
        
        R,Z,Phi,Rorig,Zorig,phiorig = LOScoords
        
        R = np.linspace(1.7,1.8, 1000)
        Z = np.linspace(-1., 1.,1000)
        
        #print tin.shape, R.shape, Z.shape
        #print tin
        
        rhoLOS = self.kkrzptfn(shot,tin,R,Z,exp=exp,diag=diag,ed=ed,rho_lbl=rho_lbl,integ_time=0,dR=dR,dZ=dZ)
        
        
        
        
        dL = np.hypot(np.diff(R), np.diff(Z))
        rhoLOS_b = (rhoLOS[:,1:]+rhoLOS[:,:-1])/2
        
        dL_drho = dL/np.diff(rhoLOS,axis=1)
        #dL_drho[rhoLOS_b] = 0
        
        sigrho = (np.sign(dL_drho)*rhoLOS_b)[0,:]
        (np.sign(dL_drho)*rhoLOS_b)[0,:] = argsort(sigrho)
        
        #LOSweight=interp1d(,dL_drho,'linear',1,False,False,0)(rhoin)

        

        testdata = np.eye(len(rhoin))
        dL = np.hypot(np.diff(R), np.diff(Z))
        dL = (np.r_[0,dL]+np.r_[dL,0])/2
        
        
        LOSweight=interp1d(rhoin,testdata,'linear',1,False,False,0)(rhoLOS)
        LOSweight = np.sum(LOSweight*dL,2)
        
        
        
        

        
        
            
        
        
        #LOScoords  
        
    


def main():
    kk = eqi_map()
 
    
    
    #from matplotlib.pylab import *
    
    #import gc
    #gc.enable()
    #gc.set_debug(gc.DEBUG_LEAK)

    ###27242 - blbě
    ###26919 - blbě
    nr = 100


    shot=30132
    r_grid=np.linspace(-0.5,0.6,nr+1)   +1.65
    z_grid=np.linspace(-1,1,nr*2)
    
    rho = np.linspace(0,1.1,20)
    data_rho = (1-rho)**2
    #time = 2
    #map_rho2RZ(rho,data_rho,r_grid,z_grid,shot,time)
    import matplotlib.pylab as plt

    
    #shot = 28053
    ##ed=0
    diag='EQH'
    #exp='AUGD'
    rin=np.arange(-nr,nr)*.5/nr+1.778
    zin=np.linspace(-0.5,0.6,nr*2)
    #rin = zin*0+1.781 
    tin=np.linspace(1,6,20)
    tin = 5
    #tin = np.atleast_1d(tin)

    X,Y = np.meshgrid(rin,zin)
    rhopol = np.linspace(0,1,10)
    #print X.ravel().shape
    #print Y.ravel().shape¨
    
    rhopol = np.linspace(0.01,0.99,100)
    
    theta = np.linspace(0,2*np.pi,20)
    data_theta = np.sin(3*theta)
    
    
    #M1 = kk.map_theta2RZ_matrix(shot,tin,theta,r_grid,z_grid)
        
    
    
    
    nLOS = 24
    Z = [6.7715e-02,6.6442e-02,6.5170e-02,6.3897e-02,6.2625e-02,6.1352e-02,6.0080e-02,
    5.8807e-02,5.7535e-02,5.6263e-02,5.4990e-02,5.3718e-02,5.2445e-02,5.1173e-02,
    4.9900e-02,4.8628e-02,4.6083e-02,4.4811e-02,4.3538e-02,4.2266e-02,3.9721e-02,3.8448e-02,3.7176e-02,3.5903e-02]

    R = [1.7123e+00,1.7292e+00,1.7454e+00,1.7627e+00,1.7769e+00,1.7927e+00,1.8066e+00,1.8223e+00,1.8398e+00,
    1.8557e+00,1.8680e+00,1.8832e+00,1.8970e+00,1.9125e+00,1.9266e+00,1.9403e+00,1.9698e+00,1.9839e+00,
    1.9980e+00,2.0126e+00,2.0400e+00,2.0559e+00,2.0706e+00,2.0848e+00]

    Phi = [4.7431e+01,4.7058e+01,4.6712e+01,4.6348e+01,4.6055e+01,4.5737e+01,4.5463e+01,4.5159e+01,
    4.4828e+01,4.4532e+01,4.4309e+01,4.4035e+01,4.3793e+01,4.3525e+01,4.3287e+01,4.3057e+01,
    4.2578e+01,4.2354e+01,4.2135e+01,4.1911e+01,4.1500e+01,4.1268e+01,4.1057e+01,4.0856e+01]

    Rorig = [2.8640,]*nLOS
    Zorig = [2.6950e-01,]*nLOS
    phiorig = [3.4369e+02/180*np.pi,]*nLOS
    R,Z,Phi,Rorig,Zorig,phiorig = np.array(R),np.array(Z),np.array(Phi),np.array(Rorig),np.array(Zorig),np.array(phiorig)
    
    kk.kkeqrinv(shot,tin,rhopol,rho_lbl='rho_pol')
    print 'end'
    exit()
    
    kk.MapLOStoRho(shot,rhopol, tin,(R,Z,Phi,Rorig,Zorig,phiorig), diag='EQH', exp='AUGD', ed=0,rho_lbl='rho_pol')
    
    #M = kk.map_RZ2rho_matrix(shot,tin,rhopol,r_grid,z_grid)
            
    
    

    #rho = kk.kkrzptfn(shot,tin,rin,zin,rho_lbl='rho_tor',integ_time=0,dR = 0.01, dZ = 0.01)
    ##print rho.shape
    #print rho
    #print#rho = kk.kkrhopto( rhopol,tin,shot,diag='EQH', exp='AUGD', ed=0,integ_time=0)
    
    #plot(rho.T)
    #show()
    #M2 = kk.mdap_RZ2rho_matrix(shot,tin,rhopol,r_grid,z_grid,rho_lbl='rho_pol')

    #M1 = kk.map_rho2RZ_matrix(shot,tin,rhopol,r_grid,z_grid,rho_lbl='rho_pol')
    exit()

    

    rho = kk.kkrhoPVo(shot,tin,rho,diag='EQH', exp='AUGD', ed=0,rho_lbl='rho_tor')
    plot(rho.T)
    show()
    #exit()
    

    #Br,Bz = kkrzbrzt(shot,tin,X.ravel(),Y.ravel())
    #Br    = Br.reshape(len(tin),len(zin),len(rin))
    #Bz    = Bz.reshape(len(tin),len(zin),len(rin))
    #print Br

    ##plt.quiver( X, Y, Br[0,...], Bz[0,...], units='width')
    
    #plt.pcolor(X,Y, np.hypot(Bz,Br)[0,...])
    #plt.show()
    #plt.imshow( Br[0,...])
    #plt.show()
    #plt.imshow( Bz[0,...])
    #plt.show()
    
    #exit()
    
    ####print tin.shape
    #rin = np.load('rin.npy')
    #zin = np.load('zin.npy')
    ####tin = np.load('tin.npy')
    ####print rin
    ####print zin
    ####print tin
    ####print rin.shape, zin.shape,tin.shape
 
    exit()
    map_RZ2rho_matrix

    exit()
    
    kkeqpsp(shot,np.linspace(3,4,10),rhopol,  diag='EQH', exp='AUGD', ed=0,rho_lbl='rho_pol')
    data_rho = (1-rhopol**2)**.5+tin[:,None]/10-0.5
    #data_rho = data_rho.T
    
    #print data_rho.shape
    #exit( )
    #plt.plot(rhopol,data_rho)
    
    ##exit()
    ##R,Z,F = map_rho2RZ(shot,tin,rhopol,data_rho,r_grid,z_grid,rho_lbl='rho_pol')
    ###print F.shape, F[10,...].shape
    ##plt.imshow(F[10,...] )
    ##plt.show()
    

    ###print F.shape
    ###exit()
    ##print tin.shape, r_grid.shape, z_grid.shape, F.shape
    
    R,Z = np.meshgrid(r_grid,z_grid)
    print 'plot'
    #print np.dot(M2[0].todense(),R.ravel()).shape, rhopol.shape
    #R02 = np.dot(M2[0].todense(),R.ravel()).T
    R02 = np.dot((M1[0]*M2[0]).todense(),R.ravel()**2).reshape(2*nr, 2*nr)
    R02[R02==0]= np.nan
    contours = kkeqpsp(shot,tin,np.linspace(0,1,10))
    extent = [R.min(),R.max() ,Z.min(),Z.max()]
    #im = plt.imshow(np.squeeze(R-R0)**2,extent=extent,aspect='auto')
    im = plt.imshow((np.array(R**2-R02))[::-1,:]/1.65**2,extent=extent,aspect='auto',interpolation='nearest')

    plt.colorbar(im)
    
    #print R0.shape, (R-R0).shape, 

    for c in contours:
        for cr in c:
            #print np.size(cr), np.shape(cr)
            if np.size(cr)>1:
                plt.plot(cr[:,0],cr[:,1],c='w')
    
 
    plt.show()

    ##print M2[0].shape,M1[0].shape 
    #plt.imshow((M2[0]*M1[0]).todense(),aspect='auto',interpolation='nearest')
    #plt.show()

    #print data_rho.shape, FSA.shape, data_rho.shape
    #plt.plot(tin,FSA,'r')
    #plt.plot(tin,data_rho,'k')

    #plt.show()
    #exit()
    
    
    #####print rhopol.shape
    #kkrhopto(rhopol,tin,shot,diag=diag)
    
    #for shot in range(30200,30000,-1):
        
        #if not dd.Open('EQH',shot):
        
            ##print 'no eqh'
            #continue 
        
        #PFM = dd.GetSignalGroup('PFM')
        #tvec_=dd.GetSignal('time')#nedávat edition of KK vždy na nulu 
        #ixti = np.int_(dd.GetSignal('ixti')-1)

        #dd.Close()

        #PFM = PFM.astype(np.float32, copy=False)
        #PFM = PFM[:,:,ixti]
        #f = interp1d(tvec_,PFM,copy=False,kind='nearest')
        #PFM  = f(tvec_)
        #del f
        
     #kkrzptfn(shot,tin,rin,zin,exp='AUGD',diag='EQH',  ed=0,rho_lbl='rho_pol'):
    ###plt.plot(rho)
    ###plt.show()
    #dump_garbage()
    #PlasmaVolume(shot,rho, tin, diag='EQH', exp='AUGD', ed=0,rho_lbl='rho_pol')
    #kkrhotpq(shot,tin,rho,diag='EQH', exp='AUGD', ed=0,rho_lbl='rho_tor')
    #dump_garbage()


if __name__ == "__main__":
    main()
#ikCAT test

#time grid not equally spaced 

