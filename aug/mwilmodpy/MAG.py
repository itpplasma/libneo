import numpy as np
import scipy.interpolate
from IPython import embed
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import pickle
from matplotlib import gridspec
from matplotlib import colorbar as cb
from matplotlib import colors

import tools
import LSQFFT


class MAGhelp:
    status = False

class MAG:
	def __init__( self ,  shot = None,Experiment = 'AUGD' ):
		self.Status = False
                self.isCorrected = False
		if shot != None:
			self.LoadAUG( shot )

	def __del__( self ):
		self.Unload( )
		self.Status = False
	

	def Unload( self ):
		if self.Status:
			del self.time
	

        def LoadAUG( self ,  shot, Experiment='AUGD', Diagnostic='MHA', Edition = 0L, tBegin=1.0, tEnd=11.0,usedArray='toroidal',calibDistance=False,calibN=2):
            		
                import dd
            
                self.Unload()


                if shot<33739:
                    Diags=['MHA','MHB','MHC','MHD']
                else:
                    Diags=['MHI']
                    

                if shot > 21496:
                    calib_shot =dd.getLastShotNumber('CMH', shot)  
                    sfCalib = dd.shotfile('CMH',calib_shot)
                else:
                    sfCalib = dd.shotfile('MHA',shot)  
               
                
                if usedArray == 'toroidal':
                    phase_balooning_coils = [1,2,3,12,13,14,40]
                    sig_names17 = ['CB17-02']
                elif usedArray == 'poloidal':
                    phase_balooning_coils= [5,6,7,8,9,10,2]
                else:  
                    print 'no choice'

                sig_names = ['CB31-%.2d'%c for c in phase_balooning_coils]
                if usedArray== 'toroidal':
                    sig_names.append ('CB17-02')

                self.phi = [sfCalib.getParameter(n,'phi').data for n in sig_names]   
                self.R = [sfCalib.getParameter(n,'R').data for n in sig_names]   
                self.z = [sfCalib.getParameter(n,'z').data for n in sig_names]  
                c_names = sfCalib.getObjectNames()  


                sfCalib.close()   
                #embed()
                useIdx=np.where(np.array(self.phi)!=0)[0]
                self.idxPhiSort=useIdx[np.argsort(np.array(self.phi)[useIdx])]
                self.phi = np.rad2deg(np.squeeze(self.phi)[self.idxPhiSort])+22.5 #correction due to different coor
                self.R = np.squeeze(self.R)[self.idxPhiSort]
                self.z = np.squeeze(self.z)[self.idxPhiSort]

                self.sig_names =  np.squeeze([n[1:] for n in sig_names if n[:2] == 'CB' and n[-3] == '-'])[self.idxPhiSort]
            
                
                
                
                intSig,sig,lbl=[],[],[]
                for di in Diags:
                    try:
                        sfMag = dd.shotfile( di,shot)

                        sf_names = sfMag.getSignalNames()
                        sf_names.sort()
                        self.time = sfMag.getTimeBase(sf_names[0],tBegin=tBegin,tEnd=tEnd)
   #                 nbeg, nend = np.array(time).searchsorted((tBegin,tEnd))

                        for n in self.sig_names:
                            if n in  sf_names:
                                sig.append(sfMag.getSignalCalibrated(n,tBegin=tBegin,tEnd=tEnd)[0])
                                lbl.append(n)

                        sfMag.close()
                    except:
                        print di,' not available'

                for n in self.sig_names:
                    if not n in lbl:
                        idx=np.where(self.sig_names==n)[0]
                        self.phi = np.delete(self.phi ,idx)
                        self.R = np.delete(self.R ,idx)
                        self.z = np.delete(self.z ,idx)
                        self.sig_names = np.delete(self.sig_names ,idx, axis=0)


                    
                #embed()
                
                if calibDistance:
                    import kk_mwillens as kk
                    useN=calibN
                    rhopIn=np.linspace(np.sqrt(0.95),1.0,30)
                    outq=-kk.KK().kkrhopfq(shot,self.time.mean(),rhopIn).q
                    useQ = np.ceil(outq[0])
                    rhopQrat = np.interp(useQ,outq,rhopIn)
                    fluxsurf=kk.KK().kkrhorz(shot,self.time.mean(),[np.squeeze(rhopQrat)],angle=np.linspace(0.,360.,100))
                    dist=np.zeros_like(self.R)
                    for i in np.arange(dist.size):
                        dist[i] =np.sqrt((fluxsurf.r-self.R[i])**2.+(fluxsurf.z-self.z[i])**2.).min()
                    
                    calib = (dist/dist.min())**(useQ*useN+1.)
                    print 'apply corrections following for corrrections for n=%d'%useN+' at q=%d : '%useQ,calib
                else:
                    calib=np.ones_like(self.R)

                
                sfMAC = dd.shotfile( 'MAC' , shot)
                MAC_Ipol = -sfMAC('Ipolsola', tBegin=tBegin, tEnd=tEnd).data
                MAC_time = sfMAC('Ipolsola', tBegin=tBegin, tEnd=tEnd).time
                sfMAC.close()
                self.MAC_time,self.MAC_Ipol = MAC_time,MAC_Ipol


                self.sig=(np.squeeze(sig).T*calib).T
                self.intSig=np.zeros_like(self.sig[:,1:])
                for s in np.arange(self.sig[:,0].size):
                    self.intSig[s]=integrate.cumtrapz(sig[s],self.time)

                self.phiRad=np.deg2rad(self.phi-360.)
                self.machine = 'AUG'
                self.shot = shot
                self.tmin,self.tmax = tBegin,tEnd 
                self.Status = True


        def LoadD3d( self ,  shot, path='/ptmp1/work/mwillens/dataDIII-D/BetaScanN1_135__/' , tBegin=1.0, tEnd=7.0,usedArray='toroidal',option=1,integratedSig=False):

            if option == 1:
                if usedArray == 'toroidal':
                    file="TorMag_%d.pkl"
          
                torMag = pickle.load( open( path+file%shot, "rb" ) )

                sigIn=[]
                phi=[]
                time=[]

                for i in np.arange(len(torMag)):
                    sigIn.append(torMag[i]['Bdot'])
                    phi.append(torMag[i]['coor'][2])
                    time.append(torMag[i]['time'])

            elif option == 2:
                if usedArray == 'toroidal':
                    if  integratedSig:
                        file="MDSTorMagF_%d.pkl"
                    else:
                        file="MDSTorMag_%d.pkl"

                torMag = pickle.load( open( path+file%shot, "rb" ) )

                sigIn=[]
                phi=[]
                time=[]
                
                for k in torMag.keys():
                    if ((np.size(torMag[k]['data'])>0) & (np.std(torMag[k]['data'])>0.) ):
                        sigIn.append(torMag[k]['data'])
                        phi.append(np.float(k[6:-1]))
                        time.append(torMag[k]['time']/1.e3)
               

            phi = np.squeeze(phi)
            phi[phi>180.] -= 360.
            idxPhi=np.argsort(phi)
            phi = np.squeeze(phi)[idxPhi]
            lbl = np.squeeze(torMag.keys())[idxPhi]

            tmin,tmax=time[0].min(),time[0].max()
            useIdx=0
            for i in np.arange(len(phi)):
                if (tmin < time[i].min()) & (tmax > time[i].max()) :
                    tmin = time[i].min()
                    tmax = time[i].max()

            useIdx = i

            newTime  = time[useIdx]    
            sig=[]

            for i in np.arange(len(phi)):
                sig.append(np.interp(newTime,time[i],sigIn[i]))

            sig=np.squeeze(sig)[idxPhi]

            intSig = np.zeros_like(sig)[:,:-1]

            
            if  integratedSig:
                intSig = sig 
            else:
                for i in np.arange(len(phi)):
                    intSig[i] = integrate.cumtrapz(sig[i],newTime)
           
            idxWeak = []    
            #check for weak probes and remove weak [robes                
            if ((np.floor(shot/100) == 1357)  ):
                idxWeak=np.where(phi == (277-360))[0]
                print 'discard ', lbl[idxWeak]

           # std=np.std(sig,axis=1)
           # idxWeak=np.where(std<(std.mean()/3.*2.))[0]
            #embed()   
           # if np.size(idxWeak)>0:
           #     print 'discard ', lbl[idxWeak]

            self.intSig=np.delete(np.squeeze(intSig[:,1:]),idxWeak,axis=0)
            if  integratedSig == False:
                self.intSig=np.delete(np.squeeze(intSig),idxWeak,axis=0)
                self.sig=np.delete(np.squeeze(sig),idxWeak,axis=0)
            self.phi = np.delete(np.squeeze(phi),idxWeak)
            self.lbl = np.delete(np.squeeze(lbl),idxWeak)
            self.shot=np.squeeze(shot)
            self.time=newTime
            self.phiRad=np.deg2rad(self.phi)
            self.machine = 'D3D'
            self.tmin,self.tmax=tmin,tmax
            self.Status = True


        def getModeAmplitudes(self, baseFreq=1.,intFreq=10.,sign=1):
            if self.Status:
                self.doBaselineCorrection(baseFreq=baseFreq,intFreq=intFreq,sign=sign)

                
                self.LS = LSQFFT.LSQFFT()
                self.LS.initializeData(self.phiRad,self.corrected, freq_in = 1/(2.*np.pi),order=3,poly=0,negSign=False)
                self.LS.leastSquareFit()

                

                embed()
        
               
        def doBaselineCorrection(self, baseFreq=1.,intFreq=20.,sign=1,doPlot=False):
            if self.Status:
                timeBase, intSigBase = tools.dataBinning( self.time[:-1], self.intSig, samplefreq = baseFreq )
                slowtime,  intSig = tools.dataBinning( self.time[:-1], self.intSig, samplefreq = intFreq )
                
                correct=[]
                for i in np.arange(np.size(intSigBase[:,0])):
                    correct.append( ( intSig[i] - np.interp(slowtime,timeBase,  intSigBase[i] ) )* sign)
#                    
                self.corrected=np.squeeze(correct)   
                self.timecorr=np.squeeze(slowtime) 
                self.isCorrected = True
                if doPlot:
                    #embed()

                    binnedtime,  binnedSig = tools.dataBinning( self.time, self.sig, samplefreq = intFreq )
                    
                    idx=np.abs(intSigBase).max(axis=1).argmin() 

                    labelsize=22
                    ticksize=18
                    plt.figure(figsize=(9, 9))
                    ax1=plt.subplot(311)
                    plt.title('%s'%self.sig_names[idx]+', filtered to %.0f kHz'%intFreq+', %d'%self.shot,fontsize=labelsize)
                    plt.plot(binnedtime, binnedSig[idx] )
                    plt.ylabel(r'$\rm{dB_r/dt\ [T/s]}$',fontsize=labelsize)
                    plt.xticks(fontsize=ticksize)
                    plt.yticks(fontsize=ticksize)
                    plt.xlim([1.,8])
                    ax2=plt.subplot(312,sharex=ax1)
                    plt.plot(slowtime, intSig[idx]*1.e3 )
                    plt.plot(timeBase, intSigBase[idx]*1.e3,lw=2 )
                    plt.ylabel(r'$\rm{B_r\ [mT]}$',fontsize=labelsize)
                    plt.xticks(fontsize=ticksize)
                    plt.yticks(fontsize=ticksize)
                    plt.xlim([1.,8])
                    ax3=plt.subplot(313,sharex=ax1)
                    plt.plot(self.timecorr, self.corrected[idx]*1.e3 )
                    plt.xlabel(r'$\rm{time\ [s]}$',fontsize=labelsize)
                    plt.ylabel(r'$\rm{\delta B_r\ [mT]}$',fontsize=labelsize)
                    plt.xticks(fontsize=ticksize)
                    plt.yticks(fontsize=ticksize)
                    plt.xlim([1.,8])
                    plt.tight_layout()
                    plt.show()

                    labelsize=22
                    ticksize=18
                    plt.figure(figsize=(9, 6))
                    ax1=plt.subplot(211)
                    plt.title('%s'%self.sig_names[idx]+', filtered to %.0f kHz'%intFreq+', %d'%self.shot,fontsize=labelsize)
                    plt.plot(binnedtime, binnedSig[idx] )
                    plt.ylabel(r'$\rm{dB_r/dt\ [T/s]}$',fontsize=labelsize)
                    plt.xticks(fontsize=ticksize)
                    plt.yticks(fontsize=ticksize)
                    plt.xlim([1.,8])
                    ax3=plt.subplot(212,sharex=ax1)
                    plt.plot(self.timecorr, self.corrected[idx]*1.e3 )
                    plt.xlabel(r'$\rm{time\ [s]}$',fontsize=labelsize)
                    plt.ylabel(r'$\rm{\delta B_r\ [mT]}$',fontsize=labelsize)
                    plt.xticks(fontsize=ticksize)
                    plt.yticks(fontsize=ticksize)
                    plt.xlim([1.,8])
                    plt.tight_layout()
                    plt.show()
                    
                    #embed()
                    
        def markELMPos(self):

            if self.Status:
                if (self.machine=='AUG'):
                    import dd
                    ELM = dd.shotfile('ELM',self.shot)
                    tbeg = ELM('t_begELM')
                    tend = ELM('t_endELM').data

                if self.isCorrected==False:
                    self.doBaselineCorrection()

                idxoELM = np.zeros_like(self.timecorr,dtype='bool')
                
            
            #embed()
                timeoELMmax = []
                sigoELMmax = []
                timeoELMmin = []
                sigoELMmin = []
                sigoELMmaxy = []
                sigoELMminy = []
                timeoELMmaxy = []
                timeoELMminy = []
                for i in np.arange(tbeg.size):
                    idxt = np.where((self.timecorr > tbeg[i]-0.001) & (self.timecorr < tend[i]))[0]
                    #              idxoELM[ (self.timecorr > tbeg[i]) & (self.timecorr < tend[i])   ] = 1
                #self.corrected[self.idxPhiSort],self.corrected[self.idxThetaSort]
                    if idxt.size > 0:
                        idxArrMax = np.unravel_index(self.corrected[self.idxPhiSort][:,idxt].argmax(), self.corrected[self.idxPhiSort][:,idxt].shape)
                        idxArrMin = np.unravel_index(self.corrected[self.idxPhiSort][:,idxt].argmin(), self.corrected[self.idxPhiSort][:,idxt].shape)
                        timeoELMmax.append(self.timecorr[idxt[idxArrMax[1]]])
                        sigoELMmax.append( self.phi[self.idxPhiSort][idxArrMax[0]])
                        
                        timeoELMmin.append(self.timecorr[idxt[idxArrMin[1]]])
                        sigoELMmin.append( self.phi[self.idxPhiSort][idxArrMin[0]])

                timeoELMmax = np.squeeze(timeoELMmax)
                sigoELMmax  = np.squeeze(sigoELMmax)
                timeoELMmin =np.squeeze( timeoELMmin)
                sigoELMmin  =np.squeeze( sigoELMmin)
                



            return timeoELMmax ,sigoELMmax, timeoELMmin, sigoELMmin





        def plotIntWiBase(self, baseFreq=1.,intFreq=10.,doPlot=False,plotTime=None,vmin=None,vmax=None,cmap='gnuplot',sign=1, addDiv=False,markELM=False):

            self.doBaselineCorrection(baseFreq=baseFreq,intFreq=intFreq,sign=sign)
            
            if doPlot:
                
                if plotTime != None:
                    idxt=np.where((self.timecorr>plotTime[0])&(self.timecorr<plotTime[1]))[0]
                else:
                    idxt=np.where((self.timecorr>self.tmin)&(self.timecorr<self.tmax))[0]

                ticksize=18
                labelsize=22

                plt.rcParams['axes.formatter.useoffset'] = False
                
                fig = plt.figure(figsize=(9,7))
   #             plt.title('%s'%self.machine+', LFS midplane, %d'%self.shot)
                ny=10
                gs = gridspec.GridSpec(4,ny) 

                if (addDiv) & (self.machine=='AUG'):    
                    yrange=1
                else:
                    yrange=0
                    
                ax1=fig.add_subplot(gs[yrange:,:-1])

               
                plt.xlim([self.timecorr[idxt].min(),self.timecorr[idxt].max()])
                extendCorr = np.vstack([self.corrected[0],self.corrected,self.corrected[-1]])

                if self.machine== 'D3D':
                    extendPhi = -1.*np.hstack([self.phi[0]-np.diff(self.phi)[0]/4.,self.phi, self.phi[-1]+np.diff(self.phi)[-1]/4.])
                else:
                    extendPhi = np.hstack([self.phi[0]-np.diff(self.phi)[0]/4.,self.phi, self.phi[-1]+np.diff(self.phi)[-1]/4.])

                if (vmin != None) & (vmax != None):
                    v=np.linspace(vmin,vmax,300)*1.e3
                else:
                    vmax=np.abs([extendCorr.min(),extendCorr.max()]).min()
                    vmin=-vmax 
                    v=np.linspace(vmin,vmax,300)*1.e3
                    
                plt.contourf(self.timecorr[idxt],extendPhi,extendCorr[:,idxt]*1.e3,v,cmap=cmap,extend='both')

                self.plotx,self.ploty,self.plotz,self.plotv= self.timecorr[idxt],extendPhi,extendCorr[:,idxt]*1.e3,v
                
                for p in self.phi:
                    if self.machine== 'D3D':
                        plt.hlines(-p,self.timecorr[idxt].min(),self.timecorr[idxt].max(),lw=0.5,linestyle='--',color='gray')
                    else:
                        plt.hlines(p,self.timecorr[idxt].min(),self.timecorr[idxt].max(),lw=0.5,linestyle='--',color='gray')


                if markELM:
                    timeoELMmax ,sigoELMmax, timeoELMmin, sigoELMmin = self.markELMPos()
                    plt.plot(timeoELMmax ,sigoELMmax,'ko',markersize=8)

                plt.xlabel(r'$\rm{time\ [s]}$',fontsize=labelsize)
                if self.machine== 'D3D':
                    plt.ylabel(r'$\rm{-\ toroidal\ \phi\ [deg]}$',fontsize=labelsize)
                else:
                    plt.ylabel(r'$\rm{toroidal\ \phi\ [deg]}$',fontsize=labelsize)

                plt.xticks(fontsize=ticksize)
                plt.yticks(fontsize=ticksize)


                ax1b=fig.add_subplot(gs[yrange:,-1])   
                norm = colors.Normalize(vmin=vmin*1e3,vmax=vmax*1e3)
                bounds = np.arange(np.ceil(vmin*1e3), np.floor(vmax*1.e3)+1.,1.)
                cbar = cb.ColorbarBase(ax1b,norm=norm,extend='both',ticks=bounds,cmap=cmap)

               # cbar.set_ticklabels(cbarlabels)
                cbar.set_label(r'$\delta\ B\ \rm{[mT]}$',fontsize=labelsize)
                for t in  cbar.ax.get_yticklabels():
                    t.set_fontsize(ticksize)


                if (addDiv) & (self.machine=='AUG'):
                    ax0 = fig.add_subplot(gs[0,:-1],sharex=ax1)
                    #                   ax0 = plt.axes([0.0, 0.8, 1.0, 0.2],sharex=ax0)
                    ax0.set_ylabel(r'$\rm{[kA]}$',fontsize=labelsize)
                    ax0.set_title('%s'%self.machine+', LFS midplane, %d'%self.shot)
                    idxtMAC = np.where((self.MAC_time>self.timecorr[idxt].min())&(self.MAC_time<self.timecorr[idxt].max()))[0]
                    MAC1_time, MAC1_Ipol = tools.dataBinning( self.MAC_time[idxtMAC], self.MAC_Ipol[idxtMAC]/1.e3, samplefreq = intFreq )
                    ax0.plot(MAC1_time, MAC1_Ipol,'k-',linewidth=1.5)
                    plt.yticks(np.arange(np.ceil(MAC1_Ipol.min()/5.)*5.,MAC1_Ipol.max(),5),fontsize=ticksize)
                    plt.ylim([MAC1_Ipol.min()*0.9,MAC1_Ipol.max()*1.1])
                    for tt in  ax0.get_xticklabels():
                        tt.set_fontsize(fontsize=0)

                    ax0.set_xlabel('')
                    plt.xlim([self.timecorr[idxt].min(),self.timecorr[idxt].max()])
                else:
                    ax1.set_title('%s'%self.machine+', LFS midplane, %d'%self.shot)
 #               gs.update(hspace=0.0)
                      
                

                #plt.rc_context(rc={'axes.formatter.offset_threshold' : 2}) 
                plt.tight_layout()
                plt.show()



                #embed()
        
        def MakeHistogram(self, baseFreq=1.,intFreq=20.,doPlot=True,timeRange=None,shot=34424,MPtime=3.050,pathNem='/tokp/work/mwillens/VMEC/AUGD/34424/micdu_eqb_2_3.0/nemec/',filename='wout.3d_rmp_n2_z64_p107',sign=1):
            if self.Status:

                self.doBaselineCorrection(baseFreq=baseFreq,intFreq=intFreq,sign=sign)
            
      
                if timeRange != None:
                    tBegin = timeRange[0]
                    tEnd = timeRange[1]
                else:
                    tBegin= self.tmin
                    tEnd = self.tmax
       
                import nemec_slim as nemec
                nem = nemec.nemec(folder=pathNem,filename=filename)
                nem.read()
                phix=np.linspace(-np.pi/2,np.pi/2,361,endpoint=True)
                Rsym,zsym,vecNorm,corrLen = nem.getSurfaceCorrugation(rhoIn=[1.0],phiIn=phix,nTheta=256,exact=True)
                idxtheta=  np.argmin( np.sqrt((np.squeeze(Rsym)-self.R[0])**2.+(np.squeeze(zsym)-self.z[0])**2.))
                Rsymp,zsymp,vecNormp,corrLenp = nem.getSurfaceCorrugation(rhoIn=[1.0],phiIn=phix*3.,nTheta=256,exact=True)

                clen=corrLen[:,idxtheta][:,0]
                idxpos=np.where(np.gradient(clen)>0)[0]
                phiZero=phix[idxpos[np.argmin(np.abs(clen[idxpos]))]]
                
                timeoELMmax ,sigoELMmax, timeoELMmin, sigoELMmin = self.markELMPos()
                
                import MAWnew as MAW
                MW=MAW.MAW(shot)
                
                phi3D=[]
#for tBegin1& tEnd1:
                for p in self.phi:
                    idxp=np.where((sigoELMmax == p))[0]
                    idxt=idxp[np.where((timeoELMmax[idxp]> tBegin) & (timeoELMmax[idxp] < tEnd))[0]]
                    phi3D=np.append(phi3D,MW.getSynPhase(timeoELMmax[idxt],tref=MPtime,phiDiag=p,usePSL=True))

                phi3Dshift= phi3D- phiZero
                phi3Dtrunc = np.copy(phi3Dshift)
                phi3Dtrunc[phi3Dshift<-np.pi/2.] += np.pi
                phi3Dtrunc[phi3Dshift>np.pi/2.] -= np.pi
                
                if doPlot:
                    
                    plt.figure(figsize=(9, 8))
                    ax1=plt.subplot(211)
                    plt.xlim([-90,90])
                    plt.plot(np.rad2deg(phix*3.-phiZero),corrLenp[:,idxtheta]*1000,'b',lw=2)
                    plt.ylim([-9,9])
                    plt.axhline([0.0],color='black',linestyle='dashed')
                    plt.ylabel(r'$\rm{[mm]}$',fontsize=18)
                    #plt.xlabel(r'$\rm{toroidal\ coordinate\ \phi\  [deg]}$',fontsize=18)
                    plt.title('Corrugation',fontsize=19)
                    plt.xticks([-90,-60,-30,0,30,60,90],fontsize=16)
                    plt.yticks(fontsize=16)
                    ax2=plt.subplot(212,sharex=ax1)
                    a,bins,patches=plt.hist(np.rad2deg(phi3Dtrunc), bins=15, normed= True,linewidth=1,range=[-90,90],facecolor='green',histtype='bar')
                    #(n, s) =np.rad2deg( norm.fit(phi3Dtrunc))
                    #h = mlab.normpdf( bins, n, s)
                    #plt.plot(bins, h, 'r--', linewidth=2)
                    plt.xlim([-90,90])
                    plt.title('probability density of ELMs, %d'%shot,fontsize=14,x=0.25,y=0.9)
                    plt.ylabel(r'$\rm{pdf}$',fontsize=18)
                    plt.xlabel(r'$\rm{toroidal\ coordinate\ \phi\  [deg]}$',fontsize=18)
                    plt.xticks([-90,-60,-30,0,30,60,90],fontsize=16)
                    plt.yticks(fontsize=16)
                    plt.tight_layout()
                    plt.savefig('%d'%shot+'_%.3f__sign+1.png'%MPtime)
                    plt.clf()
                    plt.close()
                    #plt.show()
                    
                    #embed()
                    

shotList=[34424,34852,34427,34852,34424]
timeList=[[2.5,4.5],[4.7,6.5],[2.5,4.5],[2.5,4.5],[4.9,6.9]]
i=0

path=['/tokp/work/mwillens/VMEC/AUGD/34424/micdu_eqb_2_3.0/nemec/',
      '/tokp/work/mwillens/VMEC/AUGD/34424/micdu_eqb_2_3.0/nemec/',
      '/tokp/work/mwillens/VMEC/AUGD/34424/micdu_eqb_2_3.0/nemec/',
      '/tokp/work/mwillens/VMEC/AUGD/34424/micdu_eqb_2_3.0/nemec/',
      '/tokp/work/mwillens/VMEC/AUGD/34424/micdu_eqb_2_5.1/nemec/']
MPtime=[ 3.050, 5.689, 3.092, 3.482,  5.432 ]
shift =[ 0.004,   0.001,   0.004,  0.001,  0.001  ]
pan = ['(a)','(b)','(c)','(d)','(e)']
filename=['wout.3d_rmp_n2_z64_p107',
          'wout.3d_rmp_n2_z64_p167',
          'wout.3d_rmp_n2_z64_m163',
          'wout.3d_rmp_n2_z64_m103',
          'wout.3d_rmp_n2_z64_p18']