#!/usr/bin/env python

import numpy as np
from scipy.optimize import leastsq
from IPython import embed
import matplotlib.pyplot as plt
from scipy.interpolate import splev,splrep,interp1d
from scipy import stats
from scipy.interpolate import UnivariateSpline

import dd
import MAWnew as MAW
import kk_mwillens as kk
import LSQFFT
import ECE
import REF
import CXRS
import LIN
import ElmSync
import pedestalFits
import HEB

def shiftFuncResi( shift, x,  f, ferror,coeff ):
    return  (f - splev(x-shift,coeff))/ferror


class shift:
    ###methods = 2line,shift,interpol for 2 line and shift a sepvalues is needed
    def __init__( self, shot=None,timeRef=None,timeWindow=[1.,5.],method='interpol',diag='ECE',diagExp='AUGD',useELMs=True,ELMperiod=[0.1,0.8],ELMExp='AUGD',fitNeParams= [0.8e19,1.8e19,0.05],fitTiParams= [100.,1200.,3.0e16,0.1],fitTeParams= [100.,1200.,0.1],neTrace= [1.0e19],neTraceRPS_LFS= None, TiTrace= [200.],TeTrace= [200.], CXRSuse = 'inte', doPlot=False):

        self.Status = False  
        self.diagsRead = 0
        if (np.array(method).shape <= 1) |  (np.array(method).shape == ()):
            if (method == 'spline') | (method == '2line') | (method == 'interpol')| (method == 'shine'):
                self.method = np.array([method])
            elif (method == 'all'):
                self.method  = np.array(['interpol','spline','2line'])
            else:
                self.method = ['interpol']
        else:
            self.method = np.array(method)

        #first is lower threshold for density, second width
        if np.size(fitNeParams) == 3:
            self.fitNeParams= fitNeParams
        else:
            print 'no valid input for fitNeParams'
            return
        #first is lower threshold for CXRS intensity, second width
        if np.size(fitTiParams) == 4:
            self.fitTiParams= fitTiParams
        else:
            print 'no valid input for fitTiParams'
            return

        if np.size(fitTeParams) == 3:
            self.fitTeParams= fitTeParams
        else:
            print 'no valid input for fitTeParams'
            return

        ##c
        self.neTrace= neTrace
        self.TiTrace= TiTrace
        self.TeTrace= TeTrace

        if (np.array(diag).shape <= 1) |  (np.array(diag).shape == ()):
            if (diag == 'ECE') | (diag == 'FRS') | (diag == 'RPS_LFS')| (diag == 'RPS_HFS')  | (diag == 'LIN') | (diag == 'CPZ') |  (diag == 'CMZ')|  (diag[:3] == 'RIC')  :
                self.diag = np.array([diag])
            elif (diag == 'all'):
                self.diag = np.array(['RMD','LIN','CPZ','CMZ','FRS','RPS_LFS','RPS_HFS','RIC1','RIC4','RIC8'])
            else:
                self.diag = np.array(['RMD'])
        else:
            self.diag = np.array(diag)
            
        self.ndiag = np.size(self.diag)
      
  
        if ((CXRSuse == 'Ti') | (CXRSuse == 'vrot') | (CXRSuse == 'inte')):
            self.CXRSuse = CXRSuse

        if np.all(timeRef)==None:
            self.timeRef = np.array([np.min(timeWindow),np.min(timeWindow)+0.05])
        else:
            self.timeRef = np.array([np.min(timeRef),np.max(timeRef)])

        self.timeWindow =  np.array([np.min(timeWindow),np.max(timeWindow)])
        self.useELMs = useELMs
        if self.useELMs:
            self.ELMperiod = ELMperiod

        #embed()

        if np.all(shot)==None:
            print 'no shot number'
            self.Status = False
        else:
            self.shot = shot
            self.Status = True
            self.loadData()
            if self.diagsRead > 0:
                self.calcDisplacement(plot=doPlot)
                self.evalDisplacement(plot=doPlot)
            else:
                print 'no Diagnostic read'
        return



    def __del__( self ):

        if self.Status:

            del self.method
            del self.diag
            del self.timeRef
            del self.timeWindow   
              
        self.Status= False

    #def 


#stdOutlier = 2.5, ROutlier = 0.025
    def evalDisplacement( self , stdOutlier = 3.0, ROutlier = 0.03, plot=False):
        if self.Status:
 
            if len(self.displ) < 1:
                print 'nothing to evaluate'
                return

            MW = MAW.MAW() 
            MW.Load(self.shot)
            idx = np.where((MW.PSLtime>self.timeWindow[0]) & (MW.PSLtime<self.timeWindow[1]))[0]
            MW_Icoil = MW.Icoils[idx,0,0]
            MW_time = MW.time[idx]
            idxNew = np.where((MW.PSLMPtime>self.timeWindow[0]) & (MW.PSLMPtime<self.timeWindow[1]))[0]
            dPhi = np.unwrap(MW.PSLdPhase[idxNew]).mean()
            print 'dPhi: '
            print dPhi
            if dPhi < -np.pi: 
                dPhi += 2.*np.pi
            if dPhi > np.pi: 
                dPhi -= 2.*np.pi

            nNumber = MW.n.mean()
            

            sf=dd.shotfile('TOT',self.shot)
            #exception due to NTM
            if ((self.shot==34852) & (self.timeWindow[1]<6.5)):
                betaN = (sf('beta_N',tBegin=2.5, tEnd=3.1).data).mean()
            else:       
                betaN = (sf('beta_N',tBegin=self.timeWindow[0], tEnd=self.timeWindow[1]).data).mean()
            
            ##read equilibrium and get edge flux surface. Not exactly separatrix to avoid x-point
            print 'read equilibrium to get angles'
            RSep,zSep = [], []
            RMag,zMag = [], []
            
            for t in np.arange(self.timeWindow[0],self.timeWindow[1],0.1):
                out = kk.KK().kkrhorz(self.shot,t,[0.999,0.0],angle=np.linspace(0.,360.,200,endpoint=True),diag='EQI',ed=1)
                RSep.append(np.squeeze(out.r[:,0]))
                zSep.append(np.squeeze(out.z[:,0]))
                RMag.append(np.squeeze(out.r[:,1]))
                zMag.append(np.squeeze(out.z[:,1]))
            ### get values and define zero matrices and gradient    
            RSep = np.array(RSep)
            zSep = np.array(zSep)

            RMag = np.nanmean(RMag)
            zMag = np.nanmean(zMag)
            nuller = np.zeros_like(RSep)
            einser = np.ones_like(RSep)
            gradRSep=np.gradient(RSep)[1]
            gradzSep=np.gradient(zSep)[1]
            #get normalization
            InI = np.sqrt(gradRSep*gradRSep+gradzSep*gradzSep)
            InI[InI==0.0]=1.0
            #define vector for corss product
            vecTang=np.array([gradRSep/InI,gradzSep/InI,nuller])
            vecPhi=np.array([nuller,nuller,-einser])
            #calculate vector normal to the surface
            vecNorm = np.cross(vecPhi, vecTang,axis=0)[:2]
            RSepMean,zSepMean,vecNormMean = np.nanmean(RSep,axis=0),np.nanmean(zSep,axis=0), np.nanmean(vecNorm,axis=1)
            #angle 
            alphaEq = np.arctan(vecNormMean[1]/vecNormMean[0])
            #plt.quiver(RSepMean,zSepMean,vecNormMean[0],vecNormMean[1] , headwidth=4, headlength=6)
            ##equlibirum
           
 
            twoSigma = 0.682 
            self.displAmp=[]
            for i in np.arange(len(self.displ)):
                method = self.displ[i]['method']
                
                if (method == 'interpol') | (method == 'spline') | (method == '2line')| (method == 'shine'):
                    # go through each diagnostic
                    

                    label = self.displ[i]['lbl']
                    timeTrace = self.displ[i]['time']
                    xTrace = self.displ[i]['xTrace']
                    RTrace = self.displ[i]['RTrace']
                    zTrace = self.displ[i]['zTrace']
                        

                        #embed()
                        ## remove outlier
                    idxUse1 = np.abs(xTrace-np.nanmean(xTrace)) < np.std(xTrace)*stdOutlier
                        
                        #do fit
                    LSQshift = LSQFFT.LSQFFT()
                    LSQshift.initializeData(timeTrace[idxUse1], xTrace[idxUse1] ,time_ref=MW_time,signal_ref=MW_Icoil,order=4,poly=3,negSign=True)
                    useFreq = np.round(np.copy(LSQshift.freq0))
                    LSQshift.initializeData(timeTrace[idxUse1], xTrace[idxUse1], freq_in=useFreq,order=4,poly=3,negSign=True)

                    try:
                        print 'frequency: '+ LSQshift.freq0
                    except:
                        pass
                        #embed()

                    try:
                        LSQshift.leastSquareFit()
                    except:
                        print 'error'
                        embed()
                        ### remove datapoints which are more than 2 cms away from the fit
                    fitData = np.squeeze(LSQshift.getFunc( timeTrace ))
                    idxUse2 = np.abs(xTrace - fitData) < ROutlier
                    idxUse = idxUse1 & idxUse2 & np.isfinite(xTrace)
                    del LSQshift

                        ## use only these datapoints
                    label = self.displ[i]['lbl']
                    timeTrace = self.displ[i]['time'][idxUse]
                    xTrace = self.displ[i]['xTrace'][idxUse]
                    RTrace = self.displ[i]['RTrace'][idxUse]
                    zTrace = self.displ[i]['zTrace'][idxUse]
                    
                    idxSort = np.arange(RTrace.size)[np.isfinite(RTrace) & np.isfinite(zTrace)]
                    idxXMin,idxXMax = idxSort[np.nanargmin(xTrace[idxSort])],idxSort[np.nanargmax(xTrace[idxSort])]
                        #angle of diagnostic
                    if (self.displ[i]['lbl'] == 'CMZ') | (self.displ[i]['lbl'] == 'CPZ'):
                        alphaDiag = 0.0
                    else:
                        alphaDiag = np.arctan((zTrace[idxXMax]-zTrace[idxXMin])/(RTrace[idxXMax]-RTrace[idxXMin] ))
                        
                    thetaDiag = np.arctan2(zTrace.mean()-zMag,RTrace.mean()-RMag)

                        #closest point of mean Equilibrium to diagnostic
                    idxClose = np.argmin(np.sqrt((RSepMean-RTrace.mean())**2.+(zSepMean-zTrace.mean())**2.))
                        #calculate projection on the normalsurface, add +1.e-7 to avoid infinity
                    facPerp = 1./(np.abs(np.cos(alphaEq[idxClose]-alphaDiag))+1.e-7)
                    #embed()     
                    
                        #calculate fit using selected timepoints
                    LSQshift = LSQFFT.LSQFFT()
                    LSQshift.initializeData(timeTrace, xTrace, freq_in = useFreq, order=4, poly=3, negSign=True)
                    LSQshift.leastSquareFit()

                        ### calculate uncetainty
                    fitData = np.squeeze( LSQshift.getFunc( timeTrace ))
                    std = np.std(xTrace-fitData)
                        
                ### get uncertainty, ia 68.2% of all data points
                        
                    #make histogram
                    hist = np.histogram((xTrace-fitData)[:], bins=200,density =True)
                    xhist = hist[1][:-1]+np.diff(hist[1])
                    #integrate histogram to get surface
                    intHist = np.cumsum(hist[0])
                    #68.2% of the integrations defines the uncertainty
                    UncLo,UncHi = np.interp([(1.-twoSigma)/2.,1.-(1.-twoSigma)/2.],intHist/intHist.max(),xhist)
                        #embed()

                    self.displAmp.append( {'method':method,'lbl':label ,'amp':LSQshift.amp[0]*facPerp, 'ampStd':std*facPerp,'ampUncLo':UncLo*facPerp,'ampUncHi':UncHi*facPerp, 'ampUnc':LSQshift.uncAmp[0]*facPerp, 'phase':LSQshift.phase[0], 'phaseUnc':LSQshift.uncPhase[0],'real':LSQshift.real[0]*facPerp,'imag':LSQshift.imag[0]*facPerp, 'allAmp':LSQshift.amp*facPerp, 'allAmpUnc':LSQshift.uncAmp*facPerp ,'allPhase':LSQshift.phase,'allPhaseUnc':LSQshift.phase,'allReal':LSQshift.real*facPerp, 'allRealUnc':LSQshift.uncReal*facPerp,'allImag':LSQshift.imag*facPerp, 'allImagUnc':LSQshift.uncImag*facPerp,'facPerp':facPerp,'theta':thetaDiag,'dPhi':dPhi,'shot':self.shot, 'timeWindow': [self.timeWindow[0],self.timeWindow[1]],'betaN':betaN ,'freq':LSQshift.freq0, 'n':nNumber,'time':timeTrace,'xTrace':xTrace,'RTrace':RTrace,'zTrace': zTrace,'fitData':fitData } )                          
                   

                        
                    if plot:
                        plt.title('%s: '%method+'%s'%label+' factor %.3f'%facPerp)
                        plt.plot( self.displ[i]['time'][idxUse], self.displ[i]['xTrace'][idxUse])
                        plt.plot( self.displ[i]['time'][idxUse], fitData[:], linewidth = 3)
                        plt.xlabel('time [s]', fontsize=22)
                        plt.tight_layout()
                        plt.show()
                           

                    del LSQshift

   

            #embed()
            if plot:
                fmt = ['>','o','<','^']
                color = ['r','g','b','m','k']
                if len(self.displAmp)>0:
                    for i in np.arange(len(self.displAmp)):  
                        x=float( np.where(self.diag == self.displAmp[i]['lbl'])[0])
                        j=(np.where(self.method == self.displAmp[i]['method'])[0])
                        plt.errorbar(np.array([x+float(j)/10.]),[self.displAmp[i]['amp']],label=self.displAmp[i]['lbl']+' '+self.displAmp[i]['method'],fmt=fmt[j], yerr=[self.displAmp[i]['ampUnc']])
 
                    plt.xlim([-1,self.diag.size+3])
                    plt.legend()
                    plt.show()
                #embed()
#"""

##########################################################################################
##########################################################################################
####################      calculate Displacement           
##########################################################################################
##########################################################################################

    def calcDisplacement( self, timeIn=None, xIn=None, yIn=None, plot=True , debug=False):
        if self.Status:
            
            
            self.displ=[]
            self.interDispl=[]
            self.splineDispl=[]
            self.twolineDispl=[]
######## REF value 
            for method in self.method:
                for i in np.arange(len(self.refData)):
                        
#get Reference data
                    
                    xRef = self.refData[i]['x']
                    RRef = self.refData[i]['R']
                    zRef = self.refData[i]['z']
                    timeRef = self.refData[i]['time']
                    dataRef = self.refData[i]['data']
                    dataMean = np.nanmean(dataRef)

                    try:
                        intensRef = self.refData[i]['intens']
                    except:
                        pass

                #################################################################################
                ####################   INTERPOL for calculating displacement           
                #################################################################################

                    if method == 'interpol':
                        
                        
                        timeTrace = self.winData[i]['time']
                        
                        if (self.winData[i]['lbl'] == 'LIN') | (self.winData[i]['lbl'] == 'FRS') | (self.winData[i]['lbl'] == 'RPS_LFS') | (self.winData[i]['lbl'] == 'RPS_HFS') | (self.winData[i]['lbl'][:3] == 'RIC') | (self.winData[i]['lbl'] == 'HEB_ne') :
                            fac = np.linspace(self.fitNeParams[0],self.fitNeParams[1],30,endpoint=True)/dataMean
                            ## divide through mean to avoid numerical difficutlies
                            dataTrace = self.neTrace/dataMean
                            xIso = np.zeros((timeTrace.size,fac.size))
                        elif (self.winData[i]['lbl'] == 'RMD') | (self.winData[i]['lbl'] == 'CEC')| (self.winData[i]['lbl'] == 'HEB_Te') :
                            fac = np.linspace(self.fitTeParams[0],self.fitTeParams[1],30,endpoint=True)/dataMean
                            xIso = np.zeros((timeTrace.size,fac.size))
                            dataTrace = self.TeTrace/dataMean
                            lowThres = self.fitTeParams[0]/dataMean
                            upThres = self.fitTeParams[1]/dataMean
                        elif (self.winData[i]['lbl'] == 'CMZ') | (self.winData[i]['lbl'] == 'CPZ') :
                            fac = np.linspace(self.fitTiParams[0],self.fitTiParams[1],30,endpoint=True)/dataMean
                            xIso = np.zeros((timeTrace.size,fac.size))
                            dataTrace = self.TiTrace/dataMean
                            lowThres = self.fitTiParams[0]/dataMean
                            upThres = self.fitTiParams[1]/dataMean

                        #idxTrace = np.argmin(np.abs(fac - dataTrace))
                        xTrace = np.zeros((timeTrace.size))
                        dataWin = self.winData[i]['data']/dataMean
                        dataWinMean =  dataWin.mean(axis=0)
                        #embed()
                        if (self.winData[i]['lbl'] == 'RMD') | (self.winData[i]['lbl'] == 'CEC') :
                            xWin = self.winData[i]['x']
                            for j in np.arange(timeTrace.size):
                                idxSelect = np.where((dataWin[j]>lowThres) & (dataWin[j]<upThres))[0]
                                idxUse = idxSelect[np.argsort(xWin[j][idxSelect])]
                                xIn, dataIn = xWin[j][idxUse],dataWin[j][idxUse]
                                #search for the steepest gradient
                                dTrad_dr = -np.gradient(dataIn)/ np.gradient(xIn)
                                d2Trad_d2r = -np.gradient(dTrad_dr)/ np.gradient(xIn)
                                ### test this
                                idxShineMin = np.argmax( d2Trad_d2r )-1
                                shineMin = xIn[idxShineMin]
                                idxSort = np.arange(idxShineMin,dataIn.size)[np.argsort(dataIn[idxShineMin:])]
                                xIso[j,:] = np.interp(fac, dataIn[idxSort], xIn[idxSort])
                                xTrace[j] = np.interp(dataTrace, dataIn[idxSort], xIn[idxSort])

                           
                        else:

 #                           for j in np.arange(timeTrace.size):
 #                               idxFin = np.isfinite(dataWin[j]) & np.isfinite((self.winData[i]['x'])[j]) 
 #                               idxArg=idxFin[np.argsort( (dataWin[j])[idxFin] )]
 #                               xIso[j,:] = np.interp(fac, (dataWin[j])[idxArg], ((self.winData[i]['x'])[j])[idxArg])
 #                               xTrace[j] = np.interp(dataTrace, (dataWin[j])[idxArg], ((self.winData[i]['x'])[j])[idxArg])

                        ###interpolate for given denistz 
                            if dataWinMean[0:dataWinMean.size/2].mean()>dataWinMean[dataWinMean.size/2:-1].mean():
                                for j in np.arange(timeTrace.size):
                                    idxFin = np.isfinite(dataWin[j]) & np.isfinite((self.winData[i]['x'])[j]) 
                                    xIso[j,:] = np.interp(fac, (dataWin[j])[idxFin][::-1], ((self.winData[i]['x'])[j])[idxFin][::-1])
                                    xTrace[j] = np.interp(dataTrace, (dataWin[j])[idxFin][::-1], ((self.winData[i]['x'])[j])[idxFin][::-1])
                            else:
                                for j in np.arange(timeTrace.size):
                                    idxFin = np.isfinite(dataWin[j]) & np.isfinite((self.winData[i]['x'])[j])
                                    xIso[j,:] = np.interp(fac, dataWin[j][idxFin], (self.winData[i]['x'])[j][idxFin])
                                    xTrace[j] = np.interp(dataTrace, dataWin[j][idxFin], (self.winData[i]['x'])[j][idxFin])

                        try:
                            
                            self.displ.append( {'method': method,'time':timeTrace,'xTrace':xTrace,'RTrace':(self.winData[i]['Rfunc'])(xTrace),'zTrace':(self.winData[i]['zfunc'])(xTrace),'dataTrace':dataTrace*dataMean, 'lbl':self.winData[i]['lbl'], 'data':fac*dataMean,'xIso':xIso} )
                            
                        except:
                            print 'error writing displacement'
                            embed()

                    elif method == 'shine':
                        
                        if (self.winData[i]['lbl'] == 'RMD') | (self.winData[i]['lbl'] == 'CEC') :
                        ## divide through mean to avoid numerical difficutlies
                            timeTrace = self.winData[i]['time']
                            lowThres = self.fitTeParams[0]
                            upThres = self.fitTeParams[1]
 
                        #idxTrace = np.argmin(np.abs(fac - dataTrace))
                            xTrace = np.zeros((timeTrace.size))
                            xMin = np.zeros((timeTrace.size))
                            dataWin = self.winData[i]['data']
                            xWin = self.winData[i]['x']
                            
                            for j in np.arange(timeTrace.size):                       
                                idxSelect = np.where((dataWin[j]>lowThres) & (dataWin[j]<upThres))[0]
                                idxUse = idxSelect[np.argsort(xWin[j][idxSelect])]
                                xIn, dataIn = xWin[j][idxUse],dataWin[j][idxUse]
                                #search for the steepest gradient
                                dTrad_dr = -np.gradient(dataIn)/ np.gradient(xIn)
                                d2Trad_d2r = -np.gradient(dTrad_dr)/ np.gradient(xIn)
                                #skip first and lst point
                                idxShineMin = np.argmax( d2Trad_d2r[2:-1] )
                                ##data points around the minimu
                                xSmall,dataSmall = xIn[idxShineMin-2:idxShineMin+5],dataIn[idxShineMin-2:idxShineMin+5] 
                                zPoly = np.polyfit(xSmall,dataSmall,3)
                                p = np.poly1d(zPoly)
                                xfunc = np.linspace(xSmall.min(),xSmall.max(),200)
                                xTrace[j] = xfunc[np.argmin(p(xfunc))]
                                
                            try:
                                self.displ.append( {'method': method,'time':timeTrace,'xTrace':xTrace,'RTrace':(self.winData[i]['Rfunc'])(xTrace),'zTrace':(self.winData[i]['zfunc'])(xTrace), 'lbl':self.winData[i]['lbl'] } )
                            except:
                                print 'error writing displacement'
                                if debug:
                                    embed()
                        else:
                            print 'only ECE can use this'

        #################################################################################
        ####################   SPLINE for calculating displacement           
        #################################################################################
                        #embed() 
                    elif method == 'spline':

                        #to normalize data
                        dataMean = np.nanmean(dataRef)

                        #########################
                        ####################    preSelect data, remove low densities, shine through, etc...            
                        #########################
                   
                        #### remove SOL densities and temperature channels!
                        if (self.winData[i]['lbl'] == 'HEB_ne') | (self.winData[i]['lbl'] == 'LIN') | (self.winData[i]['lbl'] == 'FRS') | (self.winData[i]['lbl'] == 'RPS_LFS') | (self.winData[i]['lbl'] == 'RPS_HFS')| (self.winData[i]['lbl'][:3] == 'RIC'):
                            lowThres = self.fitNeParams[0]
                            upThres = self.fitNeParams[1]
                            dX = self.fitNeParams[2]
                            res = 0.0025
                            ### get R for ne
                            xMean = np.nanmedian(xRef[:,dataRef.mean(axis=0) > lowThres],axis=0)
                            #print 'hallo i bins'
                            #embed()
                            xMin = xMean.min()
                            if(self.winData[i]['lbl'] == 'HEB_ne') | (self.winData[i]['lbl'] == 'LIN'):
                                idxChan = (xRef >= xMin) & (xRef <= (xMin+dX) ) & (dataRef > lowThres) & (dataRef < upThres)  
                            else:
                                idxChan = np.where((np.nanmedian(xRef,axis=0) >= xMin) & (np.nanmedian(xRef,axis=0) <= (xMin+dX) )& (dataRef.mean(axis=0) > lowThres) & (dataRef.mean(axis=0) < upThres) )[0]
                            
                            #embed()
                                #sort argument
                            #idxSort = np.argsort(xRef[idxChan].ravel())
                        elif  (self.winData[i]['lbl'] == 'CPZ') | (self.winData[i]['lbl'] == 'CMZ'):
                            lowThres = self.fitTiParams[0]
                            upThres = self.fitTiParams[1]
                            ## set threshold for intensity
                            inteThres = self.fitTiParams[-2]
                            dX = self.fitTiParams[-1]
                            #chMedian = np.nanmedian(intensRef,axis=0)
                            #idxChan = chMedian>inteThres
                            idxChan = (intensRef >= inteThres ) & (dataRef > lowThres) & (dataRef < upThres)  
                            #embed()
                            
                            #Rmax = RRef[:,Rallidx].max()
                            #idxChan = np.where((np.nanmedian(RRef,axis=0) <= Rmax) & (np.nanmedian(RRef,axis=0) >= (Rmax-dR) ) & (dataRef.mean(axis=0) > lowThres) & (dataRef.mean(axis=0) < upThres)  )[0]

                        elif  (self.winData[i]['lbl'] == 'HEB_Te') | (self.winData[i]['lbl'] == 'CEC') | (self.winData[i]['lbl'] == 'RMD'):
                            lowThres = self.fitTeParams[0]
                            upThres = self.fitTeParams[1]
                            dX = self.fitTeParams[2]
                          
                       #########################
                       ####################       get spline of Reference!!                  
                       ####################       average the data LIN, CXRS must be averaged differently than ReF
                       #########################
                        # These diagnostics have the same R but different values
                        if (self.winData[i]['lbl'] == 'LIN') | (self.winData[i]['lbl'] == 'CPZ')| (self.winData[i]['lbl'] == 'CMZ'):
                            xIn,dataIn = xRef[idxChan],dataRef[idxChan]
                            idxSelect = np.where(xIn<(xIn.min()+dX))[0]
                            #xSpline,dataSpline = np.nanmean(xIn[idxSelect],axis=0),np.nanmean(dataIn[idxSelect]/dataMean,axis=0)xBins,dataBins,dummy = stats.binned_statistic(dataIn[idxSelect].ravel()/dataMean,xIn[idxSelect].ravel(),statistic='mean', bins=10)
                            dataNewBins,xNewBins,dummy = stats.binned_statistic(xIn[idxSelect].ravel(),dataIn[idxSelect].ravel()/dataMean,statistic='mean', bins=10)
                            
                            idxUse = np.isfinite(dataNewBins)
                            dataBinsIn = dataNewBins[idxUse]
                            xBinsIn = (np.diff(xNewBins).mean()+xNewBins[:-1])[idxUse]
                            
                            xSpline,dataSpline = xBinsIn,dataBinsIn

                            try:
                                splCoeffData = splrep(xSpline,dataSpline)
                            except:
                                print 'error in spline'
                                embed()

                            xLast = xBinsIn.min()
                        # Reflectometry is different because R has to be averaged
                        elif (self.winData[i]['lbl'] == 'RPS_LFS') | (self.winData[i]['lbl'] == 'RPS_HFS') | (self.winData[i]['lbl'] == 'FRS')| (self.winData[i]['lbl'][:3] == 'RIC') :
                        
                            
                            ### bin first the data and then 
                            bin_means, bin_edges, binnumber = stats.binned_statistic(dataRef[:,idxChan].ravel()/dataMean,xRef[:,idxChan].ravel(),statistic='mean', bins=20)
                            idxUse = np.isfinite(bin_means)
                            xBins = bin_means[idxUse]
                            dataBins = (np.diff(bin_edges).mean()+bin_edges[:-1])[idxUse]
                            idxSort = np.argsort(xBins)

                            dataNewBins,xNewBins,dummy = stats.binned_statistic( xBins, dataBins ,statistic='mean', bins=np.arange(RBins.min(),RBins.max(),res)  )
                             
                            idxUse = np.isfinite(dataNewBins)
                            dataBinsIn = dataNewBins[idxUse]
                            xBinsIn = (np.diff(xNewBins).mean()+xNewBins[:-1])[idxUse]
                            idxSort = np.argsort(xBinsIn)

                            xSpline,dataSpline = xBinsIn[idxSort],dataBinsIn[idxSort]
                            splCoeffData =splrep(RSpline,dataSpline)
                            xLast = xSpline.min()

                            #if plot:
                               # plt.plot(RSteps[idxSort],neSteps[idxSort],splev(x-shift,coeff))
                                
                                #splev(x-shift,splCoeffData)
                        elif (self.winData[i]['lbl'] == 'CEC')| (self.winData[i]['lbl'] == 'RMD'):                      
                            ## average data
                            xIn,DataIn =  np.nanmedian(xRef[:,:],axis=0)[::-1], np.nanmedian(dataRef[:,:],axis=0)[::-1]
                            idxSelect = np.where((DataIn>lowThres) & (DataIn<upThres))[0]
                            idxUse = idxSelect[np.argsort(xIn[idxSelect])]
                            xIn, DataIn = xIn[idxUse],DataIn[idxUse]
                            dTrad_dr = -np.gradient(DataIn)/ np.gradient(xIn)
                            d2Trad_d2r = - np.gradient(dTrad_dr)/ np.gradient(xIn)
                                #shinetrough position minus 1 to be sure
                            shinePos = np.argmax(d2Trad_d2r)+1
                            idxPed = np.where((np.arange(xIn.size)>=shinePos) & (xIn[:] <= (xIn[shinePos]+dX)) )[0]
                            xSpline,dataSpline = xIn[idxPed],DataIn[idxPed]/dataMean
                            try:
                                splCoeffData =splrep(xSpline,dataSpline)
                            except:
                                print 'stopped at spline function'
                                embed()
                            xLast = xIn[shinePos]
                            #print 'ECE spline'
                            #embed()
                        else:

                            idxSort = np.argsort(xRef[:,idxChan].ravel())
                            xSpline,dataSpline = xRef[:,idxChan].ravel()[idxSort],dataRef[:,idxChan].ravel()[idxSort]/dataMean
                            splCoeffData = splrep(xSpline,dataSpline)
                            xLast = xSpline.min()
                        if plot:
                            plt.plot(xSpline,dataSpline*dataMean,'ro',label=self.winData[i]['lbl'])
                            plt.plot(xSpline,splev(xSpline,splCoeffData)*dataMean,'b-',linewidth=2 )
                            plt.show()
                            #embed()
                            #plt.plot(RRef[:,idxChan].ravel(),dataRef[:,idxChan].ravel(),'ro')
                             #   plt.plot(RIn,splev(RIn,splCoeffData)*dataMean,'b-',linewidth=2 )
                              #  plt.show()

                       #########################
                       ####################    got throug every time points                   
                       #########################

                        ## get through every time points
                        xTmpLi = []
                        useChan = []
                        xShift=np.zeros_like(self.winData[i]['time'])
                        for t in np.arange(xShift.size):
                            if (self.winData[i]['lbl'] == 'LIN') | (self.winData[i]['lbl'] == 'FRS') | (self.winData[i]['lbl'] == 'RPS_LFS')| (self.winData[i]['lbl'][:3] == 'RIC')| (self.winData[i]['lbl'] == 'RPS_HFS')| (self.winData[i]['lbl'] == 'HEB_ne')| (self.winData[i]['lbl'] == 'HEB_Te'):

                                xTmp = self.winData[i]['x'][t,(self.winData[i]['data'][t]>lowThres)].min()
                                useChan =  np.where((self.winData[i]['x'][t,:] >= xTmp) & (self.winData[i]['x'][t,:] <= (xTmp+dX) ))[0]
                                useChan = useChan[np.argsort(self.winData[i]['x'][t,useChan])]
                            
                            # for CXRS the intensity sets the Threshold
                            elif (self.winData[i]['lbl'] == 'CMZ') | (self.winData[i]['lbl'] == 'CPZ') :
                                xTmp = self.winData[i]['x'][t,(self.winData[i]['intens'][t]>inteThres)].max()
                                useChan =  np.where((self.winData[i]['x'][t,:] >= xTmp) & (self.winData[i]['x'][t,:] <= (xTmp+dX) ) & np.isfinite(self.winData[i]['data'][t]) & (self.winData[i]['intens'][t]>inteThres) )[0]
                                useChan = useChan[np.argsort(self.winData[i]['x'][t,useChan])]               
  
                            ### use the selected data to evaluate shift between the splines
                            elif (self.winData[i]['lbl'] == 'CEC')| (self.winData[i]['lbl'] == 'RMD'):
                                data = self.winData[i]['data'][t]
                                xx = self.winData[i]['x'][t]
                                idxSelect = np.where((data>lowThres) & (data<upThres))[0]
                                idxUse = idxSelect[np.argsort(xx[idxSelect])]
                                dTrad_dr = -np.gradient(data[idxUse])/ np.gradient(xx[idxUse])
                                d2Trad_d2r = -np.gradient(dTrad_dr)/ np.gradient(xx[idxUse])
                                ## shinePos
                                xTmp= xx[idxUse[np.argmax(d2Trad_d2r)-1]]
                                useChan =  np.where((self.winData[i]['x'][t,:] >= xTmp) & (self.winData[i]['x'][t,:] <= (xTmp+dX) ) & np.isfinite(self.winData[i]['data'][t]) & (self.winData[i]['data'][t]>lowThres)&(self.winData[i]['data'][t]<upThres) )[0]
                                useChan = useChan[np.argsort(self.winData[i]['x'][t,useChan])]
                            
                            try:
                                xShift[t],xTraceCov = leastsq(shiftFuncResi,[0.0] ,args=(self.winData[i]['x'][t,useChan],self.winData[i]['data'][t,useChan]/dataMean,np.ones_like( self.winData[i]['x'][t,useChan] ),splCoeffData),full_output=1)[0:2]
                            except:
                                print 'error in fit'
                                embed()

                        xTrace = xShift+xLast

                        self.displ.append( {'method':method, 'time':self.winData[i]['time'],'xTrace':xTrace, 'lbl':self.winData[i]['lbl'],'RTrace':(self.winData[i]['Rfunc'])(xTrace),'zTrace':(self.winData[i]['zfunc'])(xTrace)}  )
                        
                       # embed()

        #################################################################################
        ####################  END  SPLINE for calculating displacement           
        #################################################################################

                        
                    elif method == '2line':
                                               
                        if (self.winData[i]['lbl'] == 'LIN') | (self.winData[i]['lbl'] == 'FRS') | (self.winData[i]['lbl'] == 'RPS_LFS') | (self.winData[i]['lbl'] == 'RPS_HFS')| (self.winData[i]['lbl'][:3] == 'RIC'):
                            lowThres = self.fitNeParams[0]
                            upThres = self.fitNeParams[1]
                            dX = self.fitNeParams[2]
                            Trace = self.neTrace
                        elif  (self.winData[i]['lbl'] == 'CPZ') | (self.winData[i]['lbl'] == 'CMZ'):
                            ## set threshold for intensity
                            inteThres = self.fitTiParams[-2]
                            dX = self.fitTiParams[-1]
                            Trace = self.TiTrace
                        elif  (self.winData[i]['lbl'] == 'CEC') | (self.winData[i]['lbl'] == 'RMD'):
                            lowThres = self.fitTeParams[0]
                            upThres = self.fitTeParams[1]
                            dX = self.fitTeParams[2]
                            Trace = self.TeTrace 

                        dataMean = dataRef.mean() 
                        xTrace=np.zeros_like(self.winData[i]['time'])
                        xIn = np.arange(xRef.min()-0.1,xRef.max()+0.1,0.01)
                        
                        for t in np.arange(self.winData[i]['time'].size):
                            try:
                                if (self.winData[i]['lbl'] == 'LIN') | (self.winData[i]['lbl'] == 'FRS') | (self.winData[i]['lbl'] == 'RPS_LFS')| (self.winData[i]['lbl'][:3] == 'RIC')|(self.winData[i]['lbl'] == 'RPS_HFS'):
                                    xTmp = self.winData[i]['x'][t,(self.winData[i]['data'][t] > lowThres) & (self.winData[i]['data'][t] < upThres)].min()
                                    useChan =  np.where((self.winData[i]['x'][t,:] >= xTmp) & (self.winData[i]['x'][t,:] <= (xTmp+dX))&(self.winData[i]['data'][t] > lowThres) & (self.winData[i]['data'][t] < upThres))[0]
                                    useChan = useChan[np.argsort(self.winData[i]['x'][t,useChan])]
                                elif (self.winData[i]['lbl'] == 'CMZ') | (self.winData[i]['lbl'] == 'CPZ') :
                                    xTmp = self.winData[i]['x'][t,(self.winData[i]['intens'][t]>inteThres)].max()
                                    useChan =  np.where((self.winData[i]['x'][t,:] >= xTmp) & (self.winData[i]['x'][t,:] <= (xTmp-dX) ) & np.isfinite(self.winData[i]['data'][t]) & (self.winData[i]['intens'][t]>inteThres) )[0]
                                    useChan = useChan[np.argsort(self.winData[i]['x'][t,useChan])]  
                                elif (self.winData[i]['lbl'] == 'CEC')| (self.winData[i]['lbl'] == 'RMD'):
                                    data = self.winData[i]['data'][t]
                                    XX = self.winData[i]['x'][t]
                                    idxSelect = np.where((data>lowThres) & (data<upThres))[0]
                                    idxUse = idxSelect[np.argsort(XX[idxSelect])]
                                    dTrad_dr = -np.gradient(data[idxUse])/ np.gradient(XX[idxUse])
                                    d2Trad_d2r = -np.gradient(dTrad_dr)/ np.gradient(XX[idxUse])
                                ## shinePos
                                    xTmp= XX[idxUse[np.argmax(d2Trad_d2r)]]
                                    useChan =  np.where((self.winData[i]['x'][t,:] >= xTmp) & (self.winData[i]['x'][t,:] <= (xTmp+dX) ) & np.isfinite(self.winData[i]['data'][t]) & (self.winData[i]['data'][t]>lowThres)&(self.winData[i]['data'][t]<upThres) )[0]
                                    useChan = useChan[np.argsort(self.winData[i]['x'][t,useChan])]
                                else:
                                    print hallo                                             

                           
                                twoLine= pedestalFits.Fit(self.winData[i]['x'][t,useChan],self.winData[i]['data'][t,useChan]/dataMean,twoLine=True)
                                twoLine.doFit()
                                FitIn = twoLine.getFunc(xIn)
                                idxSort = np.argsort(FitIn)
                                xTrace[t] = np.interp(Trace/dataMean,FitIn[idxSort],xIn[idxSort] )
                                del twoLine
                            except:
                                print 'twoline failed'
                                #embed()
                                continue

                        self.displ.append( {'method':method, 'time':self.winData[i]['time'],'xTrace':xTrace, 'lbl':self.winData[i]['lbl'],'RTrace':(self.winData[i]['Rfunc'])(xTrace),'zTrace':(self.winData[i]['zfunc'])(xTrace)}  )
                        
                            

##########################################################################################
##########################################################################################
####################       reading data                 
##########################################################################################
##########################################################################################

    def loadData( self,  readTime = None, Rlimit = [1.95,2.2], debug=True ):

        if self.Status:
            if ( np.all(readTime) == None ):
                readTime = np.array([np.array([self.timeRef[0],self.timeWindow[0]]).min(),np.array([self.timeRef[1],self.timeWindow[1]]).max()])
             

            self.refData = []
            self.winData = []
            timeTmp,RTmp,zTmp,dataTmp,intensTmp = 0.0,0.0,0.0,0.0,0.0
            
            for diag in self.diag:
                didRead = False
                if (diag == 'RMD') | (diag == 'CEC'):
                    try:
                        RMD=ECE.ECE()
                        RMD.Load(self.shot,tBegin=readTime[0], tEnd=readTime[1],Experiment='AUGD',Diagnostic=diag)
                        timeTmp,RTmp,zTmp,phiTmp,dataTmp = RMD.time.copy(),RMD.Rall.copy(),RMD.zall.copy(),RMD.phi,RMD.Te.copy()
                        ##assuming 

                        ##assuming 
                        didRead = True
                        del RMD
                    except:
                        print 'error in reading ECE'
                        didRead = False
                        if debug:
                            embed()

                elif (diag == 'LIN'):
                    try:
                        LIZ=LIN.LIN()
                        if LIZ.Load(self.shot, tBegin=readTime[0], tEnd=readTime[1], Experiment='AUGE') == False:
                            LIZ.Load(self.shot, tBegin=readTime[0], tEnd=readTime[1], Experiment='AUGD')
                        timeTmp,RTmp,zTmp,xTmp,phiTmp,dataTmp = LIZ.time.copy(), LIZ.Rall.copy(), LIZ.zall.copy(), LIZ.xall.copy()/100., LIZ.phi, LIZ.ne.copy()
                        Rfunc = interp1d(LIZ.x/100.,LIZ.R,bounds_error=False)
                        zfunc = interp1d(LIZ.x/100.,LIZ.z,bounds_error=False)
                        didRead = True
                        del LIZ
                    except:
                        print 'error in reading LIN'
                        didRead = False
                        if debug:
                            embed()



                elif (diag[:3] == 'HEB'):
                    try:
                        HelB=HEB.HEB()
                        if diag[-2:]=='Te':
                            HelB.Load(self.shot, tBegin=readTime[0], tEnd=readTime[1], Experiment='MGRIEN',select=[0,-4],Edition=0) 
                            timeTmp,RTmp,zTmp,phiTmp,dataTmp = HelB.time.copy(), HelB.Rall.copy()[:,::-1], HelB.zall.copy()[:,::-1], HelB.phi.mean(), HelB.Te.copy()[:,::-1]
                        else:
                            HelB.Load(self.shot, tBegin=readTime[0], tEnd=readTime[1], Experiment='MGRIEN',select=[2,-1],Edition=0) 
                            timeTmp,RTmp,zTmp,phiTmp,dataTmp = HelB.time.copy(), HelB.Rall.copy()[:,::-1], HelB.zall.copy()[:,::-1], HelB.phi.mean(), HelB.ne.copy()[:,::-1]
                        idxArr= (RTmp!=0.0) & (zTmp!=0.0)
                        poly = np.polyfit(RTmp[idxArr], zTmp[idxArr], 1)
                        Rfit = np.array([RTmp[idxArr].max(),RTmp[idxArr].min()])
                        zfit = poly[0]*Rfit+poly[1]
                        xfit = np.sqrt((Rfit-RTmp[idxArr].max())**2.+(zfit-zTmp[idxArr].max())**2.)
                        xTmp = np.interp(RTmp,Rfit[::-1],xfit[::-1])
                        xTmp[RTmp==0.0] = 0.0
                        RTmp[RTmp==0.0] = float('NaN')
                        zTmp[zTmp==0.0] = float('NaN')
                        dataTmp[dataTmp==0.0] = float('NaN')
                        idxSort = np.argsort(xTmp[idxArr])
                        Rfunc =interp1d(xTmp[idxArr][idxSort], RTmp[idxArr][idxSort],bounds_error=False)
                        zfunc =interp1d(xTmp[idxArr][idxSort], zTmp[idxArr][idxSort],bounds_error=False)
                        didRead = True
                        del HelB
                    except:
                        print 'error in reading HEB'
                        didRead = False
                        if debug:
                            embed()



                elif (diag == 'FRS'):
                    try:
                        FRS=REF.REF()
                        FRS.Load(self.shot, tBegin=readTime[0], tEnd=readTime[1], Experiment='AUGD', Diagnostic='FRS')
                        timeTmp,RTmp,zTmp,dataTmp = FRS.time.copy(), FRS.R.copy(), FRS.z.copy(), FRS.ne.copy()
                        didRead = True
                        del FRS
                    except:
                        print 'error in reading FRS'
                        didRead = False
                        if debug:
                            embed()

                elif (diag == 'RPS_LFS'):
                    try:
                        RPSLFS=REF.REF()
                        RPSLFS.Load(self.shot,tBegin=readTime[0], tEnd=readTime[1],Experiment='AUGD',Diagnostic='RPS',LFSorHFS='LFS')
                        timeTmp,RTmp,zTmp,dataTmp = RPSLFS.time.copy(), RPSLFS.R.copy(), RPSLFS.z.copy(), RPSLFS.ne.copy()
                        didRead = True
                        del RPSLFS
                    except:
                        print 'error in reading RPS LFS'
                        didRead = False
                        if debug:
                            embed()

                elif (diag == 'RPS_HFS'):
                    try:
                        RPSHFS=REF.REF()
                        RPSHFS.Load(self.shot,tBegin=readTime[0], tEnd=readTime[1],Experiment='AUGD',Diagnostic='RPS',LFSorHFS='HFS')
                        timeTmp,RTmp,zTmp,dataTmp = RPSHFS.time.copy(), RPSHFS.R.copy(), RPSHFS.z.copy(), RPSHFS.ne.copy()
                        didRead = True
                        del RPSHFS
                    except:
                        print 'error in reading RPS HFS'
                        didRead = False
                        if debug:
                            embed()

                elif (diag[:3] == 'RIC'):
                    try:
                        RIC=REF.REF()
                        if  RIC.Load(self.shot,tBegin=readTime[0], tEnd=readTime[1],Experiment='AUGD',Diagnostic='RIC',RICAntenna=int(diag[3])) == False:
                            RIC.Load(self.shot,tBegin=readTime[0], tEnd=readTime[1],Experiment='RICG',Diagnostic='RIC',RICAntenna=int(diag[3]))
                        timeTmp,RTmp,zTmp,dataTmp = RIC.time.copy(), RIC.R.copy(), RIC.z.copy(), RIC.ne.copy()
                        didRead = True
                        del RIC
                    except:
                        print 'error in reading %s'%diag
                        didRead = False
                        if debug:
                            embed()

                elif (diag == 'CMZ'):
                    try:
                        CMZ=CXRS.CXRS()
                        if CMZ.Load(self.shot, tBegin=readTime[0], tEnd=readTime[1], Experiment='mcavedon', Diagnostic='CMZ') == False:
                            CMZ.Load(self.shot, tBegin=readTime[0], tEnd=readTime[1], Experiment='AUGD', Diagnostic='CMZ')
                        
                        timeTmp,RTmp,zTmp = CMZ.time.copy(), CMZ.Rall.copy(), CMZ.zall.copy()
                        
                        intensTmp =  CMZ.inte.copy()
                        if self.CXRSuse == 'Ti':
                            dataTmp =  CMZ.Ti.copy()
                        elif self.CXRSuse == 'inte':
                            dataTmp =  intensTmp
                        elif self.CXRSuse == 'vrot':
                            dataTmp =  CMZ.vrot.copy()

                        if self.shot==34852:
                            idx =np.argmin(np.abs(RTmp.mean(axis=0)-2.132)) 
                            RTmp[:,idx] = 0.0
                            dataTmp[:,idx] = 0.0
                        
                        if self.shot==34634:
                            idx1 =np.argmin(np.abs(RTmp.mean(axis=0)-2.132)) 
                            RTmp[:,idx1] = 0.0
                            dataTmp[:,idx1] = 0.0
  
                            idx2 =np.argmin(np.abs(RTmp.mean(axis=0)-2.1246)) 
                            RTmp[:,idx2] = 0.0
                            dataTmp[:,idx2] = 0.0
                      
                        ## always the same array
                        idxArr= (RTmp!=0.0) & (zTmp!=0.0)
                        poly = np.polyfit(RTmp[idxArr], zTmp[idxArr], 1)
                        Rfit = np.array([RTmp[idxArr].max(),RTmp[idxArr].min()])
                        zfit = poly[0]*Rfit+poly[1]
                        xfit = np.array([0.0,np.sqrt(np.diff(Rfit)**2.+np.diff(zfit)**2.)])
                        xTmp = np.interp(RTmp,Rfit[::-1],xfit[::-1])
                        # ein Kanal weg in CMZ


                        xTmp[RTmp==0.0] = 0.0
                        RTmp[RTmp==0.0] = float('NaN')
                        zTmp[zTmp==0.0] = float('NaN')
                        dataTmp[dataTmp==0.0] = float('NaN')
                        idxSort = np.argsort(xTmp[idxArr])
                        Rfunc =interp1d(xTmp[idxArr][idxSort], RTmp[idxArr][idxSort],bounds_error=False)
                        zfunc =interp1d(xTmp[idxArr][idxSort], zTmp[idxArr][idxSort],bounds_error=False)
                        didRead = True
                        del CMZ
                    except:
                        print 'error in reading CMZ'
                        didRead = False
                        if debug:
                            embed()

                elif (diag == 'CPZ'):
                    try:
                        CPZ=CXRS.CXRS()
                        if CPZ.Load(self.shot, tBegin=readTime[0], tEnd=readTime[1], Experiment='mcavedon', Diagnostic='CPZ') == False: 
                            CPZ.Load(self.shot, tBegin=readTime[0], tEnd=readTime[1], Experiment='AUGD', Diagnostic='CPZ')
                        timeTmp,RTmp,zTmp = CPZ.time.copy(), CPZ.Rall.copy(), CPZ.zall.copy()
                    
                        intensTmp =  CPZ.inte.copy()
                        if self.CXRSuse == 'Ti':
                            dataTmp =  CPZ.Ti.copy()
                        elif self.CXRSuse == 'inte':
                            if ((self.shot == 34424) | (self.shot == 34427 )| (self.shot == 34852 )| (self.shot == 34622 )| (self.shot == 34634 )| (self.shot == 34673 )| (self.shot == 34672 )):
                                dataTmp =  intensTmp.copy()*5. 
                            else:
                                dataTmp =  intensTmp.copy()
                        elif self.CXRSuse == 'vrot':
                            dataTmp =  CPZ.vrot.copy()

                        if self.shot==34634:
                            idx1 =np.argmin(np.abs(RTmp.mean(axis=0)-2.131)) 
                            RTmp[:,idx1] = 0.0
                            dataTmp[:,idx1] = 0.0
  
                            idx2 =np.argmin(np.abs(RTmp.mean(axis=0)-2.1246)) 
                            RTmp[:,idx2] = 0.0
                            dataTmp[:,idx2] = 0.0
                            
                            idx3 =np.argmin(np.abs(RTmp.mean(axis=0)-2.128)) 
                            RTmp[:,idx3] = 0.0
                            dataTmp[:,idx3] = 0.0



                        idxArr= (RTmp!=0.0) & (zTmp!=0.0)
                        poly = np.polyfit(RTmp[idxArr], zTmp[idxArr], 1)
                        Rfit = np.array([RTmp[idxArr].max(),RTmp[idxArr].min()])
                        zfit = poly[0]*Rfit+poly[1]
                        xfit = np.array([0.0,np.sqrt(np.diff(Rfit)**2.+np.diff(zfit)**2.)])
                        xTmp = np.interp(RTmp,Rfit[::-1],xfit)
                        xTmp[RTmp==0.0] = 0.0
                        RTmp[RTmp==0.0] = float('NaN')
                        zTmp[zTmp==0.0] = float('NaN')
                        dataTmp[dataTmp==0.0] = float('NaN')
                        idxSort = np.argsort(xTmp[idxArr])
                        Rfunc =interp1d(xTmp[idxArr][idxSort], RTmp[idxArr][idxSort],bounds_error=False)
                        zfunc =interp1d(xTmp[idxArr][idxSort], zTmp[idxArr][idxSort],bounds_error=False)
                        didRead = True
                        del CPZ
                    except:
                        print 'error in reading CPZ'
                        didRead = False
                        if debug:
                            embed()
                            
                else:
                    print 'no appropriate Diagnostic found for '+diag
                    

                if didRead & (np.size(timeTmp)>1):
                    
                    try:
                    ###get LOS of diagnostics(diag == 'RPS_HFS') 
                        if (diag == 'RPS_HFS') | (diag == 'RPS_LFS') | (diag[:3] == 'RIC') | (diag == 'RMD') | (diag == 'CEC') | (diag == 'FRS'):
                            if (diag == 'FRS'):
                                idxArr= (RTmp!=0.0)
                            else:
                                idxArr= (RTmp!=0.0) & (zTmp!=0.0)
                        ### for the HFS the LOS start from minimum R
                            if (diag == 'RPS_HFS'):
                            ## search for the first R position
                                argMax = RTmp[idxArr].argmin()
                            else:
                                argMax = RTmp[idxArr].argmax()

                            Rmax = RTmp[idxArr][argMax]
                            zmax = zTmp[idxArr][argMax]
                            xTmp = np.sqrt((RTmp-Rmax)**2.+(zTmp-zmax)**2.)
                            xTmp[idxArr==False] = 0.0
                            thetaTmp = np.arctan2((zTmp-zmax),(RTmp-Rmax))
                            idxSort = np.argsort(xTmp[idxArr])
                            Rfunc =interp1d(xTmp[idxArr][idxSort], RTmp[idxArr][idxSort],bounds_error=False)
                            zfunc =interp1d(xTmp[idxArr][idxSort], zTmp[idxArr][idxSort],bounds_error=False)
                            #poly = np.polyfit(RTmp[idxArr][idxSort],zTmp[idxArr][idxSort])
                    except:
                        print 'error in getting x axis for  ' + diag 
                        if debug:
                            embed()  

                    try:
                        if self.useELMs:
                            idxELM = ElmSync.ElmExtract(timeTmp,self.shot,plot=False,preFac = self.ELMperiod[0], postFac = self.ELMperiod[1], Experiment='AUGD')
                        else:
                            idxELM = np.arange(timeTmp.size)
                        
                        idxRef = idxELM[np.where( (timeTmp[idxELM] > self.timeRef[0]) & (timeTmp[idxELM] < self.timeRef[1]))[0]]
                        idxWin = idxELM[np.where( (timeTmp[idxELM] > self.timeWindow[0]) & (timeTmp[idxELM] < self.timeWindow[1]))[0]]

                        if len(idxWin) > 1 :

                        
                            if (diag == 'RPS_HFS'):
                                useChan =  (dataTmp.mean(axis=0) > 0.0) 
                            else:
                                useChan = ( np.nanmean(RTmp[idxRef], axis=0) > Rlimit[0]) & (np.nanmean(dataTmp,axis=0) > 0.0) & ( np.nanmean(RTmp[idxRef],axis=0) < Rlimit[1]) 
                            
                            if ( (diag == 'CPZ') | (diag == 'CMZ') | (diag == 'CNZ') ):

                                self.refData.append( {'time':timeTmp[idxRef],'R':RTmp[idxRef][:,useChan],'z':zTmp[idxRef][:,useChan],'x':xTmp[idxRef][:,useChan], 'data':dataTmp[idxRef][:,useChan],'intens':intensTmp[idxRef][:,useChan],'lbl':diag} )
                                self.winData.append( {'time':timeTmp[idxWin],'R':RTmp[idxWin][:,useChan],'z':zTmp[idxWin][:,useChan],'x':xTmp[idxWin][:,useChan], 'data':dataTmp[idxWin][:,useChan],'intens':intensTmp[idxWin][:,useChan],'Rfunc':Rfunc,'zfunc':zfunc ,'lbl':diag} )

                            else:
                            
                                self.refData.append( {'time':timeTmp[idxRef],'R':RTmp[idxRef][:,useChan],'z':zTmp[idxRef][:,useChan],'x':xTmp[idxRef][:,useChan], 'data':dataTmp[idxRef][:,useChan],'lbl':diag} )
                                self.winData.append( {'time':timeTmp[idxWin],'R':RTmp[idxWin][:,useChan],'z':zTmp[idxWin][:,useChan],'x':xTmp[idxWin][:,useChan], 'data':dataTmp[idxWin][:,useChan],'Rfunc':Rfunc,'zfunc':zfunc,'lbl':diag} )                            
                        
                            del timeTmp
                            del RTmp
                            del zTmp
                            del xTmp
                            del dataTmp
                        else:
                            print 'no Data in analysis window'

                    except:
                        print 'error in reading  ' + diag 
                        if debug:
                            embed()              

            print 'everything read'
            #embed()
            self.diagsRead = len(self.winData)
            #embed()
  





