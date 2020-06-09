import numpy as np
from IPython import embed
import matplotlib.pylab as plt
from scipy.interpolate import interp1d

import LSQFFT
import EQU as eq
import MAWnew as MAW
import ElmSync
import kk_mwillens as kk

class sepPer:

    def __init__( self ,  Shotnumber = None, Experiment1='AUGD', Experiment2='AUGD', Diagnostic='IDE', Edition1 = 0L, Edition2 = 0L,tBegin=0.0,tEnd=14.0 , ELMperiod=None, freqIn=None ):
        ### 1 fuer alter Kranz bei 157.5 und 2 fuer Kranz bei 112.5, Edition1 fuer alter Kranz bei 157.5, n...toroidal mode number, negative for negative sense of rotation
        self.IDEread=False

        self.Shotnumber = Shotnumber
        self.Experiment1 =  Experiment1
        self.Experiment2 =  Experiment2
        self.Diagnostic = Diagnostic
        self.Edition1 = Edition1 
        self.Edition2 = Edition2 

        self.IDE1 = eq.EQU()
        self.IDE1.Load( Shotnumber, Experiment1, Diagnostic, Edition1,tBegin=tBegin,tEnd=tEnd,Diagnostic2=Diagnostic )
        
        self.IDE2 = eq.EQU()
        self.IDE2.Load( Shotnumber, Experiment2, Diagnostic, Edition2,tBegin=tBegin,tEnd=tEnd,Diagnostic2=Diagnostic )
        
        MW=MAW.MAW()
        if Shotnumber==32406:
            MW.Load(32415,tBegin=tBegin,tEnd=tEnd )
        else:
            MW.Load(Shotnumber,tBegin=tBegin,tEnd=tEnd )

        if (MW.maxFreqCoils[0].mean() > 0.1):
            if((MW.rotVelUp > 0.1) & (MW.rotVelUp > 0.1)):
                self.rotDir=1.
            elif((MW.rotVelUp < -0.1) & (MW.rotVelUp < -0.1)):
                self.rotDir=-1.
            else:
                print 'something went wrong'
                return
        elif (MW.maxFreqCoils[1].mean() > 0.1):
            if((MW.rotVelLo > 0.1) & (MW.rotVelLo > 0.1)):
                self.rotDir=1.
            elif((MW.rotVelLo < -0.1) & (MW.rotVelLo < -0.1)):
                self.rotDir=-1.
            else:
                print 'something went wrong'
                return
        else:
            print 'no rotation'
        
        #self.rotDir=self.rotDir*-1.
        #embed()
        if np.all(ELMperiod) != None:
            print 'filter ELMs'
            self.idxELM1 = ElmSync.ElmExtract(self.IDE1.time,Shotnumber,plot=False,preFac = ELMperiod[0], postFac = ELMperiod[1], Experiment='AUGD') 
           # idxELM2 = ElmSync.ElmExtract(self.IDE2.time,Shotnumber,plot=False,preFac = ELMperiod[0], postFac = ELMperiod[1], Experiment='AUGD')  
        else:
            self.idxELM1 = np.arange(0,self.IDE1.time.size)
            #idxELM2 = np.arange(0,self.IDE2.time.size)

        if np.all(freqIn) != None:
            self.f = freqIn
        else:
            self.f = np.round(np.mean(MW.maxFreqCoils[MW.maxFreqCoils!=0])) #frequency

        print 'n Number is: ' ,MW.meanN
 
        #embed()
        self.n = np.round(MW.meanN)
        #lobal omega, ntor 
        self.om = 2. * np.pi * self.f
        #tor = n #toroidal mode number
        self.IDEread = True
        
    def analyseDiag( self, RIn, zIn, xIn = None,  timeDiag=None,signalDiag=None,plot=False, extendOut=False, specialPoints=None, timeTraceControl = False ,polyIn=3,orderSet=3 ):
        
        if self.IDEread:
                #use Timebase of IDE
            timeIn = self.IDE1.time[self.idxELM1] 
            self.RIn,self.zIn = RIn, zIn

            #specialPoints = np.array(specialPoints)
            if (specialPoints == None):

                if (np.size(xIn)==np.size(RIn)):
                    idxSort=np.argsort(xIn)
                    xIn = xIn[idxSort]
                    RIn = RIn[idxSort]
                    zIn = zIn[idxSort]
                    self.RIn,self.zIn = RIn, zIn
                    self.xIn = xIn
                else:
                ##calculate the most outer point
                    dist = np.sqrt((RIn-self.IDE1.Rmag.mean(axis=0))**2.+(zIn-self.IDE1.zmag.mean(axis=0))**2.)
                    idxOuter = np.argmax(dist)
                    xIn = np.sqrt((RIn-RIn[idxOuter])**2.+(zIn-zIn[idxOuter])**2.)
                    self.xIn = xIn

                rhopOut1 = self.IDE1.getRhop(timeIn,RIn,zIn) 
                rhopOut2 = self.IDE2.getRhop(timeIn,RIn,zIn)

                xSep1 = np.zeros_like(timeIn)
                xSep2 = np.zeros_like(timeIn)
                for i in np.arange(xSep1.size):
                    idxSort1=np.argsort(rhopOut1[i])
                    xSep1[i] = np.interp([1.0],rhopOut1[i,idxSort1],xIn[idxSort1])

                    idxSort2= np.argsort(rhopOut2[i])
                    xSep2[i] = np.interp([1.0],rhopOut2[i,idxSort1],xIn[idxSort1])
            
            #array 1
                    
                xSep1_mean = xSep1.mean()
                xSep2_mean = xSep2.mean()

                xSep1 = xSep1 - xSep1.mean()
                xSep2 = xSep2 - xSep2.mean()
            
            elif np.all(specialPoints == 'Raus'):
                Raus1= self.IDE1.Raus[self.idxELM1]
                Raus2= self.IDE2.Raus[self.idxELM1]
                xSep1 = Raus1 - Raus1.mean()
                xSep2 = Raus2 - Raus2.mean()
            elif np.all(specialPoints =='zaus') | np.all(specialPoints =='zmag' ):
                zmag1= self.IDE1.zmag[self.idxELM1]
                zmag2= self.IDE2.zmag[self.idxELM1]
                xSep1 = zmag1 - zmag1.mean()
                xSep2 = zmag2 - zmag2.mean()
            elif np.all(specialPoints =='Rmag'):
                Rmag1= self.IDE1.Rmag[self.idxELM1]
                Rmag2= self.IDE2.Rmag[self.idxELM1]
                xSep1 = Rmag1 - Rmag1.mean()
                xSep2 = Rmag2 - Rmag2.mean()
            elif np.all(specialPoints =='Zoben'):
                zoben1= self.IDE1.zoben[self.idxELM1]
                zoben2= self.IDE2.zoben[self.idxELM1]
                xSep1 = zoben1 - zoben1.mean()
                xSep2 = zoben2 - zoben2.mean()
            elif np.all(specialPoints =='Rzoben'):
                Rzoben1= self.IDE1.Rzoben[self.idxELM1]
                Rzoben2= self.IDE2.Rzoben[self.idxELM1]
                xSep1 = Rzoben1 - Rzoben1.mean()
                xSep2 = Rzoben2 - Rzoben2.mean()
             # calculate Eq 
           # 2*xi_eq*np.sin(-ntor*np.pi/8)*np.cos(omega*x - phi_eq - ntor*np.pi/8) + c_3
           # =>2*xi_eq*np.sin(ntor*np.pi/8)*np.sin(omega*x - phi_eq + ntor*np.pi/8-np.pi/4.) + c_
            # xiStar = xi_eq*2*np.sin(r ntor*np.pi/8)
            
            RSep1_mean = np.interp(xSep1_mean,xIn,RIn)
            RSep2_mean = np.interp(xSep2_mean,xIn,RIn)
            zSep1_mean = np.interp(xSep1_mean,xIn,zIn)
            zSep2_mean = np.interp(xSep2_mean,xIn,zIn)

            self.RSep_mean = np.mean([RSep1_mean,RSep2_mean])
            self.zSep_mean = np.mean([zSep1_mean,zSep2_mean])
            self.xSep_mean = np.mean([xSep1.mean(),xSep2.mean()])

            LSQeq_comp = []
            LSQeq_uncComp = []

            # needed because this applies only to one
            for o in np.arange(orderSet)+1:
                facAmp = (2.*np.sin(self.rotDir*self.n*o*np.pi/8.))
            # xiStar = xi_eq *facAmp          
            # phiStar = phi + np.pi*(1/2.-r*ntor/8.) 
                timeTrans =  np.pi/self.om*( self.rotDir*self.n*o/8. - 1/2.) 
            
                LSQeqTmp=LSQFFT.LSQFFT()
                LSQeqTmp.initializeData(timeIn+timeTrans,(xSep1-xSep2)/facAmp ,freq_in=self.f*o,order=1,poly=polyIn,negSign=True)
                LSQeqTmp.leastSquareFit()
   
                LSQeq_comp.append(LSQeqTmp.comp[:])
                LSQeq_uncComp.append(LSQeqTmp.uncComp[:])

            LSQeq_comp = np.array(LSQeq_comp)
            LSQeq_uncComp = np.array(LSQeq_uncComp)

            LSQ1=LSQFFT.LSQFFT()
            LSQ1.initializeData(timeIn,xSep1,freq_in=self.f,order=orderSet,poly=polyIn,negSign=True)
            LSQ1.leastSquareFit()

            LSQ2=LSQFFT.LSQFFT()
            LSQ2.initializeData(timeIn,xSep2,freq_in=self.f,order=orderSet,poly=polyIn,negSign=True)
            LSQ2.leastSquareFit()

            #complex number of the control System...
            compCl = (LSQ1.comp - LSQeq_comp)

            if (np.all(timeDiag) != None) & (np.all(signalDiag) != None):

                LSQDiag=LSQFFT.LSQFFT()
                LSQDiag.initializeData(timeDiag,signalDiag,freq_in=self.f,order=orderSet,poly=polyIn,negSign=True)
                LSQDiag.leastSquareFit()
                
                compAll = np.squeeze([LSQDiag.comp, - LSQ1.comp, LSQeq_comp])
                uncCompAll = np.squeeze([LSQDiag.uncComp, - LSQ1.uncComp, LSQeq_uncComp])
                
                fac = 1./np.abs(np.sum(compAll,axis=0))
                under = np.sqrt(np.sum(np.real(compAll),axis=0)**2.*np.sum(np.real(uncCompAll)**2.) +np.sum(np.imag(compAll),axis=0)**2.*np.sum(np.imag(uncCompAll)**2.) )

                uncAmp= fac *under  
                
                print 'Amplitudes:'
                print 'All,Diag, CL '
                print np.abs(np.sum(compAll,axis=0))*1.e3,np.abs(compAll[0])*1.e3#,np.abs(np.sum(compAll[1:],axis=0)[0])*1.e3
                


                #return np.sum(compAll,axis=0), uncAmp
            if plot:
                f, (ax1, ax2, ax3 ) = plt.subplots(3, 1,  sharex=True)                
                ax1.plot(timeIn,xSep1) 
                ax1.plot(timeIn,LSQ1.getFunc(timeIn))
                ax2.plot(timeIn,xSep2)
                ax2.plot(timeIn,LSQ2.getFunc(timeIn))
                #ax3.plot(timeDiag,signalDiag-signalDiag.mean()) 
                ax3.plot(timeIn,np.abs(compCl[0])*np.sin(self.om*timeIn-np.angle(compCl[0])))
                plt.show()

                #embed()
            if (timeTraceControl & (np.all(timeDiag) == None)):
                return timeIn,np.abs(compCl[0])*np.sin(self.om*timeIn-np.angle(compCl[0]))
            elif (timeTraceControl & (np.all(timeDiag) != None)):
                return timeDiag,np.abs(compCl[0])*np.sin(self.om*timeDiag-np.angle(compCl[0]))
            if extendOut == False:
                return np.sum(compAll,axis=0), uncAmp
            elif (np.all(timeDiag) != None) & (np.all(signalDiag) != None): 
                return timeIn,xSep1,xSep2,LSQDiag.comp,LSQ1.comp,LSQ2.comp,LSQeq_comp,compCl,self.om,uncAmp
            else:
                return timeIn,xSep1,xSep2,LSQ1,LSQ2,LSQeq_comp,compCl,self.om
                #embed()
                #return amp.
                            
            
    def getAngleCorrection( self ):

        if self.IDEread:
            
            RSep = []
            zSep = []
            RMag = []
            zMag = []

            for t in self.IDE1.time[::10]:
                out = kk.KK().kkrhorz(self.Shotnumber,t,[0.999,0.0],angle=np.linspace(0.,360.,200,endpoint=True),diag = self.Diagnostic, ed = self.Edition1, exp=self.Experiment1 )
                RSep.append(np.squeeze(out.r[:,0]))
                zSep.append(np.squeeze(out.z[:,0]))
                RMag.append(np.squeeze(out.r[:,1]))
                zMag.append(np.squeeze(out.z[:,1]))

            for t in self.IDE2.time:
                out = kk.KK().kkrhorz(self.Shotnumber,t,[0.999,0.0],angle=np.linspace(0.,360.,200,endpoint=True),diag = self.Diagnostic, ed = self.Edition2, exp=self.Experiment2 )
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
 
            #embed()
            idxXMin,idxXMax = np.nanargmin(self.xIn),np.nanargmax(self.xIn)
            alphaDiag = np.arctan((self.zIn[idxXMax]-self.zIn[idxXMin])/(self.RIn[idxXMax]-self.RIn[idxXMin] ))       
            thetaDiag = np.arctan2(self.zSep_mean-zMag,self.RSep_mean-RMag)

                        #closest point of mean Equilibrium to diagnostic
            idxClose = np.argmin(np.sqrt((RSepMean-self.RSep_mean)**2.+(zSepMean-self.zSep_mean)**2.))
                        #calculate projection on the normalsurface, add +1.e-7 to avoid infinity
            facPerp = 1./(np.abs(np.cos(alphaEq[idxClose]-alphaDiag))+1.e-7)
                    #embed()     
    
            return facPerp