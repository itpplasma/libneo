
'''
Program to do ECEI crosscalibration
'''
'''
usage:
import ECEIcalib
calECI = ECEIcalib.ECEIcalib(30839)
absShot = 30839
twinBegin = 1.7
twinEnd = 2.0
calECI.relative(absShot, tBegin = twinBegin, tEnd = twinEnd,ELMsync = False , absolute = True)
calECI.write()
calECI.write_avail()
'''

__author__='Matthias Willensdorfer (mwillens@ipp.mpg.de)'
__date__='Feb 2015'
__version__='1.0'

import numpy as np
import matplotlib.pylab as plt
import IDA
import ECE
import ECEI
import dd
import IPython
import scipy.interpolate
import ElmSync
import pickle


path = '/afs/ipp-garching.mpg.de/home/m/mwillens/Programme/python/data/ECEI/'

class ECEIcalibhelp:
    status = False


class ECEIcalib:

    def __init__( self ,  Shotnumber  ):
        #check if is already done...
        self.Shotnumber = Shotnumber  
        self.relStatus = False
        self.absStatus = False
        self.avail = np.ones((16,8),dtype=bool)
  #  def calibrate( self, relShot = None, reltBegin=-1.0, reltEnd=12.0, absShot = None, abstBegin=None, abstEnd=None, calibDiag='IDA', calibExp='AUGD', calibEd=0, ELMsync = False):
   #     int relShot

# only IDA is working at the moment.
    def relative( self, Shot = None, tBegin=-1.0, tEnd=12.0, ELMsync = False, crossDiag = 'IDA', calibExp='AUGD', calibEd=0 ,absPar='Te', absolute=False, rzExp = 'ECEI', rzDiag = 'RZO', rzEd = 0):
        
        if Shot == None:
            Shot =  self.Shotnumber

        cECI=ECEI.ECEI()
        if not cECI.Load(Shot,tBegin=0.0,tEnd=tEnd,raw=True,rzExp=rzExp, rzDiag=rzDiag, rzEd=rzEd):
            return
        cECI.Binning(samplefreq=10.0)
        idxBeg = np.argmin(np.abs(cECI.rztime-tBegin))
        idxEnd = np.argmin(np.abs(cECI.rztime-tEnd))  
        calib = np.zeros((idxEnd-idxBeg, cECI.nRows, cECI.nCols) )
        #IPython.embed()
        #calib_lin = np.squeeze( np.ones_like(cECI.phi))
        #calib_sim = np.squeeze( np.ones_like(cECI.phi))
        dt = np.mean(np.diff(cECI.rztime)/2.)
        
        if crossDiag == 'IDA':

            if absPar=='Te':
                useIDATe = True
            elif absPar=='Trad':
                useIDATe = False
            else:
                useIDATe = True
            
            cIDA=IDA.IDA()
            cIDA.Load(Shot,tBegin=tBegin,tEnd=tEnd,Experiment=calibExp, Edition = calibEd ,ECEraw=True )
            count = 0        
            for i in np.arange(idxBeg,idxEnd):
                print i
                t = cECI.rztime[i]
                tIdxIDA = np.squeeze( np.where( (cIDA.time[:] > t-dt) & (cIDA.time[:] < t+dt) ) )
                tIdxECI = np.squeeze( np.where( (cECI.time[:] > t-dt) & (cECI.time[:] < t+dt) ) )
 
                if ELMsync:
                    IDAIdxTmp = ElmSync.ElmExtract(cIDA.time[tIdxIDA],Shot,plot=False,preFac = 0.2, postFac = 0.4, Experiment='AUGD')
                    ECIIdxTmp = ElmSync.ElmExtract(cECI.time[tIdxECI],Shot,plot=False,preFac = 0.2, postFac = 0.4, Experiment='AUGD')
                    tIdxIDA = tIdxIDA[IDAIdxTmp]
                    tIdxECI = tIdxECI[ECIIdxTmp]

               
                for k in np.arange(0,cECI.nRows):
                    for l in np.arange(0,cECI.nCols):
                        #print k,l
                        if useIDATe:
                            rhoIdx = np.argmin(np.abs( cIDA.rhop[tIdxIDA[0],:] - cECI.rhop[i,k,l]) ) 
                        else:
                            rhoIdx = np.argmin(np.abs( cIDA.ece_rhop[tIdxIDA[0],:] - cECI.rhop[i,k,l]) ) 
                
                        if (rhoIdx == 0) | (rhoIdx == np.size(cIDA.rhop[tIdxIDA[0],:])-1):
                            if useIDATe:
                                IDApar = np.mean( cIDA.Te[tIdxIDA,rhoIdx])
                            else:
                                IDApar = np.mean( cIDA.ece_mod[tIdxIDA,rhoIdx])
     
                            calib[count,k,l] = IDApar / ( np.mean( cECI.data[tIdxECI,k,l] )  )
                    
                        else:
                            if useIDATe:
                                f = scipy.interpolate.interp1d(  np.mean( cIDA.rhop[tIdxIDA, ::-1 ] , axis=0) ,  np.mean( cIDA.Te[tIdxIDA, ::-1 ] , axis=0) )
                            else:
                                f = scipy.interpolate.interp1d(  np.mean( cIDA.ece_rhop[tIdxIDA, ::-1 ] , axis=0) ,  np.mean( cIDA.ece_mod[tIdxIDA, ::-1 ] , axis=0) )
                            calib[count,k,l] = f ( cECI.rhop[i,k,l] )/ ( np.mean( cECI.data[tIdxECI,k,l] ) )
                            #except:
                            #IPython.embed()
                count = count + 1
            #cECI.Te = cECI.data*np.mean(calib_lin,axis=0)
            
        elif crossDiag == 'RMD':
                
           
            cRMD =  ECE.ECE()
            cRMD.Load(Shot,Diagnostic = 'RMD', Experiment=calibExp, Edition = calibEd, loadAllRhop = True )
            
            count = 0
            for i in np.arange(idxBeg,idxEnd):
                print i
                t = cECI.rztime[i]
            
                tIdxECI = np.squeeze( np.where( (cECI.time[:] > t-dt) & (cECI.time[:] < t+dt) ) )
                tIdxRMD = np.squeeze( np.where( (cRMD.time[:] > t-dt) & (cRMD.time[:] < t+dt) ) )
               

                if ELMsync:
                    RMDIdxTmp = ElmSync.ElmExtract(cRMD.time[tIdxRMD],Shot,plot=False,preFac = 0.2, postFac = 0.4, Experiment='AUGD')
                    ECIIdxTmp = ElmSync.ElmExtract(cECI.time[tIdxECI],Shot,plot=False,preFac = 0.2, postFac = 0.4, Experiment='AUGD')
                    tIdxRMD = tIdxRMD[RMDIdxTmp]
                    tIdxECI = tIdxECI[ECIIdxTmp]
  
                minRhop =  np.min(cECI.rhop[i,:,:])*0.9 
                maxRhop =  np.max(cECI.rhop[i,:,:])*1.1  

                
                Trad = cRMD.Te[tIdxRMD,:].ravel()
                rhop = cRMD.rhop[tIdxRMD,:].ravel()
                idx = np.squeeze(np.where( (rhop>minRhop)&(rhop<maxRhop)))
                Trad = Trad[idx]
                rhop = rhop[idx]
                Trad = Trad[rhop.argsort()]
                rhop = np.sort(rhop)

                if (minRhop < 1.0) & (maxRhop > 1.0) :
                    knots = np.linspace( rhop.min()*1.05, 0.99, 3 )
                    knots = np.append(knots, np.linspace( 1.005,rhop.max()*0.95, 5 ))
                else:
                    knots = np.linspace( rhop.min()*1.05, rhop.max()*0.95, 5 )
    
                LSQspline = scipy.interpolate.LSQUnivariateSpline( rhop,Trad,knots )              

                for k in np.arange(0,cECI.nRows):
                    for l in np.arange(0,cECI.nCols):

                        calib[count,k,l] = LSQspline ( cECI.rhop[i,k,l] )/ ( np.mean( cECI.data[tIdxECI,k,l] ) )


                count = count + 1 
            
##use frequency from ECE to compare, does not strongly  depend on R sep, but on shape
## not valid 
        elif crossDiag == 'freq':
                #compare channel
           
            cRMD =  ECE.ECE()
            cRMD.Load(Shot,Diagnostic = 'RMD', Experiment=calibExp, Edition = calibEd )
            
            count = 0
            for i in np.arange(idxBeg,idxEnd):
                print i
                t = cECI.rztime[i]
               
                tIdxECI = np.squeeze( np.where( (cECI.time[:] > t-dt) & (cECI.time[:] < t+dt) ) )
                tIdxRMD = np.squeeze( np.where( (cRMD.time[:] > t-dt) & (cRMD.time[:] < t+dt) ) )
                tIdxRMDrz = np.squeeze( np.where( (cRMD.rztime[:] > t-dt) & (cRMD.rztime[:] < t+dt) ) )
                #tIdxIDA = np.squeeze( np.where( (cIDA.time[:] > t-dt) & (cIDA.time[:] < t+dt) ) )

                if ELMsync:
                    RMDIdxTmp = ElmSync.ElmExtract(cRMD.time[tIdxRMD],Shot,plot=False,preFac = 0.2, postFac = 0.4, Experiment='AUGD')
                    ECIIdxTmp = ElmSync.ElmExtract(cECI.time[tIdxECI],Shot,plot=False,preFac = 0.2, postFac = 0.4, Experiment='AUGD')
                    rzRMDIdxTmp = ElmSync.ElmExtract(cRMD.rztime[tIdxRMDrz],Shot,plot=False,preFac = 0.2, postFac = 0.4, Experiment='AUGD')
                    tIdxRMD = tIdxRMD[RMDIdxTmp]
                    tIdxRMDrz = tIdxRMDrz[rzRMDIdxTmp]
                    tIdxECI = tIdxECI[ECIIdxTmp]
                    

                #first calcule row with the closest z
          ##choose channel to compare

                chfreq = np.squeeze( np.where( (cRMD.Freq >= (cECI.Freq).min()*1.e9 ) & (cRMD.Freq <= (cECI.Freq).max()*1.e9 )) )
                if ( np.size( cRMD.z[tIdxRMDrz,0]) > 1):
                    zAvgECE = np.mean(cRMD.z[:,chfreq])
                else:
                    zAvgECE = np.mean(np.squeeze(cRMD.z)[chfreq])

                zAvgECI = np.mean(cECI.z[i],axis=1)
                chooseRow = np.argmin( np.abs( zAvgECI - zAvgECE ))
                
                rhop = cECI.rhop[i,chooseRow]
    
                freq = cECI.Freq
                quasiFreqFunc = np.poly1d( np.polyfit(rhop, freq, 3))
                ##get quasi frequencies for all channels 
                quasiFreq=quasiFreqFunc(cECI.rhop[i])
                

                Trad= np.mean(cRMD.Te[tIdxRMD],axis=0) 
                idx=np.argsort(cRMD.Freq)
                x = cRMD.Freq[idx]/1.e9
                y = Trad[idx]
                
                fcalib = scipy.interpolate.interp1d(x,y ,  bounds_error=False)	

                calib[count,:,:] = fcalib(quasiFreq)/ ( np.mean( cECI.data[tIdxECI],axis=0) )

                count = count + 1


        else:
            print 'no suitable diagnostic found'    #IPython.embed()






        absCalib = np.mean(calib,axis=0)
        self.avail[np.isinf(absCalib)]=False
        if absolute == True:
            self.Calibration= absCalib
            Te = np.mean(cECI.data,axis=0) *self.Calibration
            self.absStatus = True

        self.relCalib = np.multiply(absCalib.T,1./absCalib[:,-1]).T
        self.relStatus = True

            #Te = np.mean(cECI.data,axis=0) * absCalib
            #IPython.embed()

    def absolute( self, Shot = None, tBegin=-1.0, tEnd=12.0, ELMsync = False, absDiag = 'IDA', absExp='AUGD', absPar='Te', absEd=0 ):
        print Shot

        if Shot == None:
            Shot =  self.Shotnumber

        if not self.relStatus:
            return

        cECI=ECEI.ECEI()
        if not cECI.Load(Shot,tBegin=0.0,tEnd=tEnd,raw=True):
            return
        cECI.Binning(samplefreq=10.0)

        idxBeg = np.argmin(np.abs(cECI.rztime-tBegin))
        idxEnd = np.argmin(np.abs(cECI.rztime-tEnd))
        absCalib = np.zeros((idxEnd-idxBeg, cECI.nRows, cECI.nCols) )
        dt = np.mean(np.diff(cECI.rztime)/2.)
        
        if absDiag == 'IDA':

            if absPar=='Te':
                useIDATe = True
            elif absPar=='Trad':
                useIDATe = False
            else:
                useIDATe = True
            

            cIDA=IDA.IDA()
            cIDA.Load(Shot,tBegin=tBegin,tEnd=tEnd,Experiment=absExp, Edition = absEd,ECEraw=True )
            count=0
            for i in np.arange(idxBeg,idxEnd):
                print i
                t = cECI.rztime[i]
                tIdxIDA = np.squeeze( np.where( (cIDA.time[:] > t-dt) & (cIDA.time[:] < t+dt) ) )
                tIdxECI = np.squeeze( np.where( (cECI.time[:] > t-dt) & (cECI.time[:] < t+dt) ) )

                if ELMsync:
                    IDAIdxTmp = ElmSync.ElmExtract(cIDA.time[tIdxIDA],Shot,plot=False,preFac = 0.2, postFac = 0.4, Experiment='AUGD')
                    ECIIdxTmp = ElmSync.ElmExtract(cECI.time[tIdxECI],Shot,plot=False,preFac = 0.2, postFac = 0.4, Experiment='AUGD')
                    tIdxIDA = tIdxIDA[IDAIdxTmp]
                    tIdxECI = tIdxECI[ECIIdxTmp]

                for k in np.arange(0,cECI.nRows):
                    for l in np.arange(0,cECI.nCols):
                        #print k,l
                        if useIDATe:
                            rhoIdx = np.argmin(np.abs( cIDA.rhop[tIdxIDA[0],:] - cECI.rhop[i,k,l]) ) 
                        else:
                            rhoIdx = np.argmin(np.abs( cIDA.ece_rhop[tIdxIDA[0],:] - cECI.rhop[i,k,l]) ) 

                        if (rhoIdx == 0) | (rhoIdx == np.size(cIDA.rhop[tIdxIDA[0],:]) - 1):
                            if useIDATe:
                                IDApar = np.mean( cIDA.Te[tIdxIDA,rhoIdx])
                            else:
                                IDApar = np.mean( cIDA.ece_mod[tIdxIDA,rhoIdx])
     
                            absCalib[count,k,l] = IDApar / ( np.mean( cECI.data[tIdxECI,k,l] )  )
                        else:
                            if useIDATe:
                                f = scipy.interpolate.interp1d(  np.mean( cIDA.rhop[tIdxIDA, ::-1 ] , axis=0) ,  np.mean( cIDA.Te[tIdxIDA, ::-1 ] , axis=0) )
                            else:
                                f = scipy.interpolate.interp1d(  np.mean( cIDA.ece_rhop[tIdxIDA, ::-1 ] , axis=0) ,  np.mean( cIDA.ece_mod[tIdxIDA, ::-1 ] , axis=0) )
                            
                            absCalib[count,k,l] = f ( cECI.rhop[i,k,l] )/ ( np.mean( cECI.data[tIdxECI,k,l] )  )                       
                count = count + 1
           

            self.absCalib = np.mean(absCalib,axis=0)
            IPython.embed()
            self.Calibration=  np.multiply(self.relCalib.T,self.absCalib[:,-1]).T
            Te = np.mean(cECI.data,axis=0) *self.Calibration
            self.absStatus = True


#writing data
    def write( self ):

        if self.absStatus:
                       
            file_Name = path+'%s.dat'%self.Shotnumber
            fileObject = open(file_Name,'wb') 
            pickle.dump(self.Calibration,fileObject)   
            fileObject.close()
            
    def write_avail( self ):
        

        
        if self.Shotnumber == 30839:

            self.avail[0,3] =  False
            self.avail[0,5] =  False
            self.avail[3,4] =  False
            self.avail[4,4] =  False ## but uncertain
            self.avail[5,1] =  False
         #   avail[5,7] =  False
            self.avail[7,4] =  False
            self.avail[11,0] = False
            self.avail[11,1] = False
            self.avail[11,2] = False
            self.avail[11,3] = False
            self.avail[11,4] = False ##
            self.avail[11,5] = False ##
            self.avail[11,6] = False ##
            self.avail[11,7] = False ##
            self.avail[13,6] = False ##
            self.avail[14,6] = False

        file_avail = path+'avail_%s.dat'%self.Shotnumber
        fileObject = open(file_avail, 'wb') 
        pickle.dump(self.avail, fileObject)  
        fileObject.close()
        