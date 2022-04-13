from numpy import *
from scipy.interpolate import interp1d
import pickle

## for testing
from IPython import embed
import matplotlib.pylab as plt

##modules from mwillens which are needed
import MAW
import ECEI
import ElmSync
import thetaStar
import LSQFFT
import kk_mwillens as kk
#[4.4,6.4]

class mData:



    def __init__( self , shotnumber = None, equTime = 3.2, equDiag = 'EQI', equExp = 'AUGD',rzEd=1,trange = [2.3,4.3] ):

        self.status = False
        self.statusLoad = False
        self.unload()

        if shotnumber == None:
            print 'no shotnumber'
            return

      
##equdata:
        self.shot = shotnumber
        self.equTime = mean(trange)
        self.equDiag = equDiag
        self.equExp = equExp 
###trange, first rotation, only first NBI step
        self.trange = trange

        ## edition 2 is AUGE (BTFABB) wih 0.5% correction for ECEI
        self.rzEd=rzEd
        self.status = True
        self.statusLoad = False
 
    def __del__( self):
        self.unload()
        del self.status
        del self.statusLoad


    def unload(self):
        if self.statusLoad:
            del self.mData
            self.statusLoad = False
        if self.status:
            del    self.shot
            del    self.equTime 
            del    self.equDiag
            del    self.equExp 
            del    self.trange 
            del    self.rzEd
            self.status = False


    def load(self):
        
        if self.status:
##########################################
###First load MAW data to get n number and rotation vecolity
##########################################
            MW = MAW.MAW() 
            MW.Load(self.shot, tBegin=self.trange[0], tEnd=self.trange[1])
            MW.Binning(samplefreq = 1.0)
            MW.torAngle(mode = 'rigid')
            nNumber = mean(MW.n)
            #define function to phi rotation of upper coil set
            fphi=interp1d(MW.MPtime,MW.phase[:,0])
            MW.phi0 = unwrap( ( MW.phase[:,0]- MW.phase[0,0] )* nNumber)/nNumber


############################################
## get imaging ECE
##########################################
            ECI = ECEI.ECEI()
            ECI.Load(self.shot, raw=False,rzEd=self.rzEd)
            ### bin data to reduce amount, 2kHz
            ECI.Binning(samplefreq = 2.0)

##########################################
##ELM synchronisation
##########################################
            #to check, press plot=False
            idxELM = ElmSync.ElmExtract(ECI.time,self.shot,plot=False,preFac = 0.15, postFac = 0.8, Experiment='AUGD')
            ## index for ECEI data, tidx are indices, which are inbetween ELMs
            tidx =idxELM [squeeze(where( (ECI.time[idxELM] >= self.trange[0] ) & (ECI.time[idxELM] <= self.trange[1] )))]
            rzidx = argmin(abs(ECI.rztime - self.equTime ))
            


            meanData=(ECI.data[tidx]).mean(axis=0)
            stdData= (ECI.data[tidx]).std(axis=0



           
##########################################
### get phi and delta as new timebase
##########################################
            ECI.timePhi = fphi(ECI.time[tidx])
            ECI.timePhi0 = unwrap( (ECI.timePhi-ECI.timePhi[0])*nNumber )/nNumber
            ### define new MW phase base

##########################################
### map time on toroidal phoi
##########################################
            magPhi0 = deg2rad(mean(ECI.phi))  #toroidal position 
            ECI.spacePhi = MW.getSynPhase( ECI.time[tidx],self.equTime,magPhi0,onlyPos=True)  #get  toriodal 



            #### get shape
            origShape=shape(ECI.R[rzidx])
##########################################
### get equilibrium values.q and rhop for each ECI channel
##########################################
            output = kk.KK().kkrzq( self.shot, self.equTime, ECI.R[rzidx].ravel(), ECI.z[rzidx].ravel(), exp= self.equExp, diag=self.equDiag, ed=0)
            ### get q and reshape output
            ECI.q = reshape(-output.q,origShape)
            ECI.rhopEqu = reshape(output.rho_p,origShape)


##########################################
###calculate theta star for each Channel
##########################################
            tS=thetaStar.thetaStar(shot=self.shot,time=self.equTime,Experiment = self.equExp,Diagnostic = self.equDiag)
            #define the function  for theta star grid, thetastar is set, origin is magnetic axis and              theta=thetaStar at zero
            ## calculate big grid to get accurate thetaStar values,plt.plot(tS.R,tS.z) shows grid
            tS.define(nrhop=1000,nPol=3601,origin='LFS')
            ##get theta and thetaStar for ECI
            ECI.theta = reshape( tS.get_theta(R=ECI.R[rzidx].ravel(),z=ECI.z[rzidx].ravel()) ,origShape )              
            ## this routine need rhop value to get thetaS for a given theta
            ECI.thetaS = reshape( tS.getf_theta2thetaStar(ECI.rhopEqu.ravel(), ECI.theta.ravel())  ,origShape )              

##########################################
### simple least square fit method, if the frequency is known or extracted from a reference Signal. In this case it is MW
##########################################
            LS = LSQFFT.LSQFFT()
            ### get frequency from MW measurements and calculate 3 orders
            LS.initializeData(ECI.timePhi0,ECI.data[tidx],time_ref=MW.phi0,signal_ref=MW.I_DC1,order=3,negSign=True)
            LS.leastSquareFit()    ### result is LS.amp[0] for first harmonic




##########################################
### darken the bad channels
##########################################
            
            notAvail = ( ECI.avail==False) 
            avail = logical_not(notAvail)
           
            ECI.available = avail 
 


            embed()
##########################################
### write everything into one variable to pickle data
##########################################
            self.mData = {'n': nNumber,'time': ECI.time[tidx],'phi': ECI.timePhi,'phi0': ECI.timePhi0,'spacePhi': ECI.spacePhi,'R':ECI.R[rzidx,:],'z':ECI.z[rzidx,:] ,'data' :ECI.Trad[tidx], 'rhop':ECI.rhopEqu, 'thetaS':ECI.thetaS, 'theta':ECI.theta, 'q':ECI.q, 'ampFit':LS.amp, 'phaseFit': LS.phase, 'uncPhaseFit': LS.uncPhase, 'nPhaseFit': LS.nPhase, 'meanFit': LS.mean, 'offFit': LS.off  }
            
            self.statusLoad = True

        else: 
            print 'no successful initialisation....'

# pm.dumpData(file="mData_EQI_1.9s.dat" )
    def dumpData(self, path="/ptmp1/work/mwillens/MP/", diag='ECEI',file="mData.dat"):
        
        if self.status & self.statusLoad:
            file = path+"%s/"%self.shot+diag+"/"+file
            filehandler = open(file,"wb")
            pickle.dump(self.mData,filehandler)
            filehandler.close()
            print 'file written to ' + file
