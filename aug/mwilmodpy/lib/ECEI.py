##import ECEI
##ECI=ECEI.ECEI()
##ECI.Load(30839)



import numpy as np
import dd as dd
#import dd_mwillens as dd2
import kk as kk
import ECE
import matplotlib.pylab as plt
import scipy.interpolate
import IPython
#import eqi_map as fastkk
import pickle
import tools
import PHI
import eqi_map as fastkk

path = '/afs/ipp-garching.mpg.de/home/m/mwillens/Programme/python/data/ECEI/'

class ECEIhelp:
    status = False


class ECEI:
    def __init__( self ,  Experiment = 'AUGD', Diagnostic = 'ECI', Shotnumber = None ):
        self.Status = False
        self.TradStatus = False
        self.StatusRz = False
        self.StatusECFM = False
        self.StatusTrad = False
        self.Status_allrhop = False
        self.StatusRzWarm = False
        if Shotnumber != None :
            self.Load( Shotnumber )
		

    def __del__( self ):
        self.Unload( )
        del self.Status

    def Load( self ,  Shotnumber, Experiment='AUGD', Diagnostic='ECI', Edition = 0L, tBegin=0.0, tEnd=11.0, loadAllRhop=False, rzExp = 'ECEI', rzDiag = 'RZO', rzEd = 0, raw=True , eqExp = 'AUGD', eqDiag = 'EQH', fileEndECFM='_ECFM_all.dat'):
        self.Unload()
        if Diagnostic == 'ECI':
            try:
                sf = dd.shotfile( Diagnostic, Shotnumber, Experiment, Edition)
                self.Shotnumber = Shotnumber
            except:
                print "Error reading shotfile" 
                return False

            try:
                #			       print "get timebase"
                self.tBegin = tBegin
                self.tEnd = tEnd
                            
                self.nRows = 16
                self.time = sf.getTimeBase( 'LOS1', tBegin=tBegin, tEnd=tEnd )
                self.ntime = np.size(self.time)
			#				print "get Trad-A"
                print "Reading ECEI"
                print "Reading LOS" 
                #IPython.embed()
                LOS = sf.getSignalGroup( 'LOS1', tBegin=tBegin, tEnd=tEnd )
                self.nCols = np.shape(LOS)[1]
                for i in range(2,self.nRows+1):
                    LOS = np.concatenate((LOS, sf.getSignalGroup( 'LOS%s'%i, tBegin=tBegin, tEnd=tEnd )), axis=1)
                    print "#%s"%i
                print "LOS read"
                LOS = np.reshape(LOS,(self.ntime,self.nRows,self.nCols))*10./(2.**12.) #-5. #+/- 5V / 12bit   
                self.avail = np.ones((self.nRows,self.nCols))
          
                sf.close()
                           
            except:
                print "Error reading ECI data"
                IPython.embed()
                sf.close()
                self.Status = False
                return False

            try:
                rzsf = dd.shotfile( rzDiag, Shotnumber, rzExp, rzEd )
                self.rztime = rzsf.getTimeBase( 'R', tBegin=tBegin, tEnd=tEnd )
                self.R = rzsf.getSignalGroup( 'R', tBegin=tBegin, tEnd=tEnd )
                self.z = rzsf.getSignalGroup( 'z', tBegin=tBegin, tEnd=tEnd )
                self.rhop = rzsf.getSignalGroup( 'rho', tBegin=tBegin, tEnd=tEnd )
                # 
                self.phi = PHI.ECE()+rzsf.getSignalGroup( 'phi', tBegin=tBegin, tEnd=tEnd )
                self.dphi = rzsf.getSignalGroup( 'phi', tBegin=tBegin, tEnd=tEnd )
                self.offset = rzsf.getSignalGroup( 'offset', tBegin=tBegin, tEnd=tEnd )
                


                if tBegin==0.0:
                    idx = np.squeeze(np.where( (self.time> 0.0) & (self.time< 0.02)  ))
                    if np.size(idx) != 0:
                        self.data = LOS - np.mean(LOS[idx],axis=0)  #( np.multiply(self.time,ones.T).T*self.poly_off[:,:,0]+ones*self.poly_off[:,:,1])
                    else:
                        self.data = LOS - np.squeeze(self.offset[0])
                
                else:
                    self.data = LOS - np.squeeze(self.offset[0])

           #phi = 0, sector=9
                self.Freq  = rzsf.getParameter('PAR', 'freq').data   
                rzsf.close()

                self.StatusRz = True
                #IPython.embed()
                
            except:
                print "Error reading RZO data"
                IPython.embed()
                rzsf.close()
                self.StatusRz = False
                self.Status = False
                return False 

            
            ##read pickle
            if not raw:
                #First try to read ECFM file
                try:
                    print "use ECFM data"
                    file_ECFM = path+'%s'%self.Shotnumber+fileEndECFM
                    ##f [Ghz]  c [keV / Vs] rel. std. dev [%]  R_kin [m]    z_kin [m]        tau
                    #IPython.embed()
                    data = np.genfromtxt(file_ECFM)
                    #self.ECFMtime=np.array([2.0,1.8,1.7])
                    #idxRz = np.zeros_like(self.ECFMtime,dtype='int')
                    #for i in np.arange(idxRz.size):
                    #    idxRz[i] = np.argmin(np.abs(self.rztime-self.ECFMtime[i]))               
                    Calibration = np.reshape( np.squeeze(data[:,2]) ,(data[:,0].size/(self.nRows*self.nCols),self.nRows,self.nCols))
                    self.Calibration = (Calibration.mean(axis=0))*1.e3/(10./(2.**12.))
                    self.calStd = (Calibration.std(axis=0))*1.e3/(10./(2.**12.))/self.Calibration
                    self.Trad = self.data*self.Calibration
                    Freq =  np.reshape( np.squeeze(data[:,0]) ,np.shape(Calibration))
                    Rcold =  np.reshape( np.squeeze(data[:,3]) ,np.shape(Calibration))
                    zcold =  np.reshape( np.squeeze(data[:,4]) ,np.shape(Calibration))
                    Rkin =  np.reshape( np.squeeze(data[:,5]) ,np.shape(Calibration))
                    zkin =  np.reshape( np.squeeze(data[:,6]) ,np.shape(Calibration))                    
                    tau =  np.reshape( np.squeeze(data[:,7]) ,np.shape(Calibration))
                    Rshift = Rkin - Rcold 
                    zshift = zkin - zcold
                    self.Rshift = Rshift.mean(axis=0)
                    self.zshift = zshift.mean(axis=0)
                    self.Rcorr = self.R + self.Rshift
                    self.zcorr = self.z + self.zshift
                    self.tau = tau.mean(axis=0)
                    self.StatusTrad = True
                    self.StatusECFM = True
                    
                except:
                    print "No ECFM file found"

                if self.StatusTrad == False:

                    try:
                        print "Reading calibration"
                        file_Name = path+'%s.dat'%self.Shotnumber+''
                        fileObject = open(file_Name,'r') 
                        self.Calibration = pickle.load(fileObject)  
                        fileObject.close()
                        self.Trad = self.data*self.Calibration
                        self.StatusTrad = True
                    #IPython.embed()

                    except:
                        print "No calibration file found"
                        self.StatusTrad = False

                try:
                    file_avail = path+'avail_%s.dat'%self.Shotnumber
                    fileObject = open(file_avail, 'r') 
                    self.avail = pickle.load(fileObject)  
                    fileObject.close()
                except:
                    print "No avail file to read"

                    
            

            self.Status = True

            
            if loadAllRhop:
                self.funcLoadAllRhop(eqExp = eqExp, eqDiag = eqDiag )


            return True
        else:
            print 'Diagnotic should be ECI, but is:' ,Diagnostic
            return False

    def Binning( self, samplefreq = 10.0 ):				   
        if self.Status:				  
          
            newtime, self.data = tools.dataBinning( self.time, self.data, samplefreq = samplefreq )
        
            if  self.StatusTrad:
                self.time,self.Trad = tools.dataBinning(  self.time, self.Trad, samplefreq = samplefreq )
            else:
                self.time = newtime
            
            self.ntime = np.size(self.time)


    def funcLoadAllRhop( self, eqExp = 'AUGD', eqDiag = 'EQH' ):
		
        if self.Status:

            self.eqExp = eqExp
            self.eqDiag = eqDiag
            
            self.Rall = np.zeros_like(self.data)
            self.zall = np.zeros_like(self.data)
		   
            for i in np.arange(self.ntime):
                idx = np.argmin(np.abs(self.rztime-self.time[i]))
                self.Rall[i] = self.R[idx]
                self.zall[i] = self.z[idx]		
		
            dataShape =np.array(np.shape(self.data))  
            nChannels = np.size(self.data[0])
            R =np.reshape(self.Rall,(self.ntime,nChannels))  # numpy.zeros_like(self.Te)
            z =np.reshape(self.zall,(self.ntime,nChannels))  # numpy.zeros_like(self.Te)
            #IPython.embed()
 
            
            rhopall = fastkk.eqi_map().kkrzptfn(self.Shotnumber, self.time, R, z, exp=eqExp, diag=eqDiag )
            
            #IPython.embed()
            self.rhopall = np.reshape(rhopall,dataShape)

            self.Status_allrhop = True
        else:
            
            print 'It is not possible to read Tomas kk with limited timerange, but time is not eqidistant'
            self.Status_allrhop = False
			


    def LoadRz( self ,  Shotnumber, rzExp = 'ECEI', rzDiag = 'RZO', rzEd = 0, fileEndECFM='_ECFM_all.dat',debug=False):
        try:
            rzsf = dd.shotfile( rzDiag, Shotnumber, rzExp, rzEd )
            self.rztime = rzsf.getTimeBase( 'R' )
            self.R = rzsf.getSignalGroup( 'R' )
            self.z = rzsf.getSignalGroup( 'z' )
            self.rhop = rzsf.getSignalGroup( 'rho' )
            self.nCols = np.shape(self.R)[2]
            self.nRows = np.shape(self.R)[1]
                # 
            self.phi = rzsf.getSignalGroup( 'phi' )+PHI.ECE()
            self.dphi = rzsf.getSignalGroup( 'phi' )
 
            rzsf.close()

            try:
                print "use ECFM data"
                file_ECFM = path+'%s'%Shotnumber+fileEndECFM
                    ##f [Ghz]  c [keV / Vs] rel. std. dev [%]  R_kin [m]    z_kin [m]        tau
                data = np.genfromtxt(file_ECFM)
                self.ECFMtime=np.array([2.0,1.8,1.7])
                idxRz = np.zeros_like(self.ECFMtime,dtype='int')
                for i in np.arange(idxRz.size):
                    idxRz[i] = np.argmin(np.abs(self.rztime-self.ECFMtime[i]))               
                Calibration = np.reshape( np.squeeze(data[:,2]) ,(data[:,0].size/(self.nRows*self.nCols),self.nRows,self.nCols))

                Freq =  np.reshape( np.squeeze(data[:,0]) ,np.shape(Calibration))
                Rcold =  np.reshape( np.squeeze(data[:,3]) ,np.shape(Calibration))
                zcold =  np.reshape( np.squeeze(data[:,4]) ,np.shape(Calibration))
                Rkin =  np.reshape( np.squeeze(data[:,5]) ,np.shape(Calibration))
                zkin =  np.reshape( np.squeeze(data[:,6]) ,np.shape(Calibration))                    
                tau =  np.reshape( np.squeeze(data[:,7]) ,np.shape(Calibration))
                Rshift = Rkin - Rcold 
                zshift = zkin - zcold
                self.Rshift = Rshift.mean(axis=0)
                self.zshift = zshift.mean(axis=0)
                self.Rcorr = self.R + self.Rshift
                self.zcorr = self.z + self.zshift
                self.tau = tau.mean(axis=0)
                if debug==True:
                    self.Rcold = Rcold
                    self.zcold = zcold
                    self.Rkin = Rkin
                    self.zkin = zkin
                    
                self.StatusRzWarm = True

            except:
                print 'something went wrong with loading ECFM data'
                IPython.embed()
                self.StatusRzWarm = False
                

            self.StatusRz = True

                #IPython.embed()
                
        except:
            print "Error reading RZO data"
            IPython.embed()
            rzsf.close()
            self.StatusRz = False
            self.Status = False
            return False 

    def getRzphi(self, time=None,warm=False):
        if (self.StatusRz) & (time!=None ):
            idx = np.argmin(np.abs(self.rztime-time))
            if warm==False & self.StatusRzWarm:
                return self.R[idx],self.z[idx],self.phi[idx]
            else:
                return self.Rcorr[idx],self.zcorr[idx],self.phi[idx]
    def UnloadRz( self ):
        if self.StatusRz:
            self.StatusRz = False
            del self.rztime
            del self.R
            del self.z
            del self.phi 
            del self.avail

    def Unload( self ):
        if self.Status:
            self.Status = False
            del self.tBegin
            del self.tEnd
            del self.time	
            del self.ntime
            del self.data
            del self.nRows
            del self.nCols
            del self.rztime
            del self.R
            del self.z
            del self.phi 
            del self.avail
            if self.StatusTrad:
                del self.Trad
                self.StatusTrad = False
            
            if self.Status_allrhop:
                del self.Rall 
                del self.zall
                del self.rhopall
