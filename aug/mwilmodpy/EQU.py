##CALL
##zB: EQ = EQU.EQU()
##    EQ.Load(34444,Experiment='MICDU', Diagnostic='EQB')
##
import numpy as np
import dd
from IPython import embed
from scipy.ndimage.interpolation import map_coordinates

class EQU:

    def __init__( self ,  Shotnumber = None ):
        self.Status = False
        if Shotnumber != None :
            self.Load( Shotnumber )
		
    def __del__( self ):
        self.Unload( )
        del self.Status

    def Load( self , Shotnumber, Experiment='AUGD', Diagnostic='EQI', Edition = 0L, tBegin=0., tEnd=14., Diagnostic2='EQI'):

        self.Unload()
        if Diagnostic == 'EQI' or Diagnostic == 'EQH' or Diagnostic == 'EQB' or Diagnostic == 'IDE' or Diagnostic == 'FPP':
            try:
                if np.all(Diagnostic2 == Diagnostic):
                    sf = dd.shotfile( Diagnostic, Shotnumber, Experiment, Edition)
                else:
                    sf = dd.shotfile( Diagnostic, Shotnumber, Experiment, Edition,Diagnostic2)
                self.Shotnumber = Shotnumber
            except:
                print "Error reading shotfile"
                embed()
                return False

            self.diag = Diagnostic
            self.exp = Experiment
            self.ed = Edition
            self.Nz = sf.getParameter('PARMV','N').data+1
            self.NR = sf.getParameter('PARMV','M').data+1
            Ntime = sf.getParameter('PARMV','NTIME').data
            time = (sf.getSignal("time"))[0:Ntime]
            idxTime = np.where( (time>=tBegin) & (time<=tEnd) )[0]
            self.time = time[idxTime]
            
            self.R = (sf.getSignalGroup("Ri"))[0:Ntime,0:self.NR][idxTime]
            self.z = (sf.getSignalGroup("Zj"))[0:Ntime,0:self.Nz][idxTime]
            self.PsiOrigin = sf.getSignalGroup("PFM")[0:Ntime,0:self.Nz,0:self.NR][idxTime]
            #self.Rinv = (sf.getSignalGroup("Rinv"))[0:self.Ntime]
            self.PFL = (sf.getSignalGroup("PFL"))[0:Ntime][idxTime]
            ###magnetic axis,sxm
            self.PsiSpecial = sf.getSignalGroup("PFxx")[0:Ntime][idxTime]   
            ##time, R, z
            self.Psi = np.swapaxes(self.PsiOrigin,1,2)
            self.PsiAxis = self.PsiSpecial[:,0]
            self.PsiSep = self.PsiSpecial[:,1]
            #self.PsiAxisM = self.make_RZ_matrix(self.PsiAxis) #Zur Berechnung von rhopMatrix, selbe Werte auf alle R und z schreiben
            #self.PsiSepM = self.make_RZ_matrix(self.PsiSep)
            
            #self.rhopM = np.sqrt((self.Psi-self.PsiAxisM)/(self.PsiSepM-self.PsiAxisM)) #rhopMatrix
            self.rhopM = np.sqrt(np.abs((self.Psi.T-self.PsiAxis)/(self.PsiSep-self.PsiAxis))).T
            #IPython.embed()
            sf.close()
            self.Ntime = self.time.size
            self.Status = True
            
            try:
                if Diagnostic == 'IDE':
                    DiagnosticG = 'IDG'
                elif Diagnostic == 'EQH':
                    DiagnosticG = 'GQH'
                elif Diagnostic == 'EQI':
                    DiagnosticG = 'GQI'
                elif Diagnostic == 'FPP':
                    DiagnosticG = 'FPG'
                else:
                    return

                sfG = dd.shotfile( DiagnosticG, Shotnumber, Experiment, Edition)

                self.Rmag = sfG.getSignal("Rmag")[idxTime]
                self.zmag = sfG.getSignal("Zmag")[idxTime]
                self.Raus = sfG.getSignal("Raus")[idxTime]
                self.zaus = self.zmag
                self.zoben = sfG.getSignal("Zoben")[idxTime]
                self.Rzoben = sfG.getSignal("Rzoben")[idxTime]
                self.zunten = sfG.getSignal("Zunt")[idxTime]
                self.Rzunten = sfG.getSignal("Rzunt")[idxTime]              
                self.Zsquad = sfG.getSignal("Zsquad")[idxTime]
                sfG.close()

            except:
                print "Error reading G shotfile"
                embed()
                


    def Unload( self ):
        if self.Status:
            self.Status = False
            del self.Nz
            del self.NR
            del self.Ntime
            del self.R
            del self.z
            del self.time
            del self.Psi

    def getRhop( self, time_in, R_in, z_in ):
        
        if (np.size(R_in) ==  np.size(z_in)):
            R_matrix = self.R[0] 
            z_matrix = self.z[0]
            time_matrix = self.time
            rhopM = self.rhopM 
            
        ### calculate rhop at input coordinates
            ntime = time_in.size
            rhop = np.zeros((ntime, R_in.size))
            
            fR = np.interp(R_in,R_matrix,np.arange(R_matrix.size))
            fz = np.interp(z_in,z_matrix,np.arange(z_matrix.size))
            ftime = np.interp(time_in,time_matrix,np.arange(time_matrix.size))
        #coordinates = np.vstack((Rtmp,ztmp))
        
            RR = np.tile(fR,ntime)
            zz = np.tile(fz,ntime)
            tt = np.repeat(ftime,R_in.size)
        #rhopM  indizes time, R, z
            inputArr = np.vstack((tt,RR,zz))
            rhopOut = np.reshape(map_coordinates( rhopM, inputArr ),(ntime,R_in.size))
        #for i in range(0,IDE_in.time.size):
        #    rhop[i] = map_coordinates( rhopM[i], coordinates )
            
            return rhopOut
        else:
            print 'Rin and zIn do not have the same size'


    def getTime( self, T_begin, T_end ):
        if self.Status:
            idx_begin = np.argmin(np.abs(self.time-T_begin))
            idx_end = np.argmin(np.abs(self.time-T_end))
            return self.time[idx_begin:idx_end]

    def getTimeIdx( self, T_begin, T_end ):
        if self.Status:
            idx_begin = np.argmin(np.abs(self.time-T_begin))
            idx_end = np.argmin(np.abs(self.time-T_end))
            return idx_begin, idx_end
            
    def getRaus( self, T_begin = 0, T_end = 15 ):
        if self.Status:
            idx_begin = np.argmin(np.abs(self.time-T_begin))
            idx_end = np.argmin(np.abs(self.time-T_end))
            return self.Raus[idx_begin:idx_end]

    def getRhopM( self, timepoint ): #Rhop[R, z]
        if self.Status:
            idx = np.argmin(np.abs(self.time-timepoint))
            return self.rhopM[idx]

    def make_RZ_matrix( self, PsiIn ):
        PsiOut = np.zeros((self.Ntime, self.NR, self.Nz))
        for i in range(0,self.NR):
            for j in range(0,self.Nz):                
                PsiOut[:,i,j] = PsiIn
        return PsiOut

    def getPsi( self, timepoint,origin=False ):
        if self.Status:
            idx = np.argmin(np.abs(self.time-timepoint))
            if origin == True:
                return self.PsiOrigin[idx]
            else:
                return self.Psi[idx]
  
    def getR( self, timepoint ):
        if self.Status:
            idx = np.argmin(np.abs(self.time-timepoint))
            return self.R[idx]

    def getz( self, timepoint ):
        if self.Status:
            idx = np.argmin(np.abs(self.time-timepoint))
            return self.z[idx]

    def __call__( self , timepoint ,origin=False):
        if self.Status:
            idx = np.argmin(np.abs(self.time-timepoint))
            if origin == True:
                return self.R[idx],self.z[idx], self.PsiOrigin[idx] 
            else:               
                return self.R[idx],self.z[idx], self.Psi[idx]