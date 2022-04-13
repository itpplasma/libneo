from IPython import embed
import numpy as np
import matplotlib.pylab as plt
import pickle
from os import listdir
from scipy.ndimage.interpolation import map_coordinates
from scipy.ndimage import maximum_filter
from scipy import interpolate

import dd as dd

import LSQFFT
import MAWnew as MAW
import biotsavart
import ElmSync

class controlCoils:
	def __init__( self ,  shotnumber = None, tBegin = 0.0, tEnd = 11. ):
		
		self.statusControlRead = False	

		if shotnumber != None :
			self.load( Shotnumber, tBegin = tBegin, tEnd = tEnd )
			self.loadPsis()

        def loadCurrents( self ,  Shotnumber, tBegin=0.0, tEnd=11., useELM=False, ELMperiod=[0.1,0.6] ):

            try:
                sf = dd.shotfile( 'EQH', Shotnumber, 'AUGD', 0)	
            except:
                print "Error reading MBI shotfile" 
                return

            self.shot = Shotnumber
            self.tBegin = tBegin
            self.tEnd = tEnd 

            ## load currents
            time = sf.getSignal( 'time' )
	
	    currents = sf.getSignalGroup('CLD')
	    #currentsCLE = sf.getSignalGroup('CLE')
	    Rcl = sf.getSignalGroup('Rcl')
	    zcl = sf.getSignalGroup('Zcl')
	    sf.close()

# 1  IV1o     
#  2  IV1u    
#  3  IV2o  
#   4  IV2u  
#   5  IV3o   
#   6  IV3u     
#   7  Ipslon    
#   8  Ipslun    
#   9  ICoIo  
#  10  ICoIu      
#  11  IOH2od   
#  12  IOH2ud    
#  13  IOH     
            tidx = np.where((time >= tBegin) & (time<=tEnd))[0]
            
	    if useELM:
		    idxELM = ElmSync.ElmExtract(time[tidx],Shotnumber,plot=False,preFac = ELMperiod[0], postFac = ELMperiod[1], Experiment='AUGD')
		    tidx = tidx[idxELM]
		    
	    self.time = time[tidx]

	    fac=1.
	   
#Rcl

	    self.IV1o =  currents[tidx,0]*fac  ## to get it in kA 
	    self.RIV1o,self.zIV1o = Rcl[0],zcl[0]

            self.IV1u =  currents[tidx,1]*fac
	    self.RIV1u,self.zIV1u = Rcl[1],zcl[1]

            self.IV2o =  currents[tidx,2]*fac 
	    self.RIV2o, self.zIV2o = Rcl[2],zcl[2]

            self.IV2u =  currents[tidx,3]*fac 
	    self.RIV2u, self.zIV2u = Rcl[3],zcl[3]  
         
            self.IV3o =  currents[tidx,4]*fac  
	    self.RIV3o, self.zIV3o = Rcl[4],zcl[4]

            self.IV3u =  currents[tidx,5]*fac        
	    self.RIV3u, self.zIV3u = Rcl[5],zcl[5]

            self.Ipslon =  currents[tidx,6]*fac   
	    self.RIpslon, self.zIpslon = Rcl[6],zcl[6]

            self.Ipslun =  currents[tidx,7]*fac     
	    self.RIpslun, self.zIpslun = Rcl[7],zcl[7]

            self.ICoIo =  currents[tidx,8]*fac   
	    self.RICoIo, self.zICoIo = Rcl[8],zcl[8]  
     
            self.ICoIu =  currents[tidx,9]*fac   
	    self.RICoIu, self.zICoIu = Rcl[9],zcl[9]
   
	    self.IHO2od =  currents[tidx,10]*fac      
	    self.RIHO2od, self.zIHO2od = Rcl[10],zcl[10]
 
            self.IHO2ud =  currents[tidx,11]*fac      
	    self.RIHO2ud, self.zIHO2ud = Rcl[11],zcl[11]

            self.IOH =  currents[tidx,12]*fac
	    self.RIOH, self.zIOH = Rcl[12],zcl[12]

            ## order of list in shell important!!
            self.currents =np.array([self.ICoIo,self.ICoIu,self.IHO2od,self.IHO2ud,self.IOH,self.IV1o,self.IV1u,self.IV2o,self.IV2u,self.IV3o,self.IV3u,self.Ipslon,self.Ipslun]).T
	    self.Rcoils = np.array([self.RICoIo,self.RICoIu,self.RIHO2od,self.RIHO2ud,self.RIOH,self.RIV1o,self.RIV1u,self.RIV2o,self.RIV2u,self.RIV3o,self.RIV3u,self.RIpslon,self.RIpslun])
	    self.zcoils = np.array([self.zICoIo,self.zICoIu,self.zIHO2od,self.zIHO2ud,self.zIOH,self.zIV1o,self.zIV1u,self.zIV2o,self.zIV2u,self.zIV3o,self.zIV3u,self.zIpslon,self.zIpslun])
	    self.names = np.array(['ICoIo','ICoIu','IO2od','IO2ud','IOH','IV1o','IV1u','IV2o','IV2u','IV3o','IV3u','IPSLON','IPSLUN'])
	    self.namesErika =np.array(['ICoIo','ICoIu','IHO2od','IHO2ud','IOH','IV1o','IV1u','IV2o','IV2u','IV3o','IV3u','Ipslon','Ipslun'])


	    #Psi_1kA_ICoIo
	    #Psi_1kA_ICoIu
	    #Psi_1kA_IHO2od
	    #Psi_1kA_IHO2ud
	    #Psi_1kA_IOH
	    #Psi_1kA_IV1o
	    #Psi_1kA_IV1u
	    #Psi_1kA_IV2o
	    #Psi_1kA_IV2u
	    #Psi_1kA_IV3o
	    #Psi_1kA_IV3u
	    #Psi_1kA_Ipslon
	    #Psi_1kA_Ipslun

 #psi = biotsavart.psiFlux_Amp_Current(1.65,0.,1.e6,RR,zz)
            print 'coil currents loaded'

            self.statusCurrents = True


        ## reading Psi Matrix
        def loadPsis( self ,  path='/afs/ipp-garching.mpg.de/home/m/mwillens/Programme/python/MP/displacement/disfromMFBE/',rootName='Psi_1kA_'):

            onlyfiles = np.array([f for f in listdir(path) if f[:len(rootName)]==rootName])
	    sortedFiles = np.sort(onlyfiles)
            PsiList = []
            RCoil = []
            zCoil = []
	    
            for file in sortedFiles:
                print file
                filehandler = open(path+file,'rb')
                data=pickle.load(filehandler)
                PsiList.append(data['Psi'])
                RCoil.append(data['R'])
                zCoil.append(data['z'])
                filehandler.close()

            return  np.array(RCoil), np.array(zCoil), np.array(PsiList)*1.e-3 #1kA


        ### influence of the MP coils
        def MPinfluence(self,orderOut=1,orderFit=4):
            if self.statusCurrents:
                
                MW = MAW.MAW(self.shot)
                idx = np.where((MW.PSLtime>self.tBegin) & (MW.PSLtime<self.tEnd ))[0]
                MW_Icoil = MW.Icoils[idx,0,0]
                MW_time = MW.time[idx]

                LSQ = LSQFFT.LSQFFT()
                LSQ.initializeData(self.time, self.currents ,time_ref=MW_time,signal_ref=MW_Icoil,order=3,poly=3,negSign=True)
		LSQ.initializeData(self.time, self.currents ,freq_in=np.round(np.copy(LSQ.freq0)),order=orderFit,poly=3,negSign=True)
                LSQ.leastSquareFit()
                          
                self.currentsFiltered = LSQ.getFunc(self.time)
                currentsFund = LSQ.getFunc(self.time,order=orderOut,poly=0)
                self.currentsFund = currentsFund-currentsFund.mean(axis=0)
                self.currentsAmp = LSQ.amp[0]
		self.currentsAmpAll = LSQ.amp
		print LSQ.amp[0]
		self.currentsReal = LSQ.real[0]
		self.currentsImag = LSQ.imag[0]
		self.currentsFreq = LSQ.freq0
		self.statusControlRead = True
                #embed()

            else:
                print 'no currents loaded'
  
	def makePsiOfEachCoil(self, RIn, zIn, path='/afs/ipp-garching.mpg.de/home/m/mwillens/Programme/VMEC/coils/'):
		nCoil = self.currents[0,:].size
		RTmp,zTmp = np.meshgrid(RIn,zIn)
		RR,zz = RTmp.T, zTmp.T
		newPsiCoil = np.zeros((nCoil,RIn.size,zIn.size))
		for i in np.arange(nCoil):
			dataTmp = np.genfromtxt(path+self.names[i]).T
			nData = dataTmp[0].size
			totalTurns =  np.sum(dataTmp[2])
			PsiTmp = np.zeros((RIn.size,zIn.size))
			for j in np.arange(nData):
				PsiTmp += biotsavart.psiFlux_Amp_Current(dataTmp[0,j],dataTmp[1,j],dataTmp[2,j],RR,zz)/totalTurns
		
			newPsiCoil[i] = PsiTmp

		return newPsiCoil

	def read_Vacfield_pfcur(self,ntime=16,period=1,filename='VACFIELD_pfcur_last', path='/afs/ipp-garching.mpg.de/home/m/mwillens/Programme/VMEC/AUGD/33118/micdu_eqb_1_4.156/shot_kk/',pathOut='./',fileOutRoot='VACFIELD_pfcur_',doComplex=False):
		if self.statusControlRead:
			file = open(path+filename, 'r')
			header = file.readlines()[:5]
			file.close()
			refLine1= '            %.6s'
			refLine2= '            %.8f\n'
			data = np.genfromtxt(path+filename,skip_header=5,dtype='string').T
			current = np.array(data[1],dtype='float')
			coilLabel = np.array(data[0],dtype='string')
			#embed()
                        ## file Out
			tBegin=self.time[0]
			tEnd = self.time[0]+1./self.currentsFreq
			timeNew = np.linspace(tBegin,tEnd,ntime,endpoint=False)
			
			for t in timeNew:
				filenameOut = pathOut+fileOutRoot+'%.5fs'%t
				file = open(filenameOut, 'w')
				file.writelines(header)
				for c in np.arange(coilLabel.size):
					idxCoil = np.where(self.namesErika==coilLabel[c])[0]
					currentNew = current[c]+np.interp(t,self.time,np.squeeze(self.currentsFund[:,idxCoil]))*1.e-6
					file.writelines(refLine1%coilLabel[c]+refLine2%currentNew)
				file.close()
		

			fileRef = open(pathOut+fileOutRoot+'Ref', 'w')
			fileRef.writelines(header)

			if doComplex:
				fileRe = open(pathOut+fileOutRoot+'Real', 'w')
				fileRe.writelines(header)
				fileIm = open(pathOut+fileOutRoot+'Imag', 'w')
				fileIm.writelines(header)

			for c in np.arange(coilLabel.size):
				idxCoil = np.where(self.namesErika==coilLabel[c])
				currentRef = current[c]
				fileRef.write(refLine1%coilLabel[c]+refLine2%currentRef)
				if doComplex:
					currentReal = current[c]+self.currentsReal[idxCoil]*1.e-6 #MA
					currentImag = current[c]+self.currentsImag[idxCoil]*1.e-6 #MA
					fileRe.write(refLine1%coilLabel[c]+refLine2%currentReal)		       
					fileIm.write(refLine1%coilLabel[c]+refLine2%currentImag)	
				
			fileRef.close()	

			if doComplex:
				fileRe.close()
				fileIm.close()
			

	def calcPsiCoil(self, RIn, zIn, analytic=True):
		nCoil = self.currents[0,:].size
		newPsiCoil = np.zeros((nCoil,RIn.size,zIn.size))

		if analytic:
			RR,zz = np.meshgrid(RIn,zIn)
			for i in np.arange(nCoil):
				#1A as inout
				newPsiCoil[i] = biotsavart.psiFlux_Amp_Current(self.Rcoils[i],self.zcoils[i],1.0,RR.T,zz.T)

			return newPsiCoil
	
		else:

			RCoil,zCoil,PsiCoil = self.loadPsis()
			origShape = [RIn.size,zIn.size]

			for i in np.arange(nCoil):
				fR = np.interp(RIn,RCoil[i],np.arange(RCoil[i].size))
				fz = np.interp(zIn,zCoil[i],np.arange(zCoil[i].size))
				RR,zz = np.meshgrid(fR,fz)
				inputArr = np.vstack(((RR.T).ravel(),(zz.T).ravel()))
				newPsiCoil[i] = np.reshape(map_coordinates( PsiCoil[i],inputArr ),origShape)

			return newPsiCoil


	def newPsi(self, RIn, zIn, PsiIn, timeIn = None, unFiltered = False):
            
		if np.all(timeIn)==None:
			timeIn = self.time
			
		ntime = timeIn.size
		nCoil = self.currents[0,:].size
		selectCoils = np.ones((nCoil))
		#selectCoils[2:5] = 0.0
		unSelectCoils = np.zeros((nCoil))
		newcurrentsFund = np.zeros((timeIn.size,nCoil))
		
            
		for i in range(nCoil):
			if unFiltered:
				newcurrentsFund[:,i] = np.interp(timeIn,self.time,self.currents[:,i]-self.currents[0,i])	
			else:
				newcurrentsFund[:,i] = np.interp(timeIn,self.time,self.currentsFund[:,i])
            
		#embed()	
		newPsiCoilMod = self.makePsiOfEachCoil(RIn,zIn)
		PsiOutMod = np.einsum('nrz,tn,n->trz',newPsiCoilMod,newcurrentsFund,selectCoils) + np.array([PsiIn,]*timeIn.size)
		
#		PsiTest =  np.copy(np.array([PsiIn,]*timeIn.size) )
#		for t in range(timeIn.size):
#			for i in range(nCoil):
#				PsiTest[t] += newcurrentsFund[t,i]*newPsiCoilMod[i]

	
		RMagMod,zMagMod,PsiMagMod = getMag(RIn,zIn,PsiOutMod,resolution = 2.e-5)
		RXMod,zXMod,PsiXMod = getXpoint(RIn,zIn,PsiOutMod)
		rhopMod = getRhopM(PsiOutMod,PsiXMod,PsiMagMod)
				
#		RMagFile,zMagFile,PsiMagFile = getMag(RIn,zIn,PsiOutFile)
#		RMagTest,zMagTest,PsiMagTest = getMag(RIn,zIn,PsiTest)
		return rhopMod


#		newPsiCoilAna = self.calcPsiCoil(RIn,zIn)
#		newPsiCoilFile = self.calcPsiCoil(RIn,zIn,analytic=False)
		#unSelectCoils[4]=1.0		#embed()
#siOutFile = np.einsum('nrz,tn->trz',newPsiCoilFile,newcurrentsFund) + np.array([PsiIn,]*timeIn.size)


def getRhopM( PsiRz, PsiSep, PsiMag ):
    
    
    rhopM = np.sqrt(np.abs((PsiRz.T-PsiMag)/(PsiSep-PsiMag))).T
    return rhopM


def getXpoint( RIn, zIn, PsiRz, resolution = 2.e-4, window = 0.6 ):
    
  
    maxidx_z = PsiRz.shape[2]/2  #only first half of z indices -> only searching forXpoint at negative values of z 
    zIn = zIn[:maxidx_z]
    PsiRz = PsiRz[:,:,:maxidx_z]
    
    Rdiff = np.diff(RIn).mean()
    zdiff = np.diff(zIn).mean()
    
    ntime = PsiRz[:,0,0].size
    nArea = PsiRz[0,:,:].size
    nr = PsiRz[0,:,0].size
    nz = PsiRz[0,0,:].size
    origShape=np.shape(PsiRz[0])
    
    funcR = interpolate.interp1d(np.arange(origShape[0]),RIn)
    funcz = interpolate.interp1d(np.arange(origShape[1]),zIn)
    
    PsiRzChanged= np.reshape(PsiRz,(ntime,nArea))
    
    PsiSep, XpointR, Xpointz = np.zeros(ntime),np.zeros(ntime),np.zeros(ntime)
    #zidx, Ridx = np.zeros((ntime)), np.zeros((ntime))
    
    for t in range(ntime):
        RIdxs = PsiRz[t].argmax(axis=0) # R indexes where Psi is max (for every z row) 
        PsiMaxZ = PsiRz[t,RIdxs,np.arange(nz)] # maximum row values of Psi along z
        zIdx = np.argmin(PsiMaxZ) # z index of saddle point (the minimum of Psi onthe Psimax line)
        RIdx = RIdxs[zIdx] # R index of saddle point
        
        fR = np.arange(RIdx-window,RIdx+window,resolution/Rdiff)
        fz = np.arange(zIdx-window/2.,zIdx+window/2.,resolution/zdiff)
        
        RR,zz = np.meshgrid(fR,fz)
        inputArr = np.vstack((RR.T.ravel(),zz.T.ravel()))
        Psizoom = np.reshape(map_coordinates( PsiRz[t], inputArr ),(fR.size,fz.size))
        
        Ridxszoom = Psizoom.argmax(axis=0)
        PsiMaxZzoom = Psizoom[Ridxszoom,np.arange(fz.size)]
        zIdxzoom = np.argmin(PsiMaxZzoom) # z index of saddle point (the minimum ofPsi on the Psimax line)
        RIdxzoom = Ridxszoom[zIdxzoom] # R index of saddle point
        
        PsiSep[t] = Psizoom[RIdxzoom,zIdxzoom]
        XpointR[t], Xpointz[t] = funcR(fR[RIdxzoom]), funcz(fz[zIdxzoom])
    
    return XpointR, Xpointz, PsiSep


def getMag( RIn, zIn, PsiRz, resolution = 2.e-4, window = 0.6):

    
    Rdiff = np.diff(RIn).mean()
    zdiff = np.diff(zIn).mean()
    
    ntime = PsiRz[:,0,0].size
    nArea = PsiRz[0,:,:].size
    origShape=np.shape(PsiRz[0])
    
    funcR = interpolate.interp1d(np.arange(origShape[0]),RIn)
    funcz = interpolate.interp1d(np.arange(origShape[1]),zIn)
    
    PsiRzChanged= np.reshape(PsiRz,(ntime,nArea))
    idxMax = np.unravel_index(np.argmax(PsiRzChanged,axis=1) , origShape) 
    RIdx = idxMax[0]
    zIdx = idxMax[1]
    
    PsiMag = np.zeros(ntime)
    RMag,zMag = np.zeros_like(PsiMag),np.zeros_like(PsiMag)
    
    for t in range(ntime):
        fR = np.arange(RIdx[t]-window,RIdx[t]+window,resolution/Rdiff)
        fz = np.arange(zIdx[t]-window/2.,zIdx[t]+window/2.,resolution/zdiff)
        RR,zz = np.meshgrid(fR,fz)
        inputArr = np.vstack((RR.T.ravel(),zz.T.ravel()))
        PsiOut = np.reshape(map_coordinates( PsiRz[t], inputArr ),(fR.size,fz.size))
    
        PsiMag[t] = PsiOut.max()
        Magidxs = np.unravel_index(PsiOut.argmax(), PsiOut.shape)
        RMag[t] = funcR(fR[Magidxs[0]])
        zMag[t] = funcz(fz[Magidxs[1]])
        
    
    return RMag, zMag, PsiMag


