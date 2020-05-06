import numpy as np
from scipy.interpolate import interp1d
from scipy import stats

#from  repository
#setenv PYTHONPATH "$PYTHONPATH":/afs/ipp-garching.mpg.de/aug/ads-diags/common/python/lib
import dd

#for testing purpose
import IPython
import matplotlib.pylab as plt
#/afs/ipp-garching.mpg.de/home/m/mwillens/Programme/python/MP/IllustVideos/dphiScan
from scipy import signal

#### how to run it 
# import matplotlib.pylab as plt
# import MAWnew as MAW

#MW=MAW.MAW()
#MW.Load(32217,Binning=0.1)
#plt.title('Amplitude during 1Hz rigid rotation, 32217')
#plt.xlim([-2.,6.5])
#plt.plot(MW.MPtime,MW.amp[1],label='n=%d'%MW.nn[1])
#plt.plot(MW.PSLMPtime,MW.PSLamp[1],label='n=%d wi PSL'%MW.nn[1])
#plt.plot(MW.MPtime,MW.amp[5],label='n=%d'%MW.nn[5])
#plt.plot(MW.PSLMPtime,MW.PSLamp[5],label='n=%d wi PSL'%MW.nn[5])
#plt.xlabel('time [s]',fontsize=22)
#plt.ylabel('Amp [A]',fontsize=22)
#plt.legend(loc=2,fontsize=20)
#plt.show()

## fast switch off
#MW=MAW.MAW()
#MW.Load(32085,Binning=1.0)
#plt.title('Fast switch off, 32085')
#plt.xlim([1.,7.5])
#plt.plot(MW.time,MW.Icoils[:,0,0],label='IBu1')
#plt.plot(MW.PSLtime,MW.PSLIcoils[:,0,0],label='IBu1 wi PSL')
#plt.xlabel('time [s]',fontsize=22)
#plt.ylabel('I [A]',fontsize=22)
#plt.legend(fontsize=22)
#plt.show()

## differential phase scan
#MW = MAW.MAW() 
#MW.Load(30826,Binning=1.0)
#plt.title('Differential phase scan, 30826')
#plt.plot(MW.MPtime,np.rad2deg(MW.dPhase), label=r'$\Delta \phi$')
#plt.plot(MW.PSLMPtime,np.rad2deg(MW.PSLdPhase), label=r'$\Delta \phi$ wi PSL')
#plt.xlabel('time [s]',fontsize=22)
#plt.ylabel(r'$\Delta \phi \rm{[deg]}$',fontsize=24)
#plt.legend(fontsize=22)
#plt.show()
# print MW.n

#get effective currents for time 3.3
# MW.getIcoils([3.3],usePSL=True)
#usual currents
# MW.getIcoils([3.3],usePSL=False)


#other parameters
###  MW.time...timebase, MW.I_DC1...Current supply, MW.Icoils...currents of every coil

#evaluated parameters
### MW.dPhase.. differential Phase   MW.amp


class MPcoils:

#reading 
    def __init__( self , nRows = 2, nCoils = 8, labelRows = ['u','l']):
        self.Status = False
        self.Status_PSL=False
        ## generate Stringarray
        if nRows != np.size(labelRows):
            print 'number of Rows do not have the same'
            return False

        self.nRows =  nRows
        self.nCoils =  nCoils
        # generate Stringarray to read coil stuff#
        coilStr = np.chararray((nRows, nCoils),itemsize=2)
        for i in range(nRows):
            for j in range(nCoils):
                coilStr[i,j] = labelRows[i]+str(j+1) 
            
	self.coilStr = coilStr

    def __del__( self ):
        del self.Status

        
    def LoadCoord( self, path='/afs/ipp-garching.mpg.de/home/m/mwillens/Programme/python/data/Bcoils/', preFile='B',postFile='.asc'):
        
#generate file to read
        fileCoilStr = path + preFile +  self.coilStr + postFile
        
#X=hallo[:,0]*np.cos(hallo[:,2]/2./np.pi)
#Y=hallo[:,0]*np.sin(hallo[:,2]/2./np.pi)
        #coordPol = [] #np.array(np.shape(fileCoilStr))
        #coordCar = [] #np.array(np.shape(fileCoilStr))

        tmpRows=[]
        for i in range(self.nRows):
            tmpCoils=[]
            for j in range(self.nCoils):
                #R[m] z[m] phi[degree]
                tmpCoils.append(np.genfromtxt(fileCoilStr[i,j]))
                #ax.plot3D(coordCar[i,j,:,0],coordCar[i,j,:,1],coordCar[i,j,:,2])
            tmpRows.append(tmpCoils)

        coordPol =  np.squeeze(tmpRows)
         #X[m] Y[m] Z[m]
        coordCar = np.zeros_like(coordPol)
        coordCar[:,:,:,0] = coordPol[:,:,:,0] * np.cos(coordPol[:,:,:,2] )
        coordCar[:,:,:,1] = coordPol[:,:,:,0] * np.sin(coordPol[:,:,:,2] )
        coordCar[:,:,:,2] = coordPol[:,:,:,1] 

        self.coordPol = coordPol
        self.coordCar = coordCar

  #  ax.plot3D(coordCar[i,j,:,0],coordCar[i,j,:,1],coordCar[i,j,:,2])



## index null sind die oberen Spulen
class MAW:
	def __init__( self , Shotnumber = None, Experiment = 'AUGD', Diagnostic = 'MAW',Binning=None ,iir=False, Edition = 0L,thres=100.):
		self.Status = False
		self.Status_PSL = False
		self.Status_torangle = False
		self.Status_torangle_PSL = False
		if Shotnumber != None :
			self.Load( Shotnumber ,Binning=Binning,iir=iir,Diagnostic=Diagnostic,Edition = Edition,thres=thres)
		

	def __del__( self ):
		self.Unload( )
		del self.Status

##load data
	def Load( self, Shotnumber ,Experiment= 'AUGD', Diagnostic= 'MAW', tBegin=-1.0, tEnd=16.0, Edition = 0L, nRows = 2, nCoils = 8, labelRows = ['u','l'],Binning=None,thres=100.,iir=False):

		self.Unload()
		
		if (Diagnostic == 'MAW')  :

			##create String to read Signal
			coilStr = np.chararray((nRows, nCoils),itemsize=2)
			for i in range(nRows):
				for j in range(nCoils):
					coilStr[i,j] = labelRows[i]+str(j+1) 

			signalCoilStr = 'IB' +  coilStr 
			
			try:
				
				sf = dd.shotfile( Diagnostic , Shotnumber,Experiment, Edition)
				
				self.time = sf( signalCoilStr[0,1], tBegin=tBegin, tEnd=tEnd ).time
				self.ntimes = np.size(self.time)
				self.samplingrate = 1.0/np.mean(np.diff(self.time))

				#self.I_DC1 = sf( 'I_DC1', tBegin=tBegin, tEnd=tEnd ).data 
				#self.I_DC2 = sf( 'I_DC2', tBegin=tBegin, tEnd=tEnd ).data
				
				self.Icoils = np.zeros((self.ntimes,nRows,nCoils))
                                
                                for i in range(nRows):
                                    for j in range(nCoils):
                                        try:
					    self.Icoils[:,i,j] = sf( signalCoilStr[i,j], tBegin=tBegin, tEnd=tEnd ).data
                                        except:
                                            self.Icoils[:,i,j] = 0.0

                                             
				self.nRows = nRows
				self.nCoils = nCoils
				self.Status = True
				sf.close()

			except:
				print('Error in loading MAW')
                                IPython.embed()

                elif ((Diagnostic == 'SSV') & (Shotnumber > 0)):
                    try:
                        print('loading SSV')
                        coilStr = np.chararray((nRows, nCoils),itemsize=6)
                        id=0
                        coillabel = "Iact"
		        for i in range(nRows):
			    for j in range(nCoils):
			        coilStr[i,j] = coillabel+str(id+1)
                                id = id + 1
                            
                        sf = dd.shotfile( Diagnostic , Shotnumber,Experiment, Edition)
                        
		        self.time = sf( coilStr[0,0], tBegin=tBegin, tEnd=tEnd ).time
		        self.ntimes = np.size(self.time)
		        self.samplingrate = 1.0/np.mean(np.diff(self.time))
                    
				#self.I_DC1 = sf( 'I_DC1', tBegin=tBegin, tEnd=tEnd ).data 
				#self.I_DC2 = sf( 'I_DC2', tBegin=tBegin, tEnd=tEnd ).data
				
		        self.Icoils = np.zeros((self.ntimes,nRows,nCoils)) 
                        for i in range(nRows):
                            for j in range(nCoils):
                                try:
				    self.Icoils[:,i,j] = sf( coilStr[i,j], tBegin=tBegin, tEnd=tEnd ).data
                                except:
                                    self.Icoils[:,i,j] = 0.0

                        if  ((Shotnumber>35000) & (Shotnumber<36793)):
                            self.Icoils[:,1,4] = -self.Icoils[:,1,4]
                                                         
		        self.nRows = nRows
		        self.nCoils = nCoils
		        self.Status = True
                    except:
			print('Error in loading SSV')
                        IPython.embed()
		try:
		    self.PSLresponse(iir=iir)
		except:
		    print 'Error in calculating PSL reponse' 			

		if np.all(Binning) != None:
		    try:
			self.Binning(Binning)
		    except:
			print 'Error in Binning' 	

		try:
	       	    self.torAngle(thres=thres )
       	        except:
	       	    print 'Error in calculating toroidal angle' 

       		try:
		    self.torAngle(usePSL = True, thres=thres )
		except:
		    print 'Error in calculating toroidal angle with PSL' 

                        ###other important 
                try:
                    self.rotVelUp = np.mean(np.gradient(np.unwrap(self.meanN*self.phaseUp[int(self.meanN-1)]))/np.gradient(self.MPtime))/2./np.pi
                    self.rotVelLo = np.mean(np.gradient(np.unwrap(self.meanN*self.phaseLo[int(self.meanN-1)]))/np.gradient(self.MPtime))/2./np.pi
                except:
                    print 'no rotation'

	def Unload( self ):
		if self.Status:
			self.Status = False
			del self.time	

			del self.Icoils
			del self.nRows 
			del self.nCoils 
			if self.Status_PSL:
				del self.PSLtime
				del self.PSLIcoils





### bin data
	def Binning( self, samplefreq = 0.1 ):				   
		if self.Status:				    						   				
			
			print "MAW binning"

			self.time,self.Icoils = dataBinning(self.time,self.Icoils,samplefreq=samplefreq )
			if self.Status_PSL:
				self.PSLtime,self.PSLIcoils = dataBinning(self.PSLtime,self.PSLIcoils,samplefreq=samplefreq )
			
			self.ntimes = np.size(self.time)
			self.samplingrate = 1.0/np.mean(np.diff(self.time))


# return coil current for one timepoint
	def getIcoils(self, timepoint , usePSL = False, freq=None):

		if self.Status:
			if usePSL & self.Status_PSL :
                            idx = np.argmin( np.abs( self.PSLtime - timepoint) )
                            return np.squeeze(self.PSLIcoils[idx,:,:])
			elif np.all(freq) != None:
                            idx = np.argmin( np.abs( self.time - timepoint) )
                            Icoils = np.copy(np.squeeze(self.Icoils[idx,:,:]))
                            Icoils[0] = Icoils[0]*self.getUpperResponseAmp(freq)
                            Icoils[1] = Icoils[1]*self.getLowerResponseAmp(freq)
                            return Icoils
                        else:
                            idx = np.argmin( np.abs( self.time - timepoint) )
                            return np.squeeze(self.Icoils[idx,:,:])

### function to get toroidal number,phase, amplitude of the upper and lower coild set
	def torAngle( self, usePSL=False,useLSQ = False, order=8,thres = 100. ):
		
		if usePSL & self.Status_PSL:
			usePSL = True
			self.Status_torangle_PSL = False
		else:
			usePSL =False
			self.Status_torangle = False
		


		if self.Status:
			
			MP=MPcoils()
			MP.LoadCoord()
		
			if usePSL:
				Icoils = self.PSLIcoils
			else:
				Icoils = self.Icoils
	
			Icoilmean=(np.abs(Icoils[:,:,:]).mean(axis=2)).mean(axis=1)
			idxUse = np.where(Icoilmean >= thres)[0]
			
			if usePSL:
				MPtime = self.PSLtime[idxUse]
			else:
				MPtime = self.time[idxUse]

			
			#IPython.embed()
			#torangle = np.unwrap( MP.coordPol[:,:,:,2].mean(axis=2) )
			# indices 2 (upper/lower) 
			if useLSQ:
				nPoCoils = 11 
				MPdata = []
				minCoord = MP.coordPol[:,:,:,2].min(axis=2).ravel()
				maxCoord = MP.coordPol[:,:,:,2].max(axis=2).ravel()
				for i in np.arange(minCoord.size):
					MPdata.append(np.linspace(minCoord[i],maxCoord[i],nPoCoils,endpoint=True))


				torAngle =  np.reshape(np.array(MPdata),(2,nPoCoils*8))
				IcoilsIn = np.swapaxes(np.reshape(np.swapaxes(np.array([Icoils[idxUse,:,:].T,]*nPoCoils).T,0,1),(2,idxUse.size,torAngle[0].size)),1,2)

				idxargSort=np.argsort(torAngle[0])
				torAngleFit = torAngle[:,idxargSort]
				IcoilsFit = IcoilsIn[:,idxargSort] 
				import LSQFFT

				LSup = LSQFFT.LSQFFT()
				LSup.initializeData(torAngleFit[0],IcoilsFit[0], freq_in = 1/(2.*np.pi),order=8,poly=0,negSign=True)
				LSup.leastSquareFit()

				LSlo = LSQFFT.LSQFFT()
				LSlo.initializeData(torAngleFit[1],IcoilsFit[1], freq_in = 1/(2.*np.pi),order=8,poly=0,negSign=True)
				LSlo.leastSquareFit()

				nn = np.arange(8)+1

				nUp = np.argmax(LSup.amp,axis=0)+1.
				nLo = np.argmax(LSlo.amp,axis=0)+1.
				n = (nUp + nLo)/2.
			
				meanN= n.mean()
			
				ampUp = LSup.amp
				ampLo = LSlo.amp
				amp = (ampUp+ampLo)/2.
				
				phaseUp = (LSup.phase.T/nn).T
				phaseLo = (LSlo.phase.T/nn).T
			else:
				#indices : upper or lower, Coil,torangel 
				torCoord =  MP.coordPol[:,:,:,2]
				torCoord[torCoord>np.pi] -= 2.*np.pi
				torCoord[torCoord<-np.pi] += 2.*np.pi
				#indices : upper orlower, time,Coil,torangel 
				#IcoilsUse = np.swapaxes(np.swapaxes(np.array([Icoils[idxUse,:,:],]),0,2),2,3)
			
				torEqu = np.linspace(-np.pi,np.pi,512,endpoint=False)
				IcoilsEqu = np.zeros((2,8,512))
				
				for ul in np.arange(2):
					for c in np.arange(8):
					##sort argument after toroiddal angle
						idxArgSort = np.argsort(torCoord[ul,c])
						torValues = torCoord[ul,c][idxArgSort]
						onesKilo = np.ones_like(idxArgSort)
						#ssplit interpolation for the coil which is around _pi and +pi
						if np.abs(torCoord[ul,c].max()-torCoord[ul,c].min()) < np.pi:
							IcoilsEqu[ul,c,:] = np.interp(torEqu,torValues,onesKilo,right=0.0,left=0.0)
						else:
							#separate neg and pos values
							idxNeg = torValues<0.0
							idxPos = torValues>=0.0
							IcoilsEqu[ul,c,torEqu<0.0] = np.interp(torEqu[torEqu<0.0],torValues[idxNeg],onesKilo[idxNeg],right=0.0)
							IcoilsEqu[ul,c,torEqu>=0.0] = np.interp(torEqu[torEqu>=0.0],torValues[idxPos],onesKilo[idxPos],left=0.0)
				

				
				#Multiply the time dependent current on the normilzed one
				IcoilFit = np.einsum('kij,ijl->ikl', Icoils[idxUse],IcoilsEqu)
				IcoilFFT = np.fft.rfft(IcoilFit)
				AmpTmp = np.abs(IcoilFFT)/np.size(IcoilFFT,-1)
				# for positive signs as default for negative signs
				PhaseTmp = np.angle(IcoilFFT)*-1.
				
				nn = np.arange(8)+1

				nUp = np.argmax(AmpTmp[0],axis=-1).T
				nLo = np.argmax(AmpTmp[1],axis=-1).T 
				n = (nUp + nLo)/2.
			
				meanN= np.mean(n)
			
				ampUp = AmpTmp[0,:,1:order+1].T
				ampLo = AmpTmp[1,:,1:order+1].T
				amp = (ampUp+ampLo)/2.
				#through the n number to get geometrical phase
				phaseUp = (PhaseTmp[0,:,1:order+1]/nn).T
				phaseLo = (PhaseTmp[1,:,1:order+1]/nn).T



			
			dPhi = np.zeros_like(phaseUp)
			for k in np.arange(phaseUp[:,0].size):
				useN = float(k+1)
				# for each n the difference in phase is unwrapped
				tmpDphi=np.unwrap(phaseUp[k]*useN)/useN-np.unwrap(phaseLo[k]*useN)/useN
				tmpDphi = np.remainder(tmpDphi,-np.pi*2./useN)
				tmpDphi[tmpDphi<(-np.pi/useN)] += (2.*np.pi/useN)
				tmpDphi[tmpDphi>(np.pi/useN)] -= (2.*np.pi/useN)
				dPhi[k] = tmpDphi*useN 

			dPhase =np.zeros_like(dPhi[0])
			domN = np.zeros_like(dPhi[0])
			phase = np.zeros((n.size,2))
                        phase0 = np.zeros((n.size))
			#get the phase of the dominant toroidla number
			for i in np.arange(dPhi[0].size):
				dPhase[i] = dPhi[np.array(n[i]-1.,dtype='int'),i]
				phase[i,0] = phaseUp[np.array(nUp[i]-1.,dtype='int'),i] 
				phase[i,1] = phaseLo[np.array(nLo[i]-1.,dtype='int'),i] 
		       
				#IPython.embed()

			fphi=[interp1d(MPtime,np.unwrap(phase[:,0]*useN)/useN,assume_sorted=False,bounds_error=False,fill_value =float('nan')),interp1d(MPtime,phase[:,1],assume_sorted=False,bounds_error=False,fill_value =float('nan'))]

			if usePSL:
				##phases of the dominant currents
				self.PSLMPtime = MPtime
				self.PSLphase,self.PSLdPhase=phase,dPhase							
				self.PSLdPhi,self.PSLphaseUp,self.PSLphaseLo=dPhi,phaseUp,phaseLo
				self.PSLphaseDeg,self.PSLdPhaseDeg=np.rad2deg(phase),np.rad2deg(dPhase)
				self.PSLdPhiDeg,self.PSLphaseUpDeg,self.PSLphaseLoDeg=np.rad2deg(dPhi),np.rad2deg(phaseUp),np.rad2deg(phaseLo)
				self.PSLampUp,self.PSLampLo,self.PSLamp=ampUp,ampLo,amp 
				self.PSLnUp,self.PSLnLo,self.PSLn,self.PSLmeanN=nUp,nLo,n,meanN 
				self.PSLfphi = fphi
				self.Status_torangle_PSL = True
	      		else:
								##phases of the dominant currents
				self.MPtime = MPtime
				self.phase,self.dPhase=phase,dPhase
				self.phase0 =  (PhaseTmp[1,:,1:order+1]).T + dPhi/2.
				self.dPhi,self.phaseUp,self.phaseLo=dPhi,phaseUp,phaseLo
				self.phaseDeg,self.dPhaseDeg=np.rad2deg(phase),np.rad2deg(dPhase)
				self.dPhiDeg,self.phaseUpDeg,self.phaseLoDeg=np.rad2deg(dPhi),np.rad2deg(phaseUp),np.rad2deg(phaseLo)
				self.ampUp,self.ampLo,self.amp=ampUp,ampLo,amp 
				self.nUp,self.nLo,self.n,self.meanN=nUp,nLo,n,meanN
				self.fphi = fphi
				self.Status_torangle = True

			self.nn = nn

	def PSLresponse(self,iir=False):
		if self.Status:

                    if ((iir == True) & (np.round(self.samplingrate)==10000.)):
                        print("use iir coefficients")
                        #coefficients from scilab for 10kHy
                        ck10_upper =  [-2.056735e-07,  2.494287e-10 , 2.059229e-07 ][::-1]
                        dj10_upper =  [ 0.864943 , -1.864864 , 1.000000 ][::-1]
                        ck10_lower =  [-1.199285e-07,  1.176715e-10,  1.200461e-07][::-1]
                        dj10_lower =  [0.896064 , -1.896004 , 1.000000 ][::-1]
 
                        
                        self.PSLIcoils = np.zeros_like( self.Icoils)
                        self.PSLtime = self.time
                        
                        for i in np.arange(self.Icoils[0,:,0].size):
                            for j in np.arange(self.Icoils[0,0,:].size):
                                if i==0:
                                    ##additional coefficient since the magentic fiedl was calculated
                                    self.PSLIcoils[:,i,j] =  signal.lfilter(ck10_upper, dj10_upper, self.Icoils[:,i,j])*1000*158.41297261
                                else:
                                    self.PSLIcoils[:,i,j] =  signal.lfilter(ck10_lower, dj10_lower, self.Icoils[:,i,j])*1000*255.03395609

                        self.Status_PSL=True           
                    else:
                        print("use fft for PSL")
			##see at print("use iir coefficients")
                       # the end of the file, from http://www.aug.ipp.mpg.de/augtwiki/bin/view/CWAC/PslResponse			
			#Column 1: Frequency [Hz]
			#Column 2: Magnetic induction [T], in phase with applied coil current (In-phase or 0deg component)
			#Column 3: Magnetic induction [T], out of phase with applied coil current (quadrature or 90deg component)
			upperData = np.array(PSLupper).T
			freq = upperData[0]
			upperComplex =( upperData[1]+1j*upperData[2])/upperData[1,0]
			upperAmp,upperPhase = np.abs(upperComplex)/np.abs(upperComplex).max(),np.angle(upperComplex)-np.angle(upperComplex).max()
			lowerData = np.array(PSLlower).T
			lowerComplex = (lowerData[1]+1j*lowerData[2])/lowerData[1,0]
			lowerAmp,lowerPhase = np.abs(lowerComplex)/np.abs(lowerComplex).max(),np.angle(lowerComplex)-np.angle(lowerComplex).max()

			#is done to account for truncation
			newTime = np.fft.irfft(np.fft.rfft(self.time))

			#reserve Memory and account for truncation
			self.PSLIcoils = np.zeros_like( self.Icoils)[:newTime.size]

			#get the frequency of the this reponse
			freqIcoil = np.fft.rfftfreq(newTime.size,np.diff(newTime).mean())
			#get response function of including PSL
			upperResp = (np.interp(freqIcoil,freq,(np.real(upperComplex))))+1j*(np.interp(freqIcoil,freq,(np.imag(upperComplex))))
			lowerResp = (np.interp(freqIcoil,freq,(np.real(lowerComplex))))+1j*(np.interp(freqIcoil,freq,(np.imag(lowerComplex))))

                        #IPython.embed()
                        fupper=interp1d(freq,(upperComplex),fill_value='extrapolate',bounds_error=False,kind='cubic')
                        flower=interp1d(freq,(lowerComplex),fill_value='extrapolate',bounds_error=False,kind='cubic')

                        self.funcUpperAmp = interp1d(freqIcoil,np.abs(upperResp))
                        self.funcLowerAmp = interp1d(freqIcoil,np.abs(lowerResp))
                        self.funcUpperPhase = interp1d(freqIcoil,np.angle(upperResp))
                        self.funcLowerPhase = interp1d(freqIcoil,np.angle(lowerResp))

                        self.maxFreqCoils = np.zeros_like(self.Icoils[0])

                        for i in np.arange(self.Icoils[0,:,0].size):
                            for j in np.arange(self.Icoils[0,0,:].size):
                                FFTIn =np.fft.rfft(self.Icoils[:,i,j])
                                self.maxFreqCoils[i,j] = freqIcoil[np.abs(FFTIn).argmax()]
                                if i==0:
                                    self.PSLIcoils[:,i,j] = np.fft.irfft(FFTIn*fupper(freqIcoil))
                                else:
                                    self.PSLIcoils[:,i,j] = np.fft.irfft(FFTIn*flower(freqIcoil))

                        self.maxFreq = self.maxFreqCoils.mean()
                        
                        self.PSLtime = newTime
			self.Status_PSL=True

        ## in Hz
        def getUpperResponseAmp(self,freqIn=1.0):         
            if self.Status_PSL:
                return self.funcUpperAmp(freqIn)

        def getLowerResponseAmp(self,freqIn=1.0):         
            if self.Status_PSL:
                return self.funcLowerAmp(freqIn)


        def getPhi0(self, usePSL=True,onlyPos=False):
            if self.Status_torangle :
                if usePSL:
                    fphi = self.PSLfphi
                else:
                    fphi = self.fphi   
                    #IPython.embed()

                useN = np.round(np.mean(self.n))
                if usePSL:
                    #use upper row idx= 0, 0.615 offset and multiply by n to get 
                    Phi0= (self.PSLphase[:,0] - self.PSLdPhi[int(useN-1.)])*useN+np.deg2rad(139.25)
                else:
                    Phi0= (self.phase[:,0] - self.dPhi[int(useN-1.)])*useN+np.deg2rad(139.3)

                Phi0=np.remainder(Phi0,2.*np.pi)    
               # IPython.embed()
      #          synPhase= np.remainder(( fphi[idxCoil](tref) - (fphi[idxCoil](timeIn)) ) + np.remainder(phiDiag,2.*np.pi/np.mean(self.n) )  ,2.*np.pi/np.mean(self.n) )
                #synPhase= np.remainder( (fphi[idxCoil](tref) - fphi[idxCoil](timeIn)),2.*np.pi/np.mean(self.n)) + np.remainder(phiDiag,2.*np.pi/np.mean(self.n) )  

                #IPython.embed()
                if onlyPos:
                    Phi0[Phi0<0.0]+=(2*np.pi)
                    Phi0[Phi0>(2.*np.pi)] -= (2*np.pi)
                else:
                    Phi0[Phi0>(np.pi)]-=(2*np.pi)
                    Phi0[Phi0<(-np.pi)]+=(2*np.pi)   
                    #synPhase[synPhase<(-np.pi/np.mean(self.n))]+=(2*np.pi/np.mean(self.n))
                        
                return Phi0

            
### Transform phi values to a timebase 
        def getSynTime(self,phiInOrg,tRange, tref=2.5,phiDiag=0.0, idxCoil=0,onlyPos=False, usePSL=True,signalIn=None):
            if self.Status_torangle :
                freq = np.round(self.maxFreq)
                timeIn=np.linspace(tRange[0],tRange[0]+1./freq,1./freq*1000000)
                spacePhi = self.getSynPhase( timeIn, tref=tref,phiDiag=phiDiag, idxCoil=idxCoil,onlyPos=onlyPos, usePSL=usePSL)
                idxSort=np.argsort(spacePhi)
                phiIn = np.copy(phiInOrg)
                if onlyPos:
                    phiIn[phiIn<0.0]+=(2*np.pi/np.mean(self.n))
                    phiIn[phiIn>(2.*np.pi/np.mean(self.n))] -= (2*np.pi/np.mean(self.n))
                else:
                    phiIn[phiIn>(np.pi/np.mean(self.n))]-=(2*np.pi/np.mean(self.n))
                    phiIn[phiIn<(-np.pi/np.mean(self.n))]+=(2*np.pi/np.mean(self.n)) 

                #timeOrg = np.interp(phiIn,spacePhi[idxSort],timeIn[idxSort])
                timeOrg = np.zeros_like(phiIn)
                
                for i in np.arange(phiIn.size):
                    id=np.argmin(np.abs(phiIn[i]-spacePhi))
                    timeOrg[i] = timeIn[id]

                numRepl=np.ceil((tRange[1]-tRange[0])*freq)
                timeTmp=np.array([])

                for i in np.arange(0, numRepl):
                    timeTmp=np.append(timeTmp,timeOrg+np.float(i)/freq)
               
                useSignal=False
                   
                if np.all(signalIn) != None:
                    if np.size(phiInOrg)==np.shape(signalIn)[0]:
                        origShape = np.shape(signalIn[0])
                        origSize = np.size(signalIn[0])
                        firstDim = np.shape(signalIn)[0]
                        signalTmpX = np.tile( np.reshape(signalIn,(firstDim,origSize)).T,np.int(numRepl) ).T
                        signalTmp = np.reshape(signalTmpX,np.append(np.shape(signalTmpX)[0],origShape))
                        useSignal=True
                    else:
                        useSignal=False
                        print 'not the same size'
                #IPython.embed()

                idxUse = np.where( (timeTmp>=tRange[0]) & (timeTmp<=tRange[1] ) )[0]
                idxSortUse = idxUse[np.argsort(timeTmp[idxUse])]

                if  useSignal:
                    return timeTmp[idxSortUse],signalTmp[idxSortUse]
                else:
                    return timeTmp[idxSortUse]

              
#function to map timebase of diagnostic on geometrical phi
        def getSynPhase(self,timeIn, tref=2.5,phiDiag=0.0, idxCoil=0,onlyPos=False, usePSL=True):
        #IPython.embed()
            if self.Status_torangle :
                if usePSL:
                    fphi = self.PSLfphi
                else:
                    fphi = self.fphi   
                    #IPython.embed()
                
                if (phiDiag > np.pi) | (phiDiag < -np.pi):
                    print 'warning input phi is bigger than pi, use radians'
                
                synPhase= np.remainder(( fphi[idxCoil](tref) - (fphi[idxCoil](timeIn)) ) + np.remainder(phiDiag,2.*np.pi/np.mean(self.n) )  ,2.*np.pi/np.mean(self.n) )
                #synPhase= np.remainder( (fphi[idxCoil](tref) - fphi[idxCoil](timeIn)),2.*np.pi/np.mean(self.n)) + np.remainder(phiDiag,2.*np.pi/np.mean(self.n) )  

                #IPython.embed()
                if onlyPos:
                    synPhase[synPhase<0.0]+=(2*np.pi/np.mean(self.n))
                    synPhase[synPhase>(2.*np.pi/np.mean(self.n))] -= (2*np.pi/np.mean(self.n))
                else:
                    synPhase[synPhase>(np.pi/np.mean(self.n))]-=(2*np.pi/np.mean(self.n))
                    synPhase[synPhase<(-np.pi/np.mean(self.n))]+=(2*np.pi/np.mean(self.n))   
                    #synPhase[synPhase<(-np.pi/np.mean(self.n))]+=(2*np.pi/np.mean(self.n))
                        
                return synPhase





# return new timebase and Data
def dataBinning( time, data, samplefreq = 1.0 ):			
            						      
    print "binning with ",samplefreq," kHz"
    ntimes= np.size(time)
    samplingrate = 1.0/np.mean(np.diff(time))
    dataShape =np.array(np.shape(data))  
    #get the time index
    idxOfTime = np.squeeze(np.where(dataShape == ntimes))
    # if more index with the number of times exists, take the first one
    if np.size(idxOfTime) > 1:
        idxOfTime = idxOfTime[0]

    bins = int(ntimes*(float(samplefreq)*1.0e3/samplingrate))

    slices= np.linspace(0, ntimes, bins+1, True).astype(int)
    counts = np.diff(slices)

    #calculate new timebase
    newTime = np.add.reduceat(time, slices[:-1]) / counts
    newNtimes = np.size(newTime)

    #create new shape
    newDataShape = dataShape
    #replace old shape
    np.put(newDataShape, idxOfTime, newNtimes)
    #create new Data array
    newData = np.zeros( (newDataShape)  )

    #simplify array such as the first index is always the timebase
    newData = np.swapaxes(newData,0,idxOfTime)
    data = np.swapaxes( data,0,idxOfTime )

    storeShape = np.shape( newData )

    # rehape data to two dimensions
    data = np.reshape(data,(ntimes,np.size(data)/ntimes))
    newData = np.reshape(newData,(newNtimes,np.size(newData)/newNtimes))

    for i in np.arange(np.shape(data)[1]):
        newData[:,i] = np.add.reduceat(data[:,i], slices[:-1]) / counts

#shape back
    newData = np.reshape(newData,(storeShape))
    #swap back to original shape
    newData = np.swapaxes(newData,0,idxOfTime)

    return newTime,newData



### dataSet to get PSL response
PSLupper=[[0.0	,-0.00639344	,0],\
[0.01	,-0.00639276	,4.14117e-05],\
[0.02	,-0.00639072	,8.26923e-05],\
[0.05	,-0.00637706	,0.00020464],\
[0.1	,-0.00633366	,0.000398095],\
[0.2	,-0.00619141	,0.000744266],\
[0.5	,-0.00552656	,0.0014111],\
[1.0	,-0.0045564	,0.00161289],\
[2.0	,-0.00372656	,0.00137053],\
[3.0	,-0.00339467	,0.00119741],\
[5.0	,-0.00305732	,0.00102124],\
[6.0	,-0.00294986	,0.000965494],\
[10.0	,-0.00268725	,0.000823328],\
[20.0	,-0.002412	,0.000686552],\
[27.0	,-0.00231244	,0.000652528],\
[50.0	,-0.00212628	,0.000639539],\
[100.0	,-0.00189738	,0.000732157],\
[200.0	,-0.00153993	,0.000917016],\
[500.0	,-0.000776298	,0.000989645],\
[1000.0	,-0.000270289	,0.000735245],\
[2000.0	,-3.11631e-05	,0.000415854],\
[5000.0	,3.10123e-05	,0.00015078],\
[10000	,2.53022e-05	,6.55385e-05]]

PSLlower=[[0.00	,-0.00396941	,0],
[0.01	,-0.00396904	,2.09959e-05],
[0.02	,-0.00396794	,4.19146e-05],
[0.05	,-0.00396059	,0.000103562],
[0.1	,-0.00393763	,0.000200686],
[0.2	,-0.00386431	,0.000372881],
[0.5	,-0.00352834	,0.000700021],
[1.0	,-0.0030459	,0.000787654],
[2.0	,-0.00264888	,0.000654385],
[3.0	,-0.00249793	,0.00056721],
[5.0	,-0.00234948	,0.000486902],
[6.0	,-0.00230278	,0.000464675],
[10.0	,-0.00218773	,0.000418504],
[20.0	,-0.00205949	,0.000407054],
[27.0	,-0.00200746	,0.000425879],
[50.0	,-0.00188569	,0.000524306],
[100.0	,-0.00165041	,0.000742031],
[200.0	,-0.00117445	,0.000965846],
[500.0	,-0.000384786	,0.000809203],
[1000.0	,-8.44629e-05	,0.000482891],
[2000.0	,8.54296e-06	,0.000243242],
[5000.0	,2.06842e-05	,8.51125e-05],
[10000	,1.483e-05	,3.74586e-05]]



