import numpy as np
from scipy.optimize import curve_fit, leastsq, minimize,brute,fmin_slsqp
from scipy.interpolate import interp1d
#from  repository
#setenv PYTHONPATH "$PYTHONPATH":/afs/ipp-garching.mpg.de/aug/ads-diags/common/python/lib
import dd

#from Matthias
import tools
import MPcoils
import leastsq_bounds
#for testing purpose
import IPython
import matplotlib.pylab as plt


## Test differential n=2
# import matplotlib.pylab as plt
# import MAW
# import numpy as np
# MW = MAW.MAW() 
# MW.Load(30826, tBegin=1.3, tEnd=3.2)
# MW.Binning(samplefreq = 1.0)
# MW.torAngle(mode = 'differential')
# plt.plot(MW.MPtime,np.rad2deg(MW.dPhase))

## test Rigid, n=2
# import MAW
# MW = MAW.MAW() 
# MW.Load(30839, tBegin=2.1, tEnd=4.0)
# MW.Binning(samplefreq = 1.0)
# MW.torAngle(mode = 'rigid')
# plt.plot(MW.MPtime,np.rad2deg(MW.phase))
# print MW.n

#other parameters
###  MW.time...timebase, MW.I_DC1...Current supply, MW.Icoils...currents of every coil

## index null sind die oberen Spulen
class MAW:
	def __init__( self , Shotnumber = None, Experiment = 'AUGD', Diagnostic = 'MAW' ):
		self.Status = False
		self.Status_torangle = False
		if Shotnumber != None :
			self.Load( Shotnumber )
		

	def __del__( self ):
		self.Unload( )
		del self.Status

##load data
	def Load( self, Shotnumber ,Experiment= 'AUGD', Diagnostic= 'MAW', tBegin=-1.0, tEnd=16.0, Edition = 0L, nRows = 2, nCoils = 8, labelRows = ['u','l'] ):

		self.Unload()
		
		if (Diagnostic == 'MAW') & (Shotnumber > 0) :

			##create String to read Signal
			coilStr = np.chararray((nRows, nCoils),itemsize=2)
			for i in range(nRows):
				for j in range(nCoils):
					coilStr[i,j] = labelRows[i]+str(j+1) 

			signalCoilStr = 'IB' +  coilStr 
			
			try:
				
				sf = dd.shotfile( Diagnostic , Shotnumber,Experiment, Edition)
				
				self.time = sf( 'I_DC1', tBegin=tBegin, tEnd=tEnd ).time
				self.ntimes = np.size(self.time)
				self.samplingrate = 1.0/np.mean(np.diff(self.time))

				self.I_DC1 = sf( 'I_DC1', tBegin=tBegin, tEnd=tEnd ).data 
				self.I_DC2 = sf( 'I_DC2', tBegin=tBegin, tEnd=tEnd ).data
				
				self.Icoils = np.zeros((self.ntimes,nRows,nCoils)) 
				for i in range(nRows):
					for j in range(nCoils):
						self.Icoils[:,i,j] = sf( signalCoilStr[i,j], tBegin=tBegin, tEnd=tEnd ).data
				self.nRows = nRows
				self.nCoils = nCoils
				self.Status = True
				sf.close()
			except:
				print 'Error in loading' ,Diagnostic


	def Unload( self ):
		if self.Status:
			self.Status = False
			del self.time	
			del self.I_DC1	
			del self.I_DC2
			del self.Icoils
			del self.nRows 
			del self.nCoils 
						

### bin data
	def Binning( self, samplefreq = 0.1 ):				   
		if self.Status:				    						   				
			
			print "MAW binning"

			tmptime,self.I_DC1 = tools.dataBinning(self.time,self.I_DC1,samplefreq=samplefreq )
			tmptime,self.I_DC2 = tools.dataBinning(self.time,self.I_DC2,samplefreq=samplefreq )
			self.time,self.Icoils = tools.dataBinning(self.time,self.Icoils,samplefreq=samplefreq )
			
			self.ntimes = np.size(self.time)
			self.samplingrate = 1.0/np.mean(np.diff(self.time))


# return coil current for one timepoint
	def getIcoils(self, timepoint ):

		if self.Status:
			idx = np.argmin( np.abs( self.time - timepoint) )
			return np.squeeze(self.Icoils[idx,:,:])

### function to get toroidal number,phase, amplitude of the upper and lower coild set
	def torAngle( self, mode ='rigid' ):

		self.Status_torangle = False


		if self.Status:
			MP=MPcoils.MPcoils()
			MP.LoadCoord()
			#torangle = np.unwrap( MP.coordPol[:,:,:,2].mean(axis=2) )
			torangle =  MP.coordPol[:,:,:,2].mean(axis=2) 
			#mode unkown, could be anything and nothing
		#	if mode == 'unknown':
			Icoilmean=(np.abs(self.Icoils[:,:,:]).mean(axis=2)).mean(axis=1)
			idxTm = np.squeeze(np.where( (np.abs(Icoilmean) > 100.) | (np.abs(Icoilmean) > 100.) ))
			torTime = self.time[idxTm]
			nIdxTime = np.size(idxTm)
			torphi = np.zeros((nIdxTime,2))
			nNumber = np.zeros((nIdxTime,2),dtype=int)
			dphi = np.zeros_like(torTime)
			amp =  np.zeros_like(torphi)
			residuals =  np.zeros_like(torphi)
			## Nothing is known about the configuration
			if mode == 'unknown':
				#IPython.embed()
				p0=[np.abs(self.I_DC1).max(),0.0,1.5]
				for t in np.arange(0,nIdxTime):
					for i in np.arange(0,2):
						xdata = torangle[i]
						ydata = self.Icoils[idxTm[t],i]
						#use last result as first gues
						if t!=0:
							p0 = [amp[t-1,i],torphi[t-1,i],nNumber[t-1,i]]
						params=leastsq(sinFuncResi, p0, args=(xdata, ydata))[0]
						amp[t,i]=params[0]
						torphi[t,i] =params[1]
						nNumber[t,i] = np.round(params[2])
						residuals[t,i] = np.linalg.norm(sinFuncResi(params,xdata, ydata))
					#params=[amp[t,i],torphi[t,i],nNumber[t,i]]
			elif (mode == 'rigid') | (mode == 'differential')| (mode == 'constn') :
				#get amplitude
				maxAmp = np.max(self.Icoils)
				minAmp = np.min(self.Icoils)
				amp[:,:]= (np.abs(maxAmp)+np.abs(minAmp))*0.5
				###get nNumber
				p0=[0.0,1.5]
				torphi_test =  np.zeros_like(torphi)
				for i in np.arange(0,2):
					xdata = torangle[i]
					for t in np.arange(0,nIdxTime):
						ydata = self.Icoils[idxTm[t],i]
						if t!=0:
							p0 = [torphi[t-1,i],nNumber[t-1,i]]
						params=leastsq(sinAmpResi, [0.0,1.5],args=(amp[t,i],xdata, ydata))[0]
						torphi_test[t,i] =params[0]
						nNumber[t,i] = np.round(params[1])
						residuals[t,i]=np.linalg.norm(sinAmpResi(params,amp[t,i],xdata, ydata))
				#get finally the n number use, residuals for weigthing
				nNumber_old = np.copy(nNumber)
				nNumber[:,:]=np.round(np.sum(nNumber/residuals)/np.sum(1./residuals))
				# calculate the torphi with known amplitude and n number
				for i in np.arange(0,2):
					xdata = torangle[i]
					for t in np.arange(0,nIdxTime):
						ydata = self.Icoils[idxTm[t],i]
						A=np.array([np.sin(nNumber[t,i]*xdata),-np.cos(nNumber[t,i]*xdata), np.ones((np.size(xdata)))]).transpose()
						paramLin = np.linalg.lstsq(A, ydata)[0]
						torphi[t,i]  = np.arctan2(paramLin[1],paramLin[0])/nNumber[t,i]
						amp[t,i]=(paramLin[0] + paramLin[1])/( np.sin(nNumber[t,i]*torphi[t,i]) + np.cos(nNumber[t,i]*torphi[t,i]) )


			self.meanN= nNumber.mean()
			
			### neg sign because of the fit
			dphi[:]=np.unwrap(torphi[:,0]*self.meanN)/self.meanN-np.unwrap(torphi[:,1]*self.meanN)/self.meanN
			tmpDphi = np.remainder(dphi,-np.pi*2./self.meanN)
			tmpDphi[tmpDphi<(-np.pi/self.meanN)] += (2.*np.pi/self.meanN)
			tmpDphi[tmpDphi>(np.pi/self.meanN)] -= (2.*np.pi/self.meanN)
			

			self.dPhase = tmpDphi*self.meanN
			

			self.phase = torphi
			self.n = nNumber
			self.meanN= nNumber.mean()
			self.amp = amp
			self.MPtime = torTime

			self.fphi=[interp1d(self.MPtime,self.phase[:,0],assume_sorted=False,bounds_error=False,fill_value =float('nan')),interp1d(self.MPtime,self.phase[:,1],assume_sorted=False,bounds_error=False,fill_value =float('nan'))]


			self.Status_torangle = True

#function to map timebase of diagnostic on geometrical phi
	def getSynPhase(self,timeIn, tref=2.5,phiDiag=0.0, idxCoil=0,onlyPos=False):
		if self.Status_torangle :
			#IPython.embed()
			if (phiDiag > np.pi) | (phiDiag < -np.pi):
				print 'warning input phi is bigger than pi, use radians'
			synPhase=( self.fphi[idxCoil](tref)-self.fphi[idxCoil](timeIn) ) + np.remainder(phiDiag,2.*np.pi/np.mean(self.n) )  
			if onlyPos:
				synPhase[synPhase<0.0]+=(2*np.pi/np.mean(self.n))
			
			synPhase[synPhase>(2*np.pi/np.mean(self.n))]-=(2*np.pi/np.mean(self.n))

			#IPython.embed()

			return synPhase



### function to fit the MP coil current
def sinFuncResi(  params, x , y ):
	a = params[0]
	tphi = params[1]
	n = params[2]
	if ((a<=100.0) | (a > 1300.)):
		return np.inf
	if ((tphi < -np.pi) | (tphi > np.pi ) ):  
		return np.inf
	if ((n < 0) | (n > 4 ) ):  
		return np.inf
	return a * np.sin(n*x+tphi)-y

def sinFunc(  params, x):
	a = params[0]
	tphi = params[1]
	n = params[2]
	return a * np.sin(n*x+tphi)

def sinAmpResi(  params, amp, x , y ):
	tphi = params[0]
	n = params[1]
	if ((tphi < -np.pi) | (tphi > np.pi ) ):  
		return np.inf
	if ((n < 1) | (n > 3 ) ):  
		return np.inf
	return amp * np.sin(n*x+tphi)-y
