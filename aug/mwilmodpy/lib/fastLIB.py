import numpy
import dd 
import kk 
import scipy.interpolate
import IPython
import scipy
import matplotlib.pylab as plt

class fastLIB:
	def __init__( self ,  Experiment = 'AUGD', Diagnostic = 'LIZ', Shotnumber = None ):
	
#		self.Status_LIN = False
		self.Status_middle = False
		self.Status_upper = False
		self.Status_Rz = False
		self.Status_cut = False
		self.Status = False
		if Shotnumber != None :
			self.Load( Shotnumber )
		

	def __del__( self ):
		self.Unload( )
		del self.Status

	def Load( self , Shotnumber, Experiment = 'AUGD', Edition = 0L, tBegin=-1.0, tEnd=12.0, Optics='middle'):
		self.Unload()		
	

		try:
			sf_LIZ = dd.shotfile( 'LIZ' , Shotnumber, Experiment, Edition)       
		except:
			print "Error reading LIZ shotfile" 
			return
		

#		reading chopping signal
		self.ChoppingSignal  = sf_LIZ('ch_sig', tBegin=tBegin, tEnd=tEnd,dtype=numpy.float64).data
		self.time  = sf_LIZ('ch_sig', tBegin=tBegin, tEnd=tEnd,dtype=numpy.float64).time

		if Optics == 'middle':
			try:
				self.SIG  = sf_LIZ('MSIG', tBegin=tBegin, tEnd=tEnd, dtype=numpy.float64).data
				print "Middle signal loaded" 
				self.SIG_nchannels = numpy.size(self.SIG[0,:])
				self.Status_middle = True
			except:
				print "Error reading signals from middle optics" 

		if Optics == 'upper':
			try:
				self.SIG  = sf_LIZ('SIG', tBegin=tBegin, tEnd=tEnd, dtype=numpy.float64).data
				self.SIG_nchannels = numpy.size(self.SIG[0,:])
				self.Status_upper = True
			except:
				print "Error reading signals from upper optics" 

			
		self.ntimes = numpy.size(self.time)
		self.samplingrate = 1.0/numpy.mean(numpy.diff(self.time))
		self.Shotnumber = Shotnumber
		self.Status = True
		

	def Unload( self ):
		
		
		if self.Status:

			del self.ChoppingSignal

			if (self.Status_middle | self.Status_upper) :
				del self.SIG
				del self.SIG_nchannels

			if self.Status_middle: self.Status_middle = False
			if self.Status_upper: self.Status_upper = False
			if self.Status_cut:  self.Status_cut = False
			if self.Status_Rz:  self.Status_Rz = False
 
			del self.ntimes 
			del self.samplingrate
			self.Status = False



	def MapToRhop( self , timepoints, Experiment = 'AUGD', Diagnostic='EQH' ):
		ntimes = 0
		if self.Status:
			ntimes = numpy.size(timepoints)			
			rhop = numpy.zeros((ntimes,len(self.R[0,:])))
			try:
				for i in range(ntimes):
					
					output = kk.KK().kkrzptfn( self.Shotnumber, timepoints[i], self.R, self.z, exp= Experiment, diag=Diagnostic, ed=0)
					rhop[i,:] = output.rho_p
			except Exception:
				output = kk.KK().kkrzptfn( self.Shotnumber, timepoints, self.R, self.z, exp= Experiment, diag=Diagnostic, ed=0)
				pass
			
			if ntimes > 1: 	
				return rhop 
			else: 
				return output.rho_p

	

	def GetRz( self,Experiment = 'AUGD', Shotnumber=0 ):
		if Shotnumber == 0: Shotnumber = self.Shotnumber
		if self.Status_upper:
			try:
				sf_LIN = dd.shotfile( 'LIN' , Shotnumber, Experiment, Edition)
				ndata = numpy.size(sf_LIN('R').data)
				if ndata == self.SIG_nchannels:
					self.R = sf_LIN('R').data
					self.z = sf_LIN('z').data
					self.Status_Rz = True
				sf_LIN.close()
			except:
				print "Error reading Rz position from LIN shotfile" 
				return

		if self.Status_middle:

			if (Shotnumber < 29543):

				self.R = [2.19099,2.179097,2.1741416,2.167204,2.1612575,2.155311,2.1491663,2.1369759, 2.1305339,2.1240918,2.1179471,2.111505,2.1052612,2.0988191,2.0923771,2.0861333,2.0799885,2.0735465,2.0666089,2.0596713,2.0524364,2.0457961,2.0385612,2.0312272,2.0236949,2.0165591]	
				self.z = [320.27,319.99,319.85,319.71,319.57,319.42,319.28,318.98,318.83,318.68,318.53,318.37,318.22,318.06,317.9,317.74,317.58,317.42,317.26,317.09,316.93,316.76,316.59,316.42,316.25,316.07]
				self.Status_Rz = True

			else:

				self.R = [2.202883,2.190990,2.186035,2.179097,2.173151,2.167204,2.161059,2.148869,2.142427,2.135985,2.129840,2.123398,2.117154,2.110712,2.104270,2.09808,2.091882,2.085440,2.078502,2.071564,2.064329,2.057689,2.050454,2.043120,2.035588,2.028452]
				self.z = [0.32027,0.319993,0.319853,0.319711,0.319568,0.319423,0.319277,0.318981,0.318831,0.318680,0.318527,0.318373,0.318218,0.318061,0.317903,0.317743,0.317582,0.317420,0.317256,0.317091,0.316925,0.316757,0.316588,0.316418,0.316073]
				self.Status_Rz = True

			
	def CutBeamOff( self  ):
		
		if self.Status:
			meanChopping = numpy.mean(self.ChoppingSignal)
			indexOn = numpy.where(self.ChoppingSignal < meanChopping)[0]
			ruleOut = numpy.zeros(self.ntimes)
			ruleOut[indexOn] = 1.0
			edges = numpy.diff(self.time[indexOn])
			indexEdges = numpy.add(numpy.where(edges > numpy.mean(edges))[0],1)
			ruleOut[indexEdges] = 0.0
		#	IPython.embed()
			ruleOut[indexOn[numpy.add(indexEdges,1)]]=0.0
			ruleOut[indexOn[indexEdges]]=0.0
			ruleOut[indexOn[[0,-1]]]=0.0
			ruleOut[indexOn[numpy.add(indexEdges,-1)]]=0.0
			indexOn = numpy.where(ruleOut == 1.0)[0]
		#	IPython.embed()
			self.ChoppingSignal = self.ChoppingSignal[indexOn]
			self.time = self.time[indexOn]
			self.ntimes = numpy.size(self.time) 
			if (self.Status_upper | self.Status_middle) :
				self.SIG = numpy.reshape(self.SIG[indexOn,:],(self.ntimes,self.SIG_nchannels))
				self.Status_cut = True
		#	IPython.embed()
			  		
	def Binning( self, samplefreq = 8.0 ):	
			   
			if self.Status:	
				if samplefreq < self.samplingrate: 
					if (self.Status_cut) & (self.Status_upper | self.Status_middle):
						print "LIZ binning with ",samplefreq," kHz"	
						edges = numpy.diff(self.time)
						indexEdges = numpy.where(edges > numpy.mean(edges))[0]
						nEdges = numpy.size(indexEdges)
						newtime = []
						newSIG = []
						for i in range(nEdges-2):
							ntimes = indexEdges[i+1] - (indexEdges[i]+1)
							bins = int(ntimes*(float(samplefreq)*1.0e3/self.samplingrate))

							slices= numpy.linspace(0, ntimes, bins+1, True).astype(numpy.int)
							counts = numpy.diff(slices)
							
							tmptime =  numpy.add.reduceat(self.time[(indexEdges[i]+1):indexEdges[i+1]], slices[:-1]) / counts
							tmpSIG = numpy.zeros((numpy.size(tmptime),self.SIG_nchannels))

							for j in range(self.SIG_nchannels):
								tmpSIG[:,j] =  numpy.add.reduceat(self.SIG[(indexEdges[i]+1):indexEdges[i+1],j], slices[:-1]) / counts

							if i==0:
								newtime = tmptime
								newSIG = tmpSIG
							else:
								newtime = numpy.append( newtime, tmptime )
								newSIG = numpy.append( newSIG,tmpSIG, axis=0 )
						
						self.SIG = newSIG
						self.time = newtime
						
					else:
						print "LIZ binning with ",samplefreq," kHz"
						bins = int(self.ntimes*(float(samplefreq)*1.0e3/self.samplingrate))
						
						if (self.Status_upper | self.Status_middle) :
							tempSIG = self.SIG
	
							slices= numpy.linspace(0, self.ntimes, bins+1, True).astype(numpy.int)
							counts = numpy.diff(slices)

							self.time = numpy.add.reduceat(self.time, slices[:-1]) / counts
							self.ChoppingSignal = numpy.add.reduceat(self.ChoppingSignal, slices[:-1]) / counts
							self.ntimes = numpy.size(self.time)
							self.samplingrate = 1.0/numpy.mean(numpy.diff(self.time))

							self.SIG = numpy.zeros((self.ntimes,self.SIG_nchannels))
							for i in range(self.SIG_nchannels):
								self.SIG[:,i] = numpy.add.reduceat(tempSIG[:,i], slices[:-1]) / counts

						else:
							print "No Lithium data read"
				else:
					print "Binning frequency is higher than sampling frequency, does not make sense to bin"
				
	def TransformTimebase( self, time_in, signal_in ):
		
		if self.Status:
			signal_out = numpy.zeros(self.ntimes)
			f=scipy.interpolate.interp1d( time_in, signal_in, bounds_error=False, fill_value=0 )
			return f(self.time)
		else:
			return 0
	

		