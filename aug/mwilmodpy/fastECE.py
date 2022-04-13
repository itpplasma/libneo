##
# import fastECE
# ECE=fastECE.fastECE()
# ECE.Load(29161)
# ECE.Binning(20)
import numpy
import dd 
import kk 
import scipy.interpolate
import IPython
import scipy

class fastECE:
	def __init__( self ,  Experiment = 'AUGD', Diagnostic = 'RMD', Shotnumber = None ):
		self.Status_RMD = False
		self.Status_RMC = False
		self.Status_Magn = False
		self.Status = False
		if Shotnumber != None :
			self.Load( Shotnumber )
		

	def __del__( self ):
		self.Unload( )
		del self.Status

	def Load( self , Shotnumber, Experiment = 'AUGD', Edition = 0L, tBegin=-1.0, tEnd=12.0, Magnetics=False,Binning=20.):
		self.Unload()		
	

		try:
			sf_RMD = dd.shotfile( 'RMD' , Shotnumber,Experiment, Edition)
       
		except:
			print "Error reading RMD shotfile" 
			return
		

#			read availability and freqency first for sorting
		avail  = sf_RMD.getParameter('parms-A', 'AVAILABL',dtype=numpy.float64).data
		self.Freq = sf_RMD.getParameter('parms-A', 'f',dtype=numpy.float64).data
		
#			high frequencies first (HFS)
		SortIndex =  numpy.argsort( self.Freq )[::-1]
			# remove bad channels in the sort index
		SortIndex = SortIndex[numpy.where( avail[SortIndex] == 1 )]

		self.Freq = self.Freq[SortIndex]
		self.ChannelNr =  numpy.arange(1,numpy.size(avail)+1)[SortIndex]
	       	self.nChannels = numpy.size(SortIndex)

#		self.test_time = sf_RMD.getTimeBase( 'Trad-A', tBegin=tBegin, tEnd=tEnd)
#		self.test_Te = sf_RMD( 'Trad-A', tBegin=tBegin, tEnd=tEnd ).data[:,SortIndex]
		try:		
			self.z = sf_RMD( 'z-A', tBegin=tBegin, tEnd=tEnd ).data[:,SortIndex]
			self.R = sf_RMD( 'R-A', tBegin=tBegin, tEnd=tEnd ).data[:,SortIndex]
			self.rztime = sf_RMD.getTimeBase( 'R-A', tBegin=tBegin, tEnd=tEnd )	       	
		except:
			print "No R and z in RMD shotfile" 
			sf_RMD.close()
			return
		
		self.Btot = sf_RMD.getParameter('parms-A', 'Btot',dtype=numpy.float64).data[SortIndex]
		self.Ifgroup = sf_RMD.getParameter('parms-A', 'IFGROUP',dtype=numpy.float64).data[SortIndex]
		self.Calfact = sf_RMD.getParameter('parms-A', 'calfact',dtype=numpy.float64).data[SortIndex]
		self.Shift = sf_RMD.getParameter('parms-A', 'shift',dtype=numpy.float64).data[SortIndex]
	
		self.Multi00 = numpy.concatenate( (sf_RMD.getParameter('eCAL-A1', 'MULTIA00',dtype=numpy.float64).data,sf_RMD.getParameter('eCAL-A2', 'MULTIA00',dtype=numpy.float64 ).data),axis=0)[SortIndex]
		self.Shift00 = numpy.concatenate( (sf_RMD.getParameter('eCAL-A1', 'SHIFTB00',dtype=numpy.float64).data,sf_RMD.getParameter('eCAL-A2', 'SHIFTB00',dtype=numpy.float64).data),axis=0)[SortIndex]

		self.Sideband = sf_RMD.getParameter('METHODS', 'SIDEBAND').data
		self.Freqlo = sf_RMD.getParameter('METHODS', 'FREQLO',dtype=numpy.float64).data
		self.nIfgroups = int(sf_RMD.getParameter('METHODS', 'IFGROUPS',dtype=numpy.float64).data)
		
		self.Status_RMD = True
		self.Shotnumber = Shotnumber
				
	       	sf_RMD.close()

		#reading RMC data

		try:
			sf_RMC = dd.shotfile( 'RMC' , Shotnumber, Experiment, Edition)
       
		except:
			print "Error reading RMC shotfile" 
			return
			
		self.time = sf_RMC.getTimeBase( 'Trad-A1', tBegin=tBegin, tEnd=tEnd )
		self.ntimes = numpy.size(self.time)

		self.samplingrate = 1.0/numpy.mean(numpy.diff(self.time))
		
			# read raw rmc data Trad-A1 and Trad-A2 and sort it			
		self.Te = numpy.concatenate( (sf_RMC('Trad-A1', tBegin=tBegin, tEnd=tEnd,dtype=numpy.float64).data,sf_RMC('Trad-A2',tBegin=tBegin,tEnd=tEnd,dtype=numpy.float64).data),axis=1 )[:,SortIndex]
		if Magnetics:
			self.N1 = sf_RMC( 'N=1', tBegin=tBegin, tEnd=tEnd,dtype=numpy.float64 ).data
			self.N2 = sf_RMC( 'N=2', tBegin=tBegin, tEnd=tEnd,dtype=numpy.float64 ).data
			self.OdNEvM = sf_RMC( 'OdNEvM', tBegin=tBegin, tEnd=tEnd,dtype=numpy.float64 ).data
			self.EvNOdM = sf_RMC( 'EvNOdM', tBegin=tBegin, tEnd=tEnd,dtype=numpy.float64 ).data
			self.Status_Magn = True

		sf_RMC.close()
		print "RMC read"

		self.Te = numpy.multiply(self.Multi00,self.Te,dtype=numpy.float64)
#		self.Te = numpy.add(self.Shift00,self.Te,dtype=numpy.float64)
		if (tBegin < -0.1):
			SelectTime = numpy.where( self.time <= 0.0 ) 
			CalcShift = numpy.reshape(numpy.mean(self.Te[SelectTime,:],axis=1,dtype=numpy.float64),(self.nChannels))
			self.Te = numpy.add(-CalcShift,self.Te,dtype=numpy.float64)
		else:		
			
			self.Te = numpy.add(self.Shift00,self.Te,dtype=numpy.float64)

		self.ChannelsSorted = False
		#self.Te = numpy.add(-self.Shift,self.Te)
		#self.Te = numpy.add(-CalcShift,self.Te,dtype=numpy.float64)
		self.Status_RMC = True
		if self.Status_RMC & self.Status_RMD:
			self.Status = True
		
		if numpy.all(Binning) != None:
			self.Binning(Binning)

		self.Rall = numpy.zeros_like(self.Te)
		self.zall = numpy.zeros_like(self.Te)
		self.ntime = numpy.size(self.time)
		
		for i in numpy.arange(self.ntime):
			idx = numpy.argmin(numpy.abs(self.rztime-self.time[i]))
			self.Rall[i,:] = self.R[idx,:]
			self.zall[i,:] = self.z[idx,:]		


		#IPython.embed()

	def Unload( self ):
		
		if self.Status_RMD:

			self.Status_RMD = False
			
		#	del self.test_time	
		#	del self.test_Te	
			del self.R
			del self.z
			del self.rztime
			del self.Shotnumber

			del self.Freq
			del self.Btot
			del self.Ifgroup
			del self.Calfact
			del self.Shift
			del self.nIfgroups
			del self.Sideband
			del self.Freqlo
			del self.ChannelNr
			del self.nChannels

		if self.Status_RMC:

			self.Status_RMC = False

			del self.time
			del self.ntimes
			del self.Te
			del self.samplingrate
			if self.Status_Magn:

				self.Status_Magn = False

				del self.N1
				del self.N2
				del self.OdNEvM
				del self.EvNOdM

		if self.Status: 

			self.Status = False




        #single time point
        def getRzphi(self, time ):
            if self.StatusRz | self.Status:
                idx =  numpy.nanargmin(numpy.abs(self.rztime - time))
                R_out= numpy.squeeze(self.R[idx,:]) 
                z_out= numpy.squeeze(self.z[idx,:]) 
                phi = numpy.repeat([PHI.ECE()],numpy.size(R_out))
                return R_out, z_out, phi


	def funcLoadAllRhop( self, eqExp = 'AUGD', eqDiag = 'EQH', eqi_map=False ):
		
		if self.Status:

                    if (((self.tBegin==-1.0) & (self.tEnd == 12.0)) & eqi_map):
                        self.eqExp = eqExp
                        self.eqDiag = eqDiag
				
                        R =self.Rall  # numpy.zeros_like(self.Te)
                        z =self.zall  # numpy.zeros_like(self.Te)
                        
                        self.rhop = fastkk.eqi_map().kkrzptfn(self.Shotnumber, self.time, R, z, exp=eqExp, diag=eqDiag )
                        self.Status_allrhop = True
                    else:
                        
                        out = kk.KK().kkrzptfn(self.Shotnumber, self.time, self.Rall, self.zall, exp=eqExp, diag=eqDiag )
                        self.rhop = out.rho_p
                        del out
                        
                        self.Status_allrhop = True
			



	def getTe(self, timepoint, sorted = True):

		if self.Status:	

			if (self.ChannelsSorted) & (sorted == True):
				idx = self.SortedIndex
			else:
				idx = numpy.arange(self.nChannels)
	
			if numpy.size(timepoint)== 1:
				return (self.getData( timepoints = timepoint, value='Te'))[idx]

			elif numpy.size(timepoint) == 2:
				return  (self.getData( range = timepoint, value='Te'))[:,idx]

			return None
		else:
			print 'no Data loaded'
			return None


	def getRhop(self, timepoint, sorted = True, eqExp = 'AUGD', eqDiag='EQH'):

		if self.Status:		

			if (self.ChannelsSorted) & (sorted == True):
				idx = self.SortedIndex
			else:
				idx = numpy.arange(self.nChannels)

			if numpy.size(timepoint)== 1 :
				return (self.getData( timepoints = timepoint, value='rhop', eqExp = eqExp, eqDiag=eqDiag))[idx]

			elif numpy.size(timepoint) == 2:
				return  (self.getData( timepoints = timepoint, value='rhop', eqExp = eqExp, eqDiag=eqDiag))[:,idx]

			return None
		else:
			print 'no Data loaded'
			return None

        def MapToRhop(self, timepoint, eqExp = 'AUGD', eqDiag='EQH' ):
            if self.Status:	
                return self.getRhop( timepoint, eqExp = eqExp, eqDiag=eqDiag )



	def get(self, time, sorted = True, eqExp = 'AUGD', eqDiag='EQH'):

		if self.Status:	
			if (self.ChannelsSorted) & (sorted == True):
				idx = self.SortedIndex
			else:
				idx = numpy.arange(self.nChannels)
			
			
			if numpy.size(time)== 1 :
				output = ECEhelp()
				output.Te = (self.getData( timepoints = time, value='Te'))[idx]
				output.rhop = (self.getData( timepoints = time, value='rhop', eqExp = eqExp, eqDiag=eqDiag))[idx]
				output.time = self.getData( timepoints = time, value='time')
				return output

			elif numpy.size(time)== 2:
				output = ECEhelp()
				output.Te = (self.getData( range = time, value='Te'))[:,idx]
				output.rhop = (self.getData( range = time, value='rhop', eqExp = eqExp, eqDiag=eqDiag))[:,idx]
				output.time = self.getData( range = time, value='time')
				return output
		else:
			print 'no Data loaded'
			return None


	def getTe(self, timepoint, sorted = True):

		if self.Status:	

			if (self.ChannelsSorted) & (sorted == True):
				idx = self.SortedIndex
			else:
				idx = numpy.arange(self.nChannels)
	
			if numpy.size(timepoint)== 1:
				return (self.getData( timepoints = timepoint, value='Te'))[idx]

			elif numpy.size(timepoint) == 2:
				return  (self.getData( range = timepoint, value='Te'))[:,idx]

			return None
		else:
			print 'no Data loaded'
			return None





	#Default Te, indexbased
	def getData( self , timepoints = None, range = None, value='Te', eqExp = 'AUGD', eqDiag='EQH' ):
		ntimes = 0

		if self.Status:

			if not ( (value == 'Te') | (value == 'rhop')| (value == 'time') ):
				print 'value must be Te or ne or rhop'
				return None

			if (value == 'rhop'):
				if not self.Status_allrhop:
					print 'not all rhop values loaded, every point has to be loaded: very slow'

			if (timepoints == None) & (range == None) :
				print 'no timepoints or range are given, return None'
				return None

			if timepoints != None:
				
				timepoints = numpy.array(timepoints)
				if range != None:
					print 'either range or timepoints must be given'
					return False

				ntimes = numpy.size(timepoints)
				if ntimes > 1:
					#get min from the time matrix of the timepoints vector, 1. self.time is copied ntimes, then timepoints is subtracted row wise.
					search_array = numpy.sum([numpy.array([self.time,]*ntimes),numpy.multiply(timepoints,-1.)],axis=0)
					# the indices of the closest points are calculated
					idx = numpy.nanargmin(numpy.abs( numpy.reshape( numpy.concatenate(search_array),(ntimes,self.ntime) )),axis=1)
				
					#if (value == 'Te') | (value == 'ne') | (value == 'rhop'):
					#ne/Te/rhop should have the same area base
						#out_temp = numpy.zeros((ntimes,numpy.size(self.Te[0,:])))

				elif ntimes == 1:
					idx = numpy.argmin( numpy.abs( self.time - timepoints ) )
					if (value == 'time'): 
						return self.time[idx]

				else:
					print 'Something went wrong with time points'
					return None
				

			if range != None:
				range = numpy.array(range)

				if not numpy.size(range) == 2:
					print 'range must have two elements'
					return None
				if ( (range[0]>range[1]) & (range[0] >= self.tBegin) & (range[1] <= self.tEnd)):
					print 'second Element must be bigger than first and range must be in range of tBegin or tEnd'
					return None				
				if timepoints != None:
					print 'either range or timepoints must be given'
					return None
				idx = numpy.squeeze(numpy.where( (self.time[:] >= range[0]) & (self.time[:] <= range[1] ) ))
				
				if (value == 'time'): 
					return self.time[idx]

	
			if (value == 'Te') |  (value == 'rhop'):
					#ne/Te/rhop should have the same area base
				
				if value == 'Te':
					return  self.Te[idx,:]
				
				if (value == 'rhop'):
					if (self.Status_allrhop):
						if  (Experiment == self.eqExp) & (Diagnostic == self.eqDiag):
							return self.rhop[idx,:]
						else:
							print 'Equlibrium Experiment differ, calculate new'
							loadAllRhop( eqExp=Experiment, eqDiag=Diagnostic)
							return self.rhop[idx,:]

					else:
						nindex = numpy.size(idx)
						if nindex > 1: 
							out_temp = numpy.zeros((nindex,numpy.size(self.Te[0,:])))
							for i in numpy.arange(nindex):
								rz_idx = numpy.argmin(numpy.abs(self.rztime - self.time[idx[i]]))
								output = kk.KK().kkrzptfn( self.Shotnumber, self.time[idx[i]], self.R[rz_idx,:], self.z[rz_idx,:], exp= eqExp, diag=eqDiag, ed=0)
								out_temp[i,:] = output.rho_p	
                                                        return out_temp
						else:
							rz_idx = numpy.argmin(numpy.abs(self.rztime - self.time[idx]))
							output = kk.KK().kkrzptfn( self.Shotnumber, self.time[idx], self.R[rz_idx,:], self.z[rz_idx,:], exp= eqExp, diag=eqDiag, ed=0)
							return output.rho_p	

			
		return None
		


						  		
	def Binning( self, samplefreq = 8.0 ):				   
		if self.Status:				    						   				
			print "RMC binning with ",samplefreq," kHz"
			bins = int(self.ntimes*(float(samplefreq)*1.0e3/self.samplingrate))

			tempTe = self.Te

			slices= numpy.linspace(0, self.ntimes, bins+1, True).astype(numpy.int)
			counts = numpy.diff(slices)

			self.time = numpy.add.reduceat(self.time, slices[:-1]) / counts
			self.ntimes = numpy.size(self.time)
			self.samplingrate = 1.0/numpy.mean(numpy.diff(self.time))
			self.Te = numpy.zeros((self.ntimes,self.nChannels))
			for i in range(self.nChannels):
				self.Te[:,i] = numpy.add.reduceat(tempTe[:,i], slices[:-1]) / counts
	
			try:
				self.Rall = numpy.add.reduceat(self.Rall, slices[:-1]) / counts
				self.zall = numpy.add.reduceat(self.zall, slices[:-1]) / counts
			except:
				print 'no Rall or zall'

	# can be used to transform any Signal in to the timebase of ECE data
	def TransformTimebase( self, time_in, signal_in ):
		
		if self.Status:
			signal_out = numpy.zeros(self.ntimes)
			f=scipy.interpolate.interp1d( time_in, signal_in, bounds_error=False, fill_value=0 )
			return f(self.time)
		else:
			return 0


	def __call__( self , timepoints ):
		ntimes = 0
		if self.Status:
			ntimes = numpy.size(timepoints)			
			Te_temp = numpy.zeros((ntimes,numpy.size(self.Te[0,:])))
			try:
				for i in range(ntimes):
					idx = numpy.argmin( numpy.abs( self.time - timepoints[i] ) )
					Te_temp[i,:] = self.Te[idx,:]

			except Exception:
				print 'timepoints: ',timepoints
				idx = numpy.argmin(numpy.abs(self.time-timepoints))
				output = self.Te[idx,:]
				pass
			
			if ntimes > 1: 	
				return Te_temp
			else: 
				return output			

		
