import numpy
import dd as dd
import kk as kk
import scipy.interpolate

class ECE:
	def __init__( self ,  Experiment = 'AUGD', Diagnostic = 'MAI', Shotnumber = None ):
		self.Status = False
		self.exists_Rmaj = False
		self.exists_rhop = False
		if Shotnumber != None :
			self.Load( Shotnumber )
		

	def __del__( self ):
		self.Unload( )
		del self.Status

	def Load( self , Experiment, Diagnostic, Shotnumber, Edition = 0L ):
		self.Unload()
		if Diagnostic == 'MAI':
			sf = dd.shotfile()
			if sf.Open( Diagnostic , Shotnumber, Experiment, Edition  ):
				self.time = sf.GetTimebase( 'BTF' )
				self.BTF = sf.GetSignal( 'BTF' )
				sf.Close()
		else:
			print 'Diagnotic should be CEC or RMD, but is:' ,Diagnostic


	def Unload( self ):
		if self.Status:
			self.Status = False
			del self.time	
			del self.BTF	



	def MapToRhop_timepoints( self , timepoints, Experiment = 'AUGD', Diagnostic='EQH' ):
		ntimes = 0
		if self.Status:
			ntimes = numpy.size(timepoints)
		#	if type(timepoints) in (list, numpy.ndarray):
		#		ntimes=len(timepoints)
		#	else: 
		#	        
		#		timepoints=[timepoints]
		#		ntimes = timepoints.size()
				#ntimes=len(timepoints)
			
			rhop = numpy.zeros((ntimes,len(self.R[0,:])))
			try:
				for i in range(ntimes):
					rz_idx = numpy.argmin( numpy.abs( self.rztime - timepoints[i] ) )
					output = kk.KK().kkrzptfn( self.Shotnumber, self.rztime[rz_idx], self.R[rz_idx,:], self.z[rz_idx,:], exp= Experiment, diag=Diagnostic, ed=0)
					rhop[i,:] = output.rho_p
			except Exception:
				rz_idx = numpy.argmin(numpy.abs(self.rztime-timepoints))
				output = kk.KK().kkrzptfn( self.Shotnumber, self.rztime[rz_idx], self.R[rz_idx,:], self.z[rz_idx,:], exp= Experiment, diag=Diagnostic, ed=0)
				pass
			
			if ntimes > 1: 	
				return rhop 
			else: 
				return output.rho_p
			
	

	# Map all points to
	def MapToRhop_all( self , Experiment = 'AUGD', Diagnostic='EQH' ):
	
		if self.Status & self.exists_rhop == False:

			self.rhop = numpy.zeros_like(self.Te)
			for i in range(len(self.Te)):
							
				rz_idx = numpy.argmin( numpy.abs( self.rztime - self.time[i] ) )
				
				output = kk.KK().kkrzptfn( self.Shotnumber, self.time[i], self.R[rz_idx,:], self.z[rz_idx,:], exp= Experiment, diag=Diagnostic, ed=0)
				self.rhop[i,:] = output.rho_p
				
			self.exists_rhop = True



	def SortECEChannels( self ):

		if self.Status:
			
			sort_ECE = numpy.zeros(len(self.Freq), dtype=[('freq',float),('idx',int),('avail',float),('ifgroup',float)])
			sort_ECE['freq'] = self.Freq
			sort_ECE['idx']= numpy.arange(int(len(self.Freq)))
			sort_ECE['avail']= self.avail 
			sort_ECE['ifgroup']= self.Ifgroup

			a = numpy.sort( sort_ECE, axis=0, order='freq' ) 
	
			self.SortedIndex =  a['idx'][numpy.where(a['avail'] == 1) ]
		

			self.IfGroupIndex = []
			
			for i in range(1,self.nIfgroups+1):
				
				self.IfGroupIndex.append(a['idx'][numpy.where( (a['ifgroup'] == i) & (a['avail'] == 1) )])
			
			self.ChannelsSorted = True

#			numpy.zeros(len(numpy.where(a['avail'] == 1)))

	#		self.Ifchannels = numpy.array( self.nIfgroups )
#			self.Ifchannels[0] =  a['idx'][numpy.where(a['ifgroup'] == 1)]
#= a['idx'][numpy.where(a['avail'] == 1)]

	
#	def GetTeProfile( self, time ):					
#		if self.Status:


#	def __call__( self , time , rhop = None ):

#		if self.Status:
	#		idx = numpy.argmin( numpy.abs( self.time - time ) )
#			rz_idx = numpy.argmin( numpy.abs( self.rztime - time ) )


	#		print idx,self.time[idx],time
		#	print rz_idx,self.rztime[rz_idx]
	#		if rhop != None:
	#			if self.rhop[1,idx] > self.rhop[0,idx]:
	#				Te = scipy.interpolate.interp1d( self.rhop[:,idx] , self.ne[:,idx] )
	#			else:
	#				Te = scipy.interpolate.interp1d( self.rhop[::-1,idx] , self.ne[::-1,idx] )
	#			return Te( rhop )
	#		else:
	#			return self.rhop[:,idx] , self.ne[:,idx]

		
