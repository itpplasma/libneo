import numpy
import scipy.interpolate
import kk
import IPython
import eqi_map as fastkk
import tools

class REF:
	def __init__( self ,  Experiment = 'AUGD', Shotnumber = None ):
		self.Status = False
		if Shotnumber != None:
			self.Load( Shotnumber )
		self.Status_allrhop = False

	def __del__( self ):
		self.Unload( )
		self.Status = False
		self.Status_allrhop = False

	def Unload( self ):
		if self.Status:
			del self.time
			del self.ne
			del self.R
			del self.z
			del self.ntime

	def Load( self ,  Shotnumber, Experiment='AUGD', Diagnostic='RPS', loadAllRhop=False, loadAllRhot=False, Edition = 0L, tBegin=-1.0, tEnd=12.0, eqExp = 'AUGD', eqDiag = 'EQH',Rshift=0.0 ):

		self.Unload()
		if Shotnumber< 40000:

			path='/afs/ipp-garching.mpg.de/home/m/mwillens/Programme/python/data/REF/'+str(Shotnumber)
			filename_ne=path+'/ne.dat'
			filename_R=path+'/R.dat'
			filename_time=path+'/time.dat'

			time = numpy.mean(numpy.genfromtxt(filename_time), axis=1)
			idx = numpy.where( ( time > tBegin ) & ( time < tEnd ))
			self.time = time[idx]
			self.ntime = numpy.size(self.time)
			self.ne =numpy.genfromtxt(filename_ne)[idx]
			self.R = numpy.genfromtxt(filename_R)[idx]+Rshift
			self.z = numpy.zeros_like(self.R,dtype=float)
			self.nChannels = numpy.size(self.R[0,:])
			self.Shotnumber = Shotnumber

			if loadAllRhop:
					#check if the entire dataset is read
				if ( tBegin == -1.0 ) & ( tEnd == 12.0 ):

					self.eqExp = eqExp
					self.eqDiag = eqDiag
					self.rhop = fastkk.eqi_map().kkrzptfn(self.Shotnumber, self.time, self.R, self.z, exp=eqExp, diag=eqDiag )
				#	IPython.embed()
					self.Status_allrhop = True
				else:

					print 'It is not possible to read Tomas kk with limited timerange, but time is not eqidistant'
					self.Status_allrhop = False
					return False


			if loadAllRhot:
					#check if the entire dataset is read
				if ( tBegin == -1.0 ) & ( tEnd == 12.0 ):
					
					self.eqExp = eqExp
					self.eqDiag = eqDiag
					self.rhot = fastkk.eqi_map().kkrzptfn(self.Shotnumber, self.time, self.R, self.z, exp=eqExp, diag=eqDiag, rho_lbl='rho_tor' )
				#	IPython.embed()
					self.Status_allrhot = True
				else:

					print 'It is not possible to read Tomas kk with limited timerange, but time is not eqidistant'
					self.Status_allrhot = False
					return False



			self.Status = True
		#	/afs/ipp-garching.mpg.de/home/m/mwillens/Programme/python/data/REF/31106




	def __call__( self , timepoints , rhop = None, Experiment = 'AUGD', Diagnostic='EQH' ):
		ntimes = 0
		ne_temp = None
		if self.Status:
			ntimes = numpy.size(timepoints)			
			if ntimes > 1:
				if rhop != None:

					ne_temp = numpy.zeros((ntimes,numpy.size(rhop)))
					distance = numpy.diff(timepoints)
					stdDistance = numpy.std(distance)
					meanDistance = numpy.mean(distance)

					if (stdDistance/meanDistance > 0.25) :
						print 'Timebase not equidistant'
						return False
					
                               #get min from the time matrix of the timepoints vector, 1. self.time is copied ntimes, then timepoints is subtracted row wise.
					search_array = numpy.sum([numpy.array([self.time,]*ntimes),numpy.multiply(timepoints,-1.)],axis=0)
				# the indices of the closest points are calculated
					idx = numpy.nanargmin(numpy.abs( numpy.reshape( numpy.concatenate(search_array),(2,self.ntime) )),axis=1)
				
					rhop_out = fastkk.eqi_map().kkrzptfn(self.Shotnumber, self.time[idx], self.R[idx,:], self.z[idx,:], exp=Experiment, diag=Diagnostic )
					
					time_in =numpy.reshape(numpy.repeat(self.time[idx],numpy.size(self.ne[0,:])),(ntimes,numpy.size(self.ne[0,:])) )

					#the following should be changed, if assume sorted is available, 
					#for reflectometry take of the array the point which is most ar from the center (beg_idx), and the first index which is deepest in the plasma
					beg_idx = max(numpy.nanargmin(rhop_out,axis=1))
					end_idx = min(numpy.nanargmax(rhop_out,axis=1))

					ne = scipy.interpolate.interp2d( time_in[:,beg_idx:end_idx], rhop_out[:,beg_idx:end_idx] , self.ne[idx,beg_idx:end_idx].transpose())#,  bounds_error=True, fill_value = 0.0)
					ne_temp = ne( self.time[idx], rhop ).transpose()
					
#self.Shotnumber, self.time[rz_idx], self.R[rz_idx,:], self.z[rz_idx,:]
					#IPython.embed()
					
				else:
					ne_temp = numpy.zeros((ntimes,numpy.size(self.ne[0,:])))
					for i in range(ntimes):
						idx = numpy.argmin( numpy.abs( self.time - timepoints[i] ) )
						ne_temp[i,:] = self.ne[idx,:]
			elif ntimes == 1:
				if rhop != None:
					idx = numpy.argmin(numpy.abs(self.time-timepoints))
					rhop_temp =  self.MapToRhop( self.time[idx], Experiment = Experiment, Diagnostic=Diagnostic )
					if rhop_temp[1] > rhop_temp[0]:
						ne = scipy.interpolate.interp1d( rhop_tmp , self.ne[:,idx] )
					else:
						ne = scipy.interpolate.interp1d( rhop_tmp , self.ne[::-1,idx] )
					
					ne_temp = ne( rhop )
				else:
					idx = numpy.argmin(numpy.abs(self.time-timepoints))
					ne_temp = self.ne[idx,:]
				pass
			
			if ne_temp != None : 	
				return ne_temp
			else: 
				print "error in reading data"
				return False



	def MapToRhot( self , timepoints, Experiment = 'AUGD', Diagnostic='EQH' ):
		ntimes = 0
		if self.Status:
			ntimes = numpy.size(timepoints)
			rhot = numpy.zeros( (ntimes,numpy.size(self.R[0,:])) )
			try:
				for i in range(ntimes):
					rz_idx = numpy.argmin( numpy.abs( self.time - timepoints[i] ) )
					output = kk.KK().kkrzptfn( self.Shotnumber, self.time[rz_idx], self.R[rz_idx,:], self.z[rz_idx,:], exp= Experiment, diag=Diagnostic, ed=0)
					rhot[i,:] = output.rho_t
			except Exception:
				rz_idx = numpy.argmin(numpy.abs(self.time-timepoints))
				output = kk.KK().kkrzptfn( self.Shotnumber, self.time[rz_idx], self.R[rz_idx,:], self.z[rz_idx,:], exp= Experiment, diag=Diagnostic, ed=0)
				pass
			print rz_idx
			if ntimes > 1: 	
				
				return rhot
			else: 
				return output.rho_t
			


	def MapToRhop( self , timepoints, Experiment = 'AUGD', Diagnostic='EQH' ):
		ntimes = 0
		if self.Status:
			ntimes = numpy.size(timepoints)
			rhop = numpy.zeros( (ntimes,numpy.size(self.R[0,:])) )
			try:
				for i in range(ntimes):
					rz_idx = numpy.argmin( numpy.abs( self.time - timepoints[i] ) )
					output = kk.KK().kkrzptfn( self.Shotnumber, self.time[rz_idx], self.R[rz_idx,:], self.z[rz_idx,:], exp= Experiment, diag=Diagnostic, ed=0)
					rhop[i,:] = output.rho_p
			except Exception:
				rz_idx = numpy.argmin(numpy.abs(self.time-timepoints))
				output = kk.KK().kkrzptfn( self.Shotnumber, self.time[rz_idx], self.R[rz_idx,:], self.z[rz_idx,:], exp= Experiment, diag=Diagnostic, ed=0)
				pass
			print rz_idx
			if ntimes > 1: 	
				
				return rhop 
			else: 
				return output.rho_p
			

        def getRsep( self ):
            if self.Status_allrhop:
                ntimes = self.ntime
                data = numpy.zeros((ntimes))
                data_in = [1.0]

                x = self.rhop
                y = self.R

		if x[0,0] < x[0,-1]:
			data = tools.interpolSimpleArray(x,y,data_in)
		else:
                        data = tools.interpolSimpleArray(x[:,::-1],y[:,::-1],data_in)

		self.R_sep = self.R - numpy.reshape(numpy.repeat(data,self.nChannels),(ntimes,self.nChannels))

		return data




	#Default Te, indexbased
	def getData( self , timepoints = None, range = None, value='ne', Experiment = 'AUGD', Diagnostic='EQH' ):
		ntimes = 0

		if self.Status:

			if not ( (value == 'ne') | (value == 'rhop')| (value == 'time') ):
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

	
                        if (value == 'ne') | (value == 'rhop'):
					#ne/Te/rhop should have the same area base
								
                            if value == 'ne':
                                return  self.ne[idx,:]
				
                            if (value == 'rhop'):
                                if (self.Status_allrhop):
                                    if  (Experiment == self.eqExp) & (Diagnostic == self.eqDiag):
                                        return self.rhop[idx,:]
                                    else:
                                        print 'Equlibrium Experiment differ, calculate new'
                                        rhop_out = fastkk.eqi_map().kkrzptfn(self.Shotnumber, self.time, self.R, self.z, exp=Experiment, diag=Diagnostic )
                                        self.eqDiag = Diagnostic
                                        self.eqExp = Experiment
                                        return rhop_out[idx,:]

                                else:
                                    nindex = numpy.size(idx)
                                    if nindex > 1: 
                                        out_temp = numpy.zeros((nindex,numpy.size(self.ne[0,:])))
                                        for i in numpy.arange(nindex):	
                                            output = kk.KK().kkrzptfn( self.Shotnumber, self.time[idx[i]], self.R[idx[i],:], self.z[idx[i],:], exp= Experiment, diag=Diagnostic, ed=0)
                                            out_temp[i,:] = output.rho_p	
                                        return out_temp
                                    else:
                                        output = kk.KK().kkrzptfn( self.Shotnumber, self.time[idx], self.R[idx,:], self.z[idx,:], exp= Experiment, diag=Diagnostic, ed=0)
                                        return output.rho_p	





	def get(self, time, eqExp = 'AUGD', eqDiag='EQH'):

		if self.Status:		
			if (numpy.size(time)== 1) | ( numpy.size(time) > 2):
				output = RPShelp()
				output.ne = self.getData( timepoints = time, value='ne')
				output.rhop = self.getData( timepoints = time, value='rhop')
				output.time = self.getData( timepoints = time, value='time')
				return output

			elif numpy.size(time)== 2:
				output = RPShelp()
				output.ne = self.getData( range = time, value='ne')
				output.rhop = self.getData( range = time, value='rhop')
				output.time = self.getData( range = time, value='time')
				return output
		else:
			print 'no Data loaded'
			return None





	def getne(self, timepoint):

		if self.Status:		
			if (numpy.size(timepoint)== 1) | ( numpy.size(timepoint) > 2):
				return self.getData( timepoints = timepoint, value='ne')

			elif numpy.size(timepoint)== 2:
				return  self.getData( range = timepoint, value='ne')

			return None
		else:
			print 'no Data loaded'
			return None



	def getRhop(self, timepoint, eqExp = 'AUGD', eqDiag='EQH'):

		if self.Status:		
			if (numpy.size(timepoint)== 1) | ( numpy.size(timepoint) > 2):
				return self.getData( timepoints = timepoint, value='rhop')

			elif numpy.size(timepoint)== 2:
				return  self.getData( range = timepoint, value='rhop')

			return None
		else:
			print 'no Data loaded'
			return None



	def getTime(self, time):

		if self.Status:		
			if (numpy.size(timepoint)== 1) | ( numpy.size(timepoint) > 2):
				return self.getData( timepoints = timepoint, value='time')

			elif numpy.size(timepoint)== 2:
				return  self.getData( range = timepoint, value='time')

			return None
		else:
			print 'no Data loaded'
			return None



	def __call__( self , time , rhop = None, eqExp = 'AUGD', eqDiag='EQH' ):

		if self.Status:
			if rhop == None:
				return self.getne( time )
			else:
				ne_out = self.getne( time )
				rhop_out = self.getRhop( time )
				ntimes = numpy.size(time)
				if ntimes > 1:
					ne_temp = numpy.zeros((ntimes,numpy.size(rhop)))
					for i in numpy.arange(ntimes) :						
						if rhop_out[i,1]<= rhop_out[i,0]:
							ne = scipy.interpolate.interp1d(rhop_out[i,::-1] , ne_out[i,::-1],  bounds_error=False, fill_value = float('nan'))
						else:
							ne = scipy.interpolate.interp1d(rhop_out[i,:] , ne_out[i,:1],  bounds_error=False, fill_value = float('nan'))
						ne_temp[i,:] = ne( rhop )
					return ne_temp
				else:
					if rhop_out[i,1] <= rhop_out[i,0]:
						ne = scipy.interpolate.interp1d(rhop_out[i,::-1] , ne_out[i,::-1],  bounds_error=False, fill_value = float('nan'))
					else:
						ne = scipy.interpolate.interp1d(rhop_out[i,:] , ne_out[i,:1],  bounds_error=False, fill_value = float('nan'))
					return ne( rhop )					

		else:
			print 'no Data loaded'
			return None



