import numpy
import dd
import scipy.interpolate
import kk_mwillens as kk
import IPython

class TShelp:
    status = False

class TS:
	def __init__( self ,  Shotnumber = None,Experiment = 'AUGD' ):
		self.Status = False
		self.Status_edge = False
		self.Status_core = False
		self.Status_allrhot = False
		self.Status_allrhop = False
		if Shotnumber != None:
			self.Load( Shotnumber )

	def __del__( self ):
		self.Unload( )
		self.Status = False
		self.Status_edge = False
		self.Status_core = False
		self.Status_allrhop = False	
                self.Status_allrhot = False	

	def Unload( self ):
		if self.Status_edge:
			del self.time_edge
			del self.ne_edge 
			del self.neupp_edge  
			del self.nelow_edge 
			del self.Te_edge
 			del self.Teupp_edge  
			del self.Telow_edge  
			del self.R_edge
			del self.z_edge
			self.Status_edge = False
		if self.Status_core:
			del self.time_core
			del self.ne_core
			del self.neupp_core  
			del self.nelow_core  
			del self.Te_core 
			del self.Teupp_core  
			del self.Telow_core  
			del self.R_core
			del self.z_core	     
			self.Status_core = False
		if self.Status_allrhop:
			del self.rhop
			del self.eqExp
			del self.eqDiag
		if self.Status:
			del self.Shotnumber
			del self.tBegin
			del self.tEnd
			del self.time
			del self.ne 
			del self.Te 
			del self.R
			del self.z	
			self.Status = False
		

	def Load( self ,  Shotnumber, Experiment='AUGD', Diagnostic='VTA', Edition = 0L, tBegin=-1.0, tEnd=12.0, CoreOrEdge='core', loadAllRhop=False ,eqExp = 'AUGD', eqDiag = 'EQH'):
		self.Unload()
		sf = dd.shotfile()
		if Diagnostic == 'VTA':
			try:
				sf = dd.shotfile( Diagnostic, Shotnumber, Experiment, Edition)
			except:
				print "Error reading shotfile" 
				return False
			self.Shotnumber = Shotnumber
			if Shotnumber > 27434 :
				if (CoreOrEdge != 'edge') & (CoreOrEdge != 'core')  :
					print "CoreOrEdge keyword is not edge or core, it is: ", CoreOrEdge 
					return False
				if (CoreOrEdge == "core") :
					try:
						self.time_core = sf.getTimeBase( 'Ne_c', tBegin=tBegin, tEnd=tEnd )
						self.ne_core = sf( 'Ne_c', tBegin=tBegin, tEnd=tEnd ).data
                                                self.neupp_core = sf( 'Neupp_c', tBegin=tBegin, tEnd=tEnd ).data
                                                self.nelow_core = sf( 'Nelow_c', tBegin=tBegin, tEnd=tEnd ).data
						self.Te_core = sf( 'Te_c', tBegin=tBegin, tEnd=tEnd ).data
                                                self.Teupp_core = sf( 'Teupp_c', tBegin=tBegin, tEnd=tEnd ).data
                                                self.Telow_core = sf( 'Telow_c', tBegin=tBegin, tEnd=tEnd ).data
						self.z_core = sf( 'Z_core', tBegin=tBegin, tEnd=tEnd ).data
						self.R_core = sf( 'R_core', tBegin=tBegin, tEnd=tEnd ).data
						self.Status_core = True
					except: 
						print "Error reading data" 
						sf.close()
						return False			
					
					self.time = self.time_core
                                        self.ntime = numpy.size(self.time)
					self.ne =self.ne_core
					self.Te =self.Te_core
					self.R = numpy.reshape(numpy.repeat(self.R_core,numpy.size(self.z_core)),(numpy.size(self.R_core),numpy.size(self.z_core)))
					self.z =numpy.reshape(numpy.tile(self.z_core,numpy.size(self.R_core)),(numpy.size(self.R_core),numpy.size(self.z_core)))
					
					self.tBegin = tBegin
					self.tEnd = tEnd

					self.Status = True

				if (CoreOrEdge == "edge")  :
					try:
						self.time_edge = sf.getTimeBase( 'Ne_e', tBegin=tBegin, tEnd=tEnd )
						self.ne_edge = sf( 'Ne_e', tBegin=tBegin, tEnd=tEnd ).data
                                                self.neupp_edge = sf( 'Neupp_e', tBegin=tBegin, tEnd=tEnd ).data
                                                self.nelow_edge = sf( 'Nelow_e', tBegin=tBegin, tEnd=tEnd ).data
						self.Te_edge = sf( 'Te_e', tBegin=tBegin, tEnd=tEnd ).data
                                                self.Teupp_edge = sf( 'Teupp_e', tBegin=tBegin, tEnd=tEnd ).data
                                                self.Telow_edge = sf( 'Telow_e', tBegin=tBegin, tEnd=tEnd ).data
						self.z_edge = sf( 'Z_edge', tBegin=tBegin, tEnd=tEnd ).data
						self.R_edge = sf( 'R_edge', tBegin=tBegin, tEnd=tEnd ).data
						self.Status_edge = True
					except:
						print "Error reading data" 
						sf.close()
						return False						
					
					self.time = self.time_edge
                                        self.ntime = numpy.size(self.time)
					self.ne =self.ne_edge
					self.Te =self.Te_edge
					self.R = numpy.reshape(numpy.repeat(self.R_edge,numpy.size(self.z_edge)),(numpy.size(self.R_edge),numpy.size(self.z_edge)))
					self.z =numpy.reshape(numpy.tile(self.z_edge,numpy.size(self.R_edge)),(numpy.size(self.R_edge),numpy.size(self.z_edge)))			
					self.Status = True


			else:
				try:
					self.time = sf.getTimeBase( 'Ne', tBegin=tBegin, tEnd=tEnd )
                                        self.ntime = numpy.size(self.time)
					self.ne = sf( 'Ne', tBegin=tBegin, tEnd=tEnd ).data
				
					self.Te = sf( 'Te', tBegin=tBegin, tEnd=tEnd ).data
					self.z_tmp = sf( 'Z', tBegin=tBegin, tEnd=tEnd ).data
					self.R_tmp = sf( 'R', tBegin=tBegin, tEnd=tEnd ).data
				except:
					print "Error reading data" 
					sf.close()
					return False

				self.R = numpy.reshape(numpy.repeat(self.R_edge,numpy.size(self.z_edge)),(numpy.size(self.R_edge),numpy.size(self.z_edge)))
				self.z =numpy.reshape(numpy.tile(self.z_edge,numpy.size(self.R_core)),(numpy.size(self.R_edge),numpy.size(self.z_edge)))	
				self.Status = True

					#calcule all rhop values
			if loadAllRhop:
					#check if the entire dataset is read
			
                            self.eqExp = eqExp
                            self.eqDiag = eqDiag
                                        #IPython.embed()
                            out = kk.KK().kkrzptfn( self.Shotnumber, self.time, self.R, self.z, exp= eqExp, diag=eqDiag, ed=0)
                            self.rhop = out.rho_p 
                            self.Status_allrhop = True

			sf.close()
			return True

				

	def get(self, time):

		if self.Status:		
			if (numpy.size(time)== 1) | ( numpy.size(time) > 2):
				output = TShelp()
				output.Te = self.getData( timepoints = time, value='Te')
				output.ne = self.getData( timepoints = time, value='ne')
				output.rhop = self.getData( timepoints = time, value='rhop')
				output.time = self.getData( timepoints = time, value='time')
				return output

			elif numpy.size(time)== 2:
				output = TShelp()
				output.Te = self.getData( range = time, value='Te')
				output.ne = self.getData( range = time, value='ne')
				output.rhop = self.getData( range = time, value='rhop')
				output.time = self.getData( range = time, value='time')
				return output
		else:
			print 'no Data loaded'
			return None


	def getTe(self, timepoint):

		if self.Status:		
			if (numpy.size(timepoint)== 1) | ( numpy.size(timepoint) > 2):
				return self.getData( timepoints = timepoint, value='Te')

			elif numpy.size(timepoint)== 2:
				return  self.getData( range = timepoint, value='Te')

			return None
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

	def getRhop(self, timepoint):
		if self.Status:		
			if (numpy.size(timepoint)== 1) | ( numpy.size(timepoint) > 2):
				return self.getData( timepoints = timepoint, value='rhop')

			elif numpy.size(timepoint)== 2:
				return  self.getData( range = timepoint, value='rhop')

			return None
		else:
			print 'no Data loaded'
			return None

	def getRhot(self, timepoint):
		if self.Status:		
			if (numpy.size(timepoint)== 1) | ( numpy.size(timepoint) > 2):
				return self.getData( timepoints = timepoint, value='rhot')

			elif numpy.size(timepoint)== 2:
				return  self.getData( range = timepoint, value='rhot')

			return None
		else:
			print 'no Data loaded'
			return None

        def MapToRhot(self, timepoint):
		if self.Status:		
                    return self.getRhot(timepoint)

        def MapToRhop(self, timepoint):
		if self.Status:		
                    return self.getRhop(timepoint)


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



	#Default Te, indexbased
	def getData( self , timepoints = None, range = None, value='Te', Experiment = 'AUGD', Diagnostic='EQH' ):
		ntimes = 0

		if self.Status:

			if not ( (value == 'Te') | (value == 'ne') | (value == 'rhop')| (value == 'time') | (value == 'rhot')):
				print 'value must be Te or ne or rhop'
				return None

			if (value == 'rhop'):
				if not self.Status_allrhop:
					print 'not all rhop values loaded, every point has to be loaded: very slow'

			if (value == 'rhot'):
				if not self.Status_allrhot:
					print 'not all rhot values loaded, every point has to be loaded: very slow'

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
                            
                            if ( (range[0]>range[1]) & (range[0] >= self.time[0]) & (range[1] <= self.time[-1])):
                                print 'second Element must be bigger than first and range must be in range of tBegin or tEnd'
                                return None				
                            
                            if timepoints != None:
                                print 'either range or timepoints must be given'
                                return None
				
                            idx = numpy.squeeze(numpy.where( (self.time[:] >= range[0]) & (self.time[:] <= range[1] ) ))
				
                            if (value == 'time'): 
                                return self.time[idx]

	
                        if (value == 'Te') | (value == 'ne') | (value == 'rhop') | (value == 'rhot'):
					#ne/Te/rhop should have the same area base
				
                            if value == 'Te':
                                return  self.Te[idx,:]
				
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
                                        out_temp = numpy.zeros((nindex,numpy.size(self.Te[0,:])))
                                        for i in numpy.arange(nindex):	
                                            output = kk.KK().kkrzptfn( self.Shotnumber, self.time[idx[i]], self.R[idx[i],:], self.z[idx[i],:], exp= Experiment, diag=Diagnostic, ed=0)
                                            out_temp[i,:] = output.rho_p	
                                        return out_temp
                                    else:
                                        output = kk.KK().kkrzptfn( self.Shotnumber, self.time[idx], self.R[idx,:], self.z[idx,:], exp= Experiment, diag=Diagnostic, ed=0)
                                        return output.rho_p

	
                            if (value == 'rhot'):
                                if (self.Status_allrhot):
                                    if  (Experiment == self.eqExp) & (Diagnostic == self.eqDiag):
                                        return self.rhot[idx,:]
                                    else:
                                        print 'Equlibrium Experiment differ, calculate new'
                                        rhot_out = fastkk.eqi_map().kkrzptfn(self.Shotnumber, self.time, self.R, self.z, exp=Experiment, diag=Diagnostic )
                                        self.eqDiag = Diagnostic
                                        self.eqExp = Experiment
                                        return rhot_out[idx,:]

                                else:
                                    nindex = numpy.size(idx)
                                    if nindex > 1: 
                                        out_temp = numpy.zeros((nindex,numpy.size(self.Te[0,:])))
                                        for i in numpy.arange(nindex):	
                                            output = kk.KK().kkrzptfn( self.Shotnumber, self.time[idx[i]], self.R[idx[i],:], self.z[idx[i],:], exp= Experiment, diag=Diagnostic, ed=0)
                                            out_temp[i,:] = output.rho_t	
                                        return out_temp
                                    else:
                                        output = kk.KK().kkrzptfn( self.Shotnumber, self.time[idx], self.R[idx,:], self.z[idx,:], exp= Experiment, diag=Diagnostic, ed=0)
                                        return output.rho_t

                           
                                        
                return None
		

