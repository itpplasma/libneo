import numpy
import dd
import kk_mwillens as kk
import scipy.interpolate
import matplotlib.pylab as plt
import IPython


class initialize_output:

	status = False

class IDA:

	def __init__( self ,  Experiment = 'AUGD', Shotnumber = None ):
		self.Status = False
		self.StatusTe = False
		self.StatusPe = False
		self.Status_mean = False
		self.StatusECE = False
		if Shotnumber != None:
			self.Load( Shotnumber )

	def __del__( self ):
		self.Unload( )
		del self.Status


	def Load( self ,   Shotnumber, Experiment = 'AUGD', Edition = 0L, tBegin=-1.0, tEnd=12.0 , ECEraw = False):

		self.Unload()
		try:
			print Shotnumber
			sf = dd.shotfile( 'IDA', Shotnumber, Experiment, Edition)
		except:
			print "Error reading shotfile" 
			return
	
		print tBegin, tEnd
		self.time = sf( 'ne' ).time
		index = numpy.squeeze( numpy.where( ( self.time > tBegin ) & ( self.time < tEnd ) ))
		self.time = self.time[index]
		self.rhop = sf( 'rhop'  ).data[index]
	
		self.ne = sf( 'ne'  ).data[index]
		self.ne_lo = sf( 'ne_lo').data[index]
		self.ne_up = sf( 'ne_up'  ).data[index]
		self.dne_dr = sf( 'dne_dr'  ).data[index]

		try:
			self.Te = sf( 'Te'  ).data[index]
			self.Te_lo = sf( 'Te_lo'  ).data[index]
			self.Te_up = sf( 'Te_up'  ).data[index]
			self.dTe_dr = sf( 'dTe_dr'  ).data[index]
			self.maxTe = numpy.max(self.Te)
			self.StatusTe = True
		except:
			print 'Could not read Te'
			self.StatusTe = False


		try:
			self.pe = sf( 'pe'  ).data[index]
			self.dpe_dr = sf( 'dpe_dr'  ).data[index]
			self.StatusPe = True
		except:
			print 'Could not read pressure'
			self.StatusPe = False


		if ECEraw:
			try:
				self.ece_rhop = sf( 'ece_rhop'  ).data[index]   
			#gives segmentation fault 
			#	self.ece_dat = sf( 'ece_dat'  ).data[index]   
			#	self.ece_unc = sf( 'ece_unc'  ).data[index]
				self.ece_mod = sf( 'ece_mod'  ).data[index]
			#	self.ece_resi = sf( 'ece_resi'  ).data[index]
				self.StatusECE = True
			except:
				print 'Could not read ECE raw data'
				self.StatusECE = False



		self.maxne = numpy.max(self.ne)
		self.Shotnumber = Shotnumber
		self.Status = True

		#self.meanProfile()
		sf.close()


	def Unload( self ):

		if self.Status:

			self.Status = False
			self.Status_mean = False
			del self.time
			del self.rhop
			del self.ne
			del self.ne_lo
			del self.ne_up
			del self.dne_dr

			if self.StatusTe:

				del self.Te
				del self.Te_lo
				del self.Te_up	
				del self.dTe_dr
				self.StatusTe = False

			if self.StatusPe:

				del self.pe	
				del self.dpe_dr
				self.Statuspe = False

			if self.StatusECE:
	
				del self.ece_rhop
			#	del self.ece_dat
			#	del self.ece_unc 
				del self.ece_mod 
			#	del self.ece_resi
				self.StatusECE = False

			del self.Shotnumber

		if self.Status_mean:

			del self.meanNe
			del self.meanTe
			del self.meanRho

	def MapToRmaj( self ,  Experiment = 'AUGD', Diagnostic='EQH', Edition=0 ):
		ntimes = 0
		if self.Status:
			ntimes = numpy.size(self.time)
			Rmaj = numpy.zeros( (ntimes,numpy.size(self.rhop[0,:])) )
                        Rsep = numpy.zeros( (ntimes)) 
			try:				
				#for i in range(ntimes):
				#	rz_idx = numpy.argmin( numpy.abs( self.time - timepoints[i] ) )
				#	output = kk.KK().kkrhorz( self.Shotnumber, self.time[rz_idx], self.rhop[rz_idx],angle=0.0,exp=Experiment, diag=Diagnostic, ed=Edition)
                                output = kk.KK().kkrhorz( self.Shotnumber, self.time, self.rhop,angle=0.0,exp=Experiment, diag=Diagnostic, ed=Edition)
				Rmaj[:] = output.r
                                for i in range(ntimes):
                                            Rsep[i]=numpy.interp([1.0],self.rhop[i],Rmaj[i])
			except Exception:
				print 'Error in calc. major Radius'
                                IPython.embed()
                                pass
                        
    			self.rmaj=Rmaj
                        self.Rsep=Rsep              
			return Rmaj
                
				

	def IDAonRz( self , time , R, z , Experiment = 'AUGD', Diagnostic='EQH' ):
		if self.Status:
			idx = numpy.argmin( numpy.abs( self.time - time ) )
			output =initialize_output()
			out = kk.KK().kkrzptfn( self.Shotnumber, self.time[idx], R, z, exp= Experiment, diag=Diagnostic, ed=0)
			f = scipy.interpolate.interp1d( self.rhop[idx,:] , self.ne[idx,:] )
			output.ne = f( out.rho_p )
			f = scipy.interpolate.interp1d( self.rhop[idx,:] , self.ne_lo[idx,:] )
			output.ne_lo = f( out.rho_p )
			f = scipy.interpolate.interp1d( self.rhop[idx,:] , self.ne_up[idx,:] )
			output.ne_up = f( out.rho_p )
			f = scipy.interpolate.interp1d( self.rhop[idx,:] , self.Te[idx,:] )
			output.Te = f( out.rho_p )
			f = scipy.interpolate.interp1d( self.rhop[idx,:] , self.Te_lo[idx,:] )
			output.Te_lo =f( out.rho_p )
			f = scipy.interpolate.interp1d( self.rhop[idx,:] , self.Te_up[idx,:] )				
			output.Te_up = f( out.rho_p )


		#	rhop = numpy.where( output.rho_p < numpy.max( self.rhop[idx,:]) )

			return output 


	def MapToRhop( self , timepoints ):
		ntimes = 0
		if self.Status:
			ntimes = numpy.size(timepoints)
			rhop = numpy.zeros( (ntimes,numpy.size(self.rhop[0,:])) )
			try:
				for i in range(ntimes):
					rz_idx = numpy.argmin( numpy.abs( self.time - timepoints[i] ) )
					rhop[i,:] = self.rhop[rz_idx,:]
			except Exception:
				rz_idx = numpy.argmin(numpy.abs(self.time-timepoints))
				rhop = self.rhop[rz_idx,:]
				pass      			
			return rhop 
	

	def MapToRhot( self , timepoints, Experiment = 'AUGD', Diagnostic='EQH' ):
		ntimes = 0
		if self.Status:
			ntimes = numpy.size(timepoints)
			rhot = numpy.zeros( (ntimes,numpy.size(self.rhop[0,:])) )
			try:
				for i in range(ntimes):
					rz_idx = numpy.argmin( numpy.abs( self.time - timepoints[i] ) )
					rhot[i,:] = kk.KK().kkrhopto(self.Shotnumber,self.time[rz_idx], self.rhop[rz_idx,:]).rho_t 
			except Exception:
				rz_idx = numpy.argmin(numpy.abs(self.time-timepoints))	
				rhot = kk.KK().kkrhopto(self.Shotnumber,self.time[rz_idx], self.rhop[rz_idx,:]).rho_t 
				pass      			
			return rhot 



	def meanProfile(self, tstart = 0.0 ,tend = 10.0 , coordinates = 'rhop' ):
		
		if self.Status:
			if tend > tstart:
				idx = numpy.where( ( self.time > tstart ) & ( self.time < tend ) )
				self.meanRho = numpy.mean(self.rhop,axis=0)
				self.meanNe = numpy.mean(self.ne,axis=0)
				if self.StatusTe:
					self.meanTe = numpy.mean(self.Te,axis=0)
				self.Status_mean=True
				
			else:
				print ''
				return None
				
		else:
			print 'no data loaded'
			return None
	



	def write_ascii( self , tstart = 0.0 ,tend = 10.0 , filename = "IDA_out.dat", key=['ne', 'Te'] ):	
		
		if self.Status:
			out_file = open(filename,"w")	
			start_idx = numpy.argmin( numpy.abs( self.time - tstart ) )
			end_idx = numpy.argmin( numpy.abs( self.time - tend ) )

			out_file.write("# time rhop ne ne_lo ne_up Te Te_lo Te_up \n" )
			for j in range(start_idx,end_idx):
				for i in range(numpy.size(self.rhop[j,:])):
					out_file.write("%.6f " % self.time[j])
					out_file.write("%e " % self.rhop[j,i])
					out_file.write("%e " % self.ne[j,i])
					out_file.write("%e " % self.ne_lo[j,i])
					out_file.write("%e " % self.ne_up[j,i])
					out_file.write("%e " % self.Te[j,i])
					out_file.write("%e " % self.Te_lo[j,i])
					out_file.write("%e " % self.Te_up[j,i])
					out_file.write("\n")
				out_file.write("\n")
				out_file.write("\n")
			out_file.close()
				

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

	def getRmaj(self, timepoint, eqExp = 'AUGD', eqDiag='EQH'):

		if self.Status:		
			if (numpy.size(timepoint)== 1) | ( numpy.size(timepoint) > 2):
				return self.getData( timepoints = timepoint, value='Rmaj',Experiment = eqExp, Diagnostic=eqDiag)

			elif numpy.size(timepoint)== 2:
				return  self.getData( range = timepoint, value='Rmaj',Experiment = eqExp, Diagnostic=eqDiag)

			return None
		else:
			print 'no Data loaded'
			return None

	def getVdia(self, timepoint):

		if self.Status:		
			if (numpy.size(timepoint)== 1) | ( numpy.size(timepoint) > 2):
				ne = self.getData( timepoints = timepoint, value='ne')
				Te = self.getData( timepoints = timepoint, value='Te')
				Rmaj = self.getData( timepoints = timepoint, value='Rmaj')
				vdia = (numpy.gradient(ne)*Te+numpy.gradient(Te)*ne)/ne/numpy.gradient(Rmaj)
				return vdia 

			elif numpy.size(timepoint)== 2:
				ne = self.getData( range = timepoint, value='ne')
				Te = self.getData( range = timepoint, value='Te')
				Rmaj = self.getData( range = timepoint, value='Rmaj')
				vdia = (numpy.gradient(ne)*Te+numpy.gradient(Te)*ne)/ne/numpy.gradient(Rmaj)		  
				return vdia

			return None
		else:
			print 'no Data loaded'
			return None


	def getTime(self, timepoint):

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

			if not ( (value == 'Te') | (value == 'ne') | (value == 'rhop')| (value == 'Rmaj')| (value == 'time') ):
				print 'value must be Te or ne or rhop'
				return None

			if (value == 'rhop'):
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


			if (value == 'Rmaj') :
				outmaj=kk.KK().kkrhorz(self.Shotnumber, self.time[idx],  self.rhop[idx,:], angle=0.0, exp='AUGD', diag='EQH', ed=0)
				return outmaj.r
			    

	
			if (value == 'Te') | (value == 'ne') | (value == 'rhop'):
					#ne/Te/rhop should have the same area base
				
                            if value == 'Te':
                                return  self.Te[idx,:]
				
                            if value == 'ne':
				    return  self.ne[idx,:]
				
                            if value == 'rhop':
				    return  self.rhop[idx,:]

			                                            
                return None
		


#	def __call__( self , time , rhop = None ):


		
