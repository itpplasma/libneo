import numpy
import dd
import scipy.interpolate
import kk
import tools
import PHI
import ElmSync
#import matplotlib.pylab as plt
from IPython import embed

class LIN:
	def __init__( self ,  Experiment = 'AUGD', Shotnumber = None ):
		self.Status = False
		if Shotnumber != None:
			self.Load( Shotnumber ,Experiment = Experiment)

	def __del__( self ):
		self.Unload( )
		del self.Status

	def Load( self ,   Shotnumber, Experiment = 'AUGD',Diagnostic = 'LIN', Edition = 0L, tBegin=-1.0, tEnd=12.0 , Rshift=0.00):	
		self.Unload()
		if Diagnostic == 'LIN':
			try:
				sf = dd.shotfile( Diagnostic, Shotnumber, Experiment, Edition)
			except:
				print "Error reading shotfile" 
				return False

			time = sf.getTimeBase( 'ne' )
			index = numpy.where( ( time > tBegin ) & ( time < tEnd ) )
			output = sf( 'ne' )
			self.Rshift = Rshift
			self.ne = output.data[index]
			self.ne_unc = sf( 'ne_unc' ).data[index]
			self.time = output.time[index]
			self.ntime = numpy.size(self.time)
			self.rhop = output.area[index]
			del output
			self.R = numpy.squeeze(sf( 'R' ).data)+Rshift
			self.nChannels = numpy.size(self.R)
			self.Rall =  numpy.reshape(numpy.tile( self.R,(self.ntime) ) ,(self.ntime,self.nChannels) )
			self.z = numpy.squeeze(sf( 'Z' ).data)
			self.zall = numpy.reshape(numpy.tile( self.z,(self.ntime) ) ,(self.ntime,self.nChannels) )
			self.x = numpy.squeeze(sf( 'x' ).data)
			self.xall = numpy.reshape(numpy.tile( self.x,(self.ntime) ) ,(self.ntime,self.nChannels) )
			self.phi = PHI.LIB() #260.64
		#	self.phiRad = self.phi*numpy.pi/180.
			#self.lib_dat=numpy.squeeze(sf( 'lib_dat' ).data)[index]
			self.Shotnumber = Shotnumber
			self.Status = True
			sf.close()
			return True

	def Unload( self ):
		if self.Status:
			self.Status = False
			del self.time
			del self.ntime
			del self.rhop
			del self.ne
			del self.ne_unc
			#del self.lib_dat
			del self.nChannels
			del self.R
			del self.Rall
			del self.z
			del self.Shotnumber

	def ElmSync( self , Experiment = 'AUGD'):
		if self.Status:
			timeidx = ElmSync.ElmExtract(self.time, self.Shotnumber, preFac=0.1,postFac=0.25)
			self.time = self.time[timeidx]
			self.ntime = numpy.size(self.time)
			self.ne = self.ne[timeidx,:]
			self.rhop = self.rhop[timeidx,:]	
			self.Rall = self.Rall[timeidx,:]	
			self.zall = self.zall[timeidx,:]

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

	def getne( self , timepoints ):
		if self.Status:    			
			return  self.getData(timepoints, selectData='ne' )

	def getR( self , timepoints ):
		if self.Status:    			
			return  self.getData(timepoints, selectData='R' )

	def getz( self , timepoints ):
		if self.Status:    			
			return  self.getData(timepoints, selectData='z' )

	def getData( self , timepoints, selectData='ne' ):
		ntimes = 0
		if self.Status:
			ntimes = numpy.size(timepoints)
			data = numpy.zeros( (ntimes,self.nChannels)) 
			if ntimes > 1:
				for i in range(ntimes):
					
					idx = numpy.argmin( numpy.abs( self.time - timepoints[i] ) )
		
					if selectData == 'ne':
						data[i,:] = self.ne[idx,:]
					if selectData == 'R':
						data[i,:] = self.Rall[idx,:]
					if selectData == 'z':
						data[i,:] = self.zall[idx,:]
				
			else:
				idx = numpy.argmin(numpy.abs(self.time-timepoints))
				if selectData == 'ne':
					data = self.ne[idx,:]
				if selectData == 'R':
					data = self.Rall[idx,:]
				if selectData == 'z':
					data = self.zall[idx,:]
			     			
			return data 



	def getRhop( self , timepoints ):
		if self.Status:
			return self.MapToRhop( timepoints )


	def MapToRhot( self , timepoints,  eqExp = 'AUGD', eqDiag='EQH' ):		
		if self.Status:
			ntimes = numpy.size(timepoints)
			rhot = numpy.zeros( (ntimes,numpy.size(self.ne[0,:])) )
			try:
				for i in range(ntimes):
					idx = numpy.argmin( numpy.abs( self.time - timepoints[i] ) )
					output = kk.KK().kkrzptfn( self.Shotnumber, self.time[idx[i]], self.R[idx[i],:], self.z[idx[i],:], exp= Experiment, diag=Diagnostic, ed=0)
                                        rhot[i,:] = output.rho_t                                       
                                        return output.rho_t				
			except Exception:
				idx = numpy.argmin(numpy.abs(self.time-timepoints))
				##R and z with zero because of constant R and z
				output = kk.KK().kkrzptfn( self.Shotnumber, self.time[idx], self.R[0,:], self.z[0,:], exp= eqExp, diag=eqDiag)
				rhot = output.rho_t
				pass      			
			return rhot 			


	def getRsep( self ):

		if self.Status:
			
			LIB_Rsep = self.getTimetrace( selectData_in = 'rhop', data_in = [1.0], selectData_out = 'R' )	
			self.R_sep = self.Rall - numpy.reshape(numpy.repeat(LIB_Rsep,self.nChannels),(self.ntime,self.nChannels))
			return LIB_Rsep


			
## only one value for ne, get R and Z
	def getTimetrace( self , time_in = None, selectData_in = 'ne', data_in = [1.5e19], selectData_out = 'R' ,kind = '1D'):

		if self.Status:

			#if no time is given, every point will be used to calc timetrace
			if (numpy.all(time_in) == None) :
				time_in = self.time
			
			if selectData_in == 'ne':
				x = self.getne(time_in)
			if (selectData_in == 'Rhop')| (selectData_in == 'rhop'):
				x = self.getRhop(time_in)
			

			if selectData_out == 'R':
				y = self.getR(time_in)
			if selectData_out == 'z':
				y = self.getz(time_in)
			if (selectData_out == 'Rhop') | (selectData_out == 'rhop'):
				y = self.getRhop(time_in)
			if selectData_out == 'ne':
				y = self.getne(time_in)

			ntimes = numpy.size(time_in)
			data = numpy.zeros((ntimes,numpy.size(data_in)))
			#embed()
			if kind == '1D':
				for i in numpy.arange(ntimes):
					idxSort=numpy.argsort(x[i])
					data[i,:] =numpy.interp(data_in,x[i][idxSort],y[i][idxSort])	

				return numpy.squeeze(data)
				
			elif kind == '2D':
				time = numpy.reshape( numpy.repeat( self.time, ( self.nChannels ) ) ,( self.ntime, self.nChannels) )
				f =  scipy.interpolate.interp2d(time, x, y, kind='linear', fill_value = float('NaN'), bounds_error=False)
				return f(time_in,data_in)



	def __call__( self , time , rhop = None, eqExp = 'AUGD', eqDiag='EQH' ):
		if self.Status:
			if rhop == None:
				return self.getne( time )
			else:
				ntimes = numpy.size(time)
				#get the time points from the get routines
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
