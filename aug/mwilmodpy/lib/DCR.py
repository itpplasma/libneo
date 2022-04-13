import numpy
import dd
import kk
import scipy.interpolate
import matplotlib.pylab as plt
import IPython


class initialize_output:
	status = False

class DCR:
	def __init__( self ,  Experiment = 'AUGD', Shotnumber = None ):
		self.Status = False
	
		if Shotnumber != None:
			self.Load( Shotnumber )

	def __del__( self ):
		self.Unload( )
		del self.Status

	def Load( self ,   Shotnumber, Experiment = 'AUGD', Edition = 0L, tBegin=-1.0, tEnd=12.0 ):
		self.Unload()
		try:
			print Shotnumber
			sf = dd.shotfile( 'DCR', Shotnumber, Experiment, Edition)
		except:
			print "Error reading shotfile" 
			return
	
		print tBegin, tEnd
		self.time = sf( 'profile' ).time
		index = numpy.where( ( self.time > tBegin ) & ( self.time < tEnd ) )
		self.time = self.time[index]
		ntimes = numpy.size(self.time)
		self.rhop = numpy.arange(0.0,1.01,0.01)
		nrhop = numpy.size(self.rhop)
		self.coeff = sf( 'profile'  ).data[index]
		ncoeff = numpy.size(self.coeff[0,:])
		coeff =  numpy.tile(self.coeff.transpose(),nrhop).reshape(5,ntimes,nrhop)
		rhop = numpy.tile(self.rhop,(ntimes)).reshape(ntimes,nrhop)
		self.ne = coeff[0]+rhop*coeff[1]+rhop*rhop*coeff[2]+rhop*rhop*rhop*coeff[3]+rhop*rhop*rhop*rhop*coeff[4]
		self.Shotnumber = Shotnumber
		self.Status = True

		sf.close()

	def Unload( self ):
		if self.Status:
			self.Status = False
			self.Status_mean = False
			del self.time
			del self.rhop
			del self.coeff
			del self.ne
			del self.Shotnumber


	def __call__( self , timepoints ):
		ntimes = 0
		if self.Status:
			ntimes = numpy.size(timepoints)			
			ne_temp = numpy.zeros((ntimes,numpy.size(self.ne[0,:])))
			try:
				for i in range(ntimes):
					idx = numpy.argmin( numpy.abs( self.time - timepoints[i] ) )
					ne_temp[i,:] = self.ne[idx,:]
			except Exception:
				idx = numpy.argmin(numpy.abs(self.time-timepoints))
				output = self.ne[idx,:]
				pass
			
			if ntimes > 1: 	
				return ne_temp
			else: 
				return output

		
				

#	def __call__( self , time , rhop = None ):


		
