
import dd as dd


class MBI:
	def __init__( self ,  Experiment = 'AUGD', Diagnostic = 'MBI', Shotnumber = None ):
		self.Status = False
		if Shotnumber != None :
			self.Load( Shotnumber )
		

	def __del__( self ):
		self.Unload( )
		del self.Status

	def Load( self ,  Shotnumber, Experiment, Diagnostic, Edition = 0L ):
		self.Unload()
		if Diagnostic == 'MBI':

			try:
				sf = dd.shotfile( Diagnostic, Shotnumber, Experiment, Edition)
			except:
				print "Error reading shotfile" 
				return

			self.time = sf.getTimeBase( 'BTF' )
			self.BTF = sf.getSignal( 'BTF' )
			self.ICoIo = sf.getSignal( 'ICoIo')
			self.ICoIu = sf.getSignal( 'ICoIu')
			self.IHO2od= sf.getSignal( 'dIOH2s' )
			self.IHO2ud= sf.getSignal( 'dIOH2u' )
			self.IOH= sf.getSignal( 'IOH' )
			self.IV1o= sf.getSignal( 'IV1o' )
			self.IV1u= sf.getSignal( 'IV1u' )
			self.IV2o= sf.getSignal( 'IV2o' )
			self.IV2u= sf.getSignal( 'IV2u' )		
			self.IV3o= sf.getSignal( 'IV3o' )
			self.IV3u= sf.getSignal( 'IV3u' )

			sf.close()
		else:
			print 'Diagnotic should be MBI, but is:' ,Diagnostic


	def Unload( self ):
		if self.Status:
			self.Status = False
			del self.time	
			del self.BTF	




