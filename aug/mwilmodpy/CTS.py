import numpy as np
import dd
import scipy.interpolate
import kk
import tools
import PHI
import ElmSync
#import matplotlib.pylab as plt
from IPython import embed

class CTS:
	def __init__( self , Shotnumber = None , Experiment = 'AUGD' ):
		self.Status = False
		if Shotnumber != None:
			self.Load( Shotnumber ,Experiment = Experiment)


	def __del__( self ):
		self.Unload( )
		del self.Status

	def Load( self ,   Shotnumber, Experiment = 'AUGD',Diagnostic = 'CTA', Edition = 0L, tBegin=-1.0, tEnd=12.0 ):	
		self.Unload()
		if Diagnostic == 'CTA':
			self.nChan = 50
			try:
				sf = dd.shotfile( Diagnostic, Shotnumber, Experiment, Edition)
			except:
				print "Error reading shotfile" 
				return False

			time = sf.getTimeBase( 'ch1' )
			index = np.where( ( time > tBegin ) & ( time < tEnd ) )
			output = []

			for i in np.arange(self.nChan):
				output.append(sf( 'ch%d'%(i+1) ).data)

			
			###remove channels
			rem=[22,23,24]

			txtinfo = np.genfromtxt('/afs/ipp-garching.mpg.de/home/m/mwillens/Programme/python/data/CTA/freq.dat').T

			self.chNr = np.delete(txtinfo[0],rem)
			self.freq = np.delete(txtinfo[1],rem)
			self.IF = np.delete(txtinfo[2],rem)
			self.data = np.array(output).T
			self.data = np.delete(((self.data[:100]).mean(axis=0)-self.data)[index],rem,axis=1)
			self.time = time[index]
			self.ntime = np.size(self.time)
			
			del output
			self.Shotnumber = Shotnumber
			self.Status = True
			sf.close()



			return True

	def Unload( self ):
		if self.Status:
			self.Status = False
			del self.time
			del self.ntime
			del self.data
			del self.Shotnumber

	def Binning( self, samplefreq = 10.0 ):				   
		if self.Status:				  
          
			newtime, self.data = tools.dataBinning( self.time, self.data, samplefreq = samplefreq )
        
			self.time = newtime
			self.ntime = np.size(self.time)
