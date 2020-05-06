import numpy as np
import kk_mwillens as kk
import IPython
import matplotlib.pylab as plt
import eqi_map as fastkk
from scipy.ndimage.interpolation import map_coordinates
from scipy import interpolate

path = '/afs/ipp-garching.mpg.de/home/m/mwillens/Programme/VACFIELD/'
path_b = '/vacfield/' 

class vacfield:
	def __init__( self ,  folder='TFripple', fileEnd='_TF' ):

		self.folder = folder
		self.fileEnd = fileEnd
		self.Status = False
		self.Status_readOutput = False

	def __del__( self ):

		self.Unload( )
		self.Status = False
		self.Status_readOutput = False

	def Unload( self ):

		if self.Status_readOutput:
			del B
			if self.nDim == 3:
				del self.Br
				del self.Bz
				del self.Bphi
				del self.Btot

# read single time point
	def readOutput( self ):

		fileOutput = path + self.folder + path_b + 'vacfield' + self.fileEnd

		file = open(fileOutput, 'r')
		arrayDimStr = (file.readline()).split()
		arrayDim = np.zeros_like(arrayDimStr,dtype=int)
		nArr = np.size(arrayDimStr )
		for i in np.arange(nArr):
			arrayDim[i] = int(float(arrayDimStr[i]))
		arrayStr = (file.readline()).split()
		nDim =  int(float(arrayStr[0]))
		nPeriods = int(float(arrayStr[1]))
		nTor = int(float(arrayStr[2]))
		arrayStr = (file.readline()).split()
		Rcenter = float(arrayStr[0])
		zcenter = float(arrayStr[1])
		#half width
		Rhw = float(arrayStr[2])
		arrayStr = (file.readline()).split()
		zhw = float(arrayStr[0])

		B = np.genfromtxt(fileOutput,skip_header = 4)
		
		file.close()

		## ph, z, R  but unsure
		self.B = np.reshape(B,np.append([arrayDim[0],arrayDim[2],arrayDim[1]],[3]))
		
		if nDim == 3:
			self.Bphi = np.squeeze(self.B[:,:,:,0])
			self.Br = np.squeeze(self.B[:,:,:,1])
			self.Bz = np.squeeze(self.B[:,:,:,2])
			self.Btot = np.sqrt( self.Br*self.Br+self.Bz*self.Bz+self.Bphi*self.Bphi)
		
		self.nDim = nDim
		self.phi = np.arange(0,360./nPeriods,(360./nPeriods)/nTor)
		self.nphi = np.size(self.phi)
		self.R = np.linspace(Rcenter-Rhw,Rcenter+Rhw,arrayDim[1])
		self.nR = np.size(self.R)
		self.z = np.linspace(zcenter-zhw,zcenter+zhw,arrayDim[2])
		self.nz = np.size(self.z)
		self.nPeriods = nPeriods

		#a = map_coordinates(Bout[:,:,:,0],[[1.2],[2.3],[41]],mode='nearest',order=1,prefilter=True)

		self.Status_readOutput = True

#	def __call__( self , timepoints , rhop = None, Experiment = 'AUGD', Diagnostic='EQH' ):



## phi in degree!!!!
	def rzphibrzphi( self, rin, zin, phiin):

		if self.Status_readOutput:

			rin = np.squeeze(np.atleast_1d(rin))
			zin = np.squeeze(np.atleast_1d(zin))
			phiin = np.squeeze(np.atleast_1d(phiin))
			
			nr =  np.size(rin)
			nz =  np.size(zin)
			nphi =  np.size(phiin)

			#get number of dimensions
			rdim = np.size(np.shape(rin))
			zdim = np.size(np.shape(zin))
			phidim = np.size(np.shape(phiin))

			## check if has one dimension
			if ((rdim != zdim) | (zdim != phidim)| (rdim != phidim)):
				print 'input must have save dimension'
				return False
			
			if ((rdim != 1) & (rdim != 3)):
				print 'input must have one dimension or 3 dimension'
				return False

			if (rdim == 1):
				nArr = np.array([nr,nz,nphi])
				idxBiggerOne = np.squeeze(np.where( nArr > 1 ))
				idxIsOne = np.squeeze(np.where( nArr == 1 ))
				#case one array is given the other one have single values
				if (( np.size(idxIsOne) == 1 ) & ( np.size(idxBiggerOne) == 2 )):
					if idxBiggerOne == 0:
						rnew = rin
						znew = np.tile(zin,(nrin))
						phinew = np.tile(phiin,(nrin))
					elif idxBiggerOne == 1:
						rnew = np.tile(rin,(nzin))
						znew = zin 
						phinew = np.tile(phiin,(nzin))
					elif idxBiggerOne == 2:
						rnew = np.tile(rin,(nphiin))
						znew = np.tile(zin,(nrin))
						phinew = phiin
					else:
						print 'something went wrong'
						return False	
				#case all have more values but the they have to have the same size
				elif np.size(idxBiggerOne) == 3:
					if (nr !=nz) | (nz != nphi) | (nr != nphi):
						print 'if all are do have single values all must heve the same size'
						return False
					rnew = rin
					znew = zin
					phinew = phiin
				elif np.size(idxIsOne) == 3:
					rnew = rin
					znew = zin
					phinew = phiin
				else:
					print 'no other cases are covered'
					return False
			
			if (rdim == 3):
				if( (np.shape(rin) != np.shape(zin)) | (np.shape(phiin) != np.shape(zin))):
					print 'shape of input is not the same'
					return False	

				#store input shape
				inputShape = np.shape(rin)
				#transform shape into one dimension
				rnew = np.ravel(rin)
				znew = np.ravel(zin)
				phinew = np.ravel(phiin)			
				
			#get functions for index in R,z and phi
			fR = interpolate.interp1d(self.R, np.arange(self.nR))
			#fz = interpolate.interp1d(self.z, np.arange(self.nz))	
			fphi = interpolate.interp1d(self.phi, np.arange(self.nphi))
		
			phinew =  np.remainder(phinew,360./self.nPeriods)

			phiMap = fphi(phinew)
			RMap = fR(rnew)
			zMap = np.interp1d(znew,self.z, np.arange(self.nz))	
			
			idxInput = np.array([phiMap,zMap,RMap]) 

			Bout = []
			#IPython.embed()
			#Bout = map_coordinates(np.squeeze(self.B[:,:,:,0],[[fR],[fz],[fPhi]],order=1,prefilter=True)
			for i in np.arange(0,self.nDim):
				Bout.append(map_coordinates(np.squeeze(self.B[:,:,:,i]),idxInput,mode='nearest',prefilter=True))

			Bout = np.squeeze(np.array(Bout)).T	

	#map_coordinates(np.squeeze(self.B[:,:,:,0]),idxInput,order=1,prefilter=True)
			if (rdim == 3):
				rnew = np.reshape(rnew,inputShape)
				znew = np.reshape(znew,inputShape)
				phinew = np.reshape(phinew,inputShape)
				Bout = np.reshape(Bout,np.append(inputShape, np.size(Bout[0,:])) )

			return rnew,znew,phinew,Bout
				
	 
		#field = kk.KK().kkrzBrzt( shot, time, surf_r_sym, surf_z_sym )


		#else:
		#	print 'tschau' 

