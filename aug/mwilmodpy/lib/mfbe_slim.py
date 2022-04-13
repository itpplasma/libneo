import numpy as np
import IPython
import matplotlib.pylab as plt
from scipy.ndimage.interpolation import map_coordinates
from scipy.interpolate import griddata
from scipy import interpolate

import BinaryReader

path = '/afs/ipp-garching.mpg.de/home/m/mwillens/Programme/VACFIELD/'
path_b = '/vacfield/' 
"""
/afs/ipp/u/ers/all/users/wls/30839/vacfield/vacfield_m_pert_n2_even
>>>
>>> Axialsymmetrisches Gleichgewichtsfeld
>>> /afs/ipp-garching.mpg.de/home/e/ers/all/users/wls/30839/mfbe/field_axi_01a 
>>>
>>> 3D Gleichgewichtsfeld
>>> /afs/ipp-garching.mpg.de/home/e/ers/all/users/wls/30839/mfbe/
>>> field_3d_rmp_n2_b_even_a
"""


class mfbe:
	def __init__( self ,  folder='/afs/ipp-garching.mpg.de/home/e/ers/all/users/wls/30839/mfbe/', filename='field_3d_rmp_n2_b_even_a' ):

		self.folder = folder
		self.filename = filename		
		self.StatusRead = False
		self.StatusBn = False


	def __del__( self ):

		self.Unload( )
		del self.StatusRead 
		self.StatusSurf = False
		del self.StatusBn
	

	def Unload( self ):

		if self.StatusRead:
			del self.B

			if self.nDim == 3:
				del self.Br
				del self.Bz
				del self.Bphi

			if self.StatusBn :
				del self.Bn
				del self.vecNorm
				self.StatusBn = False

			del self.Btot
			del self.phi
			del self.nphi
			del self.R 
			del self.nR 
			del self.z 
			del self.nz 
			del self.nPeriods
			del self.folder
			del self.filename


			self.StatusRead = False

# read single time point
	def read3D( self , binary=True):
		
		fileOutput = self.folder + self.filename
		print 'reading: ' + self.filename+' in '+self.folder
		if binary:
			br = BinaryReader.BinaryReader(fileOutput)

			bl1 = br.read('int32')
			if bl1 != 48:
				print "Hmmmh, record length %d !=48. Maybe not a binary mfbe file\n"%bl1
				br.close()
				return

			dims = br.read('double',3)
			comps = br.read('double',3)
		#enr = dims(3)
		
		#enz = br.read('double')

			bl2 = br.read('int32') 
			if bl2 != bl1: 
				print "Hmmmh, record length %d"%bl1+"!=%d. Maybe not a binary mfbe file\n"%bl2 
				br.close()
				return

			bl1 = br.read('int32') 

			if bl1!=32: 
				print "Hmmmh, record length %d1=32. Maybe not a binary mfbe file\n"%bl1
				br.close()
				return
		
			grid1 = br.read('double',4)
		
			bl2 = br.read('int32') 
			if bl2 != bl1: 
				print "Hmmmh, record length %d"%bl1+"!=%d. Maybe not a binary mfbe file\n"%bl2 
				br.close()
				return		

			rc = comps[0]*dims[1]*dims[2]
			rc8=rc*8

			y = np.zeros((int(comps[2]),int(rc)))

			for i in range(int(comps[2])):
			
				bl1= br.read('int32')
				if bl1!=rc8: 
					print "Error, record length %d"%bl1+"!=%d.\n"%rc8
					br.close()
					return
				y1 =br.read('double',rc)
				bl2= br.read('int32')
				if bl2 != bl1: 
					print "Inconsistent record length %d"%bl1+"!=<>%d.\n"%bl2 
					br.close()
					return

				y[i,:]=np.array(y1)


			br.close()
 
			enr=dims[1]
			enz=dims[2]
			nDim=comps[0]
			nPeriods=comps[1]
			nTor=comps[2]
			Rcenter = grid1[0]
			zcenter = grid1[1]
			Rhw = grid1[2]
			zhw = grid1[3]
			
		else:
# in the case of ascci files
			file = open(fileOutput, 'r')
			arrayDimStr = (file.readline()).split()
			dims = np.zeros_like(arrayDimStr,dtype=int)
			nArr = np.size(arrayDimStr )
			for i in np.arange(nArr):
				dims[i] = int(float(arrayDimStr[i]))
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
		self.B = np.reshape(y,np.append([int(dims[0]),int(dims[2]),int(dims[1])],[int(nDim)]))

		
		self.Bphi = np.squeeze(self.B[:,:,:,0])
		self.Br = np.squeeze(self.B[:,:,:,1])
		self.Bz = np.squeeze(self.B[:,:,:,2])
		self.Btot = np.sqrt( self.Br*self.Br+self.Bz*self.Bz+self.Bphi*self.Bphi)
	#	self.Vec
#calculate Bn
		#IPython.embed()
		self.nDim = nDim

		if self.nDim > 3:
			self.dBphidphi = np.squeeze(self.B[:,:,:,3])
			self.dBphidr = np.squeeze(self.B[:,:,:,4])
			self.dBphidz = np.squeeze(self.B[:,:,:,5])
			self.dBrdphi = np.squeeze(self.B[:,:,:,6])
			self.dBrdr = np.squeeze(self.B[:,:,:,7])
			self.dBrdz = np.squeeze(self.B[:,:,:,8])
			self.dBzdphi = np.squeeze(self.B[:,:,:,9])
			self.dBzdr = np.squeeze(self.B[:,:,:,10])
			self.dBzdz = np.squeeze(self.B[:,:,:,11])

		self.phi = np.arange(0,360./nPeriods,(360./nPeriods)/nTor)
		self.phirad = np.deg2rad(self.phi)
		self.nphi = np.size(self.phi)
		self.phiradRange= (np.diff(self.phirad).mean())*self.nphi
		self.phiRange= (np.diff(self.phi).mean())*self.nphi
		self.R = np.linspace(Rcenter-Rhw,Rcenter+Rhw,dims[1],endpoint=False)
		self.nR = np.size(self.R)
		self.z = np.linspace(zcenter-zhw,zcenter+zhw,dims[2],endpoint=False)
		self.nz = np.size(self.z)
		self.nPeriods = nPeriods
		print self.filename+' is read'
		self.StatusRead = True

	
	def calcBn(self, vecNormIn=None):
		if self.StatusRead:
			B3=((self.B.T)[:3]).T
			#if vecNorm not given calc from Br and Bz
			if np.all(vecNormIn) == None:
										
						#calculate Bn
				IBrzI = np.sqrt( self.Br**2. + self.Bz**2. )
				nuller = np.zeros_like(IBrzI.T)
				einser = np.ones_like(IBrzI.T)
				vecTang= np.array([nuller, (self.Br/IBrzI).T, (self.Bz/IBrzI).T]).T
				vecPhi = np.array([ einser,nuller,nuller] ).T
				self.vecNorm = np.cross(vecPhi, vecTang)
				self.Bn = np.sum(B3*self.vecNorm,axis=3) 
				self.StatusBn = True

			elif np.shape(vecNormIn) == np.shape(B3):
				self.Bn = np.sum(B3*vecNormIn,axis=3) 
				self.vecNorm = vecNormIn
				self.StatusBn = True
			else:
				print 'shape must be equal'
				self.StatusBn = False

	def addBn(self, BnIn=None):
		if np.all(BnIn) == None:
			print 'no Bn to add'
			return
		if self.StatusBn == False:
			print 'Bn not calculated'
			return
		if np.shape(self.Bn) == np.shape(BnIn):
			self.Bn = self.Bn + BnIn
		else:
			print 'shape must be equal'
			
##same as addBn but 
	def addBrzphi(self, BrIn=None, BzIn=None, BphiIn=None):
		if np.all(BrIn) == None:
			print 'no Bn to add'
			return
		if np.all(BzIn) == None:
			print 'no Bn to add'
			return
		if np.all(BphiIn) == None:
			print 'no Bn to add'
			return

		if np.shape(self.Br) != np.shape(BrIn):
			print 'Br shape must be equal'

		if np.shape(self.Bz) != np.shape(BzIn):
			print 'Bz shape must be equal'

		if np.shape(self.Bphi) != np.shape(BphiIn):
			print 'Bphi shape must be equal'

		self.Br = self.Br + BrIn
		self.Bz = self.Bz + BzIn
		self.Bphi = self.Bphi + BphiIn


		

### function to replicate symmetric equilibria
	def replicate(self, ntor=180):
		
		if self.StatusRead:
			
			if self.nphi != 1:
				print 'toroidal dimension must be one'	
				return
		
			dataShape = np.shape(self.B)
			arrDataShape = np.array( dataShape )
			mask = np.zeros_like(arrDataShape,dtype=bool)
			mask[0]=True 
			np.putmask(arrDataShape,mask,ntor)
			np.putmask(arrDataShape,np.invert(mask),1)
			self.B = np.tile(self.B,arrDataShape)
			
			self.Bphi = self.B[:,:,:,0]
			self.Br = self.B[:,:,:,1]
			self.Bz = self.B[:,:,:,2]
			self.Btot = np.sqrt( self.Br*self.Br+self.Bz*self.Bz+self.Bphi*self.Bphi)
			if self.StatusBn:
				vecNorm = np.tile(self.vecNorm,arrDataShape)
				self.Bn = np.sum(self.B*vecNorm,axis=3)
				self.vecNorm = vecNorm 

	def addField(self, Bin):
		if self.StatusRead:

			if np.shape(Bin) != np.shape(self.B):
				print 'both magnetic fields must have same shape'
				return

			self.B = self.B + Bin
			if self.nDim == 3:
				self.Bphi = np.squeeze(self.B[:,:,:,0])
				self.Br = np.squeeze(self.B[:,:,:,1])
				self.Bz = np.squeeze(self.B[:,:,:,2])
				self.Btot = np.sqrt( self.Br*self.Br+self.Bz*self.Bz+self.Bphi*self.Bphi)



	def getBrzphi(self, Rin=None, zin=None, phiIn=None,phiInRad=True):
		if (self.StatusRead ):
			if ((np.all(Rin) == None) | (np.all(zin)) == None):
				print 'no input of Rin or zin'
				return
			if (np.size(Rin) != np.size(zin) ):
				print 'Rin or zin must have same size'
				return

			Rshape = np.shape(Rin)

			Rin = Rin.ravel() 
			zin = zin.ravel()
			nRin = np.size(Rin)
			fR = interpolate.interp1d(self.R,np.arange(self.nR))
			fz = interpolate.interp1d(self.z, np.arange(self.nz))

			if ((np.all(phiIn) == None)):
				if len(np.shape(self.Br)) == 3:
				## use the given phi
					nphi = np.size( self.Br[:,0,0] )
					phimap = np.repeat( np.arange(nphi),Rin.size)
					usePhiExt=False
				else:
					nphi = 1
					phimap = np.repeat( np.arange(nphi),Rin.size)
					usePhiExt=False

			else:
				#print 'so far not available'
				#return

				phiIn = phiIn.ravel()
				#input in radien
				if phiInRad:
					phiInput = np.remainder(phiIn,self.phiradRange)
					phiInput[phiInput<0] +=self.phiradRange
					if phiInput.max() > self.phirad[-1]:
						fphi = interpolate.interp1d(np.append(self.phirad,self.phiradRange), np.arange(self.nphi+1))
						usePhiExt = True
					else:
						fphi = interpolate.interp1d(self.phirad, np.arange(self.nphi))
						usePhiExt = False

				else:
					phiInput = np.remainder(phiIn,self.phiRange)
					phiInput[phiInput<0] += self.phiRange
					if phiInput.max() > self.phi[-1]:
						fphi = interpolate.interp1d(np.append(self.phi,self.phiRange), np.arange(self.nphi+1))
						usePhiExt = True
					else:
						fphi = interpolate.interp1d(self.phi, np.arange(self.nphi))
						usePhiExt = False

				nphi = np.size( phiInput )	
				phimap = np.repeat( fphi(phiInput),Rin.size)

			RR = np.tile(fR(Rin),nphi)
			zz = np.tile(fz(zin),nphi)
			inputArr = np.transpose(np.squeeze(np.dstack((phimap, zz, RR))))

			if usePhiExt:
				BrIn = np.concatenate((self.Br,self.Br[:1]),axis=0)
				BzIn = np.concatenate((self.Bz,self.Bz[:1]),axis=0)
				BphiIn = np.concatenate((self.Bphi,self.Bphi[:1]),axis=0)
			else:
				BrIn = self.Br
				BzIn = self.Bz
				BphiIn = self.Bphi

			BrSurf = np.reshape( map_coordinates(BrIn, inputArr  ), (np.append(nphi,Rshape)) )
			BzSurf = np.reshape( map_coordinates(BzIn, inputArr  ), (np.append(nphi,Rshape)) )
			BphiSurf = np.reshape( map_coordinates(BphiIn, inputArr  ), (np.append(nphi,Rshape)) )


			return BrSurf,BzSurf,BphiSurf

##
	def BnSingle(self, Rin=None, zin=None, phiIn=None,vecIn=None, calcAll=False,phiInRad=True):
			
		if (self.StatusRead ):
			
			if ((self.StatusBn==False) & (np.all(vecIn) == None)):
				print 'vector input must be available or Bn calcuted'   
				return		   

			if ((np.all(Rin) == None) | (np.all(zin)) == None):
				print 'no input of Rin or zin'
				return
			if (np.size(Rin) != np.size(zin) ):
				print 'Rin or zin must have same size'
				return
			
			Rshape = np.shape(Rin)

			Rin = Rin.ravel() 
			zin = zin.ravel()
			nRin = np.size(Rin)

	    
			useVec = False
			if (np.all(vecIn) != None):
				if ( ( np.shape(vecIn)[0]!=2 ) | ( Rshape != np.shape(vecIn)[1:] ) ):
					print 'First dimension must have 2 or eqal RIN'
					#IPython.embed()
					return
	
				useVec = True
					

			

			fR = interpolate.interp1d(self.R,np.arange(self.nR))
			fz = interpolate.interp1d(self.z, np.arange(self.nz))
			
			if ((np.all(phiIn) == None)):
				
				if len(np.shape(self.Br)) == 3:
				## use the given phi
					nphi = np.size( self.Br[:,0,0] )
					phimap = np.repeat( np.arange(nphi),Rin.size)
					usePhiExt=False
				else:
					nphi = 1
					phimap = np.repeat( np.arange(nphi),Rin.size)
					usePhiExt=False
			else:
				#print 'so far not available'
				#return

				phiIn = phiIn.ravel()
				#input in radien
				if phiInRad:
					phiInput = np.remainder(phiIn,self.phiradRange)
					phiInput[phiInput<0] +=self.phiradRange
					if phiInput.max() > self.phirad[-1]:
						fphi = interpolate.interp1d(np.append(self.phirad,self.phiradRange), np.arange(self.nphi+1))
						usePhiExt = True
					else:
						fphi = interpolate.interp1d(self.phirad, np.arange(self.nphi))
						usePhiExt = False

				else:
					phiInput = np.remainder(phiIn,self.phiRange)
					phiInput[phiInput<0] += self.phiRange
					if phiInput.max() > self.phi[-1]:
						fphi = interpolate.interp1d(np.append(self.phi,self.phiRange), np.arange(self.nphi+1))
						usePhiExt = True
					else:
						fphi = interpolate.interp1d(self.phi, np.arange(self.nphi))
						usePhiExt = False

				nphi = np.size( phiInput )	
				phimap = np.repeat( fphi(phiInput),Rin.size)


			RR = np.tile(fR(Rin),nphi)
			zz = np.tile(fz(zin),nphi)
			inputArr = np.transpose(np.squeeze(np.dstack((phimap, zz, RR))))
			
			
			if ((useVec == False) & (self.StatusBn)):
				if usePhiExt:
					BnIn = np.concatenate((self.Bn,self.Bn[:1]),axis=0)
				else:
					BnIn = self.Bn
					
				BnSurf = np.reshape( map_coordinates(BnIn, inputArr  ), (np.append(nphi,Rshape)) )
			else:
				vecR = np.reshape(np.tile(vecIn[0].ravel(),nphi), (np.append(nphi,Rshape)) )
				vecz = np.reshape(np.tile(vecIn[1].ravel(),nphi), (np.append(nphi,Rshape)))
				if usePhiExt:
					BrIn = np.concatenate((self.Br,self.Br[:1]),axis=0)
					BzIn = np.concatenate((self.Bz,self.Bz[:1]),axis=0)
				else:
					BrIn = self.Br
					BzIn = self.Bz
				#IPython.embed()
				BrSurf = np.reshape( map_coordinates(BrIn, inputArr  ), (np.append(nphi,Rshape)) )
				BzSurf = np.reshape( map_coordinates(BzIn, inputArr  ), (np.append(nphi,Rshape)) )
				#IPython.embed()
				BnSurf = BrSurf*vecR + BzSurf*vecz
				useVec = True

			if calcAll == False:
				return BnSurf
			else:
				if (useVec == False) :
					if usePhiExt:
						vecIn = np.concatenate((self.vecNorm,self.vecNorm[:1]),axis=0)
					else:
						vecIn = self.vecNorm

					vecRSurf = np.reshape( map_coordinates(((vecIn.T)[1]).T, inputArr  ), (np.append(nphi,Rshape)) )
					veczSurf = np.reshape( map_coordinates(((vecIn.T)[2]).T, inputArr  ), (np.append(nphi,Rshape)) )
				else:
					vecRSurf = vecR
					veczSurf = vecz
				
				if usePhiExt:
					BrIn = np.concatenate((self.Br,self.Br[:1]),axis=0)
					BzIn = np.concatenate((self.Bz,self.Bz[:1]),axis=0)
					BphiIn = np.concatenate((self.Bphi,self.Bphi[:1]),axis=0)
				else:
					BrIn = self.Br
					BzIn = self.Bz
					BphiIn = self.Bphi

				BrSurf = np.reshape( map_coordinates(BrIn, inputArr  ), (np.append(nphi,Rshape)) )
				BzSurf = np.reshape( map_coordinates(BzIn, inputArr  ), (np.append(nphi,Rshape)) )
				BphiSurf = np.reshape( map_coordinates(BphiIn, inputArr  ), (np.append(nphi,Rshape)) )
				return BnSurf,BrSurf,BzSurf,BphiSurf,vecRSurf,veczSurf


	def getKappa(self, Rin=None, zin=None, phiIn=None, phiInRad=True):
		if (self.StatusRead ) & (self.nDim > 3):

			if (np.size(Rin) != np.size(zin) ) | (np.size(Rin) != np.size(phiIn) ):
				print 'Rin or zin must have same size'
				return

			Br,Bz,Bphi,dBphidphi,dBphidr,dBphidz,dBrdphi,dBrdr,dBrdz,dBzdphi,dBzdr,dBzdz =  self.getAllBn(Rin,zin,phiIn=phiIn,phiInRad=phiInRad)
			IBI = np.sqrt(Br*Br+Bz*Bz + Bphi*Bphi)

			r_vec = Br/IBI * ( dBrdr / IBI - Br * Br/IBI**3.) + Bphi/IBI / Rin *  ( dBrdphi / IBI - Br * Bphi/IBI**3.)  + Bz/IBI * ( dBrdz / IBI - Br * Bz/IBI**3.) - Bphi*Bphi/Rin /IBI/IBI
			phi_vec = Br/ IBI  * ( dBphidr / IBI - Bphi * Br/IBI**3.) + Bphi/Rin/IBI*( dBphidphi / IBI - Bphi * Bphi/IBI**3.)  +Bz*( dBphidz / IBI - Bphi * Bz/IBI**3.) + Bphi*Br/Rin/IBI/IBI
			z_vec = Br/IBI * ( dBzdr / IBI - Bz * Br/IBI**3.) + Bphi/IBI/Rin*( dBzdphi / IBI - Bz * Bphi/IBI**3.)+Bz*( dBzdz / IBI - Bz * Bz/IBI**3.)
			return np.array([r_vec,z_vec,phi_vec])
 
	def getAllBn(self, Rin=None, zin=None, phiIn=None, phiInRad=True):
			
		if (self.StatusRead ):
			
	   

			if ((np.all(Rin) == None) | (np.all(zin)) == None):
				print 'no input of Rin or zin'
				return
			if (np.size(Rin) != np.size(zin) ) | (np.size(Rin) != np.size(phiIn) ):
				print 'Rin or zin must have same size'
				return
			
			Rshape = np.shape(Rin)

			Rin = Rin.ravel() 
			zin = zin.ravel()
			phiIn = phiIn.ravel()
			

			fR = interpolate.interp1d(self.R,np.arange(self.nR))
			fz = interpolate.interp1d(self.z, np.arange(self.nz))
				#input in radien
			if phiInRad:
				phiInput = np.remainder(phiIn,self.phiradRange)
				phiInput[phiInput<0] +=self.phiradRange
				fphi = interpolate.interp1d(np.append(self.phirad,self.phiradRange), np.arange(self.nphi+1))	
			else:
				phiInput = np.remainder(phiIn,self.phiRange)
				phiInput[phiInput<0] += self.phiRange
				fphi = interpolate.interp1d(np.append(self.phi,self.phiRange), np.arange(self.nphi+1))
	
			try:

				RR = fR(Rin)
				zz = fz(zin)
				phimap = fphi(phiInput)
			except:
				print 'error in interpolation'
				IPython.embed()

			inputArr = np.transpose(np.squeeze(np.dstack((phimap, zz, RR))))

			BrIn = np.concatenate((self.Br,self.Br[:1]),axis=0)
			BzIn = np.concatenate((self.Bz,self.Bz[:1]),axis=0)
			BphiIn = np.concatenate((self.Bphi,self.Bphi[:1]),axis=0)

			if self.nDim > 3:
				dBphidphiIn = np.concatenate((self.dBphidphi,self.dBphidphi[:1]),axis=0)
				dBphidrIn = np.concatenate((self.dBphidr,self.dBphidr[:1]),axis=0)
				dBphidzIn = np.concatenate((self.dBphidz,self.dBphidz[:1]),axis=0)
				dBrdphiIn = np.concatenate((self.dBrdphi,self.dBrdphi[:1]),axis=0)
				dBrdrIn = np.concatenate((self.dBrdr,self.dBrdr[:1]),axis=0)
				dBrdzIn = np.concatenate((self.dBrdz,self.dBrdz[:1]),axis=0)
				dBzdphiIn = np.concatenate((self.dBzdphi,self.dBzdphi[:1]),axis=0)
				dBzdrIn = np.concatenate((self.dBzdr,self.dBzdr[:1]),axis=0)
				dBzdzIn = np.concatenate((self.dBzdz,self.dBzdz[:1]),axis=0)		
	

			BrSurf = np.reshape( map_coordinates(BrIn, inputArr  ), Rshape )
			BzSurf = np.reshape( map_coordinates(BzIn, inputArr  ), Rshape )
			BphiSurf = np.reshape( map_coordinates(BphiIn, inputArr  ), Rshape )
			
			if self.nDim > 3:
				dBphidphiSurf = np.reshape( map_coordinates(dBphidphiIn, inputArr  ),Rshape )
				dBphidrSurf = np.reshape( map_coordinates(dBphidrIn, inputArr  ),Rshape )
				dBphidzSurf = np.reshape( map_coordinates(dBphidzIn, inputArr  ),Rshape)
		
				dBrdphiSurf = np.reshape( map_coordinates(dBrdphiIn, inputArr  ),Rshape  )
				dBrdrSurf = np.reshape( map_coordinates(dBrdrIn, inputArr  ), Rshape )
				dBrdzSurf = np.reshape( map_coordinates(dBrdzIn, inputArr  ),Rshape  )
		
				dBzdphiSurf = np.reshape( map_coordinates(dBzdphiIn, inputArr  ), Rshape )
				dBzdrSurf = np.reshape( map_coordinates(dBzdrIn, inputArr  ),Rshape  )
				dBzdzSurf = np.reshape( map_coordinates(dBzdzIn, inputArr  ),Rshape  )

			if self.nDim > 3:
				return BrSurf,BzSurf,BphiSurf,dBphidphiSurf,dBphidrSurf,dBphidzSurf,dBrdphiSurf,dBrdrSurf,dBrdzSurf,dBzdphiSurf,dBzdrSurf,dBzdzSurf
			else:
				return BrSurf,BzSurf,BphiSurf

			#self.fPsi=interpolate.RectBivariateSpline(equR, equz, equPSI)

			"""
/ mfbe3d.sci  -- read and write 3D mfbe files
// 07-Sep-2009 Wolfgang Suttrop


// read from ascii mfbe file
function [mfbe3d, err]=read_mfbe3d_ascii(fn)
  [unit,err]=file('open',fn, 'old');
  if err<>0 then
    printf("Error opening file %s for reading\n", fn); 
    file('close',unit); mfbe3d=[]; return;
  end
  dims = read(unit,1,3); enr=dims(2); enz=dims(3);
  comps = read(unit,1,3); ena=comps(1); enfd=comps(3);
  grid1 = read(unit,1,4);
  bfield = read(unit,enr*enz*enfd,ena);
  file('close',unit);   
  mfbe3d = struct( 'bfield', bfield, ...
      'enf', dims(1),  'enr', dims(2),  'enz', dims(3), ...
      'ena', comps(1), 'enp', comps(2), 'enfd',comps(3), ...
      'r00', grid1(1), 'z00', grid1(2), 'dr0', grid1(3), 'dz0', grid1(4));  
endfunction


function [mfbe3d, err]=read_mfbe3d_binary(fn)
  [fd,err]=mopen(fn,"r");
  if err<>0 then
    printf("Error opening file %s for reading\n", fn); 
    mfbe3d=[]; return;
  end

  bl1=mgeti(1,'i',fd);
  if bl1<>48 then 
    printf("Hmmmh, record length %d<>48. Maybe not a binary mfbe file\n", bl1);
    mclose(fd); mfbe3d=[]; return;
  end
  dims = mget(3,'d',fd);  
  comps = mget(3,'d',fd);
  bl2=mgeti(1,'i',fd);
  if bl2<>bl1 then 
    printf("Hmmmh, record length %d<>%d. Maybe not a binary mfbe file\n", bl1,bl2); 
    mclose(fd); mfbe3d=[]; return;
  end
  
  bl1=mgeti(1,'i',fd);
  if bl1<>32 then 
    printf("Hmmmh, record length %d<>32. Maybe not a binary mfbe file\n", bl1);
    mclose(fd); mfbe3d=[]; return;
  end
  grid1 = mget(4,'d',fd);
  bl2=mgeti(1,'i',fd);
  if bl2<>bl1 then 
    printf("Inconsistent record length %d<>%d.\n", bl1,bl2); 
    mclose(fd); mfbe3d=[]; return;
  end

  rc = comps(1)*dims(2)*dims(3); rc8=rc*8;
  y = zeros(rc,comps(3));
  for i=1:comps(3)
    bl1=mgeti(1,'i',fd);
    if bl1<>rc8 then 
      printf("Error, record length %d<>%d.\n", bl1,rc8);
      mclose(fd); mfbe3d=[]; return;
    end
    y1 = mget(rc,'d',fd);
    bl2=mgeti(1,'i',fd);
    if bl2<>bl1 then 
      printf("Inconsistent record length %d<>%d.\n", bl1,bl2); 
      mclose(fd); mfbe3d=[]; return;
   end
   y(:,i)=y1';
  end
  mclose(fd);

  mfbe3d = struct( 'bfield', matrix(y,comps(1),dims(2)*dims(3)*comps(3))', ...
      'enf', dims(1),  'enr', dims(2),  'enz', dims(3), ...
      'ena', comps(1), 'enp', comps(2), 'enfd',comps(3), ...
      'r00', grid1(1), 'z00', grid1(2), 'dr0', grid1(3), 'dz0', grid1(4));  
endfunction

"""



		
