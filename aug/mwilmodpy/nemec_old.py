import numpy as np
import IPython
import BinaryReader
from scipy.interpolate import griddata,interp1d
##plot modules
import matplotlib.pylab as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import colors
from matplotlib import gridspec


##routines from mwillens
import kk_mwillens as kk
import thetaStar




class nemec:
	def __init__( self ,  folder='/afs/ipp-garching.mpg.de/home/e/ers/all/users/wls/30839/nemec/', filename='wout.3d_rmp_n2_b_even_a', equExp='AUGD', equDiag='EQI', equShot=0.0, equTime=2.0, shift=False):

		self.folder = folder
		self.filename = filename
		self.equExp = equExp
		self.equDiag = equDiag
		self.equShot = equShot
		self.equTime = equTime
		self.shift = shift
		self.statusRead = False
		self.statusThetaStar = False
		self.statusRhop = False
		self.corrStatus = False

	def __del__( self ):

		self.Unload( )
		del self.statusRead 
		del self.corrStatus
		del self.statusRhop

	def Unload( self ):

		if self.statusRead:
			del self.rhot 
			del self.hrhot
			del self.iota 
			del self.q
			del self.hpres
			del self.frmnc
			del self.frmns
			del self.fzmnc 
			del self.fzmns 
			
			del self.ntor 
			del self.nsin 
			del self.nflux 
			del self.mpol1 

			del self.enfp 

			del self.nNum 
			del self.mNum 
			del self.Rshift
			del self.zshift
			del self.R0
			del self.z0 
			del self.Rsep 
			del self.zsep 
			self.statusRead = False

			if self.corrStatus:
				del self.corrOut
				del  self.corrRhop 
				del self.corrThetaNemec
				del self.corrThetaStar 
				self.corrStatus = False			

			if self.statusRhop:
				del self.rhop 
				del self.hrhop
				self.statusRhop = False

			if self.statusThetaStar:
				del self.theta
				del self.thetaStar
				del self.RtS
				del self.ztS
				del self.thetaOrigin 
				self.statusThetaStar = False

			print 'hallo'
# read single time point
	def read( self ):
		
		fileOutput = self.folder + self.filename
		

		file = open(fileOutput, 'r')
	
		arrayDimStr = (file.readline()).split()
		#gam adiabatic constant,enfp number of periods,number of flux surfaces 
		gam,enfp,enrho = float(arrayDimStr[0]),float(arrayDimStr[1]),float(arrayDimStr[2])
		arrayDimStr2 = (file.readline()).split()
		# number of poloidal harmonics,nb = maximum toroidal harmonic,total number of harmonics 
		empol,entor,empnt = float(arrayDimStr2[0]),float(arrayDimStr2[1]),float(arrayDimStr2[2]) 
		arrayDimStr3 = (file.readline()).split()
		#symmetry of the equilibrium, total toroidal flux
		eiasym,phiedge = float(arrayDimStr3[0]),float(arrayDimStr3[1])

		nfp    = int(enfp)
		nrho   = int(enrho)
		mpol   = int(empol)
		ntor   = int(entor)
		mpnt   = int(empnt)
		iasym  = int(eiasym)
 
		nsin   = nrho-1
		mpol1  = mpol-1
		ds     = 1./float(nsin)

#fortran allocate(frmnc(0:mpol1,-ntor:ntor,0:nsin))
		frmnc = np.zeros( (mpol1+1,2*ntor+1,nsin+1) )
		frmns = np.zeros_like(frmnc)
		fzmnc = np.zeros_like(frmnc)
		fzmns = np.zeros_like(frmnc)

		hbumnc_up = np.zeros_like(frmnc)
		hbumns_up = np.zeros_like(frmnc)
		hbvmnc_up = np.zeros_like(frmnc)
		hbvmns_up = np.zeros_like(frmnc)

		hbsmnc_dw = np.zeros_like(frmnc)
		hbsmns_dw = np.zeros_like(frmnc)
		hbumnc_dw = np.zeros_like(frmnc)
		hbumns_dw = np.zeros_like(frmnc)
		hbvmnc_dw = np.zeros_like(frmnc)
		hbvmns_dw = np.zeros_like(frmnc)

		hlmnc = np.zeros_like(frmnc)
		hlmns = np.zeros_like(frmnc)

		#allocate(hiota(1:nsin),hpres(1:nsin),hbuco(1:nsin),hbvco(1:nsin))
		hiota = np.zeros( (nsin) )
		hpres = np.zeros_like( hiota )
		hbuco = np.zeros_like( hiota )
		hbvco = np.zeros_like( hiota )

		hmass =  np.zeros_like( hiota )
		hphip = np.zeros_like( hiota )
		hphi  = np.zeros_like( hiota )
		hvp   = np.zeros_like( hiota )
		hoverr= np.zeros_like( hiota )
		fjcuru= np.zeros_like( hiota )
		fjcurv= np.zeros_like( hiota )
		hspecw= np.zeros_like( hiota )

		fsve = np.zeros((nsin+1))
		hsve = np.zeros((nsin))
		lineCounter = 0
		for j in np.arange(nsin+1):
			for m in np.arange(mpol1+1):
				nmin0=0
				if(m == 0):
					nmin0=ntor
				for n in np.arange(nmin0,2*ntor+1):
					lineCounter = lineCounter + 1
					line= (file.readline()).split()
					buffer=np.array(line)
				
					[frmnc[m,n,j],fzmns[m,n,j],frmns[m,n,j],fzmnc[m,n,j],hbumnc_up[m,n,j],hbvmnc_up[m,n,j],hbumns_up[m,n,j],hbvmns_up[m,n,j],hlmns[m,n,j],hlmnc[m,n,j],hbumnc_dw[m,n,j],hbvmnc_dw[m,n,j],hbsmns_dw[m,n,j],hbumns_dw[m,n,j],hbvmns_dw[m,n,j],hbsmnc_dw[m,n,j]]=[float(buffer[0]),float(buffer[1]),float(buffer[2]),float(buffer[3]),float(buffer[4]),float(buffer[5]),float(buffer[6]),float(buffer[7]),float(buffer[8]),float(buffer[9]),float(buffer[10]),float(buffer[11]),float(buffer[12]),float(buffer[13]),float(buffer[14]),float(buffer[15])]

		lineCounter2=0
		data = (file.read()).split()

#for j in np.arange(nsin):                
		hiota[:] = np.fromstring(np.array((data[::12])[:nsin] ,dtype='Float64')) 
		hmass[:] =   np.fromstring(np.array((data[1::12])[:nsin] ,dtype='Float64')) 
		hpres[:] = np.fromstring(np.array((data[2::12])[:nsin] ,dtype='Float64')) 
		hphip[:] =  np.fromstring(np.array((data[3::12])[:nsin] ,dtype='Float64')) 
		hbuco[:] =  np.fromstring(np.array((data[4::12])[:nsin] ,dtype='Float64')) 
		hbvco[:] =  np.fromstring(np.array((data[5::12])[:nsin] ,dtype='Float64')) 
		hphi[:]  =   np.fromstring(np.array((data[6::12])[:nsin] ,dtype='Float64')) 
		hvp[:]   =   np.fromstring(np.array((data[7::12])[:nsin] ,dtype='Float64')) 
		hoverr[:]=   np.fromstring(np.array((data[8::12])[:nsin] ,dtype='Float64')) 
		fjcuru[:]=   np.fromstring(np.array((data[9::12])[:nsin] ,dtype='Float64')) 
		fjcurv[:]=   np.fromstring(np.array((data[10::12])[:nsin] ,dtype='Float64')) 
		hspecw[:]=   np.fromstring(np.array((data[11::12])[:nsin] ,dtype='Float64')) 

		file.close()

      
		rhot = np.sqrt(( (hphi[:])/phiedge))

#allocate(fsve(0:nsin),hsve(1:nsin))
#! --- grid in s
		fsve[0]     = 0.
		fsve[-1]  = 1.
	 #do j=1,nsin-1	
		for j in np.arange(1,nsin):
#! --- half mesh
			hsve[j-1]  = (float(j)-0.5)*ds
#! --- full mesh
			fsve[j]   = float(j)*ds
      
		hsve[-1] = (float(nsin)-0.5)*ds

		###rhot corresponding to rrnzz
		
		newPhi = np.append([0.0],hphi)
		self.rhot = np.sqrt(newPhi/newPhi[-1])#**2.
	#	IPython.embed()
		self.hrhot = np.sqrt(hphi/hphi[-1])#**2.
		self.iota = hiota
		self.hq = 1./self.iota
		self.q = np.append(self.hq[0],self.hq)
		self.hpres = hpres
		self.frmnc = frmnc
		self.frmns = frmns
		self.fzmnc = fzmnc
		self.fzmns = fzmns

		self.hbumnc_up = hbumnc_up
		self.hbumns_up = hbumns_up
		self.hbvmnc_up = hbvmnc_up
		self.hbvmns_up = hbvmns_up
		
		self.hbsmnc_dw = hbsmnc_dw
		self.hbsmns_dw = hbsmns_dw

		self.ntor = ntor
		self.nsin = nsin
		self.nflux = nsin+1
		self.mpol1 = mpol1

		self.enfp = enfp

		self.nNum = np.linspace(-self.ntor,self.ntor,2*self.ntor+1,endpoint=True)
		self.mNum = np.arange(0,self.mpol1)


		self.Rshift=0.0
		self.zshift=0.0
		## add artificial shift
		if self.shift:
			if self.equShot == 30839:
				self.Rshift = 0.003
				self.zshift = 0.0045


		##calculate origin
		theta =  np.linspace(-np.pi,np.pi,720)
		theta2 = np.array( [theta,]*self.nflux).T
		nIdx = np.where(self.nNum == 0)[0]
		
		R0 = np.zeros((self.nflux,theta.size)).T
		z0 = np.zeros((self.nflux,theta.size)).T
	

		for m in np.arange(self.mpol1+1):
			R0 += self.frmnc[m,nIdx,:]*np.cos(float(m)*theta2)+self.frmns[m,nIdx,:]*np.sin(float(m)*theta2)
			z0 += self.fzmnc[m,nIdx,:]*np.cos(float(m)*theta2)+self.fzmns[m,nIdx,:]*np.sin(float(m)*theta2)
		
#geometric axis
		self.R0 = (R0.max(axis=0)+R0.min(axis=0))/2.0 + self.Rshift
		self.z0 = (z0.max(axis=0)+z0.min(axis=0))/2.0 + self.zshift
		self.Rsep = np.squeeze(R0[:,-1]) + self.Rshift
		self.zsep = np.squeeze(z0[:,-1]) + self.zshift
		
		self.statusRead = True




##calculate unperturbed surfaces!!!
	def getUnpertSurfaces( self , thetaIn = None, rhotIn=None ):
		
		if self.statusRead == False:
			return
### if no rhot is given make all			
		if np.all(rhotIn)== None:
			idxRhot = np.arange(self.nflux)
			nflux = self.nflux
		elif np.size(rhotIn) == 1 :
			idxRhot = np.argmin(np.abs(self.rhot-rhotIn))
			nflux=1
		elif np.size(rhotIn) >= 1 :
			nflux = np.size(rhotIn)
			idxRhot = np.zeros((nflux),dtype='int')
			for i in np.arange(nflux):
				idxRhot[i] = np.argmin(np.abs(self.rhot-rhotIn[i]))

		if np.all(thetaIn)== None:
			theta = np.linspace(-np.pi,np.pi,180.)
			theta2 = np.array([theta,]*nflux).T
			nTheta = np.size(theta)
		elif np.size(np.shape(thetaIn)) == 1  :
			theta = thetaIn
			theta2 =  np.array([theta,]*nflux).T
			nTheta = np.size(theta)
		elif np.where(np.array(np.shape(thetaIn))==nflux)[0].size == 1 :	
			nTheta = np.int(np.size(thetaIn)/nflux)
			idxFlux = np.where(np.array(np.shape(thetaIn))==nflux)[0]	
			theta2 = (np.swapaxes(thetaIn,1,idxFlux) )
	

		Rsym = np.zeros((nTheta,nflux))
		zsym = np.zeros_like(Rsym)

		
		nIdx = np.where(self.nNum == 0)[0]
				
		for m in np.arange(self.mpol1+1):
			Rsym += self.frmnc[m,nIdx,idxRhot]*np.cos(float(m)*theta2)+self.frmns[m,nIdx,idxRhot]*np.sin(float(m)*theta2)
			zsym += self.fzmnc[m,nIdx,idxRhot]*np.cos(float(m)*theta2)+self.fzmns[m,nIdx,idxRhot]*np.sin(float(m)*theta2)

	### add shift
		Rsym = np.squeeze(Rsym) + self.Rshift
		zsym = np.squeeze(zsym) + self.zshift

		if nflux==1:

			nuller = np.zeros_like(Rsym)
			einser = np.ones_like(Rsym)
			gradR = (np.gradient(Rsym))
			gradz = (np.gradient(zsym))
		else:
			
			nuller = np.zeros_like(Rsym.T)
			einser = np.ones_like(Rsym.T)
			gradR = (np.gradient(Rsym)[0]).T
			gradz = (np.gradient(zsym)[0]).T
		
		InI = np.sqrt(gradR*gradR+gradz*gradz)
		InI[InI==0.0]=1.0
		vecTang=np.array([nuller,gradR/InI,gradz/InI]).T
		vecPhi=np.array([-einser,nuller,nuller]).T
		vecNorm = np.cross(vecPhi, vecTang)

		if nflux==1:
			norm = np.array([vecNorm[:,1],vecNorm[:,2]])
		else:
			norm = np.array([vecNorm[:,:,1],vecNorm[:,:,2]])

		return Rsym,zsym,norm,idxRhot
#		plt.quiver(Rsym[:,-1],zsym[:,-1],vecNorm[:,-1,1],vecNorm[:,-1,2] , headwidth=4, headlength=6)
		





### get perturbed surface, be careful, could be very processing consumiung 
	def getPertSurfaces( self , phiIn=None, thetaIn = None, rhotIn=None, nNumber =None, onlyCorrugation = False ):
		
		if np.all(thetaIn)== None:
			theta = np.linspace(-np.pi,np.pi,180)
		else:
			theta = thetaIn

		if np.all(phiIn)== None:
			phi = np.linspace(0.,(2.*np.pi)/self.enfp ,360./self.enfp)
		else:
			phi = phiIn

### if no rhot is given make all			
		if np.all(rhotIn)== None:
			idxRhot = np.arange(self.nflux)
			nflux = self.nflux
		elif np.size(rhotIn) == 1 :
			idxRhot = np.argmin(np.abs(self.rhot-rhotIn))
			nflux=1
		elif np.size(rhotIn) >= 1 :
			nflux = np.size(rhotIn)
			idxRhot = np.zeros((nflux),dtype='int')
			for i in np.arange(nflux):
				idxRhot[i] = np.argmin(np.abs(self.rhot-rhotIn[i]))
	
		if nNumber == None:
			#take all
			if onlyCorrugation:
				usedN = np.abs( self.nNum ) > 0.0
			else:
				usedN = np.abs( self.nNum ) >= 0.0

		else:
			usedN = (np.abs( self.nNum) == int(np.round(nNumber/self.enfp))) | (np.abs( self.nNum) == 0.0)

		nArr=np.where(usedN)[0]
		mArr=np.arange(self.mpol1+1)

		R = np.zeros((phi.size,theta.size,nflux))
		z = np.zeros((phi.size,theta.size,nflux)) 

		theta3,phi3,rhot3 = np.meshgrid(theta,phi,self.rhot[idxRhot]) 	
		
		print 'Calculate pertubation surfaces (can take long)'
		for nIdx in nArr:
			n = self.nNum[nIdx] 
			for m in np.arange(self.mpol1+1):
				R += self.frmnc[m,nIdx,idxRhot]*np.cos(float(m)*theta3-n*phi3*self.enfp)+self.frmns[m,nIdx,idxRhot]*np.sin(float(m)*theta3-n*phi3*self.enfp)
				z += self.fzmnc[m,nIdx,idxRhot]*np.cos(float(m)*theta3-n*phi3*self.enfp)+self.fzmns[m,nIdx,idxRhot]*np.sin(float(m)*theta3-n*phi3*self.enfp)

		if onlyCorrugation == False:
			R += self.Rshift
			z += self.zshift
		
		print 'Calculate pertubation surfaces --- done'		

		return R,z,rhot3

		
		
	def defineRhop(self,equExp=None, equDiag=None, equShot=None, equTime=None):
		
		if self.statusRead:

			if (equExp == None):
				equExp = self.equExp 
			if (equDiag == None):
				equDiag = self.equDiag 
			if (equShot == None):
				equShot = self.equShot 
			if (equTime == None):
				equTime = self.equTime 
			
			out= kk.KK().kkrhotop(equShot,equTime,self.rhot, exp= equExp, diag=equDiag, ed=0)
			self.rhop = out.rhop
			self.hrhop = out.rhop[1:]
			self.statusRhop = True


	def getSynDiag( self, RIn=None, zIn=None, nNumber = None,  rhopIn=None, thetaStarIn=None, phiIn=None, nPhi=30, rhotInterval=4, nTheta=31, useRhop=False,equExp=None, equDiag=None, equShot=None, equTime=None, plotPhi=False, Rcontrol=None, zcontrol=None, debug=False, useRhoSepMax=False,perc=0.175 ):
		if self.statusRead:
			if ( (np.all(RIn) == None) | (np.all(zIn) == None) ) :
				if ( (np.all(rhopIn) != None) & (np.all(thetaStarIn) != None) ) :
					RIn, zIn =  self.getSynthRz( rhopIn=rhopIn, thetaStarIn=thetaStarIn )
				else:
					print 'no input'
					return
			elif (np.size(RIn) != np.size(zIn)):
				print 'R and z does not have the same size'
				return
			
			if (equExp == None):
				equExp = self.equExp 
			if (equDiag == None):
				equDiag = self.equDiag 
			if (equShot == None):
				equShot = self.equShot 
			if (equTime == None):
				equTime = self.equTime 
		
			if nNumber == None:
				usedN = np.abs( self.nNum ) != 0.0
			else:
				usedN = np.abs( self.nNum) == int(np.round(nNumber/self.enfp))

			nArr=np.where(usedN)[0]

			nChannels = np.size(RIn)
	
			if  nChannels==1:
				RIn = np.array([RIn])
				zIn = np.array([zIn])
			
			oriShape = np.shape(RIn)
			RIn = RIn.ravel()
			zIn = zIn.ravel()

			if (np.all(phiIn) == None):
				phiSurf = np.linspace(0.,(2.*np.pi)/self.enfp,nPhi)
			else:
				phiSurf = phiIn

			if ((np.all(Rcontrol) != None) & (np.size(Rcontrol) == np.size(phiIn)) ) :
				Rcontrol = Rcontrol
			else:
				Rcontrol = np.zeros_like(phiIn)

			if ((np.all(zcontrol) != None) & (np.size(zcontrol) == np.size(phiIn)) ) :
				zcontrol = zcontrol
			else:
				zcontrol = np.zeros_like(phiIn)
										
			# calculate theta
			theta = np.arctan2((zIn),(RIn))
			
## get rhot and rhop for diagnostic
			if debug:
				print 'use kk '
				print 'R: ',RIn
				print 'z: ',zIn
			out=kk.KK().kkrzptfn( equShot, equTime, RIn, zIn, exp= equExp, diag=equDiag, ed=0)
			rhot = out.rho_t
			rhop = out.rho_p
			
		#increase rhot range
			perc=perc
			rhotMin = rhot.min()*(1.0 - perc)
			rhotMax = rhot.max()*(1.0 + perc)

			if rhotMax >= 1.0:
				rhotMax = 1.0 

			#to avoid conflicts with X-point geometry
			if useRhoSepMax:
				rhotMax = 1.0 

			idxMin = np.argmin(np.abs( self.rhot-rhotMin ))
			idxMax = np.argmin(np.abs( self.rhot-rhotMax ))

			#first reverse to start at the outer point
			rhotSurf = ((self.rhot[idxMin:idxMax+1])[::-rhotInterval])[::-1]
			idxRhotSurf = np.arange(self.nflux)
			idxRhotSurf =  ((idxRhotSurf[idxMin:idxMax+1])[::-rhotInterval])[::-1]
			# calculate theta using lowest rhot surface for each point
			theta = np.arctan2((zIn - self.z0[idxMin]),(RIn - self.R0[idxMin]))
					##one theta array and  to get them all
			thetaSurf = np.linspace(theta.min()-np.pi/18.,theta.max()+np.pi/18.,nTheta)
			
			
			#phi,theta, fluxes
			if debug:
				print 'get Pertubation'			
			RPer,zPer,rhotPer = self.getPertSurfaces(phiSurf, thetaSurf, rhotSurf )
			if debug:
				print 'get unperturbed surface'	
			RUnp,zUnp,vecNorm,idxRhotmp = self.getUnpertSurfaces( thetaSurf , rhotSurf)
			#synthetic rhot
			synRhot = np.zeros((phiSurf.size,nChannels))

			## rhot goues only until 1.0, to get corrugation outside LCFS, we have to extend the corrugation with 1/R decay
			if useRhop:
				### get rho values
				if self.statusRhop == False:
					self.defineRhop()
				rhopPer = np.array([[self.rhop[idxRhotSurf],]*nTheta,]*phiSurf.size)
				rhop = out.rho_p
				if debug:
					IPython.embed()
				if (rhop.max()*(1.0+perc) > 1.0) & (useRhoSepMax==False) :
				#get corrugation of the separatrix
					RCorr,zCorr,rhotCorr = self.getPertSurfaces(phiSurf, thetaSurf, np.array([1.0]), onlyCorrugation = True )
				#get separatrix position
					RSep,zSep, =RUnp[:,-1],zUnp[:,-1]
				#Problem kk routine uses different theta because origin is different
					mag = kk.KK().kkrhorz( self.equShot, self.equTime,[0.0], angle = 0.0, exp= self.equExp, diag=self.equDiag, ed=0 )
				

				# angle betwenp.arccos(np.sum(v1*v2,axis=0)/Iv1I/Iv2I)
					kkTheta = np.arctan2( zSep-mag.z,RSep-mag.r )
				#define additional rhop
					diffRhop = np.diff(self.rhop[idxRhotSurf]).mean()
					addRhop = np.arange(1.0,rhop.max()*(1.0+perc),diffRhop)

					buffer=kk.KK().kkrhorz( self.equShot, self.equTime,addRhop, angle = np.rad2deg(kkTheta), exp= self.equExp, diag=self.equDiag, ed=0 )
					RSym = buffer.r
					zSym = buffer.z
				### calculate correction between CLISTE equilibrium and erika, theta depdendebt
					RCorrect = RSep-RSym[:,0]
					zCorrect = zSep-zSym[:,0]
				
				#use this index to calculate decay, minimum index where normalize vector is close to thete
					idxTT = np.argmin( np.abs( kkTheta - np.arctan2(vecNorm[1,:,-1],vecNorm[0,:,-1])))
					RNew = np.zeros( (phiSurf.size,nTheta,addRhop.size-1) )
					zNew = np.zeros_like( RNew )

					RUnpNew = np.zeros( (nTheta,addRhop.size-1) )
					zUnpNew = np.zeros_like( RUnpNew )

					decay = np.exp((RSym[idxTT,0]-RSym[idxTT,1:])*1/0.04)
					for tt in np.arange(nTheta):
						RNew[:,tt,:] = RCorr[:,tt]*decay + RSym[tt,1:] + RCorrect[tt]
						zNew[:,tt,:] = zCorr[:,tt]*decay + zSym[tt,1:] + zCorrect[tt]
						RUnpNew[tt,:] = RSym[tt,1:] + RCorrect[tt]
						zUnpNew[tt,:] = zSym[tt,1:] + zCorrect[tt]

					RPer =np.append(RPer,RNew,axis=2)
					zPer =np.append(zPer,zNew,axis=2)

					RUnp =np.append(RUnp,RUnpNew,axis=1)
					zUnp =np.append(zUnp,zUnpNew,axis=1)
				
					rhopNew= np.array([[addRhop[1:],]*nTheta,]*phiSurf.size)
					rhoPer = np.append(rhopPer,rhopNew,axis=2)
					
					rhoUnp = np.append(rhopPer[0],rhopNew[0],axis=1)
					synRho = np.zeros((phiSurf.size,nChannels))
				else:

	#synthetic rhop
					synRho = np.zeros((phiSurf.size,nChannels))
					rhoUnp = rhopPer[0]
					rhoPer = rhopPer



			else:
				#synthetic rhot
				synRho = np.zeros((phiSurf.size,nChannels))
				rhoUnp = rhotPer[0]
				rhoPer = rhotPer

			
			points = np.squeeze(np.dstack(( RUnp.ravel(),zUnp.ravel()) ))
			synRhoUnp = griddata(points, rhoUnp.ravel(), (RIn, zIn), method='cubic')

			#run throug the different phi
			for i in np.arange(phiSurf.size):
				
				points = np.squeeze(np.dstack(( RPer[i].ravel()+Rcontrol[i],zPer[i].ravel()+zcontrol[i]) ))
				synRho[i] = griddata(points, rhoPer[i].ravel(), (RIn, zIn), method='cubic')
	

		#			plt.imshow( (LSnem.amp[j])[:,::-1],vmin=0,vmax=vmax,extent=[Rmeas.min(),Rmeas.max(),zmeas.min(),zmeas.max()],aspect='auto') 
		#			cbar=plt.colorbar()
#		
					
			if plotPhi:
#
				plt.title('synthetic Data')
				plt.plot(RIn, zIn,'ro')	
				colors=['b','r','g','y']
				for i in np.arange(phiSurf.size):
					plt.contour(RPer[i],zPer[i],rhoPer[i],colors=colors[i])				
					plt.plot(0.0,0.0,color=colors[i],label="phi=%.0f"%np.rad2deg(phiSurf[i]))

				plt.plot(0.0,0.0,color='k',linewidth=2,label="unpert.")
				CS = plt.contour(RUnp,zUnp,rhoUnp,colors='k',linewidth=4) 	
				manual_locations = [(2, -0.2), (2.025, -0.2), (2.05, -0.2), (2.075, -0.2), (2.1, -0.2), (2.125, -0.2), (2.15, -0.2)]
				plt.clabel(CS, inline=2, fontsize=14,fmt='%1.2f',manual=manual_locations)
				plt.xlim([RPer.min(),RPer.max()])
				plt.ylim([zPer.min(),zPer.max()])
				plt.xlabel("R [m]",fontsize=14)
				plt.ylabel("z [m]")
				
				plt.legend(loc=4)
				plt.show()


			
			if np.size(oriShape) == 1:
				synRhoReturn = np.reshape(synRho,(phiSurf.size,oriShape[0]))
			elif np.size(oriShape) == 2:
				synRhoReturn = np.reshape(synRho,(phiSurf.size,oriShape[0],oriShape[1]))
			elif np.size(oriShape) == 3:
				synRhoReturn = np.reshape(synRho,(phiSurf.size,oriShape[0],oriShape[1],oriShape[2]))

			return np.reshape(synRhoUnp,oriShape), phiSurf,synRhoReturn 
			
			
			
	def  getSurfaceCorrugation( self, rhoIn = 1.0, thetaIn=None,nNumber=None, useRhop = False, plot2D=False, plot3D=False, nPhi=128,nTheta=256 ):

		if self.statusRead:

			nflux = np.size(rhoIn)
			idxRhoIn = np.zeros((nflux),dtype='int')

			if useRhop:
				if self.statusRhop == False:
					self.defineRhop()

				if nflux > 1:
					for i in np.arange(nflux):
						idxRhoIn[i] = np.argmin(np.abs(self.rhop-rhoIn[i]))
				else:
					idxRhoIn[0] = np.argmin(np.abs(self.rhop-rhoIn))
			else:
				if nflux > 1:
					for i in np.arange(nflux):
						idxRhoIn[i] = np.argmin(np.abs(self.rhot-rhoIn[i]))
				else:
					idxRhoIn[0] = np.argmin(np.abs(self.rhot-rhoIn))

			
			rhotIn = self.rhot[idxRhoIn]

			phi = np.linspace(0.,2.*np.pi,nPhi)
			nPhi = phi.size


###mode selection
			if nNumber == None:
				usedN = np.abs( self.nNum ) != 0.0
			else:
				usedN = np.abs( self.nNum) == int(np.round(nNumber/self.enfp))

			nArr = np.where(usedN)[0]
			mArr = np.arange(self.mpol1+1)
			
			n2D = np.array([nArr,]*mArr.size,dtype='float')
			m2D = np.array([mArr,]*nArr.size,dtype='float').T

			use5D=False

			if np.all(thetaIn) == None:
				theta = np.linspace(-np.pi,np.pi,nTheta)
				nTheta = np.size(theta)
				theta2,phi2 = np.meshgrid( theta, phi )
				theta3 = np.array([np.squeeze(theta2),]* nflux )
				theta4D= np.swapaxes(np.reshape(np.array([theta,]*mArr.size*nPhi*nArr.size),(nArr.size,mArr.size,nPhi,nTheta) ),2,3)
				Rsym,zsym,vecNorm,idxRho = self.getUnpertSurfaces( theta , rhotIn )

			elif np.size(np.shape(thetaIn))==1 :
				theta = thetaIn
				nTheta = np.size(theta)
				theta2,phi2 = np.meshgrid( theta, phi )
				theta3 = np.array([np.squeeze(theta2),]* nflux )
				theta4D= np.swapaxes(np.reshape(np.array([theta,]*mArr.size*nPhi*nArr.size),(nArr.size,mArr.size,nPhi,nTheta) ),2,3)
				Rsym,zsym,vecNorm,idxRho = self.getUnpertSurfaces( theta , rhotIn )

			elif np.where(np.array(np.shape(thetaIn))==nflux)[0].size == 1 :	
				nTheta = np.int(np.size(thetaIn)/nflux)
				theta = np.array([thetaIn,]*nPhi)
				idxFlux = np.where(np.array(np.shape(theta))==nflux)[0]	
				theta = (np.swapaxes(theta,0,idxFlux) )
				idxPhi = np.where(np.array(np.shape(theta))==nPhi)[0]
				theta3 = (np.swapaxes(theta,1,idxPhi) )
				phi2 = np.array([phi,]*nTheta).T
				theta4D= np.swapaxes(np.reshape(np.array([theta,]*mArr.size*nPhi*nArr.size),(nArr.size,mArr.size,nPhi,nTheta) ),2,3)
				use5D=True
				Rsym,zsym,vecNorm,idxRho = self.getUnpertSurfaces( thetaIn , rhotIn )

			else:
				print 'upps, something went wrong'
				return
			
			## only one surface is allowed
	##		if np.size(rhotIn) > 1:
	##			rhotIn = rhotIn[0]
			
		
			#Rsym = np.array([np.squeeze(Rsym),]* nPhi )
			#zsym = np.array([np.squeeze(zsym),]* nPhi )
			
			VecRNorm = np.array([vecNorm[0],]*nPhi )
			VeczNorm = np.array([vecNorm[1],]*nPhi )

			corru = np.zeros_like(VecRNorm)
			
			R = np.zeros_like(VecRNorm)
			z = np.zeros_like(VecRNorm)
			if nflux > 1:
				R  = np.swapaxes(R,0,1).T
				z  = np.swapaxes(z,0,1).T
			else:
				
				R = np.array([np.squeeze(R),]* 1 )
				z = np.array([np.squeeze(z),]* 1 )			



			for i in np.arange(idxRhoIn.size):
				for nIdx in nArr:
					n = self.nNum[nIdx] 
					for m in np.arange(self.mpol1+1):	
						R[i,:,:] += self.frmnc[m,nIdx,idxRhoIn[i]]*np.cos(float(m)*theta3[i]-float(n)*phi2[:]*self.enfp)+self.frmns[m,nIdx,idxRhoIn[i]]*np.sin(float(m)*theta3[i]-float(n)*phi2[:]*self.enfp)
						z[i,:,:] += self.fzmnc[m,nIdx,idxRhoIn[i]]*np.cos(float(m)*theta3[i]-float(n)*phi2[:]*self.enfp)+self.fzmns[m,nIdx,idxRhoIn[i]]*np.sin(float(m)*theta3[i]-float(n)*phi2[:]*self.enfp)
				
			dR  = np.squeeze(np.swapaxes(R.T,0,1))
			dz  = np.squeeze(np.swapaxes(z.T,0,1))
			IVecCorrI = np.sqrt(dR*dR+dz*dz)
			VecRcorr = dR/IVecCorrI
			Veczcorr = dz/IVecCorrI

			#for the case that theta is independent of rho
			if use5D==False:
			
				ntrho=np.array([nArr,]*idxRhoIn[:].size).T
				rhotn=np.array([idxRhoIn,]*nArr.size)
				m4D = np.swapaxes(np.reshape(np.array([mArr,]*nTheta*nPhi*nArr.size),(nArr.size,nPhi,nTheta,mArr.size)),1,3)
				n4D= np.reshape(np.array([nArr,]*nTheta*nPhi*mArr.size),(nPhi,nTheta,mArr.size,nArr.size)).T
				theta4D= np.swapaxes(np.reshape(np.array([thetaSi,]*mArr.size*nPhi*nArr.size),(nArr.size,mArr.size,nPhi,nTheta) ),2,3)
				phi4D= np.reshape(np.array([phi,]*mArr.size*nTheta*nArr.size),(nArr.size,mArr.size,nTheta,nPhi) )

				Rnew = np.einsum('ijk,jilm->kml', self.frmnc[:,ntrho,rhotn],np.cos(m4D*theta4D+n4D*phi4D) )
				Rnew += np.einsum('ijk,jilm->kml', self.frmns[:,ntrho,rhotn],np.sin(m4D*theta4D+n4D*phi4D) )
				
				znew = np.einsum('ijk,jilm->kml', self.fzmnc[:,ntrho,rhotn],np.cos(m4D*theta4D+n4D*phi4D) )
				znew += np.einsum('ijk,jilm->kml', self.fzmns[:,ntrho,rhotn],np.sin(m4D*theta4D+n4D*phi4D) )
			else:
				print 'hallo'
			#plt.plot(Rsym[:,-1]+10*dR[0,:,-1],zsym[:,-1]+10*dz[0,:,-1])
			#plt.plot(Rsym[:,-1],zsym[:,-1])
			#angleVecNorm=np.arctan2(VeczNorm,VecRNorm)
			#angleVecCorr=np.arctan2(Veczcorr,VecRcorr)
	#		corru = dR * VecRNorm + dz * VeczNorm
	#		rhoN = np.zeros_like(corru[0,:,-1])
	#		idxN = np.zeros_like(corru[0,:,-1],dType='int')
			
	#		thetaTest = np.linspace(0.,2*np.pi,256)
	#		phiTest = np.linspace(0.,2*np.pi,128)
	#		m3D = np.array([m2D,]*thetaTest.size).T
			##shape n,m,theta,phi
	#		m4D = np.swapaxes(np.reshape(np.array([mArr,]*thetaTest.size*phiTest.size*nArr.size),(nArr.size,phiTest.size,thetaTest.size,mArr.size)),1,3)
	#		n4D= np.reshape(np.array([nArr,]*thetaTest.size*phiTest.size*mArr.size),(phiTest.size,thetaTest.size,mArr.size,nArr.size)).T
	#		theta4D= np.swapaxes(np.reshape(np.array([thetaTest,]*mArr.size*phiTest.size*nArr.size),(nArr.size,mArr.size,phiTest.size,thetaTest.size) ),2,3)
	#		phi4D= np.reshape(np.array([phiTest,]*mArr.size*thetaTest.size*nArr.size),(nArr.size,mArr.size,thetaTest.size,phiTest.size) )
			
	#		Rnew = np.einsum('ijk,jilm->kml', self.frmnc[:,ntrho,rhotn],np.cos(m4D*theta4D+n4D*phi4D) )
	#		Rnew += np.einsum('ijk,jilm->kml', self.frmns[:,ntrho,rhotn],np.sin(m4D*theta4D+n4D*phi4D) )

	#		znew = np.einsum('ijk,jilm->kml', self.fzmnc[:,ntrho,rhotn],np.cos(m4D*theta4D+n4D*phi4D) )
	#		znew += np.einsum('ijk,jilm->kml', self.fzmns[:,ntrho,rhotn],np.sin(m4D*theta4D+n4D*phi4D) )

	#		theta3D = np.reshape(np.array([thetaTest,]*m2D.size),(nArr.size,mArr.size,thetaTest.size))
	#		nThet=10
	#		Rtest = np.zeros_like(R)
	#		ztest = np.zeros_like(R)
	#		ntrho=np.array([nArr,]*idxRhoIn[:].size).T
	#		rhotn=np.array([idxRhoIn,]*nArr.size)

	#		np.tensordot(np.cos(m3D*theta3D),self.frmnc[:,ntrho,rhotn],axes=([1,0],[0,1])
	#		np.einsum('ijk,jil->kl', np.cos(m3D*theta3D), self.frmnc[:,ntrho,rhotn])    
	#		for i in np.arange(idxRhoIn.size):
	#			np.tensordot(self.frmnc[:,ntrho,rhotn], np.cos(m3D*theta3D))

	#		for i in np.arange(corru[0,:,-1].size):
	#			idxN[i] = np.argmin(np.abs(angleVecNorm[0,i,-1]-angleVecCorr[0,:,-1]))
	#			

	#		i=48

			IPython.embed()

			if plot2D &  nflux == 1 :
				plt.imshow(corru[:,::-1].T*1.e3,extent=[phi.min(),phi.max(),theta.min(),theta.max()])
				CB=plt.colorbar()
				CB.set_label('Corrugation [mm]')
				plt.xlabel('phi [rad]')
				plt.ylabel('theta [rad]')	
				plt.show()

			if plot3D &  nflux == 1:

				X = Rsym*np.cos(phi2)
				Y = Rsym*np.sin(phi2)
				Z = zsym
				xmin,xmax,ymin,ymax,zmin,zmax = X.min(),X.max(),Y.min(),Y.max(),Z.min(),Z.max()
				vmin,vmax = corru.min()*1.e3, corru.max()*1.e3
				fig = plt.figure()
				ax = fig.gca(projection='3d')
				ax.pbaspect = [1.0, 1.0, (zmax-zmin)/(xmax-xmin)]
				linNorm = colors.Normalize(vmin=vmin,vmax=vmax)
				surf = ax.plot_surface(X, Y, Z ,rstride=1, cstride=1,linewidth=0,facecolors=plt.cm.jet_r(linNorm(corru*1.e3)) )
				sm = plt.cm.ScalarMappable(cmap='jet_r', norm=linNorm)
# fake up the array of the scalar mappable. Urgh...
				sm._A = []
				cb=plt.colorbar(sm)
				cb.set_label('corrugation [mm]')
				ax.set_xlim3d([xmin, xmax])
				ax.set_ylim3d([ymin, ymax])
				ax.set_zlim3d([zmin, zmax])

				ax.set_xlabel('x [m]')
				ax.set_ylabel('y [m]')
				ax.set_zlabel('z [m]')

#fig.colorbar(surf, shrink=0.5, aspect=5)
				plt.show()
	
			return Rsym,zsym,corru
			
				

	def calcThetaStarCorrugation(self , thetaIn = None, rhopIn = None, nThetaStar = 256 ):
		if self.statusRead:
			if self.statusThetaStar == False:
				self.defineThetaStar()

			if self.statusRhop == False:
				self.defineRhop()
				
			if np.all(rhopIn) == None:
				rhopIn = np.linspace(0.1,1.0,180)
			
			nflux = np.size(rhopIn)
			idxRhoIn = np.zeros_like(rhopIn,dtype='int')

			if nflux > 1:
				for i in np.arange(nflux):
					idxRhoIn[i] = np.argmin(np.abs(self.rhop-rhopIn[i]))
			else:
				idxRhoIn[0] = np.argmin(np.abs(self.rhop-rhopIn))

			rhopInput = self.rhop[idxRhoIn]

			#eqidistcant in thetaStaridxRhoIn[i]
			thetaStarEqui = np.zeros((nThetaStar,nflux))
			thetaOrigEqui = np.zeros((nThetaStar,nflux))

			for i in np.arange(nflux):
				thetaStarEqui[:,i] = np.linspace(self.thetaStar[:,idxRhoIn[i]].min(),self.thetaStar[:,idxRhoIn[i]].max(),nThetaStar)
				f = interp1d(self.thetaStar[:,idxRhoIn[i]], (self.thetaOrigin[:,idxRhoIn[i]])) 	
				thetaOrigEqui[:,i] = f(thetaStarEqui[:,i] )
			
####check unwrap

			Rsym,zsym,Out = self.getSurfaceCorrugation(rhoIn=rhopInput,thetaIn = thetaOrigEqui, useRhop=True)
			#Rsym,zsym,corr = self.getSurfaceCorrugation(rhoIn=rhopInput,thetaIn = thetaOrigEqui, useRhop=True)
				
			IvecI = np.sqrt( Rsym**2. + zsym**2. )
			nuller = np.zeros_like(IvecI.T)
			einser = np.ones_like(IvecI.T)
			vecTang= np.squeeze([nuller, (Rsym/IvecI).T, (zsym/IvecI).T]).T
			vecPhi = np.squeeze([ einser,nuller,nuller] ).T
			vecNorm = np.cross(vecPhi, vecTang)


		#	IPython.embed()

			R0,z0 = Rsym[:,0].mean(),zsym[:,0].mean()
			theta = np.arctan2(zsym-z0,Rsym-R0)

			self.corrRsym = Rsym
			self.corrzsym = zsym
			self.corrq = self.q[idxRhoIn]
			self.corrOut = Out
			self.corrRhop = rhopInput
			self.corrRhot = self.rhot[idxRhoIn]
			self.corrThetaNemec = thetaOrigEqui
			self.corrThetaStar = thetaStarEqui
			self.corrTheta = theta
			self.corrStatus = True

		## from thetaStar get
	def plotsm(self , thetaIn = None, rhopIn = None, nThetaStar = 256, maxN=4,maxM=18,usedN=-2,plotq=None):

		if self.statusRead:

			if self.corrStatus == False:
				self.calcThetaStarCorrugation(thetaIn = thetaIn, rhopIn =  rhopIn, nThetaStar = nThetaStar)
				nflux = np.size(self.corrRhop)
				

			corrAmpFFT = []
			for i in np.arange(nflux):
				tmp_fft = np.absolute((np.fft.fft2(self.corrOut[:,:,i],axes=(0,1))))
				tmp1 = np.concatenate( (tmp_fft[-maxN:,:(maxM+1)],tmp_fft[:(maxN+1),:(maxM+1)]),axis=0)
				tmp2 =  np.concatenate( ( tmp_fft[-maxN:,-maxM:],tmp_fft[:(maxN+1),-maxM:]),axis=0)
				tmp_out = np.concatenate((tmp2,tmp1),axis=1)
				corrAmpFFT.append(tmp_out/tmp_fft.size)

				     
			n4fft = np.arange(-maxN,maxN+1)
			m4fft = np.arange(-maxM,maxM+1)		
			corrAmp = np.array(corrAmpFFT)
			nIdx = np.where(n4fft==usedN)[0]
			
			

			if (np.all(plotq) != None):
				nPlots=2
				fig = plt.figure(figsize=(12,6))
				gs = gridspec.GridSpec(1,nPlots)
				ax1 = fig.add_subplot(gs[0,0])

			#plt.imshow(tmp_out.T,extent=[-n4fft.min(),n4fft.max(),-m4fft.min(),m4fft.max()])

			plt.imshow(np.squeeze(corrAmp[::-1,nIdx,:]),extent=[m4fft.min()-0.5,m4fft.max()+0.5,self.corrRhop.min(),self.corrRhop.max()],aspect='auto')
		#	plt.imshow(np.squeeze(corrAmp[::-1,nIdx,:]),extent=[m4fft.min(),m4fft.max(),self.corrRhop.min(),self.corrRhop.max()],aspect='auto')
#			plt.pcolor(m4fft,self.corrRhop,np.squeeze(corrAmp[::-1,nIdx,:]))
			plt.xlim([m4fft.min(),m4fft.max()])
			plt.ylim([self.corrRhop.min(),self.corrRhop.max()])
			mq=np.arange((self.hq*self.enfp).min()+1,(self.hq*self.enfp).max(),dtype='int')
			fmq = interp1d(self.hq*self.enfp,self.hrhop) 	
			mqRhop = fmq(mq)
			plt.plot(mq,mqRhop,'wo')
			plt.plot(self.hq*self.enfp,self.hrhop,'w--')
			plt.plot(-self.hq*self.enfp,self.hrhop,'w--')
		
			plt.plot(-mq,mqRhop,'wo')
			if (np.all(plotq) != None):
				plt.plot(-plotq*self.enfp,fmq(plotq*self.enfp),'ro')
				plt.plot(plotq*self.enfp,fmq(plotq*self.enfp),'ro')
			plt.colorbar()
			#plt.pcolor(rhopInput,m4fft,np.squeeze(currAmp[:,nidx,:]))
			plt.xlabel("poloidal number m")
			plt.ylabel(r' $\rm{\rho_{pol}}$',fontsize=18)
			
			if (np.all(plotq) != None):
				ax2 = fig.add_subplot(gs[0,1])
				idxQ = np.argmin(np.abs(mq-plotq*2.))
				rhopQ = mqRhop[idxQ]
				idxRhop = np.argmin(np.abs(self.corrRhop-rhopQ))
				plt.ylabel("corrugation [m]")
				plt.xlabel(r' $ \rm{\Theta^{\star}}$',fontsize=18)
				plt.plot(self.corrThetaStar[:,idxRhop],self.corrOut[0,:,idxRhop],label=r'$\rm{\rho_{pol}}\sim$ %.3f'%rhopQ+', q=%.1f'%plotq)
				plt.legend()
				#plt.plot(self.corrThetaNemec[:,idxRhop],self.corrSurf[0,:,idxRhop])
				gs.tight_layout(fig, rect=[0, 0, 1, 0.95])

			plt.show()
			
			#IPython.embed()
			

			#plt.pcolor(rhopInput,m4fft,np.squeeze(currAmp[:,nidx,:]).T)
		
			#self.getSurfaceCorrugation(rhoIn=[0.8],thetaOrigEqui)


			
	def defineThetaStar( self  ):
	
		if self.statusRead:

			if self.statusRhop == False:
					self.defineRhop()
			
			thetaIn =  np.linspace(-np.pi,np.pi,4096)
			#get unpert 
			Rsym,zsym,norm,idxRhot=self.getUnpertSurfaces(thetaIn = thetaIn)
			thetaOrigin = np.array([thetaIn,]*Rsym[0].size).T

			R0,z0 = Rsym[:,0].mean(),zsym[:,0].mean()
			theta = np.arctan2(zsym-z0,-Rsym+R0)

			idx = (np.argsort(theta,axis=0))[::-1,:]
			Rnew,znew,thetanew = np.zeros_like(theta),np.zeros_like(theta),np.zeros_like(theta)
			for i in np.arange(idx[0].size):
				thetanew[:,i] = theta[idx[:,i],i]
				Rnew[:,i] = Rsym[idx[:,i],i]
				znew[:,i] = zsym[idx[:,i],i]

			tS=thetaStar.thetaStar()
			tS.Status = True
			thetaOut,theta_star,R0tS,z0tS = tS.tomas(Rnew,znew)
			ind1 = np.int(np.mean(np.argmin(np.abs(theta_star-np.pi), axis=1)))
			
			self.theta = (np.concatenate((thetaOut[:,ind1:]-2.*np.pi,thetaOut[:,:ind1]),axis=1)).T
			self.thetaStar = (np.concatenate((theta_star[:,ind1:]-2.*np.pi,theta_star[:,:ind1]),axis=1)).T
			self.RtS = (np.concatenate(( (Rnew.T)[:,ind1:],(Rnew.T)[:,:ind1]),axis=1)).T
			self.ztS = (np.concatenate(( (znew.T)[:,ind1:],(znew.T)[:,:ind1]),axis=1)).T
			self.thetaOrigin = thetaOrigin
			#IPython.embed()
			#		self.thetaOrig = np.arctan2(((self.ztS) - self.z0),((self.RtS) - self.R0))
			self.statusThetaStar = True

	def getSepPos(self, RIn=None, zIn=None, phiIn=None):
		if self.statusRead:
			if ( (np.all(RIn) == None) | (np.all(zIn) == None) | (np.all(phiIn) == None) ) :
				print 'no input'
				return

			Rout=np.zeros_like(phiIn)
			zout=np.zeros_like(phiIn)
			LOSr=np.squeeze([RIn.min(),RIn.max()])
			LOSz=np.squeeze([zIn.min(),zIn.max()])
			#linear interpolation
			k = (LOSz[0]-LOSz[1])/(LOSr[0]-LOSr[1])
			d = -k*LOSr[1]+LOSz[1]
			##load data
			
			thetaIn = np.arctan2((zIn-self.z0[-1]),(RIn-self.R0[-1]))

			RsepNemec,zsepNemec,rhotsep = self.getPertSurfaces(phiIn=phiIn, rhotIn= np.array([1.0]),thetaIn = np.linspace(thetaIn.min()-np.pi/4.,thetaIn.max()+np.pi/4,4000) )

			#IPython.embed()

			for i in np.arange(phiIn.size):
				idxNem = np.where((np.squeeze(RsepNemec[i,:,0])>=LOSr[0])&(np.squeeze(RsepNemec[i,:,0])<=LOSr[1]))[0]
				idxSep = np.argmin(np.abs(zsepNemec[i,idxNem]-(k*RsepNemec[i,idxNem]+d)))
				Rout[i] = RsepNemec[i,idxNem[idxSep]]
				zout[i] = zsepNemec[i,idxNem[idxSep]]
			
			return Rout, zout
				
				


## for 	synthetic diagnostics calculate R z from rhop,thetastar 		
	def getSynthRz(self, rhopIn=None, thetaStarIn=None ):

		if self.statusRead:
			
			if self.statusThetaStar == False:
				self.defineThetaStar()

			if ( (np.all(rhopIn) == None) | (np.all(thetaStarIn) == None) ) :
				print 'no input'
				return

			if np.shape(rhopIn) != np.shape(thetaStarIn) :
				print 'input not the same'
				return				
			
			originShape = np.shape(rhopIn)
			rhopIn = rhopIn.ravel()
			thetaStarIn = thetaStarIn.ravel()

			Rout=np.zeros_like(rhopIn)
			zout=np.zeros_like(rhopIn)
			idxRhop = np.zeros_like(rhopIn,dtype='int')
			
			for i in np.arange(np.size(idxRhop)):
				idxRhop = np.argmin(np.abs(self.rhop-rhopIn[i]))
				idxTheta = np.argmin(np.abs(self.thetaStar[:,idxRhop]-thetaStarIn[i]))
				Rout[i] = self.RtS[idxTheta,idxRhop]
				zout[i] = self.ztS[idxTheta,idxRhop]
			
			return np.reshape(Rout,originShape), np.reshape(zout,originShape)

	#		nIdx = np.argmin(np.abs(self.nNum-nIn))
	#		r = (self.frmnc[:,nIdx,:])+ 1j * (self.frmns[:,nIdx,:])
	#		z =  (self.fzmnc[:,nIdx,:])+ 1j * (self.fzmns[:,nIdx,:])
	#		coeff = np.sqrt(np.absolute(r)**2.+np.absolute(z)**2.)
	#		plt.imshow(Bn2DAmpPlot[::-1]*1.e3,aspect='auto',extent=[nem.m4fft.min(),nem.m4fft.max(),nem.rhopSurf[0],nem.rhopSurf[-1]])
	#		plt.xlabel(r' $ \rm{poloidal \; number \; m }$',fontsize=18)
	#		plt.ylabel(r' $ \rm{\rho_{pol}}$',fontsize=18)
	#		plt.xlim([nem.m4fft.min(),nem.m4fft.max()])
	#		plt.ylim([0.0,1.0])
			#plt.plot(qSurf*(-nem.nPeriods),nem.rhopSurf,'k-',linewidth=2,label='pitch aligned')

	#		cb = plt.colorbar()
	#		cb.set_label(r'$\rm{B_n} [mT]$',fontsize=18)
	#		plt.title(r'$\rmn = %.0f}$'%-2,fontsize=20)
	#		plt.legend(loc=4)
	#		plt.show()

				#plt.pcolor()
			#

"""

	def  getmaxCorrugation( self, Rin=None, zin=None  ):


		##(m,n,i)
		nNum = np.linspace(-self.ntor,self.ntor,2*self.ntor+1,endpoint=True)
		mNum = np.arange(0,self.mpol1)
		
		nMax=2
		usedN = np.abs( nNum) <=nMax 

		theta=np.linspace(-np.pi/2.,np.pi/2.,90.)
		phi=np.linspace(0.,np.pi,)
#theta=linspace(0.,np.pi,360.)
		thetatheta,phiphi,dummy = np.meshgrid(theta,phi,np.arange(self.nsin+1))		
		R = np.zeros((nMax+1,phi.size,theta.size,self.nsin+1))
		
		z = np.zeros_like(R)
		nRun=0
		IPython.embed()
		nArr=np.where(usedN)[0]

		for n in nNum[usedN]:
			nIdx = nArr[nRun]
			for m in np.arange(self.mpol1+1):
				R[int(np.abs(n)),:,:,:] += self.frmnc[m,nIdx,:]*np.cos(float(m)*thetatheta[:]-float(n)*phiphi[:]*self.enfp)+self.frmns[m,nIdx,:]*np.sin(float(m)*thetatheta[:]-float(n)*phiphi[:]*self.enfp)
				z[int(np.abs(n)),:,:,:] += self.fzmnc[m,nIdx,:]*np.cos(float(m)*thetatheta[:]-float(n)*phiphi[:]*self.enfp)+self.fzmns[m,nIdx,:]*np.sin(float(m)*thetatheta[:]-float(n)*phiphi[:]*self.enfp)
			nRun=nRun+1

	

		## Rz, phi
		RPsiSym = R[0].mean(axis=0)
		zPsiSym = z[0].mean(axis=0)
		
		nuller = np.zeros_like(R[0].T)
		einser = np.ones_like(R[0].T)
		gradR = (np.gradient(R[0])[0]).T
		gradz = (np.gradient(z[0])[0]).T
		InI = np.sqrt(gradR*gradR+gradz*gradz)
		InI[InI==0.0]=1.0
		vecTang=np.array([nuller,gradR/InI,gradz/InI]).T
		vecPhi=np.array([-einser,nuller,nuller]).T
		vecNorm = np.cross(vecPhi, vecTang)
		
		corru_n=  R*vecNorm[:,:,1]+z*vecNorm[:,:,2]
		
		#plt.quiver(RPsi[:,-1],zPsi[:,-1],vecNorm[:,-1,1],vecNorm[:,-1,2] , headwidth=4, headlength=6)		
		#plt.show()
		IPython.embed()
"""

	#	self.frmnc[m,nIdx,j]*np.cos(mNum[m]*theta)+self.frmns[m,nIdx,j]*np.sin(mNum[m]*tt-nNum[nIdx]*pp*self.enfp)
#for m in mNum:
#	R += self.frmnc[m,nIdx,j]*np.cos(mNum[m]*tt-nNum[nIdx]*pp*self.enfp)+self.frmns[m,nIdx,j]*np.sin(mNum[m]*tt-nNum[nIdx]*pp*self.enfp)
#	z += self.fzmnc[m,nIdx,j]*np.cos(mNum[m]*tt-nNum[nIdx]*pp*self.enfp)+self.fzmns[m,nIdx,j]*np.sin(mNum[m]*tt-nNum[nIdx]*pp*self.enfp)
		
#		Rsep = RPsi[:,-1]
#		zsep = zPsi[:,-1]

	#	plt
		#R[theta,phi,j] = self.frmnc[:,nIdx,j]*np.cos(Mnum[:]*tt-Nnum[nIdx]*pp*self.enfp)+self.frmns[:,nidx,j]*np.sin(Mnum[:]*theta-Nnum[nIdx]*phi*self.enfp)
		#z[theta,phi,j] = self.fzmnc[:,nIdx,j]*np.cos(Mnum[:]*tt-Nnum[nIdx]*pp*self.enfp)+self.fzmns[:,nidx,j]*np.sin(Mnum[:]*theta-Nnum[nIdx]*phi*self.enfp)


#            enddo
#		for j in np.arange(nsin+1):
#			for m in np.arange(mpol1+1):
#				nmin0=0
#				if(m == 0):
#					nmin0=ntor
#				for n in np.arange(nmin0,2*ntor+1):	
#					R[theta,phi,j] = frmnc[midx,nidx,j]*np.cos(m*theta-n*phi*enfp)+frmns[midx,nidx,j]*np.sin(m*theta-n*phi*enfp)
#					z[theta,phi,j] = fzmnc[midx,nidx,j]*np.cos(m*theta-n*phi*enfp)+fzmns[midx,nidx,j]*np.sin(m*theta-n*phi*enfp)
#            enddo
#          enddo
#        enddo
		#numpy.fromstring
	#	IPython.embed()
		
		#ha = np.genfromtxt(fileOutput,skip_header = 3,skip_footer = nsin,invalid_raise=False)
		


"""
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  NOTE:                                                                      !
!     the program MUST be compiled in double precision:                       !
!        LINUX:  ifort -real-size 64                                          !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      module mod_nemec
!-----------------------------------------------------------------------------!
!  NEMEC parameters and fields                                                !
!-----------------------------------------------------------------------------!
      implicit none

      integer :: nfp,nrho,mpnt,nsin
      integer :: mpol,ntor,mpol1,iasym

      real    :: gam,phiedge

      real, allocatable :: hiota(:),hpres(:),hbuco(:),hbvco(:)      
      real, allocatable :: hmass(:),hphip(:),hphi(:),hvp(:)
      real, allocatable :: hoverr(:),fjcuru(:),fjcurv(:),hspecw(:)

      real, allocatable :: frmnc(:,:,:),frmns(:,:,:)
      real, allocatable :: fzmnc(:,:,:),fzmns(:,:,:)
      real, allocatable :: hbsmnc_dw(:,:,:),hbsmns_dw(:,:,:)
      real, allocatable :: hbumnc_dw(:,:,:),hbumns_dw(:,:,:)
      real, allocatable :: hbvmnc_dw(:,:,:),hbvmns_dw(:,:,:)
      real, allocatable :: hbumnc_up(:,:,:),hbumns_up(:,:,:)
      real, allocatable :: hbvmnc_up(:,:,:),hbvmns_up(:,:,:)
      real, allocatable :: hlmnc(:,:,:),hlmns(:,:,:)
      
      real,allocatable :: fsve(:),hsve(:)    ! radial mesh

      end module mod_nemec

!-----------------------------------------------------------------------------!

      subroutine read_nemec(in_equilibrium,format_type,ok)
!-----------------------------------------------------------------------------!
! purpose: reads NEMEC output file wout....                                   !
!                                                                             !
! NOTE:                                                                       !
! The values of hbumnc_up(0:mpol1,-ntor:ntor,0)                               !
!                    :                                                        !
!               hbvmns_dw(0:mpol1,-ntor:ntor,0)                               !
! are dummies. If they are needed, than the corresponding radial functions    !
!             hbumnc_up(:,:,1:nsin), ..., hbvmns_dw(:,:,1:nsin)               !
! have to be extrapolated to the magnetic axis.                               !
!-----------------------------------------------------------------------------!
      use mod_nemec

      implicit none
      
! --- input parameters
      character*250, intent(in) :: in_equilibrium
      character*25,  intent(in) :: format_type
      integer, intent(out)      :: ok
      
! --- local parameters
      integer,parameter :: inre3=8, outp6=6
      integer           :: ierr,itype,j,m,n,nmin0
      real              :: enfp,enrho,empol,entor,empnt,eiasym
      real              :: ds
      logical           :: ex 
      
      ok = 0

      ! test input
      if(format_type == 'unformatted') then
        itype=0
      elseif (format_type == 'formatted') then
        itype=1
      else
        write(outp6,1) trim(format_type)
    1   format('********** USER error **********',/,      &
               'format_type:        ',1x,a,/,             &
               '****** no valid key-word: stop *****',/)
        ok=-1
        return
      endif
      
      ! open equilibrium file
      open(inre3,file=trim(in_equilibrium),form=trim(format_type), &
           status="old",iostat=ierr)
      if(ierr.ne.0) then
        write(outp6,2) ierr,trim(in_equilibrium)
    2   format('********** USER error **********',/, &
               'ierr = ',i3,/,                       &
               'could not open file:',/,             &
               3x,A120,/,                            &
               'STOP')
        ok = -1
        return
      endif

! --- read dimensions
      if(itype.eq.0) then
        read(inre3) gam,enfp,enrho,empol,entor,empnt,eiasym,phiedge
      else
        read(inre3,*) gam,enfp,enrho,empol,entor,empnt,eiasym,phiedge
      endif
      nfp    = nint(enfp)
      nrho   = nint(enrho)
      mpol   = nint(empol)
      ntor   = nint(entor)
      mpnt   = nint(empnt)
      iasym  = nint(eiasym)
 
      nsin   = nrho-1
      mpol1  = mpol-1
      ds     = 1./float(nsin)
      
      allocate(frmnc(0:mpol1,-ntor:ntor,0:nsin))
      allocate(frmns(0:mpol1,-ntor:ntor,0:nsin))
      allocate(fzmnc(0:mpol1,-ntor:ntor,0:nsin))
      allocate(fzmns(0:mpol1,-ntor:ntor,0:nsin))
      
      allocate(hbumnc_up(0:mpol1,-ntor:ntor,0:nsin))
      allocate(hbumns_up(0:mpol1,-ntor:ntor,0:nsin))
      allocate(hbvmnc_up(0:mpol1,-ntor:ntor,0:nsin))
      allocate(hbvmns_up(0:mpol1,-ntor:ntor,0:nsin))
      
      allocate(hbsmnc_dw(0:mpol1,-ntor:ntor,0:nsin))
      allocate(hbsmns_dw(0:mpol1,-ntor:ntor,0:nsin))
      allocate(hbumnc_dw(0:mpol1,-ntor:ntor,0:nsin))
      allocate(hbumns_dw(0:mpol1,-ntor:ntor,0:nsin))
      allocate(hbvmnc_dw(0:mpol1,-ntor:ntor,0:nsin))
      allocate(hbvmns_dw(0:mpol1,-ntor:ntor,0:nsin))
      
      allocate(hlmnc(0:mpol1,-ntor:ntor,0:nsin))
      allocate(hlmns(0:mpol1,-ntor:ntor,0:nsin))

      allocate(hiota(1:nsin),hpres(1:nsin),hbuco(1:nsin),hbvco(1:nsin))
      allocate(hmass(1:nsin),hphip(1:nsin),hphi(1:nsin),hvp(1:nsin))
      allocate(hoverr(1:nsin),fjcuru(1:nsin),fjcurv(1:nsin))         
      allocate(hspecw(1:nsin))

      allocate(fsve(0:nsin),hsve(1:nsin))

! --- read VMEC output
      if(itype.eq.0) then
       
        do j=0,nsin
          do m=0,mpol1
            nmin0=-ntor
            if(m.eq.0) nmin0=0
            do n=nmin0,ntor
! --- full mesh
              read(inre3) frmnc(m,n,j),fzmns(m,n,j),         &  
                          frmns(m,n,j),fzmnc(m,n,j),         &  
! --- half mesh
                          hbumnc_up(m,n,j),hbvmnc_up(m,n,j), &  
                          hbumns_up(m,n,j),hbvmns_up(m,n,j), &  
                          hlmns(m,n,j),hlmnc(m,n,j),         &  
                          hbumnc_dw(m,n,j),hbvmnc_dw(m,n,j), &  
                          hbsmns_dw(m,n,j),                  &  
                          hbumns_dw(m,n,j),hbvmns_dw(m,n,j), &  
                          hbsmnc_dw(m,n,j)                  
            enddo
          enddo
        enddo

! --- half mesh
        read(inre3) (hiota(j),hmass(j),hpres(j),hphip(j),hbuco(j), &  
                    hbvco(j),hphi(j),hvp(j),hoverr(j),fjcuru(j),   &
                    fjcurv(j),hspecw(j),j=1,nsin)
      else
       
        do j=0,nsin
          do m=0,mpol1
            nmin0=-ntor
            if(m.eq.0) nmin0=0
            do n=nmin0,ntor
! --- full mesh
              read(inre3,*) frmnc(m,n,j),fzmns(m,n,j),       &  
                          frmns(m,n,j),fzmnc(m,n,j),         &  
! --- half mesh
                          hbumnc_up(m,n,j),hbvmnc_up(m,n,j), &  
                          hbumns_up(m,n,j),hbvmns_up(m,n,j), &  
                          hlmns(m,n,j),hlmnc(m,n,j),         &  
                          hbumnc_dw(m,n,j),hbvmnc_dw(m,n,j), &  
                          hbsmns_dw(m,n,j),                  &  
                          hbumns_dw(m,n,j),hbvmns_dw(m,n,j), &  
                          hbsmnc_dw(m,n,j)                  
            enddo
          enddo
        enddo

! --- half mesh
        read(inre3,*) (hiota(j),hmass(j),hpres(j),hphip(j),hbuco(j),   &  
                      hbvco(j),hphi(j),hvp(j),hoverr(j),fjcuru(j),     &
                      fjcurv(j),hspecw(j),j=1,nsin)
      endif
      
      close(inre3)
      
! --- grid in s
      fsve(0)     = 0.
      fsve(nsin)  = 1.
      do j=1,nsin-1
! --- half mesh
       hsve(j)  = (j-0.5)*ds
! --- full mesh
        fsve(j)   = j*ds
      enddo
      hsve(nsin) = (nsin-0.5)*ds

      !test output
      write(6,5)     
    5 format(6x,'hsve',10x'hiota',9x,'hmass',9x,'hpres',9x,'hphip', &
             9x,'hbuco',9x,'hbvco',10x,'hphi',11x,'hvp',9x,'hoverr',&
             8x,'hspecw')
      do j=1,nsin
        write(6,3) hsve(j),hiota(j),hmass(j),hpres(j),hphip(j), &
                   hbuco(j),hbvco(j),hphi(j),hvp(j),hoverr(j),hspecw(j)
      enddo
    3 format(13(2x,e12.4))
      write(6,6)
    6 format(6x,'fsve',10x,'fjcuru',8x,'fjcurv')
      do j=1,nsin
        write(6,4) fsve(j),fjcuru(j),fjcurv(j)
      enddo
    4 format(3(2x,e12.4))
      
      end subroutine read_nemec

!-----------------------------------------------------------------------------!
      
      program readnemec
!-----------------------------------------------------------------------------!
! purpose: call of read_nemec.f90                                             !
!-----------------------------------------------------------------------------!
      implicit none
      
      character*250 :: in_equilibrium    ! filename of NEMEC output 
      character*25  :: format_type       ! file format (formatted/unformatted)
      integer       :: ok
      
      in_equilibrium = 'wout.test'
      format_type    = 'formatted'
      
      call read_nemec(in_equilibrium,format_type,ok)
      
      end program readnemec
"""
