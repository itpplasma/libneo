import numpy as np
import IPython

from scipy.integrate import cumtrapz,simps
from scipy.linalg import norm
from scipy.interpolate import griddata,interp1d
import os
##plot modules
import matplotlib.pylab as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import colors
from matplotlib import gridspec


import BinaryReader



class nemec:
	def __init__( self ,  folder='/afs/ipp-garching.mpg.de/home/e/ers/all/users/wls/30839/nemec/', filename='wout.3d_rmp_n2_b_even_a', shift=False,equShot=0,Rshift=0.0,zshift=0.0):

		self.folder = folder
		self.filename = filename
		self.Rshift = Rshift
		self.zshift = zshift
		self.shift = shift
		self.equShot = equShot
		self.statusRead = False
		self.statusThetaStar = False
		self.corrStatus = False

	def __del__( self ):

		self.Unload( )
		del self.statusRead 
		del self.corrStatus

	def Unload( self ):

		if self.statusRead:
			del self.rhot 
			del self.hrhot
			del self.hiota 
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
			del self.rhop 
			self.statusRead = False


			print 'hallo'
# read single time point
	def read( self ):
		
		fileOutput = self.folder + self.filename
		
		if os.path.islink(fileOutput):
			file = open(os.readlink(fileOutput), 'r')
		else:
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

		hbumnc_up = np.zeros( (mpol1+1,2*ntor+1,nsin+1) )
		hbumns_up = np.zeros_like(hbumnc_up)
		hbvmnc_up = np.zeros_like(hbumnc_up)
		hbvmns_up = np.zeros_like(hbumnc_up)

		hbsmnc_dw = np.zeros_like(hbumnc_up)
		hbsmns_dw = np.zeros_like(hbumnc_up)
		hbumnc_dw = np.zeros_like(hbumnc_up)
		hbumns_dw = np.zeros_like(hbumnc_up)
		hbvmnc_dw = np.zeros_like(hbumnc_up)
		hbvmns_dw = np.zeros_like(hbumnc_up)

		flmnc = np.zeros_like(frmnc)
		flmns = np.zeros_like(frmnc)

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

      
		#rhot = np.sqrt(( (hphi[:])/phiedge))

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
		self.phitor = newPhi   
		self.rhot = np.sqrt(newPhi/newPhi[-1])#**2.
		#IPython.embed()
		self.hrhot = np.sqrt((hphi[:]-np.diff(hphi).mean()/2.) / phiedge)#**2.
		self.hiota = hiota
		self.hq = 1./self.hiota
		#self.q = np.append(self.hq[0],self.hq)
		#self.hrhop = cumtrapz(dtheta_star,theta,axis=1,initial=0)
		self.hpres = hpres
		self.frmnc = frmnc
		self.frmns = frmns
		self.fzmnc = fzmnc
		self.fzmns = fzmns

		
		self.flmns = flmns
		self.flmnc = flmnc

		self.ntor = ntor
		self.nsin = nsin
		self.nflux = nsin+1
		self.mpol1 = mpol1
		
		self.mSize = np.size(frmnc[:,0,0])
		self.nSize =  np.size(frmnc[0,:,0])

		self.enfp = enfp

		self.nNum = np.linspace(-self.ntor,self.ntor,2*self.ntor+1,endpoint=True)
		self.mNum = np.arange(0,self.mpol1+1)


	#	self.Rshift=0.0
	#	self.zshift=0.0
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
	

		self.statusRead = True
		mu0=4.*np.pi/1.e7
		Rsym,zsym,norm,dummy = self.getUnpertSurfaces( thetaIn = theta2 )
#geometric axis
		self.phiEdge = phiedge
		self.R0 = (Rsym.max(axis=0)+Rsym.min(axis=0))/2.0 + self.Rshift
		self.z0 = (zsym.max(axis=0)+zsym.min(axis=0))/2.0 + self.zshift
		self.Rsep = np.squeeze(Rsym[:,-1]) + self.Rshift
		self.zsep = np.squeeze(zsym[:,-1]) + self.zshift
		#polodial flux is the integration of derivative of toroidal flux times iota
		self.hphipol = np.cumsum(hphip*ds*2.*np.pi*self.hiota) 
		self.phipol = np.append([0.0],self.hphipol)
		self.rhop = np.sqrt(self.phipol/self.phipol[-1])
		self.hrhop = np.interp(self.hrhot,self.rhot,self.rhop)
		self.hpres = hpres/mu0
		self.pres = np.interp(self.rhot,self.hrhot,self.hpres)
		self.hrhop = np.interp(self.hrhot,self.rhot,self.rhop)
		self.q = np.interp(self.rhot,self.hrhot,self.hq)
		self.hIt = -hbuco/mu0*2.*np.pi
		self.hIp = -hbvco/mu0*2.*np.pi*self.enfp
		self.It = np.interp(self.rhot,self.hrhot,self.hIt)
		self.Ip = np.interp(self.rhot,self.hrhot,self.hIp)

		self.qRat = np.arange(self.enfp,self.hq.max()*self.enfp)/self.enfp
		self.rhopRat=  np.interp(self.qRat,self.hq,self.hrhop)
		self.rhotRat = np.interp(self.qRat,self.hq,self.hrhot)
		self.statusRead = True
		#IPython.embed()




	def getFourierCoeff ( self, rhotIn = None, interpolRhot = False):
		if self.statusRead == False:
			return

### if no rhot is given make all			
		if np.all(rhotIn) == None:
			idxRhot = np.arange(self.nflux)
			nflux = self.nflux
			frmnc = self.frmnc
			frmns = self.frmns
			fzmnc = self.fzmnc
			fzmns = self.fzmns
			usedRhot = self.rhot
		elif np.size(rhotIn) >= 1 :
			#use no interpolation
			if  interpolRhot == False:
				
				if  np.size(rhotIn) == 1 :
					idx = np.argmin(np.abs(self.rhot-rhotIn))
					nflux = 1
				##to have one dimension at the end
					frmnc = (self.frmnc[:,:,idx])[:,:,None]
					frmns = (self.frmns[:,:,idx])[:,:,None]
					fzmnc = (self.fzmnc[:,:,idx])[:,:,None]
					fzmns = (self.fzmns[:,:,idx])[:,:,None]
					idxRhot = idx
				else:
					nflux = np.size(rhotIn)
					idxRhot = np.zeros((nflux),dtype='int')
					for i in np.arange(nflux):
						idxRhot[i] = np.argmin(np.abs(self.rhot-rhotIn[i]))
				
					frmnc = self.frmnc[:,:,idxRhot]
					frmns = self.frmns[:,:,idxRhot]
					fzmnc = self.fzmnc[:,:,idxRhot]
					fzmns = self.fzmns[:,:,idxRhot]
				
				usedRhot=self.rhot[idxRhot]

			else:
				nflux = np.size(rhotIn)
				frmnc = np.zeros((self.mSize,self.nSize,nflux))
				frmns = np.zeros_like(frmnc)
				fzmnc = np.zeros_like(frmnc)
				fzmns = np.zeros_like(frmnc)

				for m in np.arange(self.mSize):
					for n in np.arange(self.nSize):
						## if the rhot is not within the range 
						#Value to return for rhotIn < 0.0 default is rhot = 0.0.
						#Value to return for rhotIn > 0.0 default is rhot = 1.0.
						frmnc[m,n,:]  = np.interp(rhotIn, self.rhot, self.frmnc[m,n,:])
						frmns[m,n,:]  = np.interp(rhotIn, self.rhot, self.frmns[m,n,:])
						fzmnc[m,n,:]  = np.interp(rhotIn, self.rhot, self.fzmnc[m,n,:])
						fzmns[m,n,:]  = np.interp(rhotIn, self.rhot, self.fzmns[m,n,:])
				
				usedRhot = np.interp(rhotIn, self.rhot, self.rhot)
				
		return frmnc,frmns,fzmnc,fzmns,nflux,usedRhot




##calculate unperturbed surfaces!!!
	def getUnpertSurfaces( self , thetaIn = None, rhotIn=None,interpolRhot = True ):
		
		if self.statusRead == False:
			return


		### if no rhot is given make all
		frmnc,frmns,fzmnc,fzmns,nflux,usedRhot = self.getFourierCoeff ( rhotIn = rhotIn, interpolRhot = interpolRhot)

		if np.all(thetaIn)== None:
			theta = np.linspace(-np.pi,np.pi,90.)
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
	


		nIdx = np.where(self.nNum == 0)[0]
		mArr = np.arange(self.mpol1+1)  

		
		###always use a 2 dim theta
		mrho=np.array([mArr,]*nflux).T
		rhom=np.array([np.arange(nflux),]*mArr.size)
			##assuming theta3D m,theta,rho
		theta3D = np.array([theta2,]*mArr.size)
		m3D = np.reshape(np.array([self.mNum[mArr],]*theta2.size,dtype='float').T,np.shape(theta3D))

		cosmtheta = np.cos(m3D*theta3D)
		sinmtheta = np.sin(m3D*theta3D)
		
		rcRco = np.einsum('mr,mtr -> tr', frmnc[mrho,nIdx,rhom], cosmtheta)
		rsRsi = np.einsum('mr,mtr -> tr', frmns[mrho,nIdx,rhom], sinmtheta)
		
		zcZco = np.einsum('mr,mtr -> tr', fzmnc[mrho,nIdx,rhom], cosmtheta)
		zsZsi = np.einsum('mr,mtr -> tr', fzmns[mrho,nIdx,rhom], sinmtheta)
		
		Rsym = rcRco + rsRsi
		zsym = zcZco + zsZsi

		if self.shift:
			
			Rsym = Rsym  +	self.Rshift 
			zsym = zsym  +	self.zshift 

		
		nuller = np.zeros_like(Rsym.T)
		einser = np.ones_like(Rsym.T)
		
		if nflux == 1:
			gradR = (np.gradient(np.squeeze(Rsym)))[None,:]
			gradz = (np.gradient(np.squeeze(zsym)))[None,:]
		else:
			gradR = (np.gradient(Rsym)[0]).T
			gradz = (np.gradient(zsym)[0]).T


		
		InI = np.sqrt(gradR*gradR+gradz*gradz)
		InI[InI==0.0]=1.0
		vecTang=np.array([nuller,gradR/InI,gradz/InI]).T
		vecPhi=np.array([-einser,nuller,nuller]).T
		vecNorm = np.cross(vecPhi, vecTang)

		norm = np.array([vecNorm[:,:,1],vecNorm[:,:,2]])
			
		return Rsym,zsym,norm,usedRhot

#		plt.quiver(Rsym[:,-1],zsym[:,-1],vecNorm[:,-1,1],vecNorm[:,-1,2] , headwidth=4, headlength=6)
		





### get perturbed surface, be careful, could be very processing consumiung ,
	def getPertSurfaces( self , phiIn=None, thetaIn = None, rhotIn=None, nNumber =None, onlyCorrugation = False, interpolRhot = True ):
	
		if np.all(rhotIn) == None:
			rhotIn = np.linspace(0.1,1.0,180)
			
		nflux = np.size(rhotIn)
	
		dimTheta = 1
		if np.all(thetaIn)== None:
			theta = np.linspace(-np.pi,np.pi,90)
			nTheta = np.size(theta)
			dimTheta = 1
		elif  np.size(np.shape(thetaIn)) == 1  :
			theta = np.array(thetaIn)
			nTheta = np.size(theta)
			dimTheta = 1
		elif np.where(np.array(np.shape(thetaIn))==nflux)[0].size == 1 :	
			nTheta = np.int(np.size(thetaIn)/nflux)
			idxFlux = np.where(np.array(np.shape(thetaIn))==nflux)[0]	
			theta2 = (np.swapaxes(thetaIn,1,idxFlux) )
			dimTheta = 2


		if np.all(phiIn)== None:
			phi = np.linspace(0.,(2.*np.pi)/self.enfp ,360./self.enfp)
		else:
			phi = np.array(phiIn)

		nPhi = np.size(phi)
		frmnc,frmns,fzmnc,fzmns,nflux,usedRhot = self.getFourierCoeff (rhotIn=rhotIn, interpolRhot = interpolRhot)


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
		idxRhot = np.arange(nflux)


	
			##genetare idx for  self.frmnc[:,ntrho,rhotn]
		ntrho=np.array([nArr,]*idxRhot[:].size).T
		rhotn=np.array([idxRhot,]*nArr.size)

		n2D,phi2D  =  np.meshgrid(np.array(self.nNum[nArr],dtype='float'),phi)
		cosnphi = np.cos(n2D*phi2D*self.enfp)
		sinnphi = np.sin(n2D*phi2D*self.enfp)

		if dimTheta == 1:
##assuming theta theta,m
			m2D,theta2D =  np.meshgrid(np.array(self.mNum[mArr],dtype='float'),theta)
			cosmtheta = np.cos(m2D*theta2D)
			sinmtheta = np.sin(m2D*theta2D)

			#sin(a-b)=  sin(a)*cos(b) - cos(a) * sin(b)
			#cos(a-b) = cos(a)*cos(b) + sin(a) * sin(b)
			#
			print 'Calculate pertubation surfaces (can take long, but now faster)'
		
			##R = rc*cos(mtheta-nphi)+rs*sin(mtheta-nphi)
			#R = rc*Rcoco+rc*Rsisi + rs*Rsico - rs*Rcosi
			rcRcoco = np.einsum('mnr,tm,pn->ptr', frmnc[:,ntrho,rhotn], cosmtheta,cosnphi)
			rcRsisi = np.einsum('mnr,tm,pn->ptr', frmnc[:,ntrho,rhotn], sinmtheta,sinnphi)
			
			rsRsico = np.einsum('mnr,tm,pn->ptr', frmns[:,ntrho,rhotn], sinmtheta,cosnphi)
			nrsRcosi = -np.einsum('mnr,tm,pn->ptr', frmns[:,ntrho,rhotn], cosmtheta,sinnphi)
		
			zcRcoco = np.einsum('mnr,tm,pn->ptr', fzmnc[:,ntrho,rhotn], cosmtheta,cosnphi)
			zcRsisi = np.einsum('mnr,tm,pn->ptr', fzmnc[:,ntrho,rhotn], sinmtheta,sinnphi)
			
			zsRsico = np.einsum('mnr,tm,pn->ptr', fzmns[:,ntrho,rhotn], sinmtheta,cosnphi)
			nzsRcosi = -np.einsum('mnr,tm,pn->ptr', fzmns[:,ntrho,rhotn], cosmtheta,sinnphi)
		else :
				#('mr,mtr -> tr'
			theta3D = np.array([theta2,]*mArr.size)
			m3D = np.reshape(np.array([self.mNum[mArr],]*theta2.size,dtype='float').T,np.shape(theta3D))
			
			print 'Calculate pertubation surfaces with 2D Theta(can take long, but now faster)'
			
			cosmtheta = np.cos(m3D*theta3D)
			sinmtheta = np.sin(m3D*theta3D)

			rcRcoco = np.einsum('mnr,mtr,pn->ptr', frmnc[:,ntrho,rhotn], cosmtheta,cosnphi)
			rcRsisi = np.einsum('mnr,mtr,pn->ptr', frmnc[:,ntrho,rhotn], sinmtheta,sinnphi)
			
			rsRsico = np.einsum('mnr,mtr,pn->ptr', frmns[:,ntrho,rhotn], sinmtheta,cosnphi)
			nrsRcosi = -np.einsum('mnr,mtr,pn->ptr', frmns[:,ntrho,rhotn], cosmtheta,sinnphi)
		

			zcRcoco = np.einsum('mnr,mtr,pn->ptr', fzmnc[:,ntrho,rhotn], cosmtheta,cosnphi)
			zcRsisi = np.einsum('mnr,mtr,pn->ptr', fzmnc[:,ntrho,rhotn], sinmtheta,sinnphi)
			
			zsRsico = np.einsum('mnr,mtr,pn->ptr', fzmns[:,ntrho,rhotn], sinmtheta,cosnphi)
			nzsRcosi = -np.einsum('mnr,mtr,pn->ptr', fzmns[:,ntrho,rhotn], cosmtheta,sinnphi)


	
		Rper = rcRcoco+rcRsisi+rsRsico+nrsRcosi
		zper = zcRcoco+zcRsisi+zsRsico+nzsRcosi

		print 'done!'		

		if onlyCorrugation == False:
			Rper += self.Rshift
			zper += self.zshift
		
		return Rper,zper,usedRhot

	#sign areis important
	def getRhotFromQ( self, qIn=1.0 ):
		if self.statusRead:
			return np.interp(np.abs(qIn), np.abs(self.q), self.rhot)	

	def getRhopFromQ( self, qIn=1.0 ):
		if self.statusRead:
			return np.interp(np.abs(qIn), np.abs(self.q), self.rhop)	
		
	def getQFromRhot( self, rhotIn=1.0 ):
		if self.statusRead:
			return np.interp(np.abs(rhotIn), np.abs(self.rhot), np.abs(self.q))	
	
	def getRhotFromRhop( self, rhopIn=1.0 ):
		if self.statusRead:
			return np.interp(np.abs(rhopIn), np.abs(self.rhop), np.abs(self.rhot))	
	
	def getRhopFromRhot( self, rhotIn=1.0 ):
		if self.statusRead:
			return np.interp(np.abs(rhotIn), np.abs(self.rhot), np.abs(self.rhop))			
			
	def  getSurfaceCorrugation( self, rhoIn = 1.0, thetaIn=None,nNumber=None, phiIn=None, useRhop = False, plot2D=False, plot3D=False, nPhi=64,nTheta=512, exact=True,interpolRhot = True):

		if self.statusRead:

			nflux = np.size(rhoIn)
			idxRhoIn = np.zeros((nflux),dtype='int')

			if useRhop:

				rhotIn = self.getRhotFromRhop(rhoIn)

			else:

				rhotIn = rhoIn
			
			if np.all(thetaIn) == None:
				thetaIn = np.linspace(-np.pi,np.pi,nTheta,endpoint=False)
				
			if np.all(phiIn) == None:
				phi = np.linspace(0.,2.*np.pi,nPhi,endpoint=False)
			else:
				phi = phiIn
			nPhi = phi.size	

			Rsym,zsym,vecNorm,usedRhot = self.getUnpertSurfaces(  thetaIn = thetaIn , rhotIn = rhotIn,interpolRhot = interpolRhot)

			#use the exact way to calculate the corrugatiom
			if exact:
				print "Calculate the corrugation exactly (can also take long)"
				Rper,zper,usedRhot = self.getPertSurfaces(  thetaIn = thetaIn , rhotIn = rhotIn ,phiIn = phi,onlyCorrugation = False,interpolRhot = interpolRhot)
			#change sign alternativ to make it faster
				signVec=  (np.ones_like(vecNorm[0]))
				oneVec=  (np.ones_like(vecNorm[0]))
#[:,::4,:]*-1.)[:,2::4,:]*-1.
			## sign vec has -1, 1,-1,1
				signVec[::2,:] = oneVec[::2,:]*-1 
				ampRel = 0.2
				idxRho=0
				corrLen = np.zeros_like(Rper)

				for r in np.arange(nflux):
					print "Corrugation of Surface#: ",r
					usedAmp = (Rsym[:,r].max()-Rsym[:,r].mean())*ampRel
				#define curve to go around the pertubation, in R and z
					R_in =  (np.array([Rsym[:,r]+usedAmp*signVec[:,r]*vecNorm[0,:,r],Rsym[:,r]-usedAmp*signVec[:,r]*vecNorm[0,:,r]]).T).ravel()
					z_in =  (np.array([zsym[:,r]+usedAmp*signVec[:,r]*vecNorm[1,:,r],zsym[:,r]-usedAmp*signVec[:,r]*vecNorm[1,:,r]]).T).ravel()
				#calculate the intersection between the normal curve and the Pertubation
					for p in np.arange(nPhi):
						#add last point as first to close the loop
						RperIn=np.append(Rper[p,-1,r],Rper[p,:,r])
						zperIn=np.append(zper[p,-1,r],zper[p,:,r])
					
						res=find_intersect_vec(R_in, z_in,RperIn.ravel(), zperIn.ravel())
						if res.size != 0:
							Rres=np.squeeze(res).T[0]
							zres=np.squeeze(res).T[1]
							try:
								corrLen[p,:,r] = (Rres-Rsym[:,r])*vecNorm[0,:,r]+(zres-zsym[:,r])*vecNorm[1,:,r]
							except:
								print 'error in calculating surface'
								IPython.embed()
							
			# do it in a sloppy way by calculating only the projection
			else:
				print "Calculate the corrugation via projection (fast but sloppy)"
				VecRNorm = np.array([vecNorm[0],]*nPhi )
				VeczNorm = np.array([vecNorm[1],]*nPhi )
				dR,dz,usedRhot = self.getPertSurfaces(  thetaIn = thetaIn , rhotIn = rhotIn ,phiIn = phi,onlyCorrugation = True)
				corrLen = dR * VecRNorm + dz * VeczNorm


			if plot2D &  nflux == 1 :
				plt.imshow(np.squeeze(corrLen[:,::-1]).T*1.e3,extent=[phi.min(),phi.max(),thetaIn.min(),thetaIn.max()])
				CB=plt.colorbar()
				CB.set_label('Corrugation [mm]')
				plt.xlabel('phi [rad]')
				plt.ylabel('theta [rad]')	
				plt.show()

			if plot3D &  nflux == 1:
				IPython.embed()
				X = Rsym*np.cos(phi)
				Y = Rsym*np.sin(phi)
				Z = zsym
				xmin,xmax,ymin,ymax,zmin,zmax = X.min(),X.max(),Y.min(),Y.max(),Z.min(),Z.max()
				vmin,vmax = corrLen.min()*1.e3, corrLen.max()*1.e3
				fig = plt.figure()
				ax = fig.gca(projection='3d')
				ax.pbaspect = [1.0, 1.0, (zmax-zmin)/(xmax-xmin)]
				linNorm = colors.Normalize(vmin=vmin,vmax=vmax)
				surf = ax.plot_surface(X, Y, Z ,rstride=1, cstride=1,linewidth=0,facecolors=plt.cm.jet_r(linNorm(corrLen*1.e3)) )
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
	
			return Rsym,zsym,vecNorm,corrLen
			

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
				

	def getSepOnLOS(self, RLOS, zLOS, phiIn=None ):
		
		if self.statusRead == False:
			print 'not read'
			return
	
		if np.size(RLOS)!=np.size(zLOS):
			print 'R and z are different'
			return
		else:
			idxSort = np.argsort(RLOS)
			R_in = RLOS[idxSort]
			z_in = zLOS[idxSort]

		if np.all(phiIn)== None:
			phi = np.linspace(0.,(2.*np.pi)/self.enfp ,360./self.enfp)
		else:
			phi = np.array(phiIn)

		nPhi = np.size(phi)
		theta = np.linspace(-np.pi,np.pi,720,endpoint=False)
		
		Rper,zper,usedRhot = self.getPertSurfaces(  thetaIn = theta , rhotIn = [1.0] ,phiIn = phi,onlyCorrugation = False,interpolRhot = False)
		RCross = np.zeros(phi.size)
		zCross = np.zeros(phi.size)
				#calculate the intersection between the normal curve and the Pertubation
		for p in np.arange(nPhi):
						#add last point as first to close the loop
			RperIn=np.append(Rper[p,-1,0],Rper[p,:,0])
			zperIn=np.append(zper[p,-1,0],zper[p,:,0])

			res=find_intersect_vec(R_in, z_in,RperIn.ravel(), zperIn.ravel())
			if res.size != 0:
				Rres=np.squeeze(res).T[0]
				zres=np.squeeze(res).T[1]
				try:
					RCross[p] = np.squeeze(Rres)
					zCross[p] = np.squeeze(zres)
				except:
					print 'error in calculating surface'
					IPython.embed()

		return RCross,zCross


	def getLambda( self, thetaIn = None, rhotIn=None, phiIn=None, nNumber=None ):

		if self.statusRead == False:
			return

		if np.all(phiIn)== None:
			phi = np.linspace(0.,(2.*np.pi)/self.enfp ,360./self.enfp)
		else:
			phi = np.array(phiIn)

		if np.all(rhotIn)== None:
			rhotIn = self.rhot

		nflux = np.size(rhotIn)

		if np.all(thetaIn)== None:
			theta = np.linspace(-np.pi,np.pi,90.)
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
	

		flmns = np.zeros((self.mSize,self.nSize,nflux))
		flmnc = np.zeros_like(flmns)

		for m in np.arange(self.mSize):
			for n in np.arange(self.nSize):
						## if the rhot is not within the range 
						#Value to return for rhotIn < 0.0 default is rhot = 0.0.
						#Value to return for rhotIn > 0.0 default is rhot = 1.0.
				flmns[m,n,:]  = np.interp(rhotIn, self.rhot, self.flmns[m,n,:])
				flmnc[m,n,:]  = np.interp(rhotIn, self.rhot, self.flmnc[m,n,:])
	


		if nNumber == None:
			#take all
			usedN = np.abs( self.nNum ) >= 0.0
				
		else:
			usedN = (np.abs( self.nNum) == int(np.round(nNumber/self.enfp))) | (np.abs( self.nNum) == 0.0)

		nArr=np.where(usedN)[0]
		mArr=np.arange(self.mpol1+1)
		idxRhot = np.arange(nflux)

	
			##genetare idx for  self.frmnc[:,ntrho,rhotn]
		ntrho=np.array([nArr,]*idxRhot[:].size).T
		rhotn=np.array([idxRhot,]*nArr.size)

		n2D,phi2D  =  np.meshgrid(np.array(self.nNum[nArr],dtype='float'),phi)
		cosnphi = np.cos(n2D*phi2D*self.enfp)
		sinnphi = np.sin(n2D*phi2D*self.enfp)


		theta3D = np.array([theta2,]*mArr.size)
		m3D = np.reshape(np.array([self.mNum[mArr],]*theta2.size,dtype='float').T,np.shape(theta3D))
			
		print 'Calculate lambda with 2D Theta(can take long, but now faster)'
			
		cosmtheta = np.cos(m3D*theta3D)
		sinmtheta = np.sin(m3D*theta3D)

		lacoco = np.einsum('mnr,mtr,pn->ptr', flmnc[:,ntrho,rhotn], cosmtheta,cosnphi)
		lasisi = np.einsum('mnr,mtr,pn->ptr', flmnc[:,ntrho,rhotn], sinmtheta,sinnphi)
			
		lasico = np.einsum('mnr,mtr,pn->ptr', flmns[:,ntrho,rhotn], sinmtheta,cosnphi)
		nlacosi = -np.einsum('mnr,mtr,pn->ptr', flmns[:,ntrho,rhotn], cosmtheta,sinnphi)


		
		lamb = lacoco+lasisi+lasico+nlacosi

		return lamb



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


def find_intersect_vec(x_down, y_down, x_up, y_up):
    p = np.column_stack((x_down, y_down))
    q = np.column_stack((x_up, y_up))
    p0, p1, q0, q1 = p[:-1], p[1:], q[:-1], q[1:]
    rhs = q0 - p0[:, np.newaxis, :]
    mat = np.empty((len(p0), len(q0), 2, 2))
    mat[..., 0] = (p1 - p0)[:, np.newaxis]
    mat[..., 1] = q0 - q1
    mat_inv = -mat.copy()
    mat_inv[..., 0, 0] = mat[..., 1, 1]
    mat_inv[..., 1, 1] = mat[..., 0, 0]
    det = mat[..., 0, 0] * mat[..., 1, 1] - mat[..., 0, 1] * mat[..., 1, 0]
    mat_inv /= det[..., np.newaxis, np.newaxis]
    import numpy.core.umath_tests as ut
    params = ut.matrix_multiply(mat_inv, rhs[..., np.newaxis])
    intersection = np.all((params >= 0) & (params <= 1), axis=(-1, -2))
    p0_s = params[intersection, 0, :] * mat[intersection, :, 0]
    return p0_s + p0[np.where(intersection)[0]]





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
      real, allocatable :: flmnc(:,:,:),flmns(:,:,:)
      
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
      
      allocate(flmnc(0:mpol1,-ntor:ntor,0:nsin))
      allocate(flmns(0:mpol1,-ntor:ntor,0:nsin))

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
