### standard Python routines
import numpy as np
from IPython import embed
from scipy.interpolate import interp1d
import matplotlib.pylab as plt
from matplotlib import gridspec

import os,sys 

import thetaStar as thetaStarClass
import kk_mwillens as kk
import MAWnew as MAW
import mfbe_slim as mfbe
import dd

labelsize=26
ticksize=18
nRows = 2
nCoils = 8

class calcBres:

	def __init__( self, equShot=None, equTime=3.0, equExp='AUGD',equDiag='EQI',equEd=0,n=2, nTorAngle=128, nThetaStar = 64,nRhop=1000 ,doPlot=False):

		self.statusInit = False

		if (np.all(equShot) != None):
			self.equShot = equShot
			self.equTime = equTime
			self.equExp = equExp
			self.equDiag = equDiag
			self.equEd = equEd
			self.nRhop = nRhop
			self.statusExt = False
		else:
			self.statusExt = True

		self.n = float(n)
		self.nTorAngle = nTorAngle
		self.nThetaStar = nThetaStar	

		self.statusInit = True
		self.statusBres = False
		self.statusEqu = False
		self.statusRot = False 
		self.statusVideo = False

		if (np.all(equShot)!=None):
			self.loadEqu( doPlot=doPlot )
 
	def __del__( self ):

		self.Unload( )
		del self.statusInit 
		del self.statusEqu

	def Unload( self ):

		if self.statusInit:
			del self.equShot
			del self.equTime
			del self.equExp
			del self.equDiag
			self.statusInit =False

		if self.statusEqu:
			del self.RIn 
			del self.zIn 
			del self.rhopIn 
			del self.qIn 
			del self.qRat 
			del self.rhopRat 
			self.statusEqu = False

		if self.statusRot:		
			del self.thetaStar 
			del self.torAngle 
			del self.BnCoils 
			self.statusRot = False

		if self.statusBres:
			del self.Bres 
			del self.diffPhi
			del self.qOut
			del self.rhopout
			self.statusBres = False

		if self.statusVideo:
			del self.n4fftSav
			del self.m4fftSav
			del self.fftBnSav
			del self.BnSav 
			self.statusBres = False		

	def loadEqu( self , doPlot=False ):
            
		if self.statusInit:
		    
			nThetaStar = self.nThetaStar
			nrhop = self.nRhop
			print 'calculating thetastar'
			tS = thetaStarClass.thetaStar(shot=self.equShot,time=self.equTime,Diagnostic = self.equDiag,Experiment = self.equExp,Edition = self.equEd)
			tS.define(nPol=512,nrhop=nrhop,origin='HFS')
			nflux = tS.rhop.size
			thetaMag = np.arctan2((tS.z-tS.z[0].mean()),(tS.R-tS.R[0].mean()))
 
                
			print 'get q profile'		
			outq = kk.KK().kkrhopfq( self.equShot, self.equTime, tS.rhop, diag=  self.equDiag, exp =  self.equExp , ed =  self.equEd )
			qSurf = outq.q
			#embed()
			#qRat = np.arange(np.floir([(np.abs(qSurf*self.n)).min()/self.n]),np.abs(qSurf).max())/self.n
			qRat = np.arange(np.ceil((np.abs(qSurf)*self.n).min())/self.n,np.ceil((np.abs(qSurf)*self.n).max())/self.n,1./self.n)
			rhopRat = np.interp(qRat,np.abs(qSurf), tS.rhop )

 
			if doPlot:
				plt.plot(tS.rhop,-qSurf)
				plt.plot(rhopRat,qRat,'ro')
				plt.xlabel(r'$\rho_{pol}$',fontsize=labelsize)
				plt.ylabel(r'$-q$',fontsize=labelsize)
				plt.xticks(fontsize=ticksize)
				plt.yticks(fontsize=ticksize)
				plt.xlim([0.0,1.0])
				plt.ylim([0.0,10.])
				plt.tight_layout()
				plt.show()
				
 ##to truncate the LCFS
			trunc=2
			thetaStarEqui = np.zeros((nflux-trunc,nThetaStar))
			angleIn = np.zeros_like(thetaStarEqui)
			RIn = np.zeros_like(thetaStarEqui)
			zIn = np.zeros_like(thetaStarEqui)

			for i in np.arange(nflux-trunc):
                    ## searching for closest
                    #idxRhoIn=np.argmin(np.abs(tS.rhop[:-2]-rhopRat[i]))
				thetaStarEqui[i] = np.linspace(tS.theta_star[i].min(),tS.theta_star[i].max(),nThetaStar,endpoint=False)
				angleIn[i] =np.interp(thetaStarEqui[i],tS.theta_star[i],np.unwrap(thetaMag[i]))	
				RIn[i] =np.interp(thetaStarEqui[i],tS.theta_star[i],tS.R[i]) 
				zIn[i] =np.interp(thetaStarEqui[i],tS.theta_star[i],tS.z[i])
 

			self.rhopIn = tS.rhop[:-trunc]
			self.qIn = qSurf[:-trunc]

			self.RIn = RIn
			self.zIn = zIn
			self.thetaMag = thetaMag
			self.thetaStarIn = thetaStarEqui
			self.qRat = qRat
			self.rhopRat = rhopRat
			self.statusEqu = True

		else:

			print 'no initialisation'
			self.statusEqu = False


	def prepareRotation(self, qRange=None, highDef=False, path_field='/afs/ipp/www/aug/CWAC_Twiki/CoilsData/',RsymExt=None,zsymExt=None, vecRext=None, veczext=None):
		if self.statusEqu:
			if np.all(qRange) == None:
				qRange = self.qRat


			idxQTmp = np.zeros_like(qRange,dtype='int')
			for i in np.arange(qRange.size):
				idxQTmp[i] = np.argmin(np.abs(np.abs(self.qIn)-np.abs(qRange[i])))

			#remove surfaces with the same index
			idxQ,idx = np.unique(idxQTmp,return_index=True )
			qRange = qRange[idx]
			rhopRat = self.rhopIn[idxQ]
						#prepare analysis
			torAngle = np.arange(self.nTorAngle)*360./float(self.nTorAngle)
			#torShift =  np.int(360./8./np.diff(torAngle).mean())
	
			#get q surfaces
			surf_r_sym = self.RIn[idxQ]
			surf_z_sym = self.zIn[idxQ]
			surf_thetaStar = self.thetaStarIn[idxQ]
						
			shapeSym = np.shape(surf_r_sym)
			
			print 'get B'
			field = kk.KK().kkrzBrzt( self.equShot, self.equTime, surf_r_sym.ravel(), surf_z_sym.ravel() , diag=  self.equDiag, exp =  self.equExp, ed =  self.equEd)
			B_r =  np.reshape(field.br,shapeSym)
			B_z = np.reshape(field.bz,shapeSym)
			B_t = np.reshape(field.bt,shapeSym)

					#calculate normal vectors
			IBrzI = np.sqrt( B_r**2. + B_z**2. )
			nuller = np.zeros_like(IBrzI.T)
			einser = np.ones_like(IBrzI.T)
			vecTang= np.array([(B_r/IBrzI).T, (B_z/IBrzI).T,nuller]).T
			vecPhi = np.array([ nuller,nuller,einser] ).T
			vecNorm = np.cross(vecPhi, vecTang)
# plt.quiver(surf_r_sym[8],surf_z_sym[8],vecNorm[0,8,:,0],vecNorm[0,8,:,1])

			vec_R_sym = vecNorm[:,:,0]
			vec_z_sym = vecNorm[:,:,1]

			BnCoils = self.getBn(surf_r_sym, surf_z_sym, torAngle, np.array([vec_R_sym,vec_z_sym]), highDef=highDef)
			
			self.RSymRat =  surf_r_sym
			self.zSymRat = surf_z_sym
			self.vecNorm = vecNorm
			self.thetaStar = surf_thetaStar
			self.torAngle = torAngle
			self.BnCoils = BnCoils
                        self.Bt = B_t
                        self.IBI = np.sqrt( B_r**2. + B_z**2. + B_t**2.)
			self.rhopRat = rhopRat
			self.qRat = qRange
			self.statusRot = True

	def prepareExternal(self, surf_r_sym, surf_z_sym,vec_R_sym,vec_z_sym, thetaStar, rhopRat, qRat, highDef=False, path_field='/afs/ipp/www/aug/CWAC_Twiki/CoilsData/'):
		
		torAngle = np.arange(self.nTorAngle)*360./float(self.nTorAngle)
		
		BnCoils = self.getBn(surf_r_sym, surf_z_sym, torAngle, np.array([vec_R_sym,vec_z_sym]), highDef=highDef,path_field=path_field)

		self.thetaStar = thetaStar
		self.torAngle = torAngle
		self.BnCoils = BnCoils
		self.rhopRat = rhopRat
		self.qRat = qRat
		self.qIn = -qRat
		self.statusRot = True


	def getBn(self, surf_r_sym, surf_z_sym, torAngle, vecIn, highDef=False, path_field='/afs/ipp/www/aug/CWAC_Twiki/CoilsData/'):

		#shift of the index
		torShift =  np.int(360./8./np.diff(torAngle).mean())

		if highDef:
			fileCoilU='Bu1_1kA_large.mfbe'
			fileCoilL='Bl1_1kA_large.mfbe'
		else:
			fileCoilU='Bu1_1kA.mfbe'
			fileCoilL='Bl1_1kA.mfbe'
			
## read upper coild
		coilU = mfbe.mfbe(path_field,fileCoilU)
		coilU.read3D()
## read lower coil
		coilL = mfbe.mfbe(path_field,fileCoilL)
		coilL.read3D()


		BnCoilSingleU = coilU.BnSingle(Rin=surf_r_sym, zin=surf_z_sym, phiIn = np.deg2rad(torAngle), vecIn=vecIn )
		BnCoilSingleL = coilL.BnSingle(Rin=surf_r_sym, zin=surf_z_sym, phiIn = np.deg2rad(torAngle), vecIn=vecIn )

		shapeBnCoils =np.append(np.append(BnCoilSingleL.shape,nRows),nCoils)
		BnCoils = np.zeros((shapeBnCoils))
##every coils is shifted
		for i in np.arange(nCoils):
			BnCoils[:,:,:,0,i] = np.roll(BnCoilSingleU,torShift*i,axis=0)
			BnCoils[:,:,:,1,i] = np.roll(BnCoilSingleL,torShift*i,axis=0)
		
		return BnCoils



	def rotate(self, phaseShot = None, phaseTime=None, qRange=None,ndiffPhi=50,differential=True,forVideo=False):
		if self.statusEqu | self.statusRot:

			if ((self.statusRot == False) & (self.statusExt == False)):
				self.prepareRotation()


			nTorAngle=self.nTorAngle
			noInput = (np.all(phaseShot) == None) |  (np.all(phaseTime) == None)
			if (self.n ==2) & noInput & differential:
				phaseShot = 30680
				phaseTime=np.array([2.235,2.735]) #np.array([2.233,2.735])
				ndiffPhi=21
				MAWsampling=1.
				qRange = self.qRat
				## for the case that q is not negative
				if self.qIn.mean() < 0.0:
					signQ = -1.
				else:
					signQ = 1.
				
			elif (self.n == 1) & noInput & differential :
				phaseShot = 32941
				phaseTime=np.array([4.883,5.883])
				MAWsampling=4.
				#embed()
				#phaseShot = 32116
				#phaseTime=np.array([2.566,3.064])
				qRange = self.qRat
				

			if self.qIn.mean() < 0.0:
				signQ = -1.
			else:
				signQ = 1.
			###read MAW
			MW=MAW.MAW()
			MW.Load(phaseShot,tBegin=phaseTime[0],tEnd=phaseTime[1], Binning = MAWsampling)
		#	embed()
			PhaseMPRun = np.linspace(np.rad2deg(MW.dPhase).min(),np.rad2deg(MW.dPhase).max(),ndiffPhi)
	
			if MW.n.mean() != self.n:
				print 'calculated n from MAW %.1f differs from given '%MW.n.mean()+'%.1f'%self.n
				return
			else:
				n = self.n

		
			idxSort = np.argsort(np.rad2deg(MW.dPhase))
			## get time points
			timeMPRun = np.interp(PhaseMPRun,np.rad2deg(MW.dPhase)[idxSort], MW.MPtime[idxSort] ) 
			

			### first dimension is delta Phi, second the q surfaces
			Bres = np.zeros((timeMPRun.size,qRange.size))
			BresNeg = np.zeros((timeMPRun.size,qRange.size))
			if forVideo:
				n4fftSav,m4fftSav,fftBnSav,BnSav = [],[],[],[]	
			#BnOut = []
			print 'get MP field, take long'
			for k in np.arange(timeMPRun.size):
				### negativ because of convention
				Icoef= -MW.getIcoils(timeMPRun[k])/1000.
				BnResult = np.einsum('ijklm,lm->ijk',self.BnCoils,Icoef)
				Bn = np.swapaxes( BnResult,1,2)


				
			### fftBn surface,n,m
				n4fft,m4fft,fftBn=doFFT(Bn)
				if forVideo:
					n4fftSav.append(n4fft),m4fftSav.append(m4fft),fftBnSav.append(np.absolute(fftBn)),BnSav.append(Bn)

				nIdx = np.where(n4fft==n)[0]
				mIdx = np.zeros_like(qRange,dtype='int')
				mIdxNeg = np.zeros_like(qRange,dtype='int')
				runIdx = np.arange(qRange.size)
			
				for j in runIdx:
					mIdx[j] = np.where(m4fft == int(signQ*n*qRange[j]))[0] 
					Bres[k,j] = np.absolute(fftBn[j,nIdx,mIdx[j]])
					mIdxNeg[j] = np.where(m4fft == -signQ*n*qRange[j])[0] 
					BresNeg[k,j] = np.absolute(fftBn[j,nIdx,mIdxNeg[j]])
				
#
#plt.plot(PhaseMPRun,Bres[:])
			self.Bres =Bres
			self.diffPhi = PhaseMPRun
			self.diffTime = timeMPRun
			self.qOut = self.qRat
			self.rhopout  = self.rhopRat
			if forVideo:
				self.n4fftSav, self.m4fftSav, self.fftBnSav, self.BnSav = np.squeeze(n4fftSav), np.squeeze(m4fftSav), np.squeeze(fftBnSav), np.squeeze(BnSav)
				self.statusVideo = True

			if differential:
				self.kind  = 'diff'
			else:
				self.kind  = 'rigid'

			self.statusBres = True
	


	def save(self,show=True,rule_of_thumb=True, dt=0.06):		
		if self.statusBres:
			
			linestyle=['dotted','dashed','dashdot']
			if self.statusExt==True:
				folder='./external/'
			else:
				folder = '%d_'%self.equShot+'%.2f_'%self.equTime+self.equExp+'_'+self.equDiag+'_%d'%self.equEd+'_'+'n=%d_'%int(self.n)+self.kind
			
			os.system("mkdir "+folder) 
			np.savetxt(folder+'/'+'Bres.dat',self.Bres)
			np.savetxt(folder+'/'+'phi.dat',self.diffPhi)
			np.savetxt(folder+'/'+'surfaces.dat',np.array([self.rhopout,self.qOut]).T)


			if ((rule_of_thumb == True) & (self.statusExt==False)):
				try:
					TOT = dd.shotfile('TOT',self.equShot)
					betaN = TOT('beta_N',tBegin=self.equTime-dt, tEnd=self.equTime+dt).mean()
					TOT.close()
					GQ = dd.shotfile('G'+self.equDiag[1:],self.equShot)
					q95 = GQ('q95',tBegin=self.equTime-dt, tEnd=self.equTime+dt).mean()
					GQ.close()
					
					deltaphi_opt_tot,deltaphi_opt_vac,shift = deltaphi_opt_rule_of_thumb(self.n, abs(q95), betaN)
					
					deltaphi_opt_tot[deltaphi_opt_tot>180.] -= 360.
					deltaphi_opt_tot[deltaphi_opt_tot<-180.] += 360.
					deltaphi_opt_vac[deltaphi_opt_vac>180.] -= 360.	
					deltaphi_opt_vac[deltaphi_opt_vac<-180.] += 360.


				except:
					print 'no TOT or GQI diagnostics'
					print 'use the good old 60 deg'
					shift = 60.

				print 'shift: ', shift
				thumbPhi = self.diffPhi+shift
				thumbPhi[thumbPhi>180.] -= 360.
				thumbPhi[thumbPhi<-180.] += 360.
				idxSort = np.argsort(thumbPhi)
				thumbPhi = thumbPhi[idxSort]
				thumbBres = (self.Bres[:,-2]*1.e3)[idxSort]

			fig = plt.figure(figsize=(16,11))
			if (self.statusExt==False):
				plt.title("shot: %d"%self.equShot+', time: %.2f'%self.equTime,fontsize=24)
				#plt.title("q~6 idx: %d,"%k+" time: %.3f s,"%timeMPRun[k]+" phase: %d"%PhaseMPRun[k],fontsize=28)
			for i in np.arange(self.qOut.size):
				plt.plot(self.diffPhi,self.Bres[:,i]*1.e3,label='q=%.2f'%self.qOut[i],linewidth=3,linestyle=linestyle[i/7])
				#plt.colorbar()
#			plt.axvline(deltaphi_opt_tot)
#			plt.axvline(deltaphi_opt_vac)
			if ((rule_of_thumb == True) & (self.statusExt==False)): 
				plt.plot(thumbPhi,thumbBres,'k-',linewidth=3,label='thumb')
			plt.xlabel(r'$\Delta\phi\rm{_{UL}}\ \rm{[deg]}$',fontsize=32)
			plt.xlim([-180.,180.])
			plt.xticks(fontsize=24)
			plt.ylabel(r'$B_{res}\ \rm{[mT]}$',fontsize=32)
			plt.yticks(fontsize=24)
			plt.legend(fontsize=24,loc = 'center left',bbox_to_anchor=(1.0, 0.5))
			plt.tight_layout(rect=[0.,0.,0.84,1.])
			fig.savefig(folder+'/Bres.png')
			if show:
				plt.show()
				
	

			print 'file saved to: '+folder
		else:
			print 'no Bres calculation'



        def __call__( self , qIn=None,  BrTime = None, BrShot = None, BrEd = 0, nTorAngle=128, doPlot=False , doOutput=True, getSpectra=False, getSingleSpectra=False,getBresProfile = False,usePSL=False,relative=False,relativeB=False,MAWDiag='MAW'):
            
		if self.statusRot == False:
			self.prepareRotation()
		
			
		if (np.all(BrShot) == None):
			BrShot=self.equShot
		if (np.all(BrTime) == None):
			BrTime=self.equTime

		if (np.all(qIn) == None):
			idxQ = -1
			idxQRat = -1
		else:
			idxQ = np.argmin(np.abs(np.abs(self.qIn)-np.abs(qIn)))
			idxQRat = np.argmin(np.abs(np.abs(self.qRat)-np.abs(qIn)))

		qIn = self.qRat[idxQRat]
				

		MW=MAW.MAW()
		MW.Load(BrShot,Binning=0.5, Edition=BrEd, Diagnostic=MAWDiag )
		
		#embed()
		try:
			Icoef= -MW.getIcoils(BrTime, usePSL= usePSL)/1000.
			BnResult = np.einsum('ijklm,lm->ijk',self.BnCoils,Icoef)
                        
		except:
			print 'hello'
			embed()

                
                if relative:
                        Bn = np.squeeze(np.swapaxes( BnResult,1,2))/self.Bt.T
                elif relativeB:
                        Bn = np.squeeze(np.swapaxes( BnResult,1,2))/self.IBI.T
                else:
                        Bn = np.squeeze(np.swapaxes( BnResult,1,2))
                
                
		#embed()
		labesize = 26
		ticksize = 22

				##save
		if ((doPlot) | (doOutput)) :
			fig = plt.figure(figsize=(13,10))
			plt.title("q~%.1f, "%qIn+" time: %.3f s,"%BrTime,fontsize=28)
			plt.contourf(self.torAngle,np.rad2deg(self.thetaStar[idxQRat]),Bn[:,:,idxQRat].T*2.e3,30)
			CB=plt.colorbar()
			CB.set_label(r'$B\rm{_{n}} [mT]$',fontsize=labelsize)
			plt.xlabel(r' $ \rm{toroidal} \ \phi $',fontsize=labelsize) 
			plt.ylabel(r' $ \rm{poloidal} \ \Theta^{\star} $',fontsize=labelsize) 
			plt.yticks(fontsize=ticksize)
			plt.xticks(fontsize=ticksize)
			if doOutput:
				plt.savefig("Surf_"+"%d"%BrShot+"q~%.1f_"%qIn+"t%.3fs.pdf"%BrTime)	
			if doPlot:
				plt.show()
			
			plt.clf()
                        plt.close()

		mMax=30
		nMax=4
		n4fft,m4fft,fftBn=doFFT(Bn, maxM=mMax, maxN=nMax)
		nIdx = np.where(n4fft==-self.n)[0]
                
		mIdx = np.where(m4fft== np.round(self.n*np.abs(qIn)))[0] 
	
	
                ### plot poloidal mode spectrum
		if ((doPlot) | (doOutput)) :
			try:
				#mIdx = np.where(m4fft == np.round(self.n*self.qIn[idxQ]))[0]
				plt.title("Surface q=%.1f"%self.qRat[idxQRat],fontsize=labelsize)
				plt.xlim([-mMax,mMax])
				plt.xlabel(r' $ \rm{poloidal \ number} \ m $',fontsize=labelsize) 
				plt.ylabel(r'$B\rm{_{n}} [mT]$',fontsize=labelsize)
				plt.bar(m4fft-0.4,np.squeeze(np.absolute(fftBn[idxQRat,nIdx,:]).T*2.e3),color='blue')
				plt.bar(m4fft[mIdx]-0.4,np.squeeze(np.absolute(fftBn[idxQRat,nIdx,mIdx]))*2.e3,color='r')
				plt.xticks( fontsize=ticksize)
				plt.yticks( fontsize=ticksize)
				plt.tight_layout()
				if doOutput:
					plt.savefig("Bar_"+"%d"%BrShot+"q~%.1f_"%qIn+"t%.3fs.pdf"%BrTime)	
				if doPlot:
					plt.show()
				plt.clf()
				plt.close()
			except:
				print 'error in bar Plot'
				embed()


		if ((doPlot) | (doOutput)) :	
			try:
				plt.title(r"$\rm{vacuum\ field,}\ n=-%.1f$"%self.n,fontsize=labelsize)
				plt.contourf(m4fft,self.rhopRat,np.squeeze(np.absolute(fftBn[:,nIdx,:]))*2.e3,100)

				plt.xlim([-mMax,mMax])
				plt.ylim([self.rhopRat.min(),self.rhopRat.max()])
				plt.xlabel(r' $ \rm{poloidal \ number }\ m $',fontsize=labelsize)
				plt.ylabel(r' $ \rm{\rho_{pol}}$',fontsize=labelsize)
				plt.plot(-self.qIn*float(self.n),self.rhopIn,'w--',label='pitch aligned',linewidth=2)
				plt.plot(self.qRat*float(self.n),self.rhopRat,'wo')
				plt.xticks( fontsize=ticksize)
				plt.yticks( fontsize=ticksize)
				leg=plt.legend(loc=3,fontsize=ticksize,frameon=False)
				for text in leg.get_texts():
					plt.setp(text, color = 'w')

				cb = plt.colorbar()
				cb.set_label(r'$B\rm{_{n}}\  \rm{[mT]}$',fontsize=labelsize)
				for t in cb.ax.get_yticklabels():
					t.set_fontsize(ticksize)
					
				plt.tight_layout()
				
				if doOutput:
					plt.savefig("Ice_"+"%d"%BrShot+"q~%.1f_"%qIn+"t%.3fs.pdf"%BrTime)	
				if doPlot:
					plt.show()

				plt.clf()
                                plt.close()

			except:
				print 'error in cornetto plot'
				embed()
				
		if (getSpectra):		
			return m4fft,self.rhopRat,np.squeeze(np.absolute(fftBn[:,nIdx,:]))*2.

		if (getSingleSpectra):		
			return m4fft,self.rhopRat[idxQRat],np.squeeze(np.absolute(fftBn[:,nIdx,:]))[idxQRat]*2.

                
                
		if (getBresProfile):
			
			fftBnOut = np.zeros_like(self.rhopRat)
			for i in np.arange(fftBnOut.size):
				idxUse = np.argmin(np.abs(m4fft-np.abs(self.qRat[i]*self.n)))
				fftBnOut[i] = np.squeeze(np.absolute(fftBn[i,nIdx,idxUse]))*2.
			return self.qRat,self.qRat*self.n,self.rhopRat,fftBnOut
		
		print 'out of'



	def makeVideo(self,dPhiScan=True,qIn=4.0):
		if dPhiScan & self.statusVideo:

			folder = 'Video_'+'%d_'%self.equShot+'%.2f_'%self.equTime+self.equExp+'_'+self.equDiag+'_'+'n=%d'%int(self.n)
			os.system("mkdir "+folder)
 
			idxQRat = np.argmin(np.abs(np.abs(self.qRat)-np.abs(qIn)))
			labelsize = 20
			ticksize = 16
#\rm{_{UL}}
			
			nIdx = np.where(self.n4fftSav[0]==-self.n)[0]
			Bnmax=np.absolute(self.BnSav[:,:,:,idxQRat]*2.e3).max() 
			vBn = np.arange(-np.round(Bnmax),np.round(Bnmax)+1.,1.0) 
			absBnMax = np.max( np.absolute(self.fftBnSav[:,:,nIdx,:])*2.e3 )
			vAbsBn = np.arange(0.0,np.round(absBnMax*100.)/100.+0.01,0.01 )
			
			r=0
			#plt.ion()
			arr = np.argsort(self.diffTime)
			for i in arr:
				
				fig = plt.figure(figsize=(16,6))
				gs = gridspec.GridSpec(1,2)
				plt.suptitle('$\Delta \phi_{UL} = %.0f^{\circ}$'%self.diffPhi[i],fontsize=labelsize+4)
				ax1 = fig.add_subplot(gs[0,0])
				plt.title(r"$q=-%.1f$"%qIn,fontsize=labelsize)
				plt.contourf(self.torAngle,np.rad2deg(self.thetaStar[idxQRat]),self.BnSav[i,:,:,idxQRat].T*2.e3,vBn)
				CB=plt.colorbar()
				CB.set_label(r'$B\rm{_{r}} [mT]$',fontsize=labelsize)
				for t in CB.ax.get_yticklabels():
					t.set_fontsize(ticksize)

				plt.xlabel(r' $ \rm{toroidal} \ \phi $',fontsize=labelsize) 
				plt.ylabel(r' $ \rm{poloidal} \ \Theta^{\star} $',fontsize=labelsize)
				plt.plot(self.torAngle,self.torAngle/qIn-35.,'k-',linewidth=3)
				plt.yticks(fontsize=ticksize)
				plt.ylim([np.rad2deg(self.thetaStar[idxQRat]).min(),np.rad2deg(self.thetaStar[idxQRat]).max()])
				plt.xticks(fontsize=ticksize)

				ax2 = fig.add_subplot(gs[0,1])
				plt.title(r'$n = -%.0f$'%self.n,fontsize=labelsize)	
				plt.contourf(self.m4fftSav[i],self.rhopRat,np.squeeze(np.absolute(self.fftBnSav[i][:,nIdx,:]))*2.e3,vAbsBn)
#plt.xlim([m4fft.min(),m4fft.max()])
				plt.xlim([-20.,20.])
				plt.ylim([self.rhopRat.min(),self.rhopRat.max()])
				plt.xlabel(r' $ \rm{poloidal \ number }\ m $',fontsize=labelsize)
				plt.ylabel(r' $ \rm{\rho_{pol}}$',fontsize=labelsize)
			#plt.plot(-self.qIn*float(self.n),self.rhopIn,'w--',label="pitch\naligned",linewidth=2)
				plt.plot(self.qRat*float(self.n),self.rhopRat,'w--',label="pitch\naligned",linewidth=3)
				plt.plot(self.qRat*float(self.n),self.rhopRat,'wo',markersize=8)
				plt.xticks( fontsize=ticksize)
				plt.yticks( fontsize=ticksize)
				leg=plt.legend(loc=3,fontsize=ticksize,frameon=False)
				for text in leg.get_texts():
					plt.setp(text, color = 'w')
				
				cb = plt.colorbar()
				cb.set_label(r'$|B\rm{_{r}}|\  \rm{[mT]}$',fontsize=labelsize)
				cb.set_ticks(np.arange(0.0,np.floor(absBnMax*10.)/10.+0.1,0.1 ))
				for t in cb.ax.get_yticklabels():
					t.set_fontsize(ticksize)
					
				plt.yticks(fontsize=ticksize)
				plt.xticks(fontsize=ticksize)
				gs.tight_layout(fig, rect=[-0.01, -0.01, 1.01, 0.96])
			
				fig.savefig(folder+'/dphi_%03d.png'%r)
				#fig.savefig(folder+'/dphi_%03d.fig'%r)
				r=r+1
                                plt.close()
			print 'dPhiSan'
			embed()


def doFFT(dataIn, maxN=8,maxM=32):
    FFTout = []
    for i in np.arange(dataIn[0,0,:].size):
        tmp_fft = np.fft.fft2(dataIn[:,:,i],axes=(0,1))
        tmp1 = np.concatenate( (tmp_fft[-maxN:,:(maxM+1)],tmp_fft[:(maxN+1),:(maxM+1)]),axis=0)
        tmp2 =  np.concatenate( ( tmp_fft[-maxN:,-maxM:],tmp_fft[:(maxN+1),-maxM:]),axis=0)
        tmp_out = np.concatenate((tmp2,tmp1),axis=1)
        FFTout.append(np.squeeze(tmp_out/tmp_fft.size))

    n4fft = np.arange(-maxN,maxN+1)
    m4fft = np.arange(-maxM,maxM+1)

    return np.squeeze(n4fft),np.squeeze(m4fft),np.squeeze(FFTout)

"""
##################################################################################
#####################  D Ryan, 13/04/2016  ########################################
#  Script to print an approximate value for deltaphi_opt for given 
#  betaN, q95, n. Based on a set of equilibria spanning q95 and betaN,
#  validated with deltaphi_opt MARS-F computations for several AUG equilibria 
#  taken from distinct experiments.
#  Currently only n=1,2 available. The scans have also been done for
#  n=3,4, but I don't yet have sufficient n=3,4 experimental validation
#  points to assess the scan results. Work in progress...
#  
#  Command line usage
#  python deltaphi_opt_rule_of_thumb.py <n> <q95> <betaN>
#
#  Example
#  python deltaphi_opt_rule_of_thumb.py 2 3.75 2.15
#  Inputs: q95 = 3.750, betaN = 2.150, n = 2
#  Outputs:
#    Vacuum optimum differential phase = 39.4
#    Total optimum differential phase = 103.0
#    Shift = 63.6
#
#  Deltaphi_opt computed by the scan, is closely approximated with a 
#  simple 2D quadratic.
#  Say z = deltaphi_opt, x = betan, y = q95
#  z = a(x^2y^2) + b(x^2y) + c(x^2) + d(xy^2) + e(xy) + f(x) + g(y^2) + h(y) + i
#  n=1 
#  Coefficient  vacuum                   total           
#  a      0.13897915268759675,     0.43304538551243965      
#  b      0.15841742712121487,     -5.7000097019466098     
#  c      -1.6740577104791932,     17.096809988618769        
#  d      -0.51685976163971326,    -2.7405449274386395     
#  e      -7.6899374778193206,     29.940411228776266       
#  f      18.736665061179082,      -99.266691574563851          
#  g      -1.2554970387577757,     -0.45865516031870279      
#  h      65.149596706769387,      49.966049466296141             
#  i      -312.1864777106091,      -210.17564659628684      
#  rmse   10.9379690445,            31.4492579124   
#
#  n=2
#  Coefficient  vacuum                  total           
#  a      0.14571159718066884,    0.14046757326169945    
#  b      1.7141560193471825,     1.7732027682693108    
#  c      -6.3853893303319893,    -8.5336264975127207  
#  d      -0.24970237007704377,   -0.33718746109109077
#  e      -23.718686622573557,    -22.025319804796972
#  f      56.20551093242895,      63.892089079560563
#  g      -3.1507918928682539,    -3.1757270507359081
#  h      127.83133097817949,     129.06532687833132
#  i      -327.38192639090408,    -286.33910177977776
#  rmse   20.8837413746,          23.0167719479
#
###################################################################################
###################################################################################
"""


def deltaphi_opt_rule_of_thumb(n, q95, betaN):
	n = float(n)
	q95 = float(q95)
	betaN = float(betaN)

	if (n!=1)and(n!=2):
		print('\t\tOnly n=1,2 supported currently. n=3,4 in preparation.')
		return

	y = q95
	x = betaN

  ## Define coefficients of 2D quadratic

  ## coefficients for vacuum optimum
	if n==1:
		a   =   0.13897915268759675
		b   =   0.15841742712121487
		c   =   -1.6740577104791932
		d   =   -0.51685976163971326
		e   =   -7.6899374778193206
		f   =   18.736665061179082
		g   =   -1.2554970387577757
		h   =   65.149596706769387
		i   =   -312.1864777106091
	if n==2:
		a   =   0.14571159718066884
		b   =   1.7141560193471825
		c   =   -6.3853893303319893
		d   =   -0.24970237007704377
		e   =   -23.718686622573557
		f   =   56.20551093242895
		g   =   -3.1507918928682539
		h   =   127.83133097817949
		i   =   -327.38192639090408

	deltaphi_opt_vac = a*(x**2)*(y**2) + b*(x**2)*y + c*(x**2) + d*x*(y**2) + e*x*y +f*x + g*(y**2) + h*y + i

  ## coefficients for total optimum
	if n==1:
		a   =   0.43304538551243965      
		b   =   -5.7000097019466098     
		c   =   17.096809988618769        
		d   =   -2.7405449274386395     
		e   =   29.940411228776266       
		f   =   -99.266691574563851          
		g   =   -0.45865516031870279      
		h   =   49.966049466296141             
		i   =   -210.17564659628684  
	if n==2:
		a   =   0.14046757326169945    
		b   =   1.7732027682693108    
		c   =   -8.5336264975127207  
		d   =   -0.33718746109109077
		e   =   -22.025319804796972
		f   =   63.892089079560563
		g   =   -3.1757270507359081
		h   =   129.06532687833132
		i   =   -286.33910177977776

	deltaphi_opt_tot = a*(x**2)*(y**2) + b*(x**2)*y + c*(x**2) + d*x*(y**2) + e*x*y +f*x + g*(y**2) + h*y + i

	shift = deltaphi_opt_tot - deltaphi_opt_vac
	
#	print('Inputs: q95 = %0.3f, betaN = %0.3f, n = %d'%(q95,betaN,n))
#	print('Outputs:\n\tVacuum optimum differential phase = %0.1f\n\tTotal optimum differential phase = %0.1f\n\tShift =%0.1f'%(deltaphi_opt_vac,deltaphi_opt_tot,shift))

	return np.array(deltaphi_opt_tot),np.array(deltaphi_opt_vac),np.array(shift)

#if __name__=='__main__':
#  deltaphi_opt_rule_of_thumb(n=sys.argv[1],q95=sys.argv[2],betaN=sys.argv[3])




