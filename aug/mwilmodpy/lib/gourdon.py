#!/usr/bin/env python


#import gourdon

#gd=gourdon.gourdon()
#gd.makeInput()
#gd.makeGrid()
#gd.doCalc()
import io

import numpy as np
import IPython
import os
import matplotlib.pyplot as plt
from matplotlib import gridspec

from scipy.ndimage.interpolation import map_coordinates
from scipy.interpolate import griddata
from scipy import interpolate
import kk_mwillens as kk
import EQU
import MAWnew as MAW
import fconf

#from ctypes import *
#import dd_20140403 as dd2

exe = '/afs/ipp-garching.mpg.de/home/m/mwillens/Programme/Gourdon/GRID/GRIDGENERATOR_edit/ggemc3hydra' 
#exe = '/afs/ipp-garching.mpg.de/home/d/dbrida/GRIDGENERATOR/ggemc3hydra'

class gourdon:
    def __init__( self ,  folder='/afs/ipp-garching.mpg.de/home/m/mwillens/Programme/Gourdon/GRID/',exp='AUGD', diag='EQH', shot=30839, ed =0L, time=3.2  ):

        #self.Ithres=10.0 # below this value a Imp current is assumed exactly =0

        try:
            self.Ithres=10.0# below this value a Imp current is assumed exactly =0.0
            self.resR  =' 0.8d0     2.4d0   64'
            self.resz  ='-1.5d0     1.5d0   64'
            self.resphi=' 0.0d0   360.0d0   64'

            self.exp = exp
            self.diag = diag
            self.shot = shot
            self.ed = ed
            self.time = time
            self.folder = folder
            self.StatusInit = True
            self.StatusInput = False
            self.StatusGrid = False
            self.StatusMag = False
            self.StatusBpre=True
            self.statusLCIn = False
        except:
            self.StatusInit = False

    def makeInput( self, shotCoil = None, timeCoil = None , usePSL=False, freq = None,Icoils=None,factor=None):       
        if self.StatusInit:
            self.writePsi()
            self.writeBt()
            self.writeCoilCurrent(shot = shotCoil, time = timeCoil, usePSL = usePSL, freq=freq,IcoilsIn = Icoils,factor=factor )
            self.StatusInput = True

###options trace, poincare , connectionlength
    def doCalc( self, option='poincare' ,PoincareName="poincare", traceIn=None,phi0=0.0):
        if self.StatusInit:

            self.PoincareFile1 = PoincareName+"1"+".txt"
            self.PoincareFile2 = PoincareName+"2"+".txt"

            filename = self.folder+'tmp_calc'
            f=open(filename,'w')
            f.write(header1)
            f.write("%.5f"%self.Btin+" %.3f \n"%self.Rin)
            f.write(header2)
            f.write("2 \n")
            f.write("INPUTS/screening_currents.txt\n")
            f.write(header3) 
            f.write("SCALEB\n")
            f.write("1.d0 1.d0 1.d0\n")
            if option == 'trace':
                
                f.write("TRACE\n")
                if traceIn == None:
                   # f.write("startpoints.txt\n");
                    f.write("3\n")
                   # #        R,   z  phi dphi Npr understeps  
                    f.write("200. 0.0 0.0 1 500 5\n")
                    f.write("210. 0.0 0.0 1 500 5\n")
                    f.write("220. 0.0 0.0 1 500 5\n")
                else:
                    f.write("%s"%traceIn)

                f.write("FIN\n")

                os.system("cd "+self.folder+";ln -sf F_LINES fort.8")

            if option == 'poincare':
                ##QPROFILE
                ##qprofile.txtres
                ##FIN
                ##QFINDER
                ##3.5 1E-6 10000
                ##FIN
               # qfinder=[4.5,5.0,5.5,6.0]
               # for qval in qfinder:
               #     f.write("QFINDER\n")
               #     f.write("%.1f 1E-6 10000\n"%qval)

                f.write("POINCARE\n")
                 #* Nr,Nt,psi1,psi2,phi0,deltaphi,N_F
                f.write("100 1000 0.9 0.99 -30 %.2f -360 200 \n"%phi0)
                f.write(self.PoincareFile1+"\n")
                f.write("POINCARE\n")
                f.write("100 1000 0.9 0.99 30 %.2f 360 200 \n"%phi0)
                f.write(self.PoincareFile2+"\n")
                f.write("FIN\n")

##/afs/ipp-garching.mpg.de/home/m/mwillens/Programme/Gourdon/GRID/GRIDGENERATOR_edit/ggemc3hydra < tmp_calc

            if option == 'connectionlength':
                if self.statusLCIn:
                    f.write("CONNECTIONLENGTH\n");
                    f.write("startpoints.txt\n");
                    f.write("LC_2D.txt\n");
                    f.write("FIN\n")
                else:
                    print 'no input file generated'

            f.close()

            userexe = "cd "+self.folder+"; "+exe+" < tmp_calc"
            print userexe
            os.system(userexe)   
            
    
    def makeGrid( self ):

        if self.StatusInit & self.StatusInput:

            filename = self.folder+'tmp_grid'
            f=open(filename,'w')
            f.write(header1)
            f.write("%.5f"%self.Btin+" %.3f \n"%self.Rin)
            f.write(header2)
            f.write("1 \n")
            f.write("INPUTS/coils.conf \n")
            f.write("INPUTS/screening_currents.txt\n")
            f.write(header3)          
        #    f.write(header4)
            f.close()

            os.system("cd "+self.folder+"; "+exe+" < "+filename )
            
            self.StatusGrid = True

          #  os.system("cd "+self.folder+"; ./run_grid")
  
    def readTrace( self ):
        if self.StatusInit & self.StatusInput:
            try:
                f=open(self.folder+'out', 'r')
            except:
                print 'could not open file in read Timetrace'
                return
            
            R=[]
            z=[]
            phi=[]
            fdata=np.loadtxt(f)
            Nfl = np.int(fdata[:,0].max())
            for i in np.arange(1,Nfl+1):
                idx = np.where((fdata[:,0]) == i)[0]
                R.append(fdata[idx,1])
                z.append(fdata[idx,2])
                phi.append(fdata[idx,3])

            return R,z,phi

            #IPython.embed()
 
       ## Rin, zIn and phiIn must have the same size  
    def makeStartpoints( self,RIn=None,zIn=None,phiIn=None,direction=1, Ntor = 500000):
        if self.StatusInit & self.StatusInput:
            filenameIn= self.folder+"startpoints.txt"
            filenameOut =  self.folder+"LC_2D.txt"

            if (np.all(RIn)==None) | (np.all(zIn)==None) | (np.all(phiIn)==None):
                print 'no Input or missing input'
                return

            if (np.size(RIn)!=np.size(zIn)) | (np.size(RIn)!=np.size(phiIn)) : 
                print 'Input has not the same size'
                return            


            mesh=np.array([RIn.ravel()*100.,zIn.ravel()*100.,phiIn.ravel()]).T

            N=RIn.size
            Nsteps=5
            direction=direction
            Ntor =  Ntor 
            Ntar=1
            tar1=1


            #IPython.embed()
            fw=open(filenameIn,'w')
            fw.write('   %d    2.00000    %d     %d     %d    %d \n   %d\n' % (N, Nsteps, direction, Ntor, Ntar, tar1))
            np.savetxt(fw, mesh, fmt='%.6f')
            fw.close()

            self.LC_R = RIn.ravel()
            self.LC_z = zIn.ravel()
            self.LC_phi = phiIn.ravel()
            self.LC_filenameIn = filenameIn
            self.LC_filenameOut = filenameOut
            self.statusLCIn = True

    def readConnectionLength( self ):
         if (self.statusLCIn == True):
             print 'Connection length was calculated'
             data=np.genfromtxt(self.LC_filenameOut)
             self.LCplus = data[:,0]
             self.LCminus = data[:,1]
             self.LCPsiMinPlus = data[:,2]
             self.LCPsiMinMinus = data[:,3]

    def readBprecalc( self ):
        
   #     IPython.embed()
        import BinaryReader
        
        fileOutput = self.folder + 'Bprecalc.bin'
        br = BinaryReader.BinaryReader(fileOutput)
        tmp = br.read('int32',4)
        Nfpre,Nrpre,Nzpre,Ncomppre = tmp[0], tmp[1], tmp[2], tmp[3]
        print 'grid restored from', fileOutput
        print 'Nf,Nr,Nz=',Nfpre,Nrpre,Nzpre
        Rtmp=br.read('double',2)
        Rapre,Rbpre =  Rtmp[0],Rtmp[1]
        ztmp=br.read('double',2)
        zapre,zbpre =  ztmp[0],ztmp[1]
        ftmp=br.read('double',2)
        fapre,fbpre =  ftmp[0],ftmp[1]
        print 'valid R domain ',Rapre,'..' ,Rbpre
        print 'valid z domain ',zapre,'..' ,zbpre
        print 'valid phi domain ',fapre,'..' ,fbpre
        arr=br.read('double',Ncomppre*(Nfpre-1)*Nrpre*Nzpre)
        br.close()

     #   bfield1=np.reshape(arr,(Ncomppre,Nrpre_tmp,Nzpre_tmp,(Nfpre_tmp-1)))

     #   bfield2=np.reshape(arr,((Nfpre_tmp-1),Nrpre_tmp,Nzpre_tmp,Ncomppre))

    #    bfield3=np.reshape(arr,((Nfpre_tmp-1),Nzpre_tmp,Nrpre_tmp,Ncomppre))


        bfield=np.reshape(arr,(Nzpre,Nrpre,(Nfpre-1),Ncomppre)).T
        #br=bfield[0
             
        self.Rb=np.linspace(Rapre,Rbpre,Nrpre,endpoint=True)
        self.zb=np.linspace(zapre,zbpre,Nzpre,endpoint=True)
        self.phib=np.diff(np.linspace(zapre,zbpre,Nfpre))+np.linspace(zapre,zbpre,Nfpre-1,endpoint=False)

        self.nRb = self.Rb.size
        self.nzb = self.zb.size
        self.nphib = self.phib.size
        self.phiradRange = self.phib.max()

        self.df=(fbpre-fapre)/(Nfpre-1)
        self.dR=(Rbpre-Rapre)/(Nrpre-1)
        self.dz=(zbpre-zapre)/(Nzpre-1)
        self.Borg=bfield
        self.Bphi=np.swapaxes(bfield[0],1,2)
        self.Br=np.swapaxes(bfield[1],1,2)
        self.Bz=np.swapaxes(bfield[2],1,2)
                #  nRin = np.size(Rin)
                 # fR = interpolate.interp1d(self.R,np.arange(self.nR))
			#fz = interpolate.interp1d(self.z, np.arange(self.nz))


        self.StatusBpre=True

        return


    def getBrzphi(self, Rin=None, zin=None, phiIn=None,phiInRad=True):
		if (self.StatusBpre ):
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
			fR = interpolate.interp1d(self.Rb,np.arange(self.nRb))
			fz = interpolate.interp1d(self.zb, np.arange(self.nzb))

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
					if phiInput.max() > self.phib[-1]:
						fphi = interpolate.interp1d(np.append(self.phib,self.phiradRange), np.arange(self.nphib+1))
						usePhiExt = True
					else:
						fphi = interpolate.interp1d(self.phib, np.arange(self.nphib))
						usePhiExt = False

				else:
					phiInput = np.remainder(phiIn,self.phiRange)
					phiInput[phiInput<0] += self.phiRange
					if phiInput.max() > self.phib[-1]:
						fphi = interpolate.interp1d(np.append(self.phib,self.phiRange), np.arange(self.nphib+1))
						usePhiExt = True
					else:
						fphi = interpolate.interp1d(self.phib, np.arange(self.nphib))
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


 #       f =io.open(self.folder+'Bprecalc.bin', mode='rb', buffering=-1, encoding=None, errors=None, newline=None, closefd=True)
 #       data = f.readline()
 #       f.close()

        #try:
        #        f=open(self.folder+'out', 'r')
         #   except:
          #      print 'could not open file in read Timetrace'
           #     return
  

    def plotPoincare( self , RzDiag=None):

        try:
            fname1 = self.folder+self.PoincareFile1
            #print fname1
            f1=open(fname1,'r')
            fname2 = self.folder+self.PoincareFile2
            #print fname2
            f2=open(fname2,'r')     
        except:
            return


        self.calcMagnAxis()

        

        nPlots=3
        fig = plt.figure(figsize=(9,8))
        gs = gridspec.GridSpec(1,2)
        ax1 = fig.add_subplot(gs[0,0],aspect='auto')
        ax2 = fig.add_subplot(gs[0,1],aspect='auto')
        files = [f1,f2]
        ii=0
        str=[ ".r",".r"]
        for fi in files:
            #print fi
            poinc_info=np.fromfile(fi, dtype=int, count=2, sep=' ')
            Nstpts=poinc_info[0]+1
            Nperiod=poinc_info[1]
            poinc_data=np.fromfile(fi, dtype= float, count=Nstpts*Nperiod*5, sep=' ').reshape([Nstpts,Nperiod,5])
            Nlast=Nstpts-1
    #        for stpt in range(Nlast):
            Rpoinc=np.ravel(poinc_data[:,:,0]/100.)
            zpoinc=np.ravel(poinc_data[:,:,1]/100.)
            ax1.plot(Rpoinc, zpoinc, str[ii], markersize=0.5)

            thetaMag= np.rad2deg(np.arctan2(zpoinc-self.zmag,Rpoinc-self.Rmag))
            psi=np.ravel(poinc_data[:,:,2])
            #theta=np.ravel(poinc_data[:,:,3])
            ax2.plot(thetaMag, np.sqrt(psi),str[ii], markersize=0.5)

            ii=ii+1


        if np.all(RzDiag)!=None:
             ax1.plot(RzDiag[0].ravel(),RzDiag[1].ravel(), 'bo')
             thetaDiag = np.rad2deg(np.arctan2(RzDiag[1].ravel()-self.zmag,RzDiag[0].ravel()-self.Rmag))
             out = kk.KK().kkrzptfn( self.shot, self.time, RzDiag[0].ravel(),RzDiag[1].ravel(),exp=self.exp, diag=self.diag, ed=self.ed )
             rhopDiag = out.rho_p
             ax2.plot(thetaDiag,rhopDiag, 'bo')

        #fconf.plt_eq_pol(self.shot,self.time,rhos=[1.0])

        ax1.set_xlabel('R [m]', fontsize=22)
        ax1.set_ylabel('z [m]',fontsize=22)
#ax.axis('equal')
        ax1.set_xlim([0.7,2.6])
        ax1.set_ylim([-1.7,1.5])
        fconf.plt_vessel(ax1)

        ax2 = fig.add_subplot(gs[0,1],aspect='auto')
        ax2.set_xlabel(r'$\rm{\Theta}$ [deg]')
        ax2.set_ylabel(r'$\rm{\rho_{pol}}$')
        ax2.set_xlim([-180.,180.])
        ax2.set_ylim([0.92,1.01])

        gs.tight_layout(fig, rect=[0, 0, 1, 1])

        plt.show()
        #IPython.embed()

    def writePsi( self ):

        if self.StatusInit:

           
                ### routines from mwillens
                
            EQ=EQU.EQU()
            EQ.Load(self.shot,Experiment=self.exp,Diagnostic=self.diag,Edition=self.ed)
 
            N = EQ.Nz
            M = EQ.NR
            time = EQ.time
            tshot_idx=(np.abs(EQ.time-self.time)).argmin() #time index of time neares to tshot

            Ri = EQ.R[tshot_idx]
            zj = EQ.z[tshot_idx] 
            PFL_tshot=EQ.PsiOrigin[tshot_idx].T

            del EQ
### write Psi
            #contour,pfm,Ri,zj,nlevels=100,/isotropic
            file=open(self.folder+'psi.txt','w')

            file.write('          %d         %d   ! AUG ' %(M, N) +self.diag+' #%d @   %1.3f s\n' %  (self.shot, self.time))
            file.write('     %1.6f      %1.5f     %1.5f      %1.5f\n ' % (Ri[0],Ri[M-1],zj[0],zj[N-1]))

            PFL_tshot1=PFL_tshot.T.reshape((1,-1)).T
            idim=6
            jdim=N*M/6

            PFL_tshot1a=np.reshape(PFL_tshot1[0:jdim*idim],(-1,6))
            np.savetxt(file,PFL_tshot1a,fmt='%1.8f',newline='\n ')
            np.savetxt(file,PFL_tshot1[jdim*idim::].T,fmt='%1.8f')

            file.close()


    def calcMagnAxis( self ):
        if self.StatusInit:
            output = kk.KK().kkrhorz( self.shot, self.time, [0.0], angle=0.0,exp=self.exp, diag=self.diag, ed=self.ed )
            
            self.Rmag = output.r
            self.zmag = output.z
            
           
            self.StatusMag = True

    def writeBt( self ):
        
        if self.StatusInit:
            R = 1.1
            z = 0.0

            output = kk.KK().kkrzBrzt(self.shot,self.time,[R],[z])
            print('Bt=%fT at R=%fcm' % (output.bt, R*100))
            print("don't forget to include this value in your inputfile")
            self.Btin = output.bt
            self.Rin = R*100.

    def writeCoilCurrent( self , shot=None, time=None, usePSL=False, freq=None, IcoilsIn =None,factor=None):
        if self.StatusInit:

            print('############## write coils definition file ################')

            countcoils=0L
            ul=['u','l']

            if np.all(shot) == None:
                shot = self.shot
            else:
                shot = shot

            if np.all(time) == None:
                time = self.time
            else:
                time = time

            useFactor = False 
            useFreq = False

            if np.all(freq)==None:
                usePSL = usePSL
            else:
                usePSL = False
                useFreq = True
                freq = freq

            if np.all(factor)==None:
                usePSL = usePSL
            else:
                usePSL = False
                useFactor = True
                fac=factor
                


            if (np.all(IcoilsIn) != None) & (np.shape(IcoilsIn)==(2,8)):  
                
                Ic = np.copy((-5.0)*IcoilsIn)
            else:
                MW=MAW.MAW(shot,Binning=0.5)          
                Ic = np.copy((-5.0)*MW.getIcoils(time,usePSL=usePSL))
                if useFreq:
                    print('use PSL frequency of %.f'%freq)
                    Ic[0] = Ic[0]*MW.getUpperResponseAmp(freq)
                    Ic[1] = Ic[1]*MW.getLowerResponseAmp(freq)
                elif useFactor:
                    print('use Factor of %.f'%fac)
                    Ic[0] = Ic[0]*fac 
                    Ic[1] = Ic[1]*fac            
                del MW

            countcoils = np.size(np.where(np.abs(Ic)>=self.Ithres)[0]) 
      
            f=open(self.folder+'INPUTS/coils.conf','w')
            f.write('%d  ! coils defined for AUG # %d  @  %f s\n' % (countcoils, self.shot, self.time))
            f.write(self.resR + ' ! size of the grid in R direction\n')
            f.write(self.resz + ' ! size of the grid in z direction\n')
            f.write(self.resphi+' ! size of the grid in phi direction\n')
            f.write('12 ! components of the field + derivatives\n')
            #IPython.embed()
            for iul in range(2):
                for isec in range(1,9):
                    coilname='B'+ul[iul]+str(isec)
                    
                    Icoil=Ic[iul,isec-1]
                    #if iul ==1: Icoil=-Icoil
                    if abs(Icoil) >= self.Ithres:
			f.write(str(Icoil)+'\n')
	       	 	f.write(' INPUTS/'+coilname+'.coils\n')

            
            f.close()



header1 ="*--- poloidal field ---------------------------------------------------------\n* give file defining the poloidal flux matrix\npsi.txt \n0 1 \n143.8  -95.6 \n*--- toroidal field --------------------------------------------------------- \n* give a toroidal field B0 at R0(cm), Bf(R)=B0*R0/R \n"


#-2.78512  160.000 \n

header2= '* You can also provide Bt through toro_field_user.f by putting put B0=0 here.\n*--- additional field ------------------------------------------------------- \n* give a file name of a binary file to read from / write to \nBprecalc.bin \n* give a switch parameter\n* 0 : no additional field \n* 1 : pre-calculate field defined in a config file given in the following line \n* 2 : restore the precalculated field from Bprecalc.mag \n'

header3 = '* attention! by default the scaling factors for the poloidal and toroidal \n* field contributions are set to 1 while that for the addtional field is set to 0 \n* if you want to change these factors call SCALEB \n*---Targets and wall ----------------------------------------------------- \n* Number of targets NPFC and number of wall elements NWALL followed by NPFC+NWALL file names \n2 1 \ninner_target_curve_mod.txt \nouter_target_curve_mod.txt \nvessel.txt \n* end of target specification  ------------------------------------------------- \n* Grids \n* poloidal surfaces through X-point(s) must be sepcified. Such surfaces for CORE    \n* and private flux region are defined by the line from O-point to X-point. \n* You need to sepcify  \n*phi-begin, end and toroidal cuts of the mesh  \n*-22.5  22.5  5  2 ! phi_a, phi_b (deg.), NitHalf, Nsegm  \n*0.0  22.5  8  1 ! phi_a, phi_b (deg.), NitHalf, Nsegm  \n0.0  180.0  7  8 ! phi_a, phi_b (deg.), NitHalf, Nsegm  \n*0.0  22.5  7  1 ! phi_a, phi_b (deg.), NitHalf, Nsegm  \n* Poloial grid points on the LCFS, grid refinement on targets  \n400  2  \n***********************************************************************\n'


#header4='\n/EOF \nln -sf F_LINES               fort.8 \n/afs/ipp-garching.mpg.de/home/m/mwillens/Programme/Gourdon/GRID/GRIDGENERATOR_edit/ggemc3hydra < temporary_file \n \nrm -rf temporary_file fort.*'