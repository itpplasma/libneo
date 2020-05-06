#!/usr/bin/env python


#import gourdon

#gd=gourdon.gourdon()
#gd.makeInput()
#gd.makeGrid()
#gd.doCalc()


import numpy as np
import IPython
import os
import matplotlib.pyplot as plt

import kk_mwillens as kk
import EQU

from ctypes import *
import dd_20140403 as dd2

#exe = '/afs/ipp-garching.mpg.de/home/m/mwillens/Programme/Gourdon/GRID/GRIDGENERATOR_edit/ggemc3hydra' 
exe = '/afs/ipp-garching.mpg.de/home/d/dbrida/GRIDGENERATOR/ggemc3hydra'

class gourdon:
    def __init__( self ,  folder='/afs/ipp-garching.mpg.de/home/m/mwillens/Programme/Gourdon/GRID/',exp='AUGD', diag='EQH', shot=30839, ed =0L, time=3.2  ):

        self.Ithres=100.0 # below this value a Imp current is assumed exactly =0

        try:
            self.Ithres=100.0# below this value a Imp current is assumed exactly =0.0
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

        except:
            self.StatusInit = False

    def makeInput( self ):       
        if self.StatusInit:
            self.writePsi()
            self.writeBt()
            self.writeCoilCurrent( )
            self.StatusInput = True

###options trace, poincare , connectionlength
    def doCalc( self, option='poincare' ):
        if self.StatusInit:

            filename = self.folder+'tmp_calc'
            f=open(filename,'w')
            f.write(header1)
            f.write("%.5f"%self.Btin+" %.3f \n"%self.Rin)
            f.write(header2)
            f.write("2 \n")
            f.write(header3) 
            f.write("SCALEB\n")
            f.write("1.d0 1.d0 1.d0\n")
            if option == 'trace':
                
                f.write("TRACE\n")
                f.write("3\n")
                f.write("200. 0.0 0.0 1 500 5")
                f.write("210. 0.0 0.0 1 500 5")
                f.write("220. 0.0 0.0 1 500 5")
                f.write("FIN\n")

            if option == 'poincare':
                f.write("POINCARE\n")
                f.write("100 1000 0.4 0.98 0 360 200  \n")
                f.write("poincare.txt\n")
                f.write("FIN\n")

##/afs/ipp-garching.mpg.de/home/m/mwillens/Programme/Gourdon/GRID/GRIDGENERATOR_edit/ggemc3hydra < tmp_calc

            if option == 'connectionlength':
                f.write("CONNECTIONLENGTH\n");
                f.write("startpoints_2D.txt");
                f.write("Lc_2D.txt");
                f.write("FIN\n")

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
            f.write(header3)          
        #    f.write(header4)
            f.close()

            os.system("cd "+self.folder+"; "+exe+" < "+filename )
            
            self.StatusGrid = True

          #  os.system("cd "+self.folder+"; ./run_grid")
            

    def plotPoincare( self ):
        try:
            fname = self.folder + 'poincare.txt'
            f=open(fname,'r')
        except:
            return


        poinc_info=np.fromfile(f, dtype=int, count=2, sep=' ')
        Nstpts=poinc_info[0]+1
        Nperiod=poinc_info[1]
        poinc_data=np.fromfile(f, dtype= float, count=Nstpts*Nperiod*3, sep=' ').reshape([Nstpts,Nperiod,3])
        Nlast=Nstpts-1
        for stpt in range(Nlast):
            Rpoinc=poinc_data[stpt,:,0]/100.
            zpoinc=poinc_data[stpt,:,1]/100.
            plt.plot(Rpoinc, zpoinc, ".k", markersize=0.5)
            

        plt.xlabel('R [m]', fontsize=22)
        plt.ylabel('z [m]',fontsize=22)
#ax.axis('equal')
        plt.xlim([0.7,2.6])
        plt.ylim([-1.7,1.5])

        plt.show()


    def writePsi( self ):

        if self.StatusInit:

            try:
                ### routines from mwillens
                
                EQ=EQU.EQU()
                EQ.Load(self.shot,Experiment=self.exp,Diagnostic=self.diag,Edition=self.ed)
                Ri = EQ.R.T
                zj = EQ.z.T
                N = EQ.Nz
                M = EQ.NR
                time = EQ.time
                tshot_idx=(np.abs(time-self.time)).argmin() #time index of time neares to tshot
                PFL_tshot=EQ.PsiOrigin[tshot_idx].T

            except:

                import dd_20140403 as dd2

                sf=dd2.shotfile()
                sf.Open(self.diag,self.shot)

                Ri=sf.GetSignal('Ri')
                zj=sf.GetSignal('Zj')
                M=len(Ri)
                N=len(zj)
    
                time=sf.GetSignal('time')
                tshot_idx=(np.abs(time-self.time)).argmin() #time index of time neares to tshot
                
                IPython.embed()
                PFL=sf.GetSignal('PFM')
                PFL_tshot=PFL[:,:,tshot_idx]

                sf.Close()

            IPython.embed()
### write Psi
            #contour,pfm,Ri,zj,nlevels=100,/isotropic
            file=open(self.folder+'psi.txt','w')

            file.write('          %d         %d   ! AUG ' %(M, N) +self.diag+' #%d @   %1.3f s\n' %  (self.shot, self.time))
            file.write('     %1.6f      %1.5f     %1.5f      %1.5f\n ' % (Ri[0,0],Ri[M-1,0],zj[0,0],zj[N-1,0]))

            PFL_tshot1=PFL_tshot.T.reshape((1,-1)).T
            idim=6
            jdim=N*M/6

            PFL_tshot1a=np.reshape(PFL_tshot1[0:jdim*idim],(-1,6))
            np.savetxt(file,PFL_tshot1a,fmt='%1.8f',newline='\n ')
            np.savetxt(file,PFL_tshot1[jdim*idim::].T,fmt='%1.8f')

            file.close()



    def writeBt( self ):
        
        if self.StatusInit:
            R = 1.1
            z = 0.0

            try:
                
                output = kk.KK().kkrzBrzt(self.shot,self.time,[R],[z])
                print('Bt=%fT at R=%fcm' % (output.bt, R*100))
                print("don't forget to include this value in your inputfile")
                self.Btin = output.bt
                self.Rin = R*100.

            except:

                
            
                if sizeof(c_long) == 8:
                    libkkso = '/afs/ipp/aug/ads/lib64/@sys/libkk.so'
                else:
                    libkkso = '/afs/ipp/aug/ads/lib/@sys/libkk.so'

                libkk = cdll.LoadLibrary(libkkso)

                print('############## read toroidal magnetic field on axis ################')
                R=c_float(R)
                z=c_float(z)
#nshot_c=c_int(nshot)
                tshot_c=c_float(self.time)
                error=c_int(0); Br=c_float(0.); Bz=c_float(0.); Bt=c_float(0.); fpf=c_float(0.);jpol=c_float(0.);edition=c_int(0); #initialize variables for kkrzbrzt
                lin=c_int(1)

                res=libkk.kkrzbrzt(byref(error), c_char_p(self.exp), c_char_p(self.diag), self.shot, byref(self.ed), byref(tshot_c), byref(R), byref(z), lin,  byref(Br), byref(Bz), byref(Bt), byref(fpf), byref(jpol))

                print('Bt=%fT at R=%fcm' % (Bt.value, R.value*100))
                print("don't forget to include this value in your inputfile")
                self.Btin = Bt.value
                self.Rin = R.value*100 

    def writeCoilCurrent( self ):
        if self.StatusInit:

            print('############## write coils definition file ################')
            sf=dd2.shotfile()
            sf.Open('MAW',self.shot)

            if not sf:
                print 'ERROR OCCURED'
                exit()
            else:
                print 'COIL DIAGNOSTICS LOADED'


            countcoils=0L
            ul=['u','l']
            t=sf.GetTimebase('T-MAW')


            for iul in range(2):
                for isec in range(1,9):
                    coilname='B'+ul[iul]+str(isec)
                    s=sf.GetSignal('I'+coilname)
                    if abs(5.0*np.interp(self.time,t,s)) >= self.Ithres: countcoils=countcoils+1 
           

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
                    s = sf.GetSignal('I'+coilname)
                    Icoil=(-5.0)*np.interp(self.time,t,s)
                    if abs(Icoil) >= self.Ithres:
			f.write(str(Icoil)+'\n')
	       	 	f.write(' INPUTS/'+coilname+'.coils\n')

            sf.Close()
            f.close()



header1 ="*--- poloidal field ---------------------------------------------------------\n* give file defining the poloidal flux matrix\npsi.txt \n0 1 \n143.8  -95.6 \n*--- toroidal field --------------------------------------------------------- \n* give a toroidal field B0 at R0(cm), Bf(R)=B0*R0/R \n"


#-2.78512  160.000 \n

header2= '* You can also provide Bt through toro_field_user.f by putting put B0=0 here.\n*--- additional field ------------------------------------------------------- \n* give a file name of a binary file to read from / write to \nBprecalc.bin \n* give a switch parameter\n* 0 : no additional field \n* 1 : pre-calculate field defined in a config file given in the following line \n* 2 : restore the precalculated field from Bprecalc.mag \n'

header3 = '* attention! by default the scaling factors for the poloidal and toroidal \n* field contributions are set to 1 while that for the addtional field is set to 0 \n* if you want to change these factors call SCALEB \n*---Targets and wall ----------------------------------------------------- \n* Number of targets NPFC and number of wall elements NWALL followed by NPFC+NWALL file names \n2 1 \ninner_target_curve_mod.txt \nouter_target_curve_mod.txt \nvessel.txt \n* end of target specification  ------------------------------------------------- \n* Grids \n* poloidal surfaces through X-point(s) must be sepcified. Such surfaces for CORE    \n* and private flux region are defined by the line from O-point to X-point. \n* You need to sepcify  \n*phi-begin, end and toroidal cuts of the mesh  \n*-22.5  22.5  5  2 ! phi_a, phi_b (deg.), NitHalf, Nsegm  \n*0.0  22.5  8  1 ! phi_a, phi_b (deg.), NitHalf, Nsegm  \n0.0  180.0  7  8 ! phi_a, phi_b (deg.), NitHalf, Nsegm  \n*0.0  22.5  7  1 ! phi_a, phi_b (deg.), NitHalf, Nsegm  \n* Poloial grid points on the LCFS, grid refinement on targets  \n400  2  \n***********************************************************************\n'


#header4='\n/EOF \nln -sf F_LINES               fort.8 \n/afs/ipp-garching.mpg.de/home/m/mwillens/Programme/Gourdon/GRID/GRIDGENERATOR_edit/ggemc3hydra < temporary_file \n \nrm -rf temporary_file fort.*'
