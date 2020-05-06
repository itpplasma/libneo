import numpy as np
import matplotlib.pylab as plt
import sys
sys.path.append('/afs/ipp/aug/ads-diags/common/python/lib')
sys.path.append('/afs/ipp/u/mcavedon/repository/python')
import os
from ddPlus import ddPlus as dd
import ddPlus.utils.Style, ddPlus.utils.Ipsh

class CXZ(dd.shotfile):
    def __init__(self, diagnostic=None, pulseNumber=None,\
                 experiment='AUGD', edition=0):
        self.shot = pulseNumber
        self.diagnostic = diagnostic
        self.pulseNumber = pulseNumber
        self.experiment = experiment
        self.edition = edition
        self.shift = 0
        self.tmid = 0
        super(CXZ, self).__init__(diagnostic=diagnostic,\
                                   pulseNumber=pulseNumber,\
                                   experiment=experiment,\
                                   edition=edition)


    def readProfiles(self,filter=None):
        self.shift = 0
        if self.diagnostic.decode('utf-8')[-1] == 'Z':
            self.profiles = ['Ti_c','err_Ti_c','vrot','err_vrot',\
                             'inte','err_inte']
            self.nch = self.getParameter('LOSInfo','Number').data
            losNames = self.getParameter('LOSInfo','Name').data[:self.nch]
            self.losNames = []
            for jch in range(self.nch):
                self.losNames.append(losNames[jch].strip())
            self.losNames = np.array(self.losNames)
        elif self.diagnostic.decode('utf-8')[-1] == 'S':
            self.profiles = ['nimp','err_nimp']
            self.losNames = self.getParameter('Info','LOSName')
            self.losNames = self.losNames.data[np.where(self.losNames.data != 'NOT-USED')]
            self.nch = self.losNames.size
        for pro in self.profiles:
            self.__dict__[pro] = self(pro)
            self.__dict__[pro].data = self.__dict__[pro].data[:,:self.nch]
            self.__dict__[pro].area.data = \
                    self.__dict__[pro].area.data[0,:self.nch]
            # Set zero to nan to do not see data in the plot
            self.__dict__[pro].data[self.__dict__[pro].data==0] = np.nan
        self.time = self.__dict__[pro].time
        self.R = self('R').data[0,:self.nch]
        self.z = self('z').data[0,:self.nch]
        # Exclude channels which are always nan
        jch2take = []
        for jch in range(self.nch):
            if not np.isnan(self.__dict__[self.profiles[0]].data[:,jch]).all():
                jch2take.append(jch)
        if len(jch2take) > 0:
            jch2take = np.array(jch2take)
            for pro in self.profiles:
                self.__dict__[pro].data = self.__dict__[pro].data[:,jch2take]
                self.__dict__[pro].area.data = \
                        self.__dict__[pro].area.data[jch2take]
            self.R = self.R[jch2take]
            self.z = self.z[jch2take]
            self.nch = jch2take.size
            self.losNames = self.losNames[jch2take]
        if self.diagnostic[-1] == 'Z':
            self.cxwavel = self.getParameter('LineInfo','cxwavel').data
            self.cxline = self.getParameter('LineInfo','cxline').data
            self.wvl_cx = self.cxwavel[self.cxline-1]
            self.atom,self.charge = self.wavelegth2atom(self.wvl_cx)
            # Read passive signal relative to active CX line
            jpassive = np.where(np.delete(self.cxwavel,self.cxline-1) == self.wvl_cx)[0] # index of the passive signal
            if jpassive.size != 1:
                print('No passive fit')
            else:
                flux = self.getSignalGroup('flux')[:,:self.nch,jpassive[0]]
                self.profiles.append('passive')
                self.passive = ddPlus.signalGroup('passive',self.inte.header,flux,\
                        time=self.time,area=self.inte.area)
            if filter is not None:
                # Set to nan is error on Ti is greater than 100 eV
                index = self.err_Ti_c.data <= 0
                if 'err_Ti' in filter.keys():
                    index = np.logical_or(index,self.err_Ti_c.data > filter['err_Ti'])
                if 'min_Ti' in filter.keys():
                    index = np.logical_or(index,self.Ti_c.data < filter['min_Ti'])
                if 'err_inte' in filter.keys():
                    index = np.logical_or(index,self.err_inte.data > filter['err_inte'])
                if 'err_vrot' in filter.keys():
                    index = np.logical_or(index,self.err_vrot.data > filter['err_vrot'])
                if 'p2a_ratio' in filter.keys():
                    if jpassive.size == 1:
                        p2a_ratio = self.passive.data/self.inte.data
                        index = np.logical_or(index,p2a_ratio > filter['p2a_ratio'])
                        index = np.logical_or(index,self.passive.data <= 1.e14)
                if 'only_good' in filter.keys():
                    if filter['only_good']:
                        fit_stat = self('fit_stat')
                        index = np.logical_or(index,fit_stat.data[:,:self.nch] != 1)
                for pro in self.profiles:
                    self.__dict__[pro].data[index] = np.nan
        else:
            self.which_active = self.getParameter('LineInfo','active').data
            jactive = np.where(self.which_active == 1)[0]
            self.wvl_cx = self.getParameter('LineInfo','cxwavel').data[jactive]
            self.atom,self.charge = self.wavelegth2atom(self.wvl_cx)

    def wavelegth2atom(self,wvl):
        if np.abs(wvl-567) < 1:
            atom = 'N'
            charge = 7
        if np.abs(wvl-494.467) < 0.1:
            atom = 'B'
            charge = 6
        if np.abs(wvl-468.5) < 0.5:
            atom = 'He'
            charge = 2
        return atom,charge

    def selectTime(self, tBegin=False, tEnd=False, index=False, timeRanges=None):
        if timeRanges != None:
            timeRanges = np.atleast_2d(timeRanges)
            output = self.__call__(timeRanges[0,0],timeRanges[0,1])
            for range in timeRanges[1:]:
                output = output.append(self.__call__(range[0],range[1]))
            return output
        if not tBegin and not tEnd:
            pass
        elif tBegin == tEnd:
            index = np.argmin(np.abs(self.time - tBegin))
        else:
            index = np.arange(self.time.size)\
                    [(self.time >= tBegin)*(self.time <= tEnd)]
        output = CXZ(self.diagnostic, self.pulseNumber,\
                     self.experiment, self.edition)
        output.profiles = self.profiles
        output.time = np.atleast_1d(self.time[index])
        output.nt = output.time.size
        output.nch = self.nch
        for pro in self.profiles:
            output.__dict__[pro] = self.__dict__[pro](index=index)
            if hasattr(output.__dict__[pro].area,'rMaj'):
                output.__dict__[pro].area.rMaj = output.__dict__[pro].area.rMaj[index]
                output.__dict__[pro].area.zMaj = output.__dict__[pro].area.zMaj[index]
        output.R = self.R
        output.z = self.z
        output.losNames = self.losNames
        output.atom = self.atom
        output.charge = self.charge
        try:
            output.cxwavel = self.cxwavel
            output.cxline = self.cxline
        except:
            pass
        output.wvl_cx = self.wvl_cx
        output.shift = output.shift
        output.tmid = self.tmid
        if hasattr(self,'Bp'):
            output.Bp = self.Bp[index]
            output.Bt = self.Bt[index]
            output.gridR = self.gridR[index]
            output.gridRho = self.gridRho[index]
            output.ngrid = self.ngrid
            output.eqm = self.eqm
        if hasattr(self,'timeBaseSync'):
            output.timeBaseSync = self.timeBaseSync
        return output

    def shiftProfiles(self, shiftR=None,ExpEq='AUGD',DiagEq='EQH',EdEq=0):
        """
        Shift profile by a certain deltaRmaj. Only usable if the original timebase is 
        conserved.
        """
        # Check if the time base has been syncronized and raise error
        if self.time[0] < 0:
            raise ValueError('Shifting can be applied only with the original time base.')

        rMaj = self.Ti_c.area.rMaj+shiftR
        rhoPol = self.eqm.rz2rho(rMaj,self.Ti_c.area.zMaj,t_in=self.time)
        for pro in self.profiles:
            self.__dict__[pro].area.rMaj = rMaj
            self.__dict__[pro].area.data = rhoPol
        self.shift = shiftR

    def normalizeIntensity(self,factor=None,normalize=False,los_name=None,plot=False):
        """
        Normalize the intensities to one
        """
        if factor != None and los_name != None:
            for jch in range(factor.size):
                jch_data = np.where(los_name[jch]==self.losNames)
                self.inte.data[:,jch_data] *= factor[jch]
                self.err_inte.data[:,jch_data] *= factor[jch]
        if normalize:
            for jt in range(self.nt):
                norm = self.inte.data[jt].max()
                if np.isnan(norm):
                    index_nan = np.where(np.isnan(self.inte.data[jt]))
                    if np.delete(self.inte.data[jt],index_nan).size == 0:
                        continue
                    norm = np.delete(self.inte.data[jt],index_nan).max()
                self.inte.data[jt] = self.inte.data[jt]/norm
                self.err_inte.data[jt] = self.err_inte.data[jt]/norm
        if plot:
            f = plt.figure()
            ax = f.add_subplot(111)
            colors = Style.colors(self.nch)
            for jch in range(self.nch):
                ax.plot(self.inte.area.data[:,jch],\
                        self.inte.data[:,jch],'o',label=self.losNames[jch])
            ax.legend()
            plt.show()

    def fluxCoordinate(self,ExpEq='AUGD',DiagEq='EQH',EdEq=0,fast_correction=False):
        """
        Map mesuraments onto flux surfaces
        v 1.7 -> switch to map_equ.py
        """
        import map_equ_test
        self.eqm = map_equ_test.equ_map()
        self.nt = self.time.size
        if not self.eqm.Open(self.shot,diag=DiagEq,exp=ExpEq,ed=EdEq):
            return
        self.Br,self.Bz,self.Bt = self.eqm.rz2brzt(self.R,self.z,t_in=self.time)
        # Take only the diag of the r,z grid
        self.Br = np.diagonal(self.Br,axis1=1,axis2=2)
        self.Bt = np.diagonal(self.Bt,axis1=1,axis2=2)
        self.Bz = np.diagonal(self.Bz,axis1=1,axis2=2)
        self.Bp = np.sqrt(self.Br**2+self.Bz**2)

        rhoPol = self.eqm.rz2rho(self.R,self.z,t_in=self.time)
        self.ngrid = 100
        # +/- 0.1 to avoid problems due to shifts in R
        rho = np.linspace(np.max([rhoPol.min()-0.1,0.001]),rhoPol.max()+0.1,num=self.ngrid)
        self.gridRho = np.tile(rho,self.nt).reshape((self.nt,self.ngrid))
        self.gridR,_ = self.eqm.rhoTheta2rz(rho,0,t_in=self.time)
        self.gridR = self.gridR[:,0,:]
        rMaj,zMaj = self.eqm.rhoTheta2rz(rhoPol,0,t_in=self.time)
        rMaj = rMaj[:,0,]
        zMaj = zMaj[:,0,]
        
        # ---- Correction for fast diagnostics --- 

        if fast_correction:
            # Read FPG:Raus
            fpg = ddPlus.shotfile('FPG',self.shot)
            raus = fpg('Raus')
            for jt in range(self.nt):
                jtEQM = np.argmin(np.abs(self.eqm.t_eq-self.time[jt]))
                jtRausAtEq = np.argmin(np.abs(self.eqm.t_eq[jtEQM]-raus.time))
                # Raus at the time point on which the eq is used
                rausAtEq = raus.data[jtRausAtEq]
                # average of raus during exposure time
                if jt == 0:
                    dtPlus = self.time[jt+1]-self.time[jt]
                    dtMinus = dtPlus
                elif jt == self.nt-1:
                    dtMinus = self.time[jt]-self.time[jt-1]
                    dtPlus = dtMinus
                else:
                    dtPlus = self.time[jt+1]-self.time[jt]
                    dtMinus = self.time[jt]-self.time[jt-1]
                tBegInt = self.time[jt] - dtMinus/2.
                tEndInt = self.time[jt] + dtPlus/2.
                indexAvRaus = np.arange(raus.time.size)[\
                        (raus.time > tBegInt)*(raus.time < tEndInt)]
                rausAvExp = raus.data[indexAvRaus].mean()
                rMaj[jt] += -rausAtEq + rausAvExp
            rhoPol = self.eqm.rz2rho(rMaj,zMaj,t_in=self.time)

        # ------------------------------------------------------------------------

        for pro in self.profiles:
            self.__dict__[pro].area.data = rhoPol
            self.__dict__[pro].area.rMaj = rMaj
            self.__dict__[pro].area.zMaj = zMaj

    def average(self,axis=0):
        for pro in self.profiles:
            self.__dict__[pro].time = np.average(\
                        self.__dict__[pro].time,axis=axis)
            self.__dict__[pro].data = np.average(\
                        self.__dict__[pro].data,axis=axis)
            try:
                self.__dict__[pro].area.data = np.average(\
                            self.__dict__[pro].area.data,axis=axis)
                self.__dict__[pro].area.rMaj = np.average(\
                            self.__dict__[pro].area.rMaj,axis=axis)
            except:
                pass
        try:
            self.Bt = np.average(self.Bt,axis=axis)
            self.Bp = np.average(self.Bp,axis=axis)
        except:
            pass

    def plotTimetraces(self,bin=None,xlim=None):
        from utils import OnClick,Style
        from matplotlib.ticker import MaxNLocator
        nrow = 3
        ncol = 1
        f = plt.figure()
        f.canvas.mpl_connect('button_press_event',OnClick.on_click)
        colors = Style.colors(self.nch)
        for jtrack in range(self.nch):
            jplt = 1
            ax = f.add_subplot(nrow,ncol,jplt)
            jplt += 1
            ax.errorbar(self.inte.time,self.inte.data[:,jtrack],\
                        self.err_inte.data[:,jtrack],\
                        label='Inte',color=colors[jtrack])
            ax.set_ylabel('$Inte [ph/..]$')
            ax.yaxis.set_major_locator(MaxNLocator(4))
            ax = f.add_subplot(nrow,ncol,jplt,sharex=ax)
            jplt += 1
            ax.errorbar(self.Ti_c.time,self.Ti_c.data[:,jtrack],\
                        self.err_Ti_c.data[:,jtrack],\
                        label='Ti',color=colors[jtrack])
            ax.set_ylabel('$T_i [eV]$')
            ax.yaxis.set_major_locator(MaxNLocator(4))
            ax = f.add_subplot(nrow,ncol,jplt,sharex=ax)
            jplt += 1
            ax.errorbar(self.vrot.time,self.vrot.data[:,jtrack]/1.e3,\
                        self.err_vrot.data[:,jtrack]/1.e3,\
                        color=colors[jtrack])
            ax.set_ylabel('$v_{rot} [km/s]$')
            ax.yaxis.set_major_locator(MaxNLocator(4))
        ax.set_xlabel('Time [s]')
        ax.xaxis.set_major_locator(MaxNLocator(4))
        plt.tight_layout(h_pad=-0.2)
        plt.show()

    def plotProfile(self,pro,x='rho',ax=None,**plotkwargs):
        show = False
        if ax is None:
            f = plt.figure()
            ax = f.add_subplot(111)
            show = True
        norm = {
            'Ti_c':1.e3,\
            'vrot':1.e3,\
            'inte':1.e16}
        ylabel = {
            'Ti_c':'$T_i$ [keV]',\
            'vrot':'$v_{rot}$ [km/s]',\
            'inte':'inte [ph/...]'}
        index = np.where(self.__dict__[pro].data.ravel() != 0)
        if x == 'rho':
            ax.errorbar(self.__dict__[pro].area.data.ravel()[index],\
                        self.__dict__[pro].data.ravel()[index]/norm[pro],\
                        self.__dict__['err_'+pro].data.ravel()[index]/norm[pro],\
                        **plotkwargs)
            ax.set_xlabel(r'$\rho_{pol}$')
        if x == 'R':
            ax.errorbar(self.__dict__[pro].area.rMaj.ravel()[index],\
                        self.__dict__[pro].data.ravel()[index]/norm[pro],\
                        self.__dict__['err_'+pro].data.ravel()[index]/norm[pro],\
                        **plotkwargs)
            ax.set_xlabel(r'$R_{maj}$')
        ax.set_ylabel(ylabel[pro])
        if show:
            plt.show()

    def excludeChannel(self,name=None,index=None):
        if name == None and index == None:
            return
        if name != None:
            losNames = self.getParameter('LOSInfo','Name').data
            jdel = np.where(name == losNames)
        if index != None:
            jdel = index
        self.nch -= len(jdel)
        for pro in self.profiles:
            self.__dict__[pro].data = np.delete(self.__dict__[pro].data,\
                                                jdel,axis=1)
            if hasattr(self.__dict__[pro].area,'rMaj'):
                self.__dict__[pro].area.data = np.delete(\
                                                self.__dict__[pro].area.data,\
                                                jdel,axis=1)
                self.__dict__[pro].area.rMaj = np.delete(\
                                                self.__dict__[pro].area.rMaj,\
                                                jdel,axis=1)
            else:
                self.__dict__[pro].area.data = np.delete(\
                                                self.__dict__[pro].area.data,\
                                                jdel,axis=0)
        self.R = np.delete(self.R,jdel,axis=0)
        self.z = np.delete(self.z,jdel,axis=0)

    def elmSync(self,timeBaseSync=np.array([]),ShotELM=None,\
                      ExpELM='AUGD',EdELM=0,full_output=False):
        from utils import Utils
        if timeBaseSync.size == 0:
            elm = ddPlus.shotfile('ELM',self.shot,ExpELM,EdELM)
            timeBaseSync = elm('t_begELM')
        if not ((timeBaseSync > self.time[0])*(timeBaseSync < self.time[-1])).any():
            print('no ELM sync possible')
            return
        for pro in self.profiles:
            self.__dict__[pro] = self.__dict__[pro].conditionalSynchronization(\
                                timeBaseSync)
            _, self.__dict__[pro].area.rMaj = Utils.ConditionalSynchronization(self.time,\
                                                       self.__dict__[pro].area.rMaj,\
                                                       timeBaseSync)
            _, self.__dict__[pro].area.zMaj = Utils.ConditionalSynchronization(self.time,\
                                                       self.__dict__[pro].area.zMaj,\
                                                       timeBaseSync)
        if hasattr(self,'Bp'):
            _,self.Bp = Utils.ConditionalSynchronization(self.time,\
                                                         self.Bp,\
                                                         timeBaseSync)
            _,self.Bt = Utils.ConditionalSynchronization(self.time,\
                                                         self.Bt,\
                                                         timeBaseSync)
            _,self.gridR = Utils.ConditionalSynchronization(self.time,\
                                                         self.gridR,\
                                                         timeBaseSync)
            _,self.gridRho = Utils.ConditionalSynchronization(self.time,\
                                                         self.gridRho,\
                                                         timeBaseSync)
        index = np.arange(timeBaseSync.size)[\
                (timeBaseSync > self.time.min())*(timeBaseSync < self.time.max())]
        self.timeBaseSync = timeBaseSync[index]
        self.time = self.Ti_c.time
        self.nt = self.time.size

    def flatten(self,time_offset=0):
        for pro in self.profiles:
            # Sort Rho_pol
            index = np.argsort(self.__dict__[pro].area.data.flatten())
            self.__dict__[pro].area.data = np.array([self.__dict__[pro].area.data.flatten()[index],])
            self.__dict__[pro].data = np.array([self.__dict__[pro].data.flatten()[index],])
            self.__dict__[pro].time = np.array([self.__dict__[pro].time.mean()\
                                                +time_offset,])
            self.__dict__[pro].area.rMaj = np.array([self.__dict__[pro].area.rMaj.flatten()[index],])
            self.__dict__[pro].area.zMaj = np.array([self.__dict__[pro].area.zMaj.flatten()[index],])
        self.Bp = np.array([self.Bp.flatten()[index],])
        self.Bt = np.array([self.Bt.flatten()[index],])
        index = np.argsort(self.gridRho.flatten())
        self.gridRho = np.array([self.gridRho.flatten()[index],])
        self.gridR = np.array([self.gridR.flatten()[index],])
        self.time = np.array([self.time.mean()+time_offset,])
        self.nt = 1
        self.nch = self.Ti_c.data.shape[1]
        self.ngrid = self.gridR.shape[1]

    def aggregate(self,dt=None,nt=None,tmid=0):
        """
        Downsample the time points and increase the number of channels. This funtion is particularly 
        usefull for ELM sync analysis. The number of channels has to be, however, constant. Hence the 
        reshaping is done on a constant number of time points.
            dt -> downsample dt window (if timebase has negative values (means ELM sync) calculate nt close to 0)
            nt -> number of frames to reform
        """
        if dt is not None:
            if self.time.min() < 0:
                j0 = np.argmin(np.abs(self.time)) # find index closest to 0
                # Calculate measurament frequency close to j0
                window = 5
                dt_sample = np.average(np.diff(self.time[j0-window:j0+window]))
                if tmid == 0:
                    tmid = self.timeBaseSync[np.int(self.timeBaseSync.size/2)]
                    print("I need to remap after aggregating. No time point given,")
                    print("I use the middle of the sync base: %.3f s"%tmid)
            else:
                dt_sample =  np.average(np.diff(self.time))
            nt = np.rint(dt/dt_sample)
            if nt == 0:
                raise ValueError("Aggregating window smaller than sample frequency")
        elif nt is None:
            return
        # I need to keep matrixes...
        to_delete = 0
        if self.nt%nt > 0:
            to_delete = self.nt%nt
        # Convert to int to avoid warnings
        to_delete = np.int(to_delete)
        self.nt = np.int(self.nt/nt)
        self.nch = np.int(nt*self.nch)
        nt = np.int(nt)
        # Reshape all the time dependent variables and increase the number of channels. Averaging over time 
        # Then sort everything
        for pro in self.profiles:
            self.__dict__[pro].area.data = self.__dict__[pro].area.data[to_delete:,:].reshape((self.nt,self.nch))
            index = np.argsort(self.__dict__[pro].area.data)
            self.__dict__[pro].area.data = self.__dict__[pro].area.data[np.arange(self.nt)[:,np.newaxis],index]
            self.__dict__[pro].data = self.__dict__[pro].data[to_delete:,:].reshape((self.nt,self.nch))
            self.__dict__[pro].data = self.__dict__[pro].data[np.arange(self.nt)[:,np.newaxis],index]
            self.__dict__[pro].time = self.__dict__[pro].time[to_delete:].reshape((self.nt,nt)).mean(axis=1)
        self.time = self.Ti_c.time
        # Remap everything starting from rhoPol
        rMaj,zMaj = self.eqm.rhoTheta2rz(self.Ti_c.area.data,0,t_in=self.time+tmid)
        rMaj = rMaj[:,0,]
        zMaj = zMaj[:,0,]
        for pro in self.profiles:
            self.__dict__[pro].area.rMaj = rMaj
            self.__dict__[pro].area.zMaj = zMaj

        self.Br,self.Bz,self.Bt = self.eqm.rz2brzt(rMaj,zMaj,t_in=self.time+tmid)
        # Take only the diag of the r,z diagonal of the grid
        self.Br = np.diagonal(self.Br,axis1=1,axis2=2)
        self.Bt = np.diagonal(self.Bt,axis1=1,axis2=2)
        self.Bz = np.diagonal(self.Bz,axis1=1,axis2=2)
        self.Bp = np.sqrt(self.Br**2+self.Bz**2)

        # +/- 0.1 to avoid problems due to shifts in R
        rho = np.linspace(self.Ti_c.area.data.min()-0.1,self.Ti_c.area.data.max()+0.1,num=self.ngrid)
        self.gridRho = np.tile(rho,self.nt).reshape((self.nt,self.ngrid))
        self.gridR,_ = self.eqm.rhoTheta2rz(rho,0,t_in=self.time+tmid)
        self.gridR = self.gridR[:,0,:]
        self.tmid = tmid


if __name__ == "__main__":
    shot = 34244
    tBeg = 3.0
    tEnd = 4.5
    shiftCPZ=-0.006
    cez = CXZ('CPZ',shot,experiment='AUGD')
    cez.readProfiles(filter={'err_Ti':50,'err_inte':1.e16,'err_vrot':10e3,\
            'p2a_ratio':0.3,'only_good':True})
    cez = cez.selectTime(tBeg,tEnd)
    cez.fluxCoordinate(fast_correction=True)
    cez.shiftProfiles(0.001)
    cez.elmSync(full_output=True)
    cez.aggregate(1e-3)
    cez.plotTimetraces()
    Ipsh.ipsh()
