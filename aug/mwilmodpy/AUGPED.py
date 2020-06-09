"""
Parse augped output files
"""
from IPython import embed
import numpy as np
import matplotlib.pylab as plt
import sys
sys.path.append('/afs/ipp/aug/ads-diags/common/python/lib')
sys.path.append('/afs/ipp/u/mwillens/Programme/python/lib')
import dd
import CXZ
import kk_mwillens as kk2

class AUGPED:
    def __init__(self,file_name):
        """
        Class to parse AUGPED output files
        """
        self.file_name = file_name
        self.lines = open(self.file_name,'r').readlines()
        self._parse_header()
        # Find how many profiles are saved into the file
        self.nprofiles = 0
        self.j_start_fit_profiles = []
        self.j_start_data_profiles = []
        self.profile_names = []
        for jline,line in enumerate(self.lines):
            if 'Fit to' in line:
                self.nprofiles += 1
                self.j_start_fit_profiles.append(jline)
                self.profile_names.append(line.split()[2])
            if 'Data used in fit' in line:
                self.j_start_data_profiles.append(jline)
        self._parse_profiles()

    def _parse_profiles(self):
        self.profiles = {}
        for jprofile in range(self.nprofiles-1):
            pro = self.profile_names[jprofile]
            self.profiles[pro] = {}
            self.profiles[pro]['description'] = \
                    self.lines[self.j_start_fit_profiles[jprofile]+1]
            
            try:
                nfit = self.j_start_data_profiles[jprofile]-1 -\
                       self.j_start_fit_profiles[jprofile]-4
            except:
                embed()
            jend_data = len(self.lines) if jprofile == self.nprofiles-1 else\
                    self.j_start_fit_profiles[jprofile+1]-1
            ndata = jend_data-self.j_start_data_profiles[jprofile]-5
            self.profiles[pro]['fit'] = {}
            for data in ['Rmaj','rho','psi','fit','dfit']:
                self.profiles[pro]['fit'][data] = np.zeros(nfit,dtype=np.float)
            jfit = 0
            offset = 4
            if self.lines[self.j_start_fit_profiles[jprofile]+1].split()[0] == 'Linear':
                offset = 6
            for jline in range(self.j_start_fit_profiles[jprofile]+offset,self.j_start_data_profiles[jprofile]-1):
                self.profiles[pro]['fit']['Rmaj'][jfit],self.profiles[pro]['fit']['rho'][jfit],\
                self.profiles[pro]['fit']['psi'][jfit],self.profiles[pro]['fit']['fit'][jfit],\
                self.profiles[pro]['fit']['dfit'][jfit] = self.lines[jline].split()
                if pro == 'vT':
                    self.profiles[pro]['fit']['fit'][jfit] *= self.profiles[pro]['fit']['Rmaj'][jfit]
                    self.profiles[pro]['fit']['dfit'][jfit] *= self.profiles[pro]['fit']['Rmaj'][jfit]
                jfit += 1
            self.profiles[pro]['data'] = {}
            for data in ['Rmaj','rho','psi','data','errb']:
                self.profiles[pro]['data'][data] = np.zeros(ndata,dtype=np.float)
            jdata = 0
            for jline in range(self.j_start_data_profiles[jprofile]+5,self.j_start_data_profiles[jprofile]+5+ndata):
                if self.lines[jline] == '\n':
                    break
                self.profiles[pro]['data']['Rmaj'][jdata],self.profiles[pro]['data']['rho'][jdata],\
                self.profiles[pro]['data']['psi'][jdata],self.profiles[pro]['data']['data'][jdata],\
                self.profiles[pro]['data']['errb'][jdata] = self.lines[jline].split()
                if pro == 'vT':
                    self.profiles[pro]['data']['data'][jdata] *= self.profiles[pro]['data']['Rmaj'][jdata]
                    self.profiles[pro]['data']['errb'][jdata] *= self.profiles[pro]['data']['Rmaj'][jdata]
                jdata += 1
            index = np.argsort(self.profiles[pro]['data']['Rmaj'])
            for data in ['Rmaj','rho','psi','data','errb']:
                self.profiles[pro]['data'][data] = self.profiles[pro]['data'][data][index]

    def _parse_header(self):
        #Parse Header
        jline = 0
        self.version = np.float(self.lines[jline].split()[6])
        jline += 1
        self.shot = np.int(self.lines[jline].split()[1])
        self.t1 = np.float(self.lines[jline].split()[3])
        self.t2 = np.float(self.lines[jline].split()[5])
        self.tmid = np.float(self.lines[jline].split()[7])
        jline += 2
        try:
            #ELM sync
            self.tELM1 = np.float(self.lines[jline].split()[4])*1.e-3
            self.tELM2 = np.float(self.lines[jline].split()[6])*1.e-3
            self.ELM_sync = True
            jline += 1
        except:
            #not ELM sync
            self.ELM_sync = False
        self.EQ_EXP = self.lines[jline].split()[3]
        self.EQ_DIAG = self.lines[jline].split()[6]
        jline += 1
        self.axis_r = np.float(self.lines[jline].split()[5][:-1])
        self.axis_z = np.float(self.lines[jline].split()[6][:-1])
        jline += 1
        self.sep_r = np.float(self.lines[jline].split()[6][:-1])
        self.sep_z = np.float(self.lines[jline].split()[7][:-1])
        self.jline_end_header = jline

    def plot(self,profile,ax=None,onlyData=False):
        show = False
        if ax is None:
            f = plt.figure()
            ax = f.add_subplot(111)
            show = True
        try:
            ax.errorbar(self.profiles[profile]['data']['rho'],\
                        self.profiles[profile]['data']['data'],\
                        self.profiles[profile]['data']['errb'],fmt='.')
            if not onlyData:
                ax.errorbar(self.profiles[profile]['fit']['rho'],\
                            self.profiles[profile]['fit']['fit'],\
                            self.profiles[profile]['fit']['dfit'])
        except:
            ax.plot(self.profiles[profile]['data']['rho'],\
                    self.profiles[profile]['data']['data'],'.')
            if not onlyData:
                ax.plot(self.profiles[profile]['fit']['rho'],\
                        self.profiles[profile]['fit']['fit'],lw=3)
            ax.set_xlim(self.profiles[profile]['data']['rho'].min(),\
                        self.profiles[profile]['data']['rho'].max())
            ax.set_ylim(self.profiles[profile]['data']['data'].min(),\
                        self.profiles[profile]['data']['data'].max())
        ax.set_xlabel(r'$\rho_{\rm pol}$')
        ax.set_ylabel(profile)
        if show:
            plt.show()

    def makeEr(self,useRhoFit='vP', useNimp = False,nImpShift=0.00,useSpline=False):
        requiredProfiles = ['ne','Ti','vP','vT','Te']
        for pro in requiredProfiles:
            if not pro in self.profile_names:
                print(pro,' not in the data, no Er calculation possible')
                return
        for pro in ['Er','vExB','ErNeo','vExBNeo']:
            self.profiles[pro] = {}
            self.profiles[pro]['data'] = {}
            self.profiles[pro]['fit'] = {}
            self.profile_names.append(pro)
            # Calculate dTi/dr and dne/ne
            if (useRhoFit == "vP") | (useRhoFit == "vT"):
                self.profiles[pro]['fit']['rho'] = self.profiles[useRhoFit]['fit']['rho']
                self.profiles[pro]['fit']['Rmaj'] = self.profiles[useRhoFit]['fit']['Rmaj']
                self.profiles[pro]['fit']['psi'] = self.profiles[useRhoFit]['fit']['psi']
            else:
                self.profiles[pro]['fit']['rho'] = self.profiles['vP']['fit']['rho']
                self.profiles[pro]['fit']['Rmaj'] = self.profiles['vP']['fit']['Rmaj']
                self.profiles[pro]['fit']['psi'] = self.profiles['vP']['fit']['psi']

            self.profiles[pro]['data']['rho'] = self.profiles['vP']['data']['rho']
            self.profiles[pro]['data']['Rmaj'] = self.profiles['vP']['data']['Rmaj']
            self.profiles[pro]['data']['psi'] = self.profiles['vP']['data']['psi']
            
        dr_Er = np.gradient(self.profiles['Er']['fit']['Rmaj'])
        dpsi_Er = np.gradient(self.profiles['Er']['fit']['psi'])
        dr_Ti = np.gradient(self.profiles['Ti']['fit']['Rmaj'])
        dpsi_Ti = np.gradient(self.profiles['Ti']['fit']['psi'])
        dTi = np.gradient(self.profiles['Ti']['fit']['fit'])
        dpsi_ne = np.gradient(self.profiles['ne']['fit']['psi'])
        dr_ne = np.gradient(self.profiles['ne']['fit']['Rmaj'])
        dne = np.gradient(self.profiles['ne']['fit']['fit'])
        dpsi_Te = np.gradient(self.profiles['Te']['fit']['psi'])
        dr_Te = np.gradient(self.profiles['Te']['fit']['Rmaj'])
        dTe = np.gradient(self.profiles['Te']['fit']['fit'])


        
        dTidr = np.interp(self.profiles['Er']['fit']['rho'],\
                          self.profiles['Ti']['fit']['rho'],dTi/dr_Ti)
        dnedr = np.interp(self.profiles['Er']['fit']['rho'],\
                          self.profiles['ne']['fit']['rho'],dne/dr_ne)
        dTedr = np.interp(self.profiles['Er']['fit']['rho'],\
                          self.profiles['Te']['fit']['rho'],dTe/dr_Te)

        dTidpsi = np.interp(self.profiles['Er']['fit']['rho'],\
                          self.profiles['Ti']['fit']['rho'],dTi/dpsi_Ti)
        dnedpsi = np.interp(self.profiles['Er']['fit']['rho'],\
                          self.profiles['ne']['fit']['rho'],dne/dpsi_ne)
        dTedpsi = np.interp(self.profiles['Er']['fit']['rho'],\
                          self.profiles['Te']['fit']['rho'],dTe/dpsi_Te)

        vpol =  np.interp(self.profiles['Er']['fit']['rho'],\
                          self.profiles['vP']['fit']['rho'],\
                          self.profiles['vP']['fit']['fit'])        
        vtor =  np.interp(self.profiles['Er']['fit']['rho'],\
                          self.profiles['vT']['fit']['rho'],\
                          self.profiles['vT']['fit']['fit'])
        ne = np.interp(self.profiles['Er']['fit']['rho'],\
                          self.profiles['ne']['fit']['rho'],\
                          self.profiles['ne']['fit']['fit'])
        Ti = np.interp(self.profiles['Er']['fit']['rho'],\
                          self.profiles['Ti']['fit']['rho'],\
                          self.profiles['Ti']['fit']['fit'])
        Te = np.interp(self.profiles['Er']['fit']['rho'],\
                          self.profiles['Te']['fit']['rho'],\
                          self.profiles['Te']['fit']['fit'])

        dTidrData = np.interp(self.profiles['Er']['data']['rho'],\
                              self.profiles['Ti']['fit']['rho'],dTi/dr_Ti)
        dnedrData = np.interp(self.profiles['Er']['data']['rho'],\
                              self.profiles['ne']['fit']['rho'],dne/dr_ne)
        dTedrData = np.interp(self.profiles['Er']['data']['rho'],\
                              self.profiles['Te']['fit']['rho'],dTe/dr_Te)
        
        dTidpsiData = np.interp(self.profiles['Er']['data']['rho'],\
                                self.profiles['Ti']['fit']['rho'],dTi/dpsi_Ti)
        dnedpsiData = np.interp(self.profiles['Er']['data']['rho'],\
                                self.profiles['ne']['fit']['rho'],dne/dpsi_ne)
        dTedpsiData = np.interp(self.profiles['Er']['data']['rho'],\
                                self.profiles['Te']['fit']['rho'],dTe/dpsi_Te)
        
        vtorData =  np.interp(self.profiles['Er']['data']['rho'],\
                          self.profiles['vT']['fit']['rho'],\
                          self.profiles['vT']['fit']['fit'])
        neData = np.interp(self.profiles['Er']['data']['rho'],\
                          self.profiles['ne']['fit']['rho'],\
                          self.profiles['ne']['fit']['fit'])
        TiData = np.interp(self.profiles['Er']['data']['rho'],\
                          self.profiles['Ti']['fit']['rho'],\
                          self.profiles['Ti']['fit']['fit'])
        TeData = np.interp(self.profiles['Er']['data']['rho'],\
                          self.profiles['Te']['fit']['rho'],\
                          self.profiles['Te']['fit']['fit'])
        # Determine Bt and Bp
        import kk_mcavedon as kk
        kk = kk.KK()
        out = kk.kkrhorz(self.shot,self.tmid,self.profiles['Er']['fit']['rho'])
        self.rMajGrid = out.r
        self.rMajData = np.interp(self.profiles['Er']['data']['rho'],\
                              self.profiles['Ti']['fit']['rho'],self.rMajGrid)
        
        outB = kk.kkrzbrzt(self.shot,self.tmid,out.r,out.z)
        self.Bt = outB.bt
        self.Bp = np.sqrt(outB.br**2+outB.bz**2)
        BtData = np.interp(self.profiles['Er']['data']['rho'],\
                          self.profiles['Er']['fit']['rho'],self.Bt)
        BpData = np.interp(self.profiles['Er']['data']['rho'],\
                          self.profiles['Er']['fit']['rho'],self.Bp)
        # Determine charge of measurements
        cpz = CXZ.CXZ(pulseNumber=self.shot,diagnostic=r'CPZ')
        cpz.readProfiles()
        if useNimp:

            if self.shot == 35712:
                useChCore=-5
            else:
                useChCore=-1

            ces=dd.shotfile('CES',self.shot,'AUGD',0)
            nimpCoreRaw_data=ces('nimp',tBegin=self.t1,tEnd=self.t2).data[:,:useChCore]
            nimpCoreRaw_time=ces('nimp',tBegin=self.t1,tEnd=self.t2).time
            nimpCoreRaw_err=ces('err_nimp',tBegin=self.t1,tEnd=self.t2).data[:,:useChCore]
            nimpCore_R=ces('R',tBegin=self.t1,tEnd=self.t2).data[:,:useChCore]
            nimpCore_z=ces('z',tBegin=self.t1,tEnd=self.t2).data[:,:useChCore]
            nimpCore_rhop= kk2.KK().kkrzptfn( self.shot, nimpCoreRaw_time, nimpCore_R, nimpCore_z, exp='AUGD', diag='EQH', ed=0).rho_p
            ces.close()
            nimpCoreRaw_fit=np.copy(nimpCoreRaw_data)
            nimpCoreRaw_errfit=np.copy(nimpCoreRaw_err)
            if self.shot == 35712:
                nimpCoreRaw_err[nimpCore_rhop>0.81] *= 2.
                nimpCoreRaw_errfit[nimpCore_rhop>0.81] *= 2.
                nimpCoreRaw_fit[nimpCore_rhop>0.81] *= 1.1
                nimpCoreRaw_fit[nimpCore_rhop>0.9] *= 1.15
                nimpCoreRaw_fit[(nimpCore_rhop>0.84) & (nimpCore_rhop<0.9)] *= 1.1
                fac=0.7
            else:
                fac=1.0


            if self.shot == 35712:
                delChannel=9
                useChEdge=7
            else:
                delChannel=[]
                useChEdge=0

            cms=dd.shotfile('CMS',self.shot,'AUGD',0)
            nimpEdgeRaw_data=np.delete(cms('nimp',tBegin=self.t1,tEnd=self.t2).data[:,useChEdge:]*fac,delChannel-useChEdge,axis=1)
            nimpEdgeRaw_time=cms('nimp',tBegin=self.t1,tEnd=self.t2).time
            nimpEdgeRaw_err=np.delete(cms('err_nimp',tBegin=self.t1,tEnd=self.t2).data[:,useChEdge:]*fac,delChannel-useChEdge,axis=1)
            nimpEdge_R=np.delete(cms('R',tBegin=self.t1,tEnd=self.t2).data[:,useChEdge:],delChannel-useChEdge,axis=1)
            nimpEdge_z=np.delete(cms('z').data[:,useChEdge:],delChannel-useChEdge,axis=1)
            nimpEdge_rhop= kk2.KK().kkrzptfn( self.shot, nimpEdgeRaw_time, nimpEdge_R, nimpEdge_z, exp='AUGD', diag='EQH', ed=0).rho_p
            cms.close()
            
            if self.shot == 35712:
                delChannel=[]
                useChEdge2=7
            else:
                delChannel=[]
                useChEdge2=0
                
            cps=dd.shotfile('CPS',self.shot,'AUGD',0)
            nimpEdge2Raw_data = cps('nimp',tBegin=self.t1,tEnd=self.t2).data[:,useChEdge2:]
            nimpEdge2Raw_time =  cps('nimp',tBegin=self.t1,tEnd=self.t2).time
            nimpEdge2Raw_err = cps('err_nimp',tBegin=self.t1,tEnd=self.t2).data[:,useChEdge2:]
            nimpEdge2_R=cps('R').data[:,useChEdge2:]
            nimpEdge2_z=cps('z').data[:,useChEdge2:]
            nimpEdge2_rhop= kk2.KK().kkrzptfn( self.shot, nimpEdge2Raw_time, nimpEdge2_R, nimpEdge2_z, exp='AUGD', diag='EQH', ed=0).rho_p
            cps.close()

            
            nimpEdgeZeros_rhop=np.arange(1.02,1.4,0.2)
            nimpEdgeZeros=np.linspace(0.05,0.0,nimpEdgeZeros_rhop.size)
            nimpEdgeOnes=np.ones_like(nimpEdgeZeros_rhop)

            #embed()
            X1=np.append(np.append(np.append(np.ravel(nimpCore_rhop),np.ravel(nimpEdge_rhop)),np.ravel(nimpEdge2_rhop)),nimpEdgeZeros_rhop)
            y1=np.append(np.append(np.append(np.ravel(nimpCoreRaw_fit),np.ravel(nimpEdgeRaw_data)),np.ravel(nimpEdge2Raw_data))/1.e17,nimpEdgeZeros)
            dy1=np.append(np.append(np.append(np.ravel(nimpCoreRaw_errfit),np.ravel(nimpEdgeRaw_err)),np.ravel(nimpEdge2Raw_err))/1.e17,nimpEdgeOnes)

            potFac=10
            xx=np.append(-(X1**potFac),X1**potFac)
            yy=np.append(y1,y1)
            dydy=np.append(dy1,dy1)
            sort_index = np.where(yy>0.005 )[0][np.argsort(xx[yy>0.005])]
            
           #embed()
            self.profiles['ni'] = {}
            self.profiles['ni']['data'] = {}
            self.profiles['ni']['data']['rho'] =  np.append(np.append(np.ravel(nimpCore_rhop),np.ravel(nimpEdge_rhop)),np.ravel(nimpEdge2_rhop))
            self.profiles['ni']['data']['data'] = np.append(np.append(np.ravel(nimpCoreRaw_data),np.ravel(nimpEdgeRaw_data)),np.ravel(nimpEdge2Raw_data))
            self.profiles['ni']['data']['errb'] = np.append(np.append(np.ravel(nimpCoreRaw_err),np.ravel(nimpEdgeRaw_err)),np.ravel(nimpEdge2Raw_err))
            self.profiles['ni']['fit'] = {}
            self.profiles['ni']['fit']['rho'] = self.profiles['Ti']['fit']['rho']
 
            if useSpline:
                
                knots= np.linspace(xx[sort_index].min()+0.05,xx[sort_index].max()-0.05,10)
                from scipy.interpolate import LSQUnivariateSpline
                #embed()
                LSSpline = LSQUnivariateSpline(xx[sort_index],yy[sort_index],knots,w=1/dydy[sort_index]**2.,k=5)
                
                self.profiles['ni']['fit']['fit'] = np.squeeze(LSSpline(self.profiles['Ti']['fit']['rho']**potFac))*1.e17
                #self.profiles['ni']['fit']['errb'] = np.squeeze(sigma1)*1.e17*2.
                
                #plt.plot(X1,y1,'ro') 
                #plt.plot(self.profiles['Ti']['fit']['rho'],LSSpline(self.profiles['Ti']['fit']['rho']**potFac))
                #plt.show()
                embed()

            else:
                
                from sklearn.gaussian_process import GaussianProcessRegressor
                from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C

                kernel = C(1.0, (5e-1, 1e1)) * RBF(1, (2.75e-1, 1e2))
                gp1 = GaussianProcessRegressor(kernel=kernel, alpha=dydy[sort_index]**2., n_restarts_optimizer=10)
                gp1.fit(xx[sort_index].reshape(-1, 1), yy[sort_index])
                y_pred1, sigma1 = gp1.predict(self.profiles['Ti']['fit']['rho'].reshape(-1, 1)**potFac, return_std=True)
                
                self.profiles['ni']['fit']['fit'] = np.squeeze(y_pred1)*1.e17
                self.profiles['ni']['fit']['errb'] = np.squeeze(sigma1)*1.e17*4
            
            nimp = np.interp(self.profiles['Er']['fit']['rho'],self.profiles['ni']['fit']['rho'],self.profiles['ni']['fit']['fit'])
            dnimpdr = np.gradient(nimp)/dr_Er
            dnimpdpsi = np.gradient(nimp)/dpsi_Er

            dnimpdrData = np.interp(self.profiles['Er']['data']['rho'],\
                              self.profiles['Er']['fit']['rho'],dnimpdr)
        
            dnimpdpsiData =np.interp(self.profiles['Er']['data']['rho'],\
                                     self.profiles['Er']['fit']['rho'],dnimpdpsi)
        
            nimpData =np.interp(self.profiles['Er']['data']['rho'],\
                                     self.profiles['Er']['fit']['rho'],nimp)
            #cms = CXZ.CXZ(pulseNum
            #cms = CXZ.CXZ(pulseNum
            #cms = CXZ.CXZ(pulseNumber=self.shces('nimp',tBegin=self.t1,tEnd=self.t2)
            #embed()
        else:
            nimp=ne
            dnimpdr = dnedr
            dnimpdpsi = dnedpsi
            dnimpdrData = dnedrData
            nimpData = neData
            dnimpdpsiData = dnedpsiData
            

######## CMZ
        BpData_2 = np.interp(self.profiles['vT']['data']['rho'],\
                          self.profiles['Er']['fit']['rho'],self.Bp)
        self.profiles['vtorBpol_2'] = {}
        self.profiles['vtorBpol_2']['data'] = {}
        self.profiles['vtorBpol_2']['data']['rho'] = self.profiles['vT']['data']['rho']
        self.profiles['vtorBpol_2']['data']['data'] = self.profiles['vT']['data']['data']*BpData_2
        self.profiles['vtorBpol_2']['data']['errb'] = self.profiles['vT']['data']['errb']*BpData_2
        self.profiles['vtorBpol_2']['fit'] = {}
        self.profiles['vtorBpol_2']['fit']['rho'] = self.profiles['vT']['fit']['rho']
        self.profiles['vtorBpol_2']['fit']['fit'] = vtor*self.Bp

############## ER
        self.profiles['Er']['fit']['fit'] = (dTidr+Ti*dnimpdr/nimp)/cpz.charge-\
                self.Bt*vpol +self.Bp*vtor 
        self.profiles['ErNeo']['fit']['fit'] = dTidr+Ti*dnedr/ne
 
        self.profiles['Er']['data']['data'] = (dTidrData+TiData*dnimpdrData/nimpData)/cpz.charge-\
                BtData*self.profiles['vP']['data']['data']+BpData*vtorData
        self.profiles['ErNeo']['data']['data'] = (dTidrData+TiData*dnedrData/neData)
        self.profiles['Er']['data']['errb'] = np.abs(BtData*self.profiles['vP']['data']['errb'])
        self.profiles['ErNeo']['data']['errb'] = BtData*0.

############## vExB        
        self.profiles['vExBNeo']['fit']['fit'] = (dTidr+Ti*dnedr/ne)/np.sqrt(self.Bt**2+self.Bp**2)
        self.profiles['vExBNeo']['data']['data'] = (dTidrData+TiData*dnedrData/neData)/np.sqrt(BtData**2+BpData**2)
        self.profiles['vExBNeo']['data']['errb'] = BtData*0.
        self.profiles['vExB']['fit']['fit'] = self.profiles['Er']['fit']['fit']/np.sqrt(self.Bt**2+self.Bp**2)
        self.profiles['vExB']['data']['data'] = self.profiles['Er']['data']['data']/np.sqrt(BtData**2+BpData**2)
        self.profiles['vExB']['data']['errb'] = self.profiles['Er']['data']['errb']/np.sqrt(BtData**2+BpData**2)

############## ER contributions



        self.profiles['vtorBpol'] = {}
        self.profiles['vtorBpol']['data'] = {}
        self.profiles['vtorBpol']['data']['rho'] = self.profiles['Er']['data']['rho']
        self.profiles['vtorBpol']['data']['data'] = vtorData*BpData
        self.profiles['vtorBpol']['fit'] = {}
        self.profiles['vtorBpol']['fit']['rho'] = self.profiles['Er']['fit']['rho']
        self.profiles['vtorBpol']['fit']['fit'] = vtor*self.Bp
        
        self.profiles['vpolBtor'] = {}
        self.profiles['vpolBtor']['data'] = {}
        self.profiles['vpolBtor']['data']['rho'] = self.profiles['Er']['data']['rho']
        self.profiles['vpolBtor']['data']['data'] = self.profiles['vP']['data']['data']*BtData
        self.profiles['vpolBtor']['data']['errb'] = self.profiles['vP']['data']['errb']*BtData
        self.profiles['vpolBtor']['fit'] = {}
        self.profiles['vpolBtor']['fit']['rho'] = self.profiles['Er']['fit']['rho']
        self.profiles['vpolBtor']['fit']['fit'] = vpol*self.Bt

        self.profiles['diagImp'] = {}
        self.profiles['diagImp']['data'] = {}
        self.profiles['diagImp']['data']['rho'] = self.profiles['Er']['data']['rho']
        self.profiles['diagImp']['data']['data'] = (dnimpdrData/nimpData*TiData+dTidrData)/cpz.charge
        self.profiles['diagImp']['data']['errb'] = BtData*0.+1.
        self.profiles['diagImp']['fit'] = {}
        self.profiles['diagImp']['fit']['rho'] = self.profiles['Er']['fit']['rho']
        self.profiles['diagImp']['fit']['fit'] = (dnimpdr/nimp*Ti+dTidr)/cpz.charge
        self.profiles['diagImp']['fit']['errb'] = self.Bt*0.+1.
        
        self.profiles['diagEl'] = {}
        self.profiles['diagEl']['data'] = {}
        self.profiles['diagEl']['data']['rho'] = self.profiles['Er']['data']['rho']
        self.profiles['diagEl']['data']['data'] = (dnedrData/neData*TeData+dTedrData)
        self.profiles['diagEl']['fit'] = {}
        self.profiles['diagEl']['fit']['rho'] = self.profiles['Er']['fit']['rho']
        self.profiles['diagEl']['fit']['fit'] = (dnedr/ne*Te+dTedr)
        
        self.profiles['vDiagEl'] = {}
        self.profiles['vDiagEl']['data'] = {}
        self.profiles['vDiagEl']['data']['rho'] = self.profiles['Er']['data']['rho']
        self.profiles['vDiagEl']['data']['data'] = (dnedrData/neData*TeData+dTedrData)/np.sqrt(BtData**2.+BpData**2.)
        self.profiles['vDiagEl']['fit'] = {}
        self.profiles['vDiagEl']['fit']['rho'] = self.profiles['Er']['fit']['rho']
        self.profiles['vDiagEl']['fit']['fit'] = (dnedr/ne*Te+dTedr)/np.sqrt(self.Bt**2.+self.Bp**2.)
         
        self.profiles['vEperp'] = {}
        self.profiles['vEperp']['data'] = {}
        self.profiles['vEperp']['data']['rho'] = self.profiles['Er']['data']['rho']
        self.profiles['vEperp']['data']['data'] = self.profiles['vExB']['data']['data'] + self.profiles['vDiagEl']['data']['data']
        
        self.profiles['vEperp']['fit'] = {}
        self.profiles['vEperp']['fit']['rho'] = self.profiles['Er']['fit']['rho']
        self.profiles['vEperp']['fit']['fit'] = self.profiles['vExB']['fit']['fit'] + self.profiles['vDiagEl']['fit']['fit']

        
        self.profiles['omegaExB'] = {}
        self.profiles['omegaExB']['data'] = {}
        self.profiles['omegaExB']['data']['rho'] = self.profiles['Er']['data']['rho']

        if self.profiles['Er']['data']['rho'][0] == 0.0:
            addfac=0.1
        else:
            addfac=0.0

            self.profiles['omegaExB']['data']['data'] = self.profiles['Er']['data']['data']/np.abs((BpData+addfac)*self.rMajData)
            self.profiles['omegaExB']['data']['errb'] = np.abs(BtData*self.profiles['Er']['data']['errb'])/np.abs((BpData+addfac)*self.rMajData)

        self.profiles['omegaExB']['fit'] = {}
        self.profiles['omegaExB']['fit']['rho'] = self.profiles['Er']['fit']['rho']
        self.profiles['omegaExB']['fit']['fit'] = self.profiles['Er']['fit']['fit']/np.abs((self.Bp+addfac)*self.rMajGrid)
        self.profiles['omegaExB']['fit']['errb'] = self.Bt*0.+1.

        self.profiles['omegaDiagEl'] = {}
        self.profiles['omegaDiagEl']['data'] = {}
        self.profiles['omegaDiagEl']['data']['rho'] = self.profiles['Er']['data']['rho']
        self.profiles['omegaDiagEl']['data']['data'] = self.profiles['vDiagEl']['data']['data']*np.sqrt(BtData**2.+BpData**2.)/np.abs(BpData*self.rMajData) #(dnedpsiData/neData*TeData+dTedpsiData)
        self.profiles['omegaDiagEl']['data']['errb'] = BtData*0.+1.
        self.profiles['omegaDiagEl']['fit'] = {}
        self.profiles['omegaDiagEl']['fit']['rho'] = self.profiles['Er']['fit']['rho']
        self.profiles['omegaDiagEl']['fit']['fit'] = self.profiles['vDiagEl']['fit']['fit']*np.sqrt(self.Bt**2.+self.Bp**2.)/np.abs(self.Bp*self.rMajGrid)  #(dnedpsi/ne*Te+dTedpsi)
        self.profiles['omegaDiagEl']['fit']['errb'] = self.Bt*0.+1.

        self.profiles['omegaDiagEl2'] = {}
        self.profiles['omegaDiagEl2']['data'] = {}
        self.profiles['omegaDiagEl2']['data']['rho'] = self.profiles['Er']['data']['rho']
        self.profiles['omegaDiagEl2']['data']['data'] = self.profiles['diagEl']['data']['data'] /np.abs(BpData*self.rMajData)
        self.profiles['omegaDiagEl2']['data']['errb'] = BtData*0.+1.
        self.profiles['omegaDiagEl2']['fit'] = {}
        self.profiles['omegaDiagEl2']['fit']['rho'] = self.profiles['Er']['fit']['rho']
        self.profiles['omegaDiagEl2']['fit']['fit'] = self.profiles['diagEl']['fit']['fit'] /np.abs(self.Bp*self.rMajGrid)
        self.profiles['omegaDiagEl2']['fit']['errb'] = self.Bt*0.+1.

        
        self.profiles['omegaDiagImp'] = {}
        self.profiles['omegaDiagImp']['data'] = {}
        self.profiles['omegaDiagImp']['data']['rho'] = self.profiles['Er']['data']['rho']
        self.profiles['omegaDiagImp']['data']['data'] = (dnimpdpsiData/nimpData*TiData+dTidpsiData)/cpz.charge
        self.profiles['omegaDiagImp']['data']['errb'] = BtData*0.+1.
        self.profiles['omegaDiagImp']['fit'] = {}
        self.profiles['omegaDiagImp']['fit']['rho'] = self.profiles['Er']['fit']['rho']
        self.profiles['omegaDiagImp']['fit']['fit'] = (dnimpdpsi/nimp*Ti+dTidpsi)/cpz.charge
        self.profiles['omegaDiagImp']['fit']['errb'] = self.Bt*0.+1.
      

        
        self.profiles['omegaEperp'] = {}
        self.profiles['omegaEperp']['data'] = {}
        self.profiles['omegaEperp']['data']['rho'] = self.profiles['Er']['data']['rho']
        self.profiles['omegaEperp']['data']['data'] = self.profiles['omegaExB']['data']['data'] + self.profiles['omegaDiagEl']['data']['data']
        self.profiles['omegaEperp']['data']['errb'] = self.profiles['omegaExB']['data']['errb']
        self.profiles['omegaEperp']['fit'] = {}
        self.profiles['omegaEperp']['fit']['rho'] = self.profiles['Er']['fit']['rho']
        self.profiles['omegaEperp']['fit']['fit'] = self.profiles['omegaExB']['fit']['fit'] + self.profiles['omegaDiagEl']['fit']['fit']
        self.profiles['omegaEperp']['fit']['errb'] = self.Bt*0.+1. 
   
        

      

if __name__ == "__main__":
    aped = AUGPED('/afs/ipp/u/mcavedon/python/er/HDL/augped_fits_34955_2.000-3.000s_pre.dat')
    aped.makeEr()
    #aped.plot('Er')

