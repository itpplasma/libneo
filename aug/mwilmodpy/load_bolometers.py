#!/usr/bin/env python
# -*- coding: utf-8 -*-


import matplotlib
from numpy import *
from matplotlib.pylab import *
import sys,os
sys.path.append('/afs/ipp/home/t/todstrci/TRANSP')
import dd_tomasO as dd 
import kk_tomasO as kk
import IPython
from scipy.stats.mstats import mquantiles

    
    
    
class loader_diod_bolometers:
    def __init__(self,shot,exp='AUGD',ed=0):
        self.shot = int(shot)
        self.exp = exp
        self.ed = ed

        self.phase_diag = True
       
        bolo_shotfiles = 'XVR', 'XVS'

        self.dd = dd.shotfile()
  #      IPython.embed()
        calib_shot = self.dd.cShotNr('AUGD', 'BLC',self.shot)
       
        self.dd.Open( 'BLC', calib_shot)
       
        names = [n.strip() for n in self.dd.GetNames()]

        cams = [b for b in names if b[0] == 'D' and len(b) == 3]
       
        active_los = {s:bool_(self.dd.GetParameter(s, 'active' )) for s in cams}
        activ_cam = [c for c in cams if any(active_los[c])]
        self.Phi = {s:self.dd.GetParameter(s, 'P_Blende' )[active_los[s]] for s in activ_cam}
        self.R_start = {s:self.dd.GetParameter(s, 'R_Blende' )[active_los[s]] for s in activ_cam}
        self.z_start = {s:self.dd.GetParameter(s, 'z_Blende' )[active_los[s]] for s in activ_cam}
        self.R_end   = {s:self.dd.GetParameter(s, 'R_end' )[active_los[s]] for s in activ_cam}
        self.z_end   = {s:self.dd.GetParameter(s, 'z_end' )[active_los[s]] for s in activ_cam}
        self.theta   = {s:arctan2(self.z_end[s]-self.z_start[s],self.R_end[s]-self.R_start[s])  for s in activ_cam}
        self.delta   = {s:self.dd.GetParameter(s, 'delta' )[active_los[s]] for s in activ_cam}


        sig_dict = {s:self.dd.GetParameter(s, 'RAW' ) for s in activ_cam}
        
        self.subcam_ind = {c:[self.delta[c]+self.R_start[c] == r for r in unique(self.delta[c]+self.R_start[c])] for c in activ_cam}
        

        #geometric corrections
        n =  self.R_end['DHC']-self.R_start['DHC'], self.z_end['DHC']-self.z_start['DHC']
        corr = 0.04
        self.R_end['DHC'][self.subcam_ind['DHC'][1]]+= -n[1][self.subcam_ind['DHC'][1]]*tan(corr)
        self.z_end['DHC'][self.subcam_ind['DHC'][1]]+= +n[0][self.subcam_ind['DHC'][1]]*tan(corr)
            
        n =  self.R_end['DVC']-self.R_start['DVC'], self.z_end['DVC']-self.z_start['DVC']
        corr = 0.02
        self.R_end['DVC'][self.subcam_ind['DVC'][0]]+= -n[1][self.subcam_ind['DVC'][0]]*tan(corr)
        self.z_end['DVC'][self.subcam_ind['DVC'][0]]+= +n[0][self.subcam_ind['DVC'][0]]*tan(corr)
        
            

        self.groups = activ_cam
       
        xv_names = []
        bolo_los = {}
        for bs in bolo_shotfiles:
            bolo_los[bs] = []
            self.dd.Open(bs, self.shot)
            names = self.dd.GetNames()
            xv_names.append([n.strip() for n in names])
            self.dd.Close()
           
        self.signals = {}
        for c in activ_cam:
            sig = sig_dict[c]
            activ = active_los[c]
            ch = 1+arange(len(activ))
            sf = []
                       
            for s in sig[activ]:
                s = s.strip()
                if str(s) in xv_names[0]:
                    sf.append(bolo_shotfiles[0])
                elif str(s) in xv_names[1]:
                    sf.append(bolo_shotfiles[1])
                else:
                    raise Exception('missing camera %s'%s)

            self.signals[c] = zip(sig[activ],ch[activ],array(sf))
           



    def get_signal_groups(self):
        return  self.groups

    def get_names(self,group):
        return sorted([ch for sig,ch,sf in  self.signals[group]])
    
    def get_rho_tan(self,group,names,shot,time):

        out = kk.KK().kkeqpfx( int(shot), time,diag='EQI')

        n =  (self.R_end[group]-self.R_start[group], self.z_end[group]-self.z_start[group])
        n_perp = c_[-n[1],n[0]]/hypot(n[0],n[1])[:,None]
        ddist = -sum(n_perp*(c_[self.R_start[group],self.z_start[group]]-r_[out.Raxis,out.zaxis]),1)

        t = linspace(0,1,50)
        
        X = array(n)[:,:,None]*t[None,None,:]+c_[self.R_start[group],self.z_start[group]].T[:,:,None]
       
        channels = sorted([ch for sig,ch,sf in  self.signals[group]])
        ind = in1d(channels, names)

        rho = kk.KK().kkrzptfn( int(shot),time, X[0][ind].ravel(),X[1][ind].ravel(), diag='EQI')
        rho_p = rho.rho_p.reshape(X[0][ind].shape)
        rho_tg = amin(rho_p,1)*sign(ddist[ind])
        
        Rlim = array([interp([-1,1], r,x) for r,x in zip(rho_p,X[0])])
        zlim = array([interp([-1,1], r,x) for r,x in zip(rho_p,X[1])])
        self.los_length = hypot(Rlim[:,1]-Rlim[:,0], zlim[:,1]-zlim[:,0])

        #import IPython
        #IPython.embed()
        
        return rho_tg
    
    

    def get_signal(self,group, names,calib=False,tmin=0.,tmax=10):
        if isinstance(names,str) or not  hasattr(names, '__iter__'):
            names = (names,)
        
        sf_open = None
        data = []
        for name in names:
                
            for sig,ch,sf in self.signals[group]:
                if ch == int(name):
                    break
            
            if sf_open!= sf:
                self.dd.Open(sf,self.shot, experiment=self.exp, edition=self.ed)
                sf_open = sf
                tvec = self.dd.GetTimebase('Dio-Time')
                offset_start = tvec.searchsorted(-.1)
                offset_end = tvec.searchsorted(0)
                ind = slice(tvec.searchsorted(tmin), tvec.searchsorted(tmax))
                tvec = tvec[ind]


            offset =    self.dd.GetSignal(sig,cal=calib,nbeg=offset_start+1,nend=offset_end).mean()
            data.append(self.dd.GetSignal(sig,cal=calib,nbeg=ind.start+1,nend=ind.stop)-offset)  

        #correstions of the diod bolometers

        if self.shot > 31591:
            pass
        elif self.shot > 30135:
            if group == 'DVC':
                ind = where(in1d(names, where(self.subcam_ind['DVC'][0])[0]))[0]
                for i in ind:  data[i]*= 1.31
        elif self.shot > 28523:
            pass      
        elif self.shot > 27352:
            if group == 'DHC':
                ind = where(in1d(names, where(self.subcam_ind['DHC'][1])[0]))[0]
                for i in ind:  data[i]/= 1.05
                ind = where(in1d(names, where(self.subcam_ind['DHC'][2])[0]))[0]
                for i in ind:  data[i]*= 1.3
            if group == 'DVC':
                ind = where(in1d(names, where(self.subcam_ind['DVC'][0])[0]))[0]
                for i in ind:  data[i]*= 1.1
    
       
        self.dd.Close()
        return tvec, data
    
    
class loader_foil_bolometers:
    def __init__(self,shot,exp='AUGD',ed=0):
        self.shot = int(shot)
        self.exp = exp
        self.ed = ed

        self.phase_diag = True
       
        self.dd = dd.shotfile()
       
        calib_shot = self.dd.cShotNr('AUGD', 'BLC',self.shot)
       
        self.dd.Open( 'BLC', calib_shot)
       
        names = [n.strip() for n in self.dd.GetNames()]

        cams = [b for b in names if b[0] == 'F' and len(b) == 3]
        #print bool_(self.dd.GetParameter('FHC', 'active' ))
        #exit()
        active_los = {s:bool_(self.dd.GetParameter(s, 'active' )) for s in cams}
        activ_cam = [c for c in cams if any(active_los[c])]
        self.channels = {s:arange(len(a))[a]  for s,a in  active_los.iteritems()}
        #print self.channels['FHC']
        #exit()
        self.Phi = {s:self.dd.GetParameter(s, 'P_Blende' )[active_los[s]] for s in activ_cam}
        self.R_start = {s:self.dd.GetParameter(s, 'R_Blende' )[active_los[s]] for s in activ_cam}
        self.z_start = {s:self.dd.GetParameter(s, 'z_Blende' )[active_los[s]] for s in activ_cam}
        self.R_end   = {s:self.dd.GetParameter(s, 'R_end' )[active_los[s]] for s in activ_cam}
        self.z_end   = {s:self.dd.GetParameter(s, 'z_end' )[active_los[s]] for s in activ_cam}
        self.delta   = {s:self.dd.GetParameter(s, 'delta' )[active_los[s]] for s in activ_cam}

        self.theta = {s:arctan2(self.z_end[s]-self.z_start[s],self.R_end[s]-self.R_start[s])  for s in activ_cam}
        self.subcam_ind = {c:[self.R_start[c] == r for r in unique(self.R_start[c])] for c in activ_cam}

        self.groups = activ_cam



    def get_signal_groups(self):
        return  self.groups

    def get_names(self,group):
        return self.channels[group]
    
    def get_rho_tan(self,group,names,shot,time):
        
        out = kk.KK().kkeqpfx(int(shot),time,diag='EQI')

        n =  (self.R_end[group]-self.R_start[group], self.z_end[group]-self.z_start[group])
        n_perp = c_[-n[1],n[0]]/hypot(n[0],n[1])[:,None]
        fdist = -sum(n_perp*(c_[self.R_start[group],self.z_start[group]]-r_[out.Raxis,out.zaxis]),1)
        
        t = linspace(0,1,50)
        X = array(n)[:,:,None]*t[None,None,:]+c_[self.R_start[group],self.z_start[group]].T[:,:,None]

        #import IPython
        #IPython.embed()
        ind = in1d(self.channels[group], names)
        rho = kk.KK().kkrzptfn( int(shot), time, X[0][ind].ravel(),X[1][ind].ravel(), diag='EQI')
        rho_p = rho.rho_p.reshape(X[0][ind].shape)
        rho_tg = amin(rho_p,1)*sign(fdist[ind])
        
        Rlim = array([interp([-1,1], r,x) for r,x in zip(rho_p,X[0])])
        zlim = array([interp([-1,1], r,x) for r,x in zip(rho_p,X[1])])
        self.los_length = hypot(Rlim[:,1]-Rlim[:,0], zlim[:,1]-zlim[:,0])

        
        return rho_tg
        
    
    def get_signal(self,group, names,calib=False,tmin=0.,tmax=10):

        self.dd.Open('BLB',self.shot, experiment=self.exp, edition=self.ed)

        tvec = self.dd.GetTimebase('pow'+group)
        sig  = self.dd.GetSignalGroupCalibrated('pow'+group)
        self.dd.Close()
               
        
        sig = sig[:, names]
        ind = slice(tvec.searchsorted(tmin), tvec.searchsorted(tmax))
        
        
        tvec, sig = tvec[ind], sig[ind]
        #BUG corrections of the foil bolometers
        if self.shot > 31591 and self.shot < 32300: #BUG find the accurate upper boundary!!
            if group == 'FVC':
                ind =  in1d(names, r_[16:names[-1]+1])
                sig[:,ind]*= 0.58
        
        if self.shot > 30135:
            if group == 'FVC':
                ind =  in1d(names, 16)
                sig[:,ind]*= 2.4
            if group == 'FHC':
                ind =  in1d(names, r_[36: names[-1]+1])
                sig[:,ind]*= 0.76
                ind =  in1d(names, 33)
                sig[:,ind]/= 1.1
                ind =  in1d(names, 34)
                sig[:,ind]/= 1.05
                ind =  in1d(names, r_[20:24])
                sig[:,ind]*= 1.2
                ind =  in1d(names, r_[0:12])
                sig[:,ind]*= 0.7
        elif self.shot > 28523:
            if group == 'FVC':
                ind =  in1d(names, r_[16: names[-1]+1])
                sig[:,ind]/= 1.75

            if group == 'FHC':
                ind =  in1d(names, r_[36: names[-1]+1])
                sig[:,ind]/= 1.1
                ind =  in1d(names, r_[33,36])
                sig[:,ind]/= 1.1
                

        elif self.shot > 27352:
            if group == 'FHC':
                ind =  in1d(names, 33)
                sig[:,ind]/= 1.2
                ind =  in1d(names, 34)
                sig[:,ind]/= 1.1

        #if self.shot == 32470:  #BUG!! zkontrolovat nové výboje!!
            #if group == 'FVC':
                #ind =  in1d(names, r_[16:names[-1]+1])
                #sig[:,ind]/= 0.58


        return  tvec, sig

    
    

def prepare_data(shot,tbeg,tend,plt_show=False):

    #cams = 'FHC','FVC','FLX','FDO','FDI','FHS'
    #cams = 'FHC', 'FVC','FHS','FLX','FDO','FDI'
    #if shot <  26503:
        #cams = 'FHC', 'FVC','FLX','FDO','FDI'
    #if shot >  31500:
    cams = 'FHC', 'FVC','FHS','FLX','FLH','FDO','FDI'

    
    fbolo = loader_foil_bolometers(shot)
    signals = []
    useful_cams = []
    indexes = {}
    ic = 0
    for c in cams:
        fnames = fbolo.get_names(c)
        tvec_,fsig = fbolo.get_signal(c, fnames,calib=True,tmin=tbeg,tmax=tend)
        if size(fsig) == 0: continue
        tvec = tvec_
        useful_cams.append(c)
        signals.append(fsig)
        indexes[c] = slice(ic,ic+fsig.shape[1])
        ic+= fsig.shape[1]
        
    sig = hstack(signals)
    weight = ones(sig.shape[1])
    noisy_cams =  'FHS','FLH','FDO','FDI','FLX'
    for nc in noisy_cams: 
        if nc in indexes:
            weight[indexes[nc]]/= 30
    
    ind = all(isfinite(sig),0)
    sig[:,~ind] = 0
    
    vmin,vmax  = mquantiles(sig,(.0001, .9999)) #remove outliers
    sig[(sig<vmin)|(sig>vmax)] = 0
    #sig*= weight[None,:]

    #sig  
    #sig[:,-16:] /= 100  #reduce weight of the extra noisy channels
    #sig[~isfinite(sig)] = 0
    
    
    #print sum(~isfinite(sig))
    #exit()
    u,s,v = linalg.svd(sig*weight[None,:], full_matrices=0)

    #plot(s);show()
    #plot(v[3:10].T);show()
    #BUG!!! this noidsy singular values must be specified by hand!!
    #4,5,6
    #noisy_sv = [5,6]
    #vmin,vmax  = mquantiles(sig,(.01, .99))
    #for i in range(10):
        #close()
        #subplot(121)
        #title(i)
        #sig_svd = outer(v[i],u[:,i]).T*s[i]
        #vmin_,vmax_  = mquantiles(sig_svd,(.01, .99))
        #imshow( sig_svd, interpolation='nearest', aspect='auto',vmin=vmin_,vmax=vmax_,origin='lower')
        #subplot(122)
        #imshow(sig, interpolation='nearest', aspect='auto',vmin=vmin,vmax=vmax,origin='lower');colorbar()
        #show()
    
    #s[noisy_sv] = 0
    n = 10
    fsig_ = dot(u[:,:n],s[:n,None]*v[:n,:])
    fsig_/= weight[None,:]
    #sig  /= weight[None,:]
    #for nc in noisy_cams: 
        #if nc in indexes:
            #fsig_[:,indexes[nc]]*= 10
    
    #sig[:,-16:] *= 100
    #fsig_[:,-16:] *= 100

    #import IPython
    #IPython.embed() 
    err = outer(std(fsig_-sig, 0),ones(len(tvec))).T+1e-6
    n = 7
    for nc in noisy_cams: 
        if nc in indexes:
            err[:,indexes[nc]]*= 30

    u,s,v = linalg.svd(sig/err, full_matrices=0)
    fsig_ = dot(u[:,:n],s[:n,None]*v[:n,:])*err
    
    if plt_show:
        subplot(131)
        vmin,vmax  = mquantiles(fsig_,(.01, .95))

        imshow(sig, interpolation='nearest', aspect='auto',vmin=vmin,vmax=vmax)#;colorbar()#;show()

        subplot(132)

        imshow(fsig_, interpolation='nearest', aspect='auto',vmin=vmin,vmax=vmax)#;show()
        subplot(133)
        vmin,vmax  = mquantiles( fsig_-sig,(.01, .99))

        imshow(fsig_-sig, interpolation='nearest', aspect='auto',vmin=min(vmin,-vmax),vmax=max(-vmin,vmax));colorbar()#;show()
        show()
    #err[:,-8:]/= 10

    #plot(std(fsig_-sig, 0))
    #plot(fsig_.mean(0))
    #show()
    
    for nc in noisy_cams: 
        if nc in indexes:
            err[:,indexes[nc]]*= 5
    

    tvec[0] = 0
    tvec[-1] = 9
    
    err[:,~ind] = infty

    savez_compressed('bolo_data/tvec_'+str(shot)+".npz",tvec=tvec)
    dets = 0

    for det, sig in zip(useful_cams,signals):
        ind = slice(dets,dets+sig.shape[1])
        savez_compressed('bolo_data/'+det+"_"+str(shot)+".npz", data = fsig_[:,ind].T, err = err[:,ind].T)
        dets+= sig.shape[1]

##savez_compressed('FHS'+"_"+str(shot)+".npz", data = fsig_[:,ind], err = err[:,ind])

##plot(s);show()
     
#imshow(fsig_,interpolation='nearest',aspect='auto',vmin=-fsig_.max()/10,origin='lower',vmax=fsig_.max()/10)
#show()                       
#import IPython
#IPython.embed()

#exit()
#imshow(v,interpolation='nearest',aspect='auto')
#show()

###shot,tbeg,tend = loadtxt('all_shots_.txt',unpack=True)
##for shot,tbeg,tend in zip(shot,tbeg,tend):
    
    ###if shot != 31605:
        ###continue
    
    ###if shot < 27900: continue
    ##shot = 31871
    ##tbeg = 2
    ##tend = 5.6
    ###if shot <= 26445:
        ###continue
    ###cams = ['VC', 'HC', 'DC']
    ##cams = ['HC', 'VC']

    ###close()
    
    ###f, axis = subplots(len(cams),2,figsize=(10,5*len(cams)))
    ###axis = atleast_2d(axis)
    ##try:
        ##fbolo = loader_foil_bolometers(shot)

    ##except Exception as e:
        ##print e
        ##continue
    ##for ic,c in enumerate(cams):
    ###for c in ['FHC','FVC']:


 
        ##try:
            ##fnames = fbolo.get_names('F'+c)
            ##print fnames
            ##_,fsig = fbolo.get_signal('F'+c, fnames,calib=True,tmin=tbeg,tmax=tend)
     
        ##except:
            ##pass

        ##fsig = mean(fsig,0)
        ##plot(fsig)
        ##show()






#shot,tbeg,tend = loadtxt('w_accumulations.txt',unpack=True)

#shot,tbeg,tend = loadtxt('all_shots_.txt',unpack=True)

#for shot,tbeg,tend in zip(shot,tbeg,tend):
    #if shot <= 26445:
        #continue
    ##close()
    ##if shot < 27900: continue
    ##shot = 30382
    ##tbeg = 5
    ##tend = 5.2
    ##plot()
    #try:
        #fbolo = loader_foil_bolometers(shot)
    #except Exception as e:
        #continue

    #for ic,c in enumerate(['FHC','FVC']):
        #try:
            #fnames = fbolo.get_names(c)
            #_,fsig = fbolo.get_signal(c, fnames,calib=True,tmin=tbeg,tmax=tend)
            #fsig = mean(fsig,0)
        #except Exception as e:
            #print e
            #break
        
        ##plot(fsig)
        ##show()

        ##TODO normalizovat délkou chordy? tak aby to por plocjé profily více méně sedělo. 
        #frho_tg = fbolo.get_rho_tan(c,fnames,shot,(tend+tbeg)/2)
        #los_length = fbolo.los_length 
        ##los_length= 1
        #fsig/= los_length
        #plot(frho_tg, fsig/1e6,'k',ls=['-','--'][ic], picker=2)
        #plot(frho_tg, fsig/1e6,'k',ls=['-','--'][ic], picker=2)

        #for fcam in fbolo.subcam_ind[c]:
            #plot(frho_tg[fcam], fsig[fcam]/1e6,'o', picker=2)
            
        #if c == 'FVC':
            #plot(-frho_tg, fsig/1e6,'k-',ls=['-','--'][ic], picker=2)

            #for fcam in fbolo.subcam_ind[c]:
                #plot(-frho_tg[fcam], fsig[fcam]/1e6,'s', picker=2)
                
    #xlim(-1,1.2)
    #ylim(0,median(fsig[abs(frho_tg)<0.4])*2/1e6)
    
    ##show()

    ##savez('./bolo_compare/compare_%d_%.2f'%(c, shot,tbeg), drho_tg=drho_tg[dind],frho_tg=frho_tg[find],dsig=dsig[dind],fsig=fsig[find], fd_ratio =  fd_ratio,dind= dind,find=find )
    #savefig('./bolo_compare/compare%d_%.2f.png'%(shot,tbeg),  bbox_inches='tight')
    #clf()
    #show()

##29639, 4,6


#exit()
#import kk
def show_bolo_data(shot,tbeg,tend):

    cams = ['VC', 'HC', 'DC']
    cams = ['VC', 'HC',]

    close()
    
    f, axis = subplots(len(cams),2,figsize=(10,5*len(cams)))
    axis = atleast_2d(axis)
    try:
        fbolo = loader_foil_bolometers(shot)
        dbolo = loader_diod_bolometers(shot)

    except Exception as e:
        print e
        return
    
    
    for ic,c in enumerate(cams):



        loaded = False
        try:
            fnames = fbolo.get_names('F'+c)
            dnames = dbolo.get_names('D'+c)
            _,dsig = dbolo.get_signal('D'+c, dnames,calib=True,tmin=tbeg,tmax=tend)
            _,fsig = fbolo.get_signal('F'+c, fnames,calib=True,tmin=tbeg,tmax=tend)
            if c == 'VC':
                
                dnames2 = dbolo.get_names('D13')
                _,dsig2 = dbolo.get_signal('D13', dnames2,calib=True,tmin=tbeg,tmax=tend)

            
        except Exception as e:
            print e
            return
       
            
        print 'loaded'
        loaded = True
    
    
        dsig = array([mean(d) for d in  dsig])
        fsig = mean(fsig,0)
        
        
        find = fsig > 1e4
        dind = dsig > 1e4
        if sum(find)<4 or sum(dind)<4:
            return 
        
        
        frho_tg = fbolo.get_rho_tan('F'+c,fnames,shot,(tend+tbeg)/2)
        drho_tg = dbolo.get_rho_tan('D'+c,dnames,shot,(tend+tbeg)/2)


    
        if c == 'VC':
            dsig2 = array([mean(d) for d in  dsig2])
            drho_tg2 = dbolo.get_rho_tan('D13',dnames2,shot,(tend+tbeg)/2)
            dind2 = dsig2 > 1e4


        sind = argsort(frho_tg[find])
        fd_ratio = dsig[dind]/interp(drho_tg[dind],frho_tg[find][sind],fsig[find][sind])
        

        if c == 'VC':
            fd_ratio2 = dsig2[dind2]/interp(drho_tg2[dind2],frho_tg[find][sind],fsig[find][sind])
            for dcam in dbolo.subcam_ind['D13']:
                axis[ic,0].plot(drho_tg2[dind2&dcam], fd_ratio2[dcam[dind2]] ,
                                '-*',markeredgewidth=0.0,markersize=10, picker=2) 

        for dcam in dbolo.subcam_ind['D'+c]:
            axis[ic,0].plot(drho_tg[dind&dcam], fd_ratio[dcam[dind]] ,'-o', picker=2)

        axis[ic,0].set_ylim(.05,2)
        axis[ic,0].set_ylabel('diod/foil bolometers, F%s and D%s'%(c,c))
        axis[ic,0].axhline(1)
        axis[ic,0].set_xlabel('rho_tan ')
        axis[ic,0].set_xlim(-1.07,1.07)
        axis[ic,0].grid('on')
        axis[ic,0].set_yscale("log", nonposx='clip')


        f.suptitle('Comparison of diod and foil bolometers %d, t: %.1f-%.1fs'%(shot,tbeg,tend))

        
        axis[ic,1].plot(drho_tg[dind], dsig[dind]/1e6,'-k',label='D'+c)
        axis[ic,1].plot(frho_tg[find], fsig[find]/1e6,'--k',label='F'+c)
        if c == 'VC':axis[ic,1].plot(drho_tg2[dind2], dsig2[dind2]/1e6,':k',label='D13')
        axis[ic,1].legend(loc='best')

        for dcam in dbolo.subcam_ind['D'+c]:
            axis[ic,1].plot(drho_tg[dind&dcam], dsig[dind&dcam]/1e6,'s', picker=2)
            

        for fcam in fbolo.subcam_ind['F'+c]:
            axis[ic,1].plot(frho_tg[find&fcam], fsig[find&fcam]/1e6,'o', picker=2)

        if c == 'VC':
            for dcam in dbolo.subcam_ind['D13']:
                axis[ic,1].plot(drho_tg2[dind2&dcam], dsig2[dind2&dcam]/1e6,'*',
                                markeredgewidth=0.0,markersize=10, picker=2) 

        axis[ic,1].set_ylim(2e-2,nanmean(fsig[find][abs(frho_tg[find])<0.9]/1e6)*5)
        axis[ic,1].set_xlim(-1.07,1.07)
        if c == 'HS':
            axis[ic,0].set_xlim(0.2,1.07) 
            axis[ic,1].set_xlim(0.2,1.07)
            axis[ic,1].set_ylim(1e-2, amax(fsig[find])/1e6*1.2)

        axis[ic,1].set_yscale("log", nonposx='clip')
        axis[ic,1].set_xlabel('rho_tan ')
        axis[ic,1].set_ylabel('Power [MW/m$^2$]')
        axis[ic,1].grid('on')
        
    
        
        savez('./bolo_compare/%s_%d_%.2f'%(c, shot,tbeg), drho_tg=drho_tg[dind],
              frho_tg=frho_tg[find],dsig=dsig[dind],fsig=fsig[find], fd_ratio =  fd_ratio,dind= dind,find=find )
        
    if not loaded:
        return
        
    f.savefig('./bolo_compare/HFS_%d_%.2f.png'%(shot,tbeg),  bbox_inches='tight')

    def onpick1(event):
        thisline = event.artist
        xdata = thisline.get_xdata()
        ydata = thisline.get_ydata()

        for c in cams:

            data = load('./bolo_compare/%s_%d_%.2f.npz'%(c, shot,tbeg))
            frho_tg = data['frho_tg']
            drho_tg = data['drho_tg']
            if any(in1d(xdata,drho_tg2)):
                dnames2 = dbolo.get_names('D13') 
                ind = where(xdata[event.ind[0]] == drho_tg2 )[0][0]
                print c+' diod2',dnames2[ind],ind
            if any(in1d(xdata,drho_tg )):
                dnames = dbolo.get_names('D'+c)
                ind = where(xdata[event.ind[0]] == drho_tg )[0][0]
                print c+' diod1', dnames[ind],ind
            if any(in1d(xdata, frho_tg)):
                fnames = fbolo.get_names('F'+c)
                ind = where(xdata[event.ind[0]] == frho_tg )[0][0]
                print c+' foil', fnames[ind],ind

    f.canvas.mpl_connect('pick_event', onpick1)

    #show()
    
    


    
#show_bolo_data(32284,1,7)
#show()

    
 
def  main():
    prepare_data(33118,1,5)
    #exit()

    #shot,tbeg,tend = loadtxt('all_shots_.txt',unpack=True)
    show_bolo_data(33118,2,3)
    show()
    for shot in range(32470,32000,-1):
        #if shot <= 26445:
            #continue
        show_bolo_data(shot,1,5)
        
    exit()
    
    for shot,tbeg,tend in zip(shot,tbeg,tend):
        if shot <= 26445:
            continue
        show_bolo_data(shot,tbeg,tend)
    
    show_bolo_data(32470,1,5)
    show()

    #prepare_data(32470,1,5)

    #exit()
    show_bolo_data(30812,4.8,5)
    show()
    
    shot,tbeg,tend = loadtxt('w_accumulations.txt',unpack=True)

    shot,tbeg,tend = loadtxt('all_shots_.txt',unpack=True)
    for shot,tbeg,tend in zip(shot,tbeg,tend):
        if shot <= 26445:
            continue
        show_bolo_data(shot,tbeg,tend)


    shot,tbeg,tend = loadtxt('w_accumulations.txt',unpack=True)
    for shot,tbeg,tend in zip(shot,tbeg,tend):
        if shot <= 26445:
            continue
        show_bolo_data(shot,tbeg,tend)







if __name__ == "__main__":
    main()








