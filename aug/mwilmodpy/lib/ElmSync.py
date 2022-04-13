"""

    Routine for elm syncronization. It uses the augped approach.
    For reference see P. Schneider, Dissertation, pp. 54

    Version 2.0: added configuration via dictionary

"""

__author__='Marco Cavedon (mcavedon@ipp.mpg.de)'
__date__='24 Jul 2014'
__version__='2.0'


import matplotlib.pylab as plt
#import dd_20180130 as dd
import dd
import numpy as np
import IPython
defaultSetting = {'dt_elm':2.2e-3,'t_start':-1.e-3, 't_end':-2.25e-3,'max_dt_elm_cycle':20.e-4}


def ElmSync(time, nshot, syncSetting=None,
            plot=False, Debug=False, Experiment='AUGD'):
    """
    ELM sync

    IF t_signal > t_begin_elm(nearest) THEN
        IF (t_signal > t_begin_elm(nearest) + dt_elm) AND
           (t_begin_elm(nearest+1) + t_end < t_signal < t_begin_elm(nearest+1) + t_end) THEN
            t_signal OK!

    IF t_signal < t_begin_elm(nearest)
        IF (t_signal > t_begin_elm(nearest-1) + dt_elm) AND
           (t_begin_elm(nearest) + t_end < t_signal < t_begin_elm(nearest) + t_end) THEN
            t_signal OK!

    IF t_signal - t_begin_elm(nearest) > max_dt_elm_cycle THEN
        t_signal OK!
    """
    if syncSetting == None:
        syncSetting = defaultSetting
    dt_elm = syncSetting['dt_elm'] #duration of an ELM = 4 ms
    #Relative time of the considered region
    t_start = syncSetting['t_start']
    t_end = syncSetting['t_end'] # = -6.5e-3 for "slow" diags
    max_dt_elm_cycle = syncSetting['max_dt_elm_cycle']
    sf = dd.shotfile('ELM',nshot,experiment=Experiment)
    t_beg_ELM =sf.getTimeBase('f_ELM') 
    jt = 0
    good_jtime = []
    for t_signal in time:
        jelm_nearest = np.argmin(np.abs(t_signal-t_beg_ELM))
        delta_t_nearest = t_signal-t_beg_ELM[jelm_nearest]
        if (jelm_nearest == 0. or jelm_nearest == 1) and \
                np.abs(delta_t_nearest) > max_dt_elm_cycle:
            if Debug:
                print ('Warning! No ELM data for t = %.3f s'%t_signal)
            good_jtime.append(jt)
            jt+=1
            continue
        if delta_t_nearest > 0.:
            if jelm_nearest == len(t_beg_ELM)-1 and \
                t_signal > (t_beg_ELM[jelm_nearest]+dt_elm):
                good_jtime.append(jt)
                jt+=1
                if Debug:
                    print ('t =%.3f  s accepted!'%t_signal)
                continue
            if t_signal>(t_beg_ELM[jelm_nearest]+dt_elm) and\
                t_signal>(t_beg_ELM[jelm_nearest+1]+t_end) and\
                t_signal<(t_beg_ELM[jelm_nearest+1]+t_start):
                good_jtime.append(jt)
                jt+=1
                if Debug:
                    print ('t = %.3f s accepted! '%t_signal)
                continue
        if delta_t_nearest < 0. and \
            t_signal>(t_beg_ELM[jelm_nearest-1]+dt_elm) and\
            t_signal>(t_beg_ELM[jelm_nearest]+t_end) and\
            t_signal<(t_beg_ELM[jelm_nearest]+t_start):
            good_jtime.append(jt)
            jt+=1
            if Debug:
                print ('t ==%.3f s accepted!'%t_signal)
            continue
        if Debug:
            print ('t =  s discard!'%t_signal)
        jt+=1
    good_jtime = np.array(good_jtime)
    if plot:
        sf = dd.shotfile('MAC',nshot)
        ipolsola =  sf.getSignal('Ipolsola')
        time_ipolsola = sf.getTimebase('Ipolsola')
        f = plt.figure()
        ax = f.add_subplot(111)
        ax.plot(time_ipolsola,ipolsola,'k',label='Ipolsola')
        ax.plot(t_beg_ELM,np.zeros(len(t_beg_ELM)),'rD',label='ELM signal')
        ax.plot(time[good_jtime],np.zeros(len(time[good_jtime])),'bo',label='Selected time points')
        dt = np.diff(time)
        dt = np.append(dt,dt[-1])
        t_start = time - dt/2.
        t_end = time + dt/2.
        ylim = ax.get_ylim()
        up = ylim[1]+np.zeros_like(time)
        down = ylim[0]+np.zeros_like(time)
        ax.set_xlim([np.min(t_beg_ELM),np.max(t_beg_ELM)])
        ax.vlines(t_start[good_jtime],down[good_jtime],up[good_jtime],linestyles='dashed',color='b')
        ax.vlines(t_end[good_jtime],down[good_jtime],up[good_jtime],linestyles='dashed',color='b')
        plt.legend()
        plt.show()
    return good_jtime


def ElmExtract(time, nshot, syncSetting=None, plot=False, preFac = 0.1, postFac = 0.2, Experiment='AUGD' ):

    #IPython.embed()
    sf = dd.shotfile('ELM',nshot,experiment=Experiment)
    t_beg_ELM =sf.getTimeBase('f_ELM')#  sf('f_ELM').time
    f_ELM = sf.getSignal('f_ELM') 
    t_end_ELM =  sf.getSignal('t_endELM')  
    sf.close()

    if syncSetting == None:
        syncSetting = defaultSetting

    #Relative time of the considered region
    t_start = syncSetting['t_start']
    t_end = syncSetting['t_end']

    good_jtime = np.arange(np.size( time ))
    deleteIdx=[]


    for i in np.arange(np.size(t_beg_ELM)):
        try:
            if (i > 1) & (i < np.size(t_beg_ELM)-2):
                t_beg_range = np.abs(t_beg_ELM[i]-t_end_ELM[i-1])*preFac
                t_end_range = np.abs(t_beg_ELM[i+1]-t_end_ELM[i])*postFac
            else:
                t_beg_range =  1./f_ELM[i]*preFac
                t_end_range =  1./f_ELM[i]*postFac
        except:
            IPython.embed()

        ind = np.squeeze(np.where( (time >= (t_beg_ELM[i] - t_beg_range) ) & (time <= (t_end_ELM[i] + t_end_range )) )[0]) 
        
        if np.size(ind) > 0:
            if np.size(deleteIdx) > 0:
                deleteIdx = np.append( deleteIdx, ind )
            else:
                deleteIdx = np.array( ind )

    #IPython.embed()
    good_jtime = np.delete(good_jtime,deleteIdx)

    if plot:
        sf = dd.shotfile('MAC',nshot)
       
        time_ipolsola = sf.getTimeBase('Ipolsola')
        ipolsola = sf.getSignal('Ipolsola')
        sf.close()
        f = plt.figure()
        ax = f.add_subplot(111)
        ax.plot(time_ipolsola,ipolsola,'k',label='Ipolsola')
        ax.plot(t_beg_ELM,np.zeros(len(t_beg_ELM)),'rD',label='ELM signal')
        ax.plot(t_end_ELM,np.zeros(len(t_end_ELM)),'gD',label='ELM signal end')
        ax.plot(time[good_jtime],np.zeros(len(time[good_jtime])),'bo',label='Selected time points')
        ylim = ax.get_ylim()
        up = ylim[1]+np.zeros_like(time)
        down = ylim[0]+np.zeros_like(time)
        ax.set_xlim([np.min(t_beg_ELM),np.max(t_beg_ELM)])
        ax.vlines(time[good_jtime],down[good_jtime],up[good_jtime],linestyles='dashed',color='b')
        plt.legend()
        plt.show()
    return good_jtime



if __name__=="__main__":
    sf = dd.shotfile('CMR',31303)
    sig = sf('CCD-Sig')
    ElmSync(sig.time,31303,plot=True)
