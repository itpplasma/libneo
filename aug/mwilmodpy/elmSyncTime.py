
'''
determines the timebase towards the next ELM and from the past ELM
'''

__author__='Florian Laggner (flag@ipp.mpg.de)'
__date__='Dec 2014'
__version__='2.0'

import numpy as np
import random as rand

import dd

def elmSyncTime(time, nshot,**kwargs):
    sf = dd.shotfile('ELM', nshot, **kwargs)
    t_begELM = sf('f_ELM').time
    t_endELM = sf('t_endELM').data
    sf.close()
    #limit to requested time interval
    ind = np.where((t_begELM >= np.min(time)) & (t_endELM <= np.max(time)))
    t_begELM = t_begELM[ind]
    t_endELM = t_endELM[ind]
    #get begin time of following ELM (corresponding to t_begELM of t_endELM)
    t_begELM_next = np.roll(t_begELM,-1)
    t_endELM_next = np.roll(t_endELM,-1)
    #define output
    Next = np.array([np.nan]*len(time))
    Last = np.array([np.nan]*len(time))
    NextEnd = np.array([np.nan]*len(time))
    LastEnd = np.array([np.nan]*len(time))
    #time before first ELM
    ind = np.where(time <= t_begELM[0])
    Next[ind] = t_begELM[0] - time[ind]
    ind = np.where(time <= t_endELM[0])
    NextEnd[ind] = (t_endELM[0]) - time[ind]
    #iterate over ELMs
    for elm in enumerate(t_begELM):
        ind = np.where((time > elm[1]) & (time <= t_begELM_next[elm[0]]))
        Next[ind] = t_begELM_next[elm[0]] - time[ind]
        Last[ind] = time[ind] - elm[1]
        ind = np.where((time > t_endELM[elm[0]]) & (time <= t_endELM_next[elm[0]]))
        NextEnd[ind] = t_endELM_next[elm[0]] - time[ind]
        LastEnd[ind] = time[ind] - t_endELM[elm[0]]
    #time after last ELM
    ind = np.where(time > t_begELM[-1])
    Last[ind] = time[ind] - t_begELM[-1]
    ind = np.where(time > t_endELM[-1])
    LastEnd[ind] = time[ind] - t_endELM[-1]
    
    #define output dictionary
    t={'next':Next,'last':Last,'nextEnd':NextEnd,'lastEnd':LastEnd}
    return t

def elmSyncTimeWOelm(time, nshot,**kwargs):
    sf = dd.shotfile('ELM', nshot, **kwargs)
    t_begELM = sf('f_ELM').time
    t_endELM = sf('t_endELM').data
    sf.close()
    #limit to requested time interval
    ind = np.where((t_begELM >= np.min(time)) & (t_endELM <= np.max(time)))
    t_begELM = t_begELM[ind]
    t_endELM = t_endELM[ind]
    #get begin time of following ELM (corresponding to t_begELM of t_endELM)
    t_begELM_next = np.roll(t_begELM,-1)
    t_endELM_next = np.roll(t_endELM,-1)
    #define output
    Next = np.array([np.nan]*len(time))
    Last = np.array([np.nan]*len(time))
    NextEnd = np.array([np.nan]*len(time))
    LastEnd = np.array([np.nan]*len(time))
    #time before first ELM
    ind = np.where(time <= t_begELM[0])
    Next[ind] = t_begELM[0] - time[ind]
    NextEnd[ind] = t_endELM[0] - time[ind]
    #iterate over ELMs
    for elm in enumerate(t_endELM):
        ind = np.where((time > elm[1]) & (time <= t_begELM_next[elm[0]]))
        Next[ind] = t_begELM_next[elm[0]] - time[ind]
        Last[ind] = time[ind] - t_begELM[elm[0]]
        NextEnd[ind] = t_endELM_next[elm[0]] - time[ind]
        LastEnd[ind] = time[ind] - t_endELM[elm[0]]
    #time after last ELM
    ind = np.where(time > t_endELM[-1])
    Last[ind] = time[ind] - t_begELM[-1]
    LastEnd[ind] = time[ind] - t_endELM[-1]
    
    #define output dictionary
    t={'next':Next,'last':Last,'nextEnd':NextEnd,'lastEnd':LastEnd}
    return t

def elmFastSlow(time, nshot,dtfastslow,begEnd=False,**kwargs):    
    #extract fast and slow ELMs (dtfastslow in [s])
    elm = dd.shotfile('ELM', nshot, **kwargs)
    t_begELM = elm('f_ELM',tBegin=np.min(time),tEnd=np.max(time)).time
    t_endELM = elm('t_endELM',tBegin=np.min(time),tEnd=np.max(time)).data
    elm.close()
    #check if determination of fast/slow ELM from begin/end of ELM
    if begEnd == False:
        indFastELM = np.squeeze(np.where(np.subtract(np.delete(np.roll(t_begELM,-1),-1),np.delete(t_begELM,-1)) <= dtfastslow))+1
        indSlowELM = np.squeeze(np.where(np.subtract(np.delete(np.roll(t_begELM,-1),-1),np.delete(t_begELM,-1)) > dtfastslow))+1
    else:
        indFastELM = np.squeeze(np.where(np.subtract(np.delete(np.roll(t_endELM,-1),-1),np.delete(t_endELM,-1)) <= dtfastslow))+1
        indSlowELM = np.squeeze(np.where(np.subtract(np.delete(np.roll(t_endELM,-1),-1),np.delete(t_endELM,-1)) > dtfastslow))+1
    #define output dictionary
    t={'fast':np.nan,'slow':np.nan,'indFast':np.nan,'indSlow':np.nan}
    #loop over all elements of time
    for timeStep in enumerate(time):
        #get indexof next ELM
        if begEnd == False:
            indNextELM = np.nanmin(np.where(t_begELM >= timeStep[1],np.arange(0,len(t_begELM)),np.nan))
        else:
            indNextELM = np.nanmin(np.where(t_endELM >= timeStep[1],np.arange(0,len(t_begELM)),np.nan))
        #check if next ELM is a fast or slow ELM
        if indNextELM in indFastELM:
            t['fast'] = np.append(t['fast'],timeStep[1])
            t['indFast'] = np.append(t['indFast'],timeStep[0])
        elif indNextELM in indSlowELM:
            t['slow'] = np.append(t['slow'],timeStep[1])
            t['indSlow'] = np.append(t['indSlow'],timeStep[0])
        #else:
        #    print 'error'

    #delete first element from output dictionary -> only for initialisation
    t['fast'] = np.delete(t['fast'],0)
    t['slow'] = np.delete(t['slow'],0)
    t['indFast'] = np.array(map(int,np.squeeze(np.delete(t['indFast'],0))))
    t['indSlow'] = np.array(map(int,np.squeeze(np.delete(t['indSlow'],0))))
    return t

def elmSyncCycle(time,nshot,begEnd=False,dtfastslow=False,remELM=False,**kwargs):
    if remELM != False:
        tbase = elmSyncTimeWOelm(time, nshot,**kwargs)
    else:
        tbase = elmSyncTime(time, nshot,**kwargs)
    if begEnd == False:               
        tsortIndNext = np.argsort(tbase['next'])
        tsortIndLast = np.argsort(tbase['last'])
        #sort arrays
        tbase['next'] = tbase['next'][tsortIndNext]
        tbase['last'] = tbase['last'][tsortIndLast]

        #   time = time[tsortIndNext]
        #   timeLast = time[tsortIndLast]

        #reverse times before Sync time
        tbase['next'] = -tbase['next'][::-1]
        tsortIndNext = tsortIndNext[::-1]
        #combine arrays for output
        tELM = np.append(tbase['next'],tbase['last'])
        tindELM = np.append(tsortIndNext,tsortIndLast)
        #define output dictionary
        t={'tELM':tELM,'tindELM':tindELM}
    else:       
        tsortIndNext = np.argsort(tbase['nextEnd'])
        tsortIndLast = np.argsort(tbase['lastEnd'])
        #sort arrays
        tbase['nextEnd'] = tbase['nextEnd'][tsortIndNext]
        tbase['lastEnd'] = tbase['lastEnd'][tsortIndLast]
        
        #   time = time[tsortIndNext]
        #   timeLast = time[tsortIndLast]
        
        #reverse times before Sync time
        tbase['nextEnd'] = -tbase['nextEnd'][::-1]
        tsortIndNext = tsortIndNext[::-1]
        #combine arrays for output
        tELM = np.append(tbase['nextEnd'],tbase['lastEnd'])
        tindELM = np.append(tsortIndNext,tsortIndLast)
        #define output dictionary
        t={'tELM':tELM,'tindELM':tindELM}

    if dtfastslow != False:
        if remELM != False:
            base = elmSyncTimeWOelm(time, nshot,**kwargs)
        else:
            base = elmSyncTime(time, nshot,**kwargs)
        FS = elmFastSlow(time, nshot,dtfastslow,begEnd=False,**kwargs)
        if begEnd == False: 
            tsortIndNext_Fast = np.argsort(base['next'][FS['indFast']])
            tsortIndLast_Fast = np.argsort(base['last'][FS['indFast']])
            tsortIndNext_Slow = np.argsort(base['next'][FS['indSlow']])
            tsortIndLast_Slow = np.argsort(base['last'][FS['indSlow']])
            #sort arrays
            nextFast = np.array(base['next'][FS['indFast']])[tsortIndNext_Fast]
            lastFast = np.array(base['last'][FS['indFast']])[tsortIndLast_Fast]
            nextSlow = np.array(base['next'][FS['indSlow']])[tsortIndNext_Slow]
            lastSlow = np.array(base['last'][FS['indSlow']])[tsortIndLast_Slow]
            #reverse times before Sync time
            nextFast = -nextFast[::-1]
            nextSlow = -nextSlow[::-1]
            tsortIndNext_Fast = tsortIndNext_Fast[::-1]
            tsortIndNext_Slow = tsortIndNext_Slow[::-1]
            #combine arrays for output
            t['tELM_fast'] = np.append(nextFast,lastFast)
            t['tELM_slow'] = np.append(nextSlow,lastSlow)
            t['tindELM_fast'] = np.append(FS['indFast'][tsortIndNext_Fast],FS['indFast'][tsortIndLast_Fast])
            t['tindELM_slow'] = np.append(FS['indSlow'][tsortIndNext_Slow],FS['indSlow'][tsortIndLast_Slow])
        else:
            tsortIndNext_Fast = np.argsort(base['nextEnd'][FS['indFast']])
            tsortIndLast_Fast = np.argsort(base['lastEnd'][FS['indFast']])
            tsortIndNext_Slow = np.argsort(base['nextEnd'][FS['indSlow']])
            tsortIndLast_Slow = np.argsort(base['lastEnd'][FS['indSlow']])
            #sort arrays
            nextFast = np.array(base['nextEnd'][FS['indFast']])[tsortIndNext_Fast]
            lastFast = np.array(base['lastEnd'][FS['indFast']])[tsortIndLast_Fast]
            nextSlow = np.array(base['nextEnd'][FS['indSlow']])[tsortIndNext_Slow]
            lastSlow = np.array(base['lastEnd'][FS['indSlow']])[tsortIndLast_Slow]
            #reverse times before Sync time
            nextFast = -nextFast[::-1]
            nextSlow = -nextSlow[::-1]
            tsortIndNext_Fast = tsortIndNext_Fast[::-1]
            tsortIndNext_Slow = tsortIndNext_Slow[::-1]
            #combine arrays for output
            t['tELM_fast'] = np.append(nextFast,lastFast)
            t['tELM_slow'] = np.append(nextSlow,lastSlow)
            t['tindELM_fast'] = np.append(FS['indFast'][tsortIndNext_Fast],FS['indFast'][tsortIndLast_Fast])
            t['tindELM_slow'] = np.append(FS['indSlow'][tsortIndNext_Slow],FS['indSlow'][tsortIndLast_Slow])
        
    return t
