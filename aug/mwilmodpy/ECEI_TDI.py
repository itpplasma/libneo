##import ECEI
##ECI=ECEI.ECEI()
##ECI.Load(30839)



import numpy as np
import dd as dd
#import dd_mwillens as dd2
import kk as kk
import ECE
import matplotlib.pylab as plt
from scipy.interpolate import interp1d
import IPython
#import eqi_map as fastkk
import pickle
import tools
import PHI


path = '/afs/ipp-garching.mpg.de/home/m/mwillens/Programme/python/data/ECEI/'
nCols = 8
nRows_old = 16
nRows_new = 20
###
###
###
class ECEIhelp:
    status = False


class ECEI:
    def __init__( self ,  Experiment = 'AUGD', Diagnostic = 'ECI', Shotnumber = None ):
        self.Status = False
        self.TradStatus = False
        self.StatusRz = False
        self.StatusECFM = False
        self.StatusTrad = False
        self.Status_allrhop = False
        self.StatusRzWarm = False
        if Shotnumber != None :
            self.Load( Shotnumber )
		

    def __del__( self ):
        self.Unload( )
        del self.Status

    def checkTimebase( self, time, timeIn, dataIn ):
        if dataIn[:,0].size != time.size:
            print "time base not the same, fill in take longer" 
            dtime=np.diff(time).mean()
     #                   from scipy.interpolate import griddata
            
            runIdx = np.arange(timeIn.size)
            corrTime=timeIn
            ## get the missing index
            missingIdx = np.remainder(np.interp(time,timeIn,runIdx),1.)>0.0
            availIdx = missingIdx==False
            dataOut=np.zeros((time.size,dataIn[0,:].size))
            dataOut[availIdx] = dataIn
            f=interp1d(timeIn, dataIn, kind='linear', axis=0, copy=True, bounds_error=False, fill_value=0.0)
            dataOut[missingIdx] = f(time[missingIdx])
 
            return time,dataOut            
        else:
            return timeIn,dataIn


    def Load( self ,  Shotnumber, Experiment='AUGD', Diagnostic='TDI', Edition = 0L, tBegin=-0.1, tEnd=11.0, loadAllRhop=False, rzExp = 'ECEI', rzDiag = 'RZN', rzEd = 0, raw=True , eqExp = 'AUGD', eqDiag = 'EQH', fileEndECFM='_ECFM_all.dat', binning=10., system='both',useEnd=False, skip=[0.0,0.0],specialOffset=False,specialdt=0.5):
        self.Unload()
        
        if (system != 'new') & (system != 'old') & (system != 'both'):
            print "System not appropriate defined, take 'both'"
            system = 'both'

        if Diagnostic == 'TDI':
            try:
                sf = dd.shotfile( Diagnostic, Shotnumber, Experiment, Edition)
                self.Shotnumber = Shotnumber
            except:
                print "Error reading shotfile" 
                return False

            try:
                #			       print "get timebase"
                self.tBegin = tBegin
                self.tEnd = tEnd
                print "Reading TDI"          
               
                if ( (np.all(skip) != 0) & (np.size(skip) == 2) & (skip[0] < skip[1]) ):
                    doSkip = True
                else:
                    doSkip = False

                self.nCols = nCols
                
                print "Reading timebase" 
                time = sf.getTimeBase( 'Sig1' )
                ntime = np.size(time)
                tidx = np.where((time>=tBegin)&(time<=tEnd))[0]
                
                discardNew=np.array([0,1,11,19])
                self.useRowsNew = np.delete(np.arange(nRows_new),discardNew)
                discardOld=[12]
                self.useRowsOld = np.delete(np.arange(nRows_old),discardOld)   
                #IPython.embed()
                self.nCols = nCols
                

                if (system=='old') |  (system=='both'):
                    print "Reading data 1" 
                    data1 = sf.getSignalGroup( 'Sig1' ).T
                    offEnd1 = np.mean(data1[time>10.],axis=0) 
                    offBeg1 = np.mean(data1[time<0.01],axis=0) 

                    newtime1, data1 = tools.dataBinning( time[tidx], data1[tidx].T, samplefreq = binning )
                    if doSkip == False:
                        if useEnd:
                            data1 = (data1.T - offEnd1).T
                        else:
                            data1 = (data1.T - offBeg1).T
                    else:
                        arr = (newtime1<=skip[0]) | (newtime1>=skip[1])
                        newtime1, data1 =newtime1[arr], data1[:,arr]
                        data1[:,newtime1<=skip[0]] = ((data1.T)[newtime1<=skip[0]] - offBeg1).T
                        if specialOffset:
                            addOffset = np.median(data1[:,(newtime1<=skip[0]) & (newtime1>=(skip[0]-specialdt))],axis=1) - np.median(data1[:,(newtime1>=skip[1]) & (newtime1<=(skip[1]+specialdt))],axis=1)
                            data1[:,newtime1>=skip[1]] = (data1[:,newtime1>=skip[1]].T + addOffset).T
                        else:
                            data1[:,newtime1>=skip[1]] = ((data1.T)[newtime1>=skip[1]] - offEnd1).T


                print "Reading data 2" 
                data2 = sf.getSignalGroup( 'Sig2' ).T
                time2 = sf.getTimeBase( 'Sig2' )
                tidx2 = np.where((time2>=tBegin)&(time2<=tEnd))[0]

                offEnd2 = np.mean(data2[time2>10.],axis=0) 
                offBeg2 = np.mean(data2[time2<0.01],axis=0) 

                if (tidx2.size!=tidx.size):
                    time2,data2 = self.checkTimebase( time[tidx], time2[tidx2], data2[tidx2] )
                    newtime2, data2 = tools.dataBinning( time2, data2.T, samplefreq = binning )
                else:
                    newtime2, data2 = tools.dataBinning( time[tidx], data2[tidx].T, samplefreq = binning )


                if doSkip == False:
                    if useEnd:
                        data2 = (data2.T - offEnd2).T
                    else:
                        data2 = (data2.T - offBeg2).T
                else: 
                    arr = (newtime2<=skip[0]) | (newtime2>=skip[1])
                    newtime2, data2 =newtime2[arr], data2[:,arr]
                    data2[:,newtime2<=skip[0]] = ((data2.T)[newtime2<=skip[0]] - offBeg2).T
                    if specialOffset:
                        addOffset = np.median(data2[:,(newtime2<=skip[0]) & (newtime2>=(skip[0]-specialdt))],axis=1) - np.median(data2[:,(newtime2>=skip[1]) & (newtime2<=(skip[1]+specialdt))],axis=1)
                        data2[:,newtime2>=skip[1]] = (data2[:,newtime2>=skip[1]].T + addOffset).T
                    else:
                        data2[:,newtime2>=skip[1]] = ((data2.T)[newtime2>=skip[1]] - offEnd2).T
                    

                
                if (system=='new') |  (system=='both'):
                    print "Reading data 3" 
                    data3 = sf.getSignalGroup( 'Sig3' ).T
                    time3 = sf.getTimeBase( 'Sig3' )
                    tidx3 = np.where((time3>=tBegin)&(time3<=tEnd))[0]
               
                    offEnd3 = np.mean(data3[time3>10.],axis=0) 
                    offBeg3 = np.mean(data3[time3<0.01],axis=0)              

   
                    if np.all(tidx3!=tidx):
                        time3,data3 = self.checkTimebase( time[tidx], time3[tidx3], data3[tidx3])
                        newtime3, data3 = tools.dataBinning( time3, data3.T, samplefreq = binning )
                    else:
                        newtime3, data3 = tools.dataBinning( time[tidx], data3[tidx].T, samplefreq = binning )

                    if doSkip == False:
                        if useEnd :
                            data3 = (data3.T - offEnd3).T
                        else:
                            data3 = (data3.T - offBeg3).T
                    else:
                        arr = (newtime3<=skip[0]) | (newtime3>=skip[1])
                        newtime3, data3 =newtime3[arr], data3[:,arr]
                        data3[:,newtime3<=skip[0]] = ((data3.T)[newtime3<=skip[0]] - offBeg3).T
                        if specialOffset:
                            addOffset = np.median(data3[:,(newtime3<=skip[0]) & (newtime3>=(skip[0]-specialdt))],axis=1) - np.median(data3[:,(newtime3>=skip[1]) & (newtime3<=(skip[1]+specialdt))],axis=1)
                            data3[:,newtime3>=skip[1]] = (data3[:,newtime3>=skip[1]].T + addOffset).T
                        else:
                            data3[:,newtime3>=skip[1]] = ((data3.T)[newtime3>=skip[1]] - offEnd3).T
   
                        


                    print "Reading data 4"
                    data4 = sf.getSignalGroup( 'Sig4' ).T
                    time4 = sf.getTimeBase( 'Sig4' )
                    tidx4 = np.where((time4>=tBegin)&(time4<=tEnd))[0]

                    offEnd4 = np.mean(data4[time4>10.],axis=0) 
                    offBeg4 = np.mean(data4[time4<0.01],axis=0)              


                    if (tidx4.size!=tidx.size):
                        time4,data4 = self.checkTimebase( time[tidx], time4[tidx4], data4[tidx4])
                        newtime4, data4 = tools.dataBinning( time4, data4.T, samplefreq = binning )
                    else:
                        newtime4, data4 = tools.dataBinning( time[tidx], data4[tidx].T, samplefreq = binning )

                    if doSkip == False:
                        if useEnd:
                            data4 = (data4.T - offEnd4).T
                        else:
                            data4 = (data4.T - offBeg4).T
                    else:
                        arr = (newtime4<=skip[0]) | (newtime4>=skip[1])
                        newtime4, data4 =newtime4[arr], data4[:,arr]
                        data4[:,newtime4<=skip[0]] = ((data4.T)[newtime4<=skip[0]] - offBeg4).T
                        if specialOffset:
                            addOffset = np.median(data4[:,(newtime4<=skip[0]) & (newtime4>=(skip[0]-specialdt))],axis=1) - np.median(data4[:,(newtime4>=skip[1]) & (newtime4<=(skip[1]+specialdt))],axis=1)
                            data4[:,newtime4>=skip[1]] = (data4[:,newtime4>=skip[1]].T + addOffset).T
                        else:
                            data4[:,newtime4>=skip[1]] = ((data4.T)[newtime4>=skip[1]] - offEnd4).T
   
                    
                        self.nRows_new = np.shape(data4)[0]/self.nCols

 
                if (system=='old'):
                    data_old = np.reshape(np.append(data1[:],data2[:7*self.nCols],axis=0),(16,8,data2[0].size))
                    del data1
                    del data2
                    print 'data of the old System is read'
                    rzDiag = 'RZO'
                    self.dataOld = np.swapaxes(data_old.T,1,2)[:,self.useRowsOld,:]

                if (system=='new'):
                    data_new =  np.reshape(np.append(np.append(data2[7*self.nCols:],data3[:],axis=0),data4[:],axis=0),(20,8,data2[0].size))
                    
                    del data2
                    del data3
                    del data4
                    print 'data of the new System is read'
                    rzDiag = 'RZN'
                    self.dataNew = np.swapaxes(data_new.T,1,2)[:,self.useRowsNew,:]

                if (system=='both'):
                    data_old = np.reshape(np.append(data1[:],data2[:7*self.nCols],axis=0),(16,8,data2[0].size))
                    data_new =  np.reshape(np.append(np.append(data2[7*self.nCols:],data3[:],axis=0),data4[:],axis=0),(20,8,data2[0].size))
                    del data1
                    del data2
                    del data3
                    del data4
                    print 'data of both Systems are read'
                    self.dataNew = np.swapaxes(data_new.T,1,2)[:,self.useRowsNew,:]
                    self.dataOld = np.swapaxes(data_old.T,1,2)[:,self.useRowsOld,:]

                self.time =  newtime2
                self.ntime = np.size(self.time)
                
                #self.avail = np.ones((self.nRows,self.nCols))
          
                sf.close()
                           
            except:
                print "Error reading ECI data"
                IPython.embed()
                sf.close()
                self.Status = False
                return False




            try:
                if ((system=='new')  |  (system=='old')):
                    if (system=='old') :
                        rzDiag = 'RZO'
                        useRows = self.useRowsOld
                    else:
                        rzDiag = 'RZN'
                        useRows = self.useRowsNew

                    rzsf = dd.shotfile( rzDiag, Shotnumber, rzExp, rzEd )
                    self.rztime = rzsf.getTimeBase( 'R', tBegin=tBegin, tEnd=tEnd )
                    self.R = rzsf.getSignalGroup( 'R', tBegin=tBegin, tEnd=tEnd )[:,useRows,:]
                    self.z = rzsf.getSignalGroup( 'z', tBegin=tBegin, tEnd=tEnd )[:,useRows,:]
                    self.rhop = rzsf.getSignalGroup( 'rho', tBegin=tBegin, tEnd=tEnd )[:,useRows,:]
                # 
                    self.phi = rzsf.getSignalGroup( 'phi', tBegin=tBegin, tEnd=tEnd )[:,useRows,:]+PHI.ECE()
                    self.dphi = rzsf.getSignalGroup( 'phi', tBegin=tBegin, tEnd=tEnd )[:,useRows,:]
                    self.offset = rzsf.getSignalGroup( 'offset', tBegin=tBegin, tEnd=tEnd )[:,useRows,:]

                    self.Freq  = rzsf.getParameter('PAR', 'freq').data 
                    rzsf.close()

                    self.StatusRz = True
                elif (system=='both'):
                    rzsfNew = dd.shotfile( 'RZN', Shotnumber, rzExp, rzEd )
                    self.rztimeNew = rzsfNew.getTimeBase( 'R', tBegin=tBegin, tEnd=tEnd )
                    self.RNew = rzsfNew.getSignalGroup( 'R', tBegin=tBegin, tEnd=tEnd )[:,self.useRowsNew,:]
                    self.zNew = rzsfNew.getSignalGroup( 'z', tBegin=tBegin, tEnd=tEnd )[:,self.useRowsNew,:]
                    self.rhopNew = rzsfNew.getSignalGroup( 'rho', tBegin=tBegin, tEnd=tEnd )[:,self.useRowsNew,:]
                    self.phiNew = rzsfNew.getSignalGroup( 'phi', tBegin=tBegin, tEnd=tEnd )[:,self.useRowsNew,:]+PHI.ECE()
                    self.dphiNew = rzsfNew.getSignalGroup( 'phi', tBegin=tBegin, tEnd=tEnd )[:,self.useRowsNew,:]
                    self.offsetNew = rzsfNew.getSignalGroup( 'offset', tBegin=tBegin, tEnd=tEnd )[:,self.useRowsNew,:]

                    self.FreqNew  = rzsfNew.getParameter('PAR', 'freq').data   
                    rzsfNew.close()
 

                    rzsfOld = dd.shotfile( 'RZO', Shotnumber, rzExp, rzEd )
                    self.rztimeOld = rzsfOld.getTimeBase( 'R', tBegin=tBegin, tEnd=tEnd )
                    self.ROld = rzsfOld.getSignalGroup( 'R', tBegin=tBegin, tEnd=tEnd )[:,self.useRowsOld,:]
                    self.zOld = rzsfOld.getSignalGroup( 'z', tBegin=tBegin, tEnd=tEnd )[:,self.useRowsOld,:]
                    self.rhopOld = rzsfOld.getSignalGroup( 'rho', tBegin=tBegin, tEnd=tEnd )[:,self.useRowsOld,:]
                    self.phiOld = rzsfOld.getSignalGroup( 'phi', tBegin=tBegin, tEnd=tEnd )[:,self.useRowsOld,:]+PHI.ECE()                   
                    self.dphiOld = rzsfOld.getSignalGroup( 'phi', tBegin=tBegin, tEnd=tEnd )[:,self.useRowsOld,:]
                    self.offsetOld = rzsfOld.getSignalGroup( 'offset', tBegin=tBegin, tEnd=tEnd )[:,self.useRowsOld,:]

                    self.FreqOld  = rzsfOld.getParameter('PAR', 'freq').data   
                    rzsfOld.close()
                   
                #IPython.embed()
                
            except:
                print "Error reading RZ data"
                #IPython.embed()
                self.StatusRz = False
                self.Status = False
                try:
                    rzsf.close()
                except:
                    pass
           
            

            self.Status = True

            return True
        else:
            print 'Diagnotic should be TDI, but is:' ,Diagnostic
            return False


		

    def Unload( self ):
        if self.Status:
            self.Status = False
            del self.tBegin
            del self.tEnd
            del self.time	
            del self.ntime
            del self.data
            del self.nRows
            del self.nCols
            del self.rztime
            del self.R
            del self.z
            del self.phi 
            del self.avail
            if self.StatusTrad:
                del self.Trad
                self.StatusTrad = False
            
            if self.Status_allrhop:
                del self.Rall 
                del self.zall
                del self.rhopall