#!/usr/bin/env python

import numpy as np
import dd
from IPython import embed
import matplotlib.pyplot as plt

class Btheta:

    def __init__( self ,  Shotnumber = None ):
        self.Status = False
        if Shotnumber != None :
            self.Load( Shotnumber )
	
    def Load( self, Shotnumber,edition=2 ):
        
        sf = dd.shotfile( 'IDF', Shotnumber, 'AUGE', edition)
        
        self.Bprobr,self.Bprobz,self.Bprobphi =  sf('Bprobr').data, sf('Bprobz').data, sf('Bprobphi').data

        #self.Bprobnam = sf('Bprobnam').data

        self.Bprobnam =sf.getSignal('Bprobnam',dtype='CHAR_8')
        #self.Bprobfla =sf.getSignal('Bprobfla',dtype='bool')
        

        dataTmp = sf('Bprobmes')
        self.Bprobmeas = dataTmp.data
        self.Bprobtime = dataTmp.time
        sf.close()
        
        sfM = dd.shotfile( 'MAY', Shotnumber, 'AUGD', edition)
        dataTmp = sfM(self.Bprobnam[0])
        self.MAYtime  = dataTmp.time
        self.MAYsignal = np.zeros((self.MAYtime.size,self.Bprobnam.size))
  
        i=0

        for sigName in self.Bprobnam:
            try:
                self.MAYsignal[:,i] = sfM.getSignal(sigName)
                
            except:
                print sigName, ' not found ,',i

            i=i+1
             
        self.MAYsignal = self.MAYsignal - self.MAYsignal[self.MAYtime<-8.].mean(axis=0)

        sfM.close()
        

    def getRzphi( self ):
        
        return self.Bprobr, self.Bprobz, self.Bprobphi

    def LoadFromMirnov( self, Shotnumber ,edition=0, n=2 ):

        embed()

        data=[]
        time=[]
        idx=[]
        signal=[]
        torpolRz = []
        torpolRzinfo = np.loadtxt('/afs/ipp/u/mrm/idl/mtr/angle.dat',skiprows=2,usecols=[1,2,3,4])
        sigName = np.loadtxt('/afs/ipp/u/mrm/idl/mtr/angle.dat',skiprows=2,usecols=[0],dtype='string')
                
        Diags=['MIR','MHA','MHB','MHC','MHD','MHE','MHF','MHG','MHH','MAN','MOD']

        if n==2:
            sigName = sigName[(np.array(sigName, dtype='|S3')=='C09') | (np.array(sigName, dtype='|S3')=='C07' )]

        for di in Diags:

            try:
                sf = dd.shotfile( di, Shotnumber, 'AUGD', edition)
                i=0
                for sig in sigName:
                    try:
                        dataTmp= sf( sig )
                        data.append(dataTmp.data)
                        time.append(dataTmp.time)
                        idx.append(i)
                        signal.append(sig)
                        torpolRz.append(torpolRzinfo[i])
                        print sig+' worked.'
                    except:
                        p=3.
 
                    i=i+1

                sf.close()
            except:
                print "Error reading shotfile" 
				
        embed()     
        R0,z0=1.65,0.03
        probeData=np.array(data)
        rzData = np.array(torpolRz)
        angle= np.array(np.arctan2(z0-rzData[:,3],R0-rzData[:,2]))
        torAngle = np.array(rzData[:,0])
        polAngle = np.array(rzData[:,1])

        phi1Data = probeData[((torAngle>4.5) &  (torAngle<4.8))] 
        phi1pol = polAngle[((torAngle>4.5) &  (torAngle<4.8))] 
        polSort = np.argsort(phi1pol)
        phitime = np.array(time)[((torAngle>4.5) &  (torAngle<4.8))].mean(axis=0) 
   
