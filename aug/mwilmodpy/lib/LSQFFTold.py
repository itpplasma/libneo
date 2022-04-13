import math
import scipy.signal
import numpy
import matplotlib.pylab as plt
import scipy.optimize as optimization
import IPython
import types


class LSQFFT:

    def __init__( self ):
       #first index is always timebase
        self.time = []
        self.data = []
        self.refSignal = []
        #self.freq = []
        
        self.status  = False
        self.statusFreq  = False
        self.statusUnc  = False
        self.statusCalculated  = False

    def Unload(self):

        self.status  = False
        self.statusFreq  = False
        self.statusUnc  = False
        self.statusCalculated  = False

    def initializeData( self, timebase, data, order=1,poly=1,unc = None, time_ref = None, signal_ref = None, freq_in = None, tBegin = None, tEnd = None, negSign=False):

        #check if timebase is available
        if (numpy.all(timebase) != None) | (numpy.all(data) != None):
            #check if timebase and data have the same size
            if (numpy.size(data,axis=0) == numpy.size(timebase)):
                self.order = order
                self.poly = poly
                #If freq is not given use external
                if (numpy.all(freq_in) != None):
                     if (freq_in > 0):
                         self.freq0 = freq_in
                         self.n = 2.*numpy.pi*self.freq0
                     else:
                         print "frequency is not positive"

                elif ( numpy.all(time_ref) != None) & (numpy.all(signal_ref) != None) :
                    dtime = numpy.mean(numpy.diff(time_ref))
                    refSigFFT = numpy.fft.fft(signal_ref-numpy.mean(signal_ref))/numpy.size(signal_ref)
                    reffreq = numpy.fft.fftfreq(numpy.size(signal_ref), dtime)
                    selectedIndex = numpy.absolute( refSigFFT[:int(numpy.size(refSigFFT)/2.0)] ).argmax()
                    self.freq0 = reffreq[selectedIndex]
                    self.n = 2.*numpy.pi*self.freq0
                    self.statusFreq = True
                else:
                    print "no way to get frequency"
            else:
                print "Timebase is data have not the same size"
            

            if (  numpy.all(tBegin) != None) & ( numpy.all(tEnd) != None) :

               idx = numpy.where( ( timebase >= tBegin ) & ( timebase <= tEnd)  )
               self.time = timebase[idx]
               self.data = data[idx,:]
               if (numpy.all(unc) != None):
                   try:
                       self.unc = unc[idx,:]
                       self.statusUnc = True
                   except:
                       self.unc = None
                       self.statusUnc = False

               self.status = True
               if negSign:
                   self.sign = -1.
               else:
                   self.sign = 1.
            else:
                
                self.time = timebase
                self.data = data
                if (numpy.all(unc) != None):
                    try:
                        self.unc = unc
                        self.statusUnc = True
                    except:
                        self.unc = None
                        self.statusUnc = False

                self.status = True
                if negSign:
                    self.sign = -1.
                else:
                    self.sign = 1.                
        else:

            print "no timebase or/and data"

        # trigometric formula to calculate A*sin(wt+phi) as A*cos(phi)*sin(wt)+A*sin(phi)*cos(wt)
 #   def func(self, x, a, b, c):
 #       return  a*numpy.sin(self.n*x)+b*numpy.cos(self.n*x) + c

#    def getFunc(self, idx1,idx2, order=1):
#        return self.amp[idx1,idx2]*numpy.sin(2.*numpy.pi*self.freq0*self.time+self.phase[idx1,idx2])+self.mean[idx1,idx2]


##take harmonics into account
#    def harmFunc(self, x, a1, a2, a3, a4, a5, a6, a0):
#        self.n=2.*numpy.pi*self.freq0
#        return  a1*numpy.sin(self.n*x)+a2*numpy.cos(self.n*x)+a3*numpy.sin(2.*self.n*x)+a4*numpy.cos(2.*self.n*x)+a5*numpy.sin(3.*self.n*x)+a6*numpy.cos(3.*self.n*x) + a0

    def getFunc(self, timeIn=None,idxArrIn = None, idx=None , order = None, poly=None):
        
        if numpy.all(order) == None:
            order = self.order

        if numpy.all(poly) == None:
            poly = self.poly

        if numpy.all(timeIn) == None:
            timeOrig = self.time
        else:
            timeOrig = timeIn
            
        idxArr = numpy.zeros(self.dataShape,dtype='bool')  
        if (numpy.all(idxArrIn) == None) & (numpy.all(idx) == None):
            idxArr[:]=True

        elif numpy.all(idx) != None:           
            if numpy.size(idx) == 1:
                numpy.ravel(idxArr)[idx]=True
            elif numpy.size(idx) == 2:
                idxArr[idx[0],idx[1]]=True
            elif numpy.size(idx) == 3:
                idxArr[idx[0],idx[1],idx[2]]=True 
            else:
                idx = [0]
                numpy.ravel(idxArr)[idx]=True


        elif ((type(idxArrIn) == types.BooleanType) & (numpy.shape(idxArrIn)==numpy.shape(idxArr) )):
            idxArr=idx
        else:
            #print 
            idxArr[:]=True

        
        nChan = numpy.where(idxArr)[0].size
        time=numpy.array([timeOrig,]*nChan).T
        #IPython.embed()

        if order == 1:
            if self.order == 1:
                f = self.amp[idxArr]*numpy.sin(self.n*time+self.phase[idxArr]*self.sign) + self.mean[idxArr] #+ self.off[idxArr]
            else:
                f = self.amp[0][idxArr]*numpy.sin(self.n*time+self.phase[0][idxArr]*self.sign)+ self.mean[idxArr]#+ self.off[idxArr]
           
        elif order >= 2:
            f = numpy.array([self.mean[idxArr],]*timeOrig.size)
            for i in numpy.arange(0,order):
                f += self.amp[i][idxArr]*numpy.sin(float(i+1)*self.n*time+self.phase[i][idxArr]*self.sign)


        #IPython.embed()
        
        if poly==0:
            f += self.off[idxArr]
        elif poly >= 1:           
            for j in numpy.arange(0,poly+1):
                f += self.polyCoeff[j][idxArr]*time**float(j)

        if numpy.all(idxArr==True):
            return numpy.reshape(f,numpy.append(timeOrig.size,self.dataShape ) )
        else:
            return f
###
###    """ 
# def getFuncArr(self, idx):
#        f = numpy.zeros_like(self.time)
#
#        if self.order == 1:
#            f = self.amp[idx]*numpy.sin(self.n*self.time+self.phase[idx]) + self.mean[idx]#
#
#        if self.order >= 2:
#            f = numpy.zeros_like(numpy.size(self.time),self.mean[idx])
#            f[:] = self.mean[idx]
#            for i in numpy.arange(0,self.order):
#                f = f + self.amp[i,idx1,idx2]*numpy.sin(float(i+1)*self.n*self.time+self.phase[i,idx1,idx2])
#
#        return f
#"""

    def leastSquareFit(self):
            
        if self.status:
            
            nChannels = numpy.size(self.data[0])
            ntime = numpy.size(self.time)

            if numpy.size(numpy.shape(self.data)) >= 2:
                dataShape = numpy.shape(self.data[0])
                data = numpy.reshape(self.data,(ntime,nChannels))
                dataMean = numpy.mean( data, axis=0 )
                factor = numpy.std( data, axis=0 )
            else:
                data = numpy.reshape(self.data,(ntime,1))
                dataShape = numpy.shape([1])
                dataMean = numpy.array([numpy.mean( data, axis=0 ),])
                factor = numpy.array([numpy.std( data, axis=0 ),])#0e19
                 

            
            #dataSigma =  numpy.std( self.data, axis=1 )
            sigma = numpy.ones((ntime))
 
            if self.order >= 1:
               
                ###Ak sin(k n x +phi) = A(k) sin(kwt) cos(phi) + A(k) cos(kwt) sin(phi) +C
                ### a(k) = A(k)cos(phi)
                ### b(k) = A(k)sin(phi)
                ### Ak sin(k n x +phi) = a(k) sin(kwt) + b(k) cos(kwt) = func(t)
                ### phi = arctan (bk/ak)
                ### A(k) = (a(k)+b(k))/(cos(phi)+sin(phi))

                ## linarized equation system A*b =y
                ## A ... sin(kwt), cos(kwt)

                A=[]
                
                # go through all harmonics until order of self.order   
                for i in numpy.arange(0,self.order):
                    ## definning colums for A
                    ## calculating  sin(kwt) , cos(kwt)...
                    A.append(numpy.sin(float(i+1)*self.n*self.time))
                    A.append(numpy.cos(float(i+1)*self.n*self.time)*self.sign)
                
                if self.poly>0:
                    for j in numpy.arange(self.poly,0,-1):
                        A.append( numpy.ones(numpy.size(self.time))*self.time**float(j) )
                
                ## adding ones for offset
                A.append( numpy.ones(numpy.size(self.time)) )
                A = numpy.vstack(A).T                                
                
                if self.poly > 0:
                    polyCoeff = numpy.zeros((self.poly+1,nChannels))

                amp = numpy.zeros((self.order,nChannels))
                phase = numpy.zeros((self.order,nChannels))
                uncPhase = numpy.zeros((self.order,nChannels))
                uncAmp = numpy.zeros((self.order,nChannels))
                
                off = numpy.zeros((nChannels))

                covA = numpy.linalg.inv(numpy.dot(A.T,A))
                n = numpy.size(self.time)
                m = A[0].size
                runM = numpy.arange(m)
               # cov = numpy.dot(A.T,A)
              #  print 'IPython'
              #  IPython.embed()
                dataIsInf = numpy.isinf(dataMean)
                for ch in numpy.squeeze(numpy.where(dataIsInf == False)):
                # substract mean to calculate y
                    y = numpy.array( (data[:,ch]-dataMean[ch])/factor[ch] )
                   
       #not used
                #param = optimization.curve_fit(self.func, self.time,y, x0, sigma)[0]
                    paramLin,resi,rank,s = numpy.linalg.lstsq(A, y)
                    a = paramLin[::2]
                    b = paramLin[1::2]
                    aSig = numpy.zeros_like(a)
                    bSig = numpy.zeros_like(b)

                    try:
                        Sig = numpy.sqrt(covA[runM,runM]*resi/float(n-m))
                        aSig = Sig[::2]
                        bSig = Sig[1::2]
                    except:
                        aSig[:] = 0.0
                        bSig[:] = 0.0
              #      if phaseSame:

               ###Ak sin(k n x +phi) = A(k) sin(kwt) cos(phi) + A(k) cos(kwt) sin(phi)
                ### a(k) = A(k)cos(phi)
                ### b(k) = A(k)sin(phi)
                ### phi = arctan (bk/ak)
                ### A(k) = (a(k)+b(k))/(cos(phi)+sin(phi))

                    for i in numpy.arange(0,self.order):
                        phase[i,ch] = numpy.arctan2(b[i],a[i])
                        amp[i,ch] = ( a[i]+b[i])/( numpy.sin(phase[i,ch]) + numpy.cos(phase[i,ch]) )*factor[ch]
                        ##d/dx arctan=1/(1+x*2)
                        ##sig = 1/(1+(b/a)*2) sqrt(bSig**2( 1/a**2 ) +aSig**2( (-b/a**2)**2 ) )
                        
                        uncPhase[i,ch] = 1./(a[i]**2.+b[i]**2.)*numpy.sqrt( aSig[i]**2.*(b[i]**2.)+ bSig[i]**2.*( a[i]**2. ) )
                        uncAmp[i,ch] = numpy.sqrt( aSig[i]**2.+bSig[i]**2.)/numpy.abs( numpy.sin(phase[i,ch]) + numpy.cos(phase[i,ch]) )*factor[ch]

                    if self.poly>0:
                        for j in numpy.arange(0,self.poly+1):
                            polyCoeff[j,ch] = paramLin[-(j+1)]*factor[ch]

                    off[ch] = a[-1]*factor[ch]
                    
                #IPython.embed()
                
                self.amp = numpy.squeeze(numpy.reshape(amp,numpy.append(self.order,dataShape) ))
                self.phase = numpy.squeeze(numpy.reshape(phase,numpy.append(self.order,dataShape) ))
                self.nPhase = numpy.squeeze(numpy.reshape(phase,numpy.append(self.order,dataShape) ))/self.n
                self.uncPhase = numpy.squeeze(numpy.reshape(uncPhase,numpy.append(self.order,dataShape) ))
                self.uncAmp = numpy.squeeze(numpy.reshape(uncAmp,numpy.append(self.order,dataShape) ))
                if self.poly>0:
                    self.polyCoeff = numpy.squeeze(numpy.reshape(polyCoeff,numpy.append(self.poly+1,dataShape) ))

                self.mean =  numpy.squeeze(numpy.reshape(dataMean,(dataShape)))
                self.factor = numpy.squeeze(numpy.reshape(factor,(dataShape)))
                self.off = numpy.squeeze(numpy.reshape(off,(dataShape)))


                #IPython.embed()

                if (numpy.squeeze(dataShape)==1).all():
                    self.amp = numpy.array([self.amp,]).T
                    self.mean = numpy.array([self.mean,]).T
                    self.off = numpy.array([self.off,]).T
                    self.phase = numpy.array([self.phase,]).T
                    self.nPhase = numpy.array([self.nPhase,]).T
                    self.uncPhase = numpy.array([self.uncPhase,]).T
                    self.uncAmp = numpy.array([self.uncAmp,]).T
                    if self.poly>0:
                        self.polyCoeff = numpy.array([self.polyCoeff,]).T
                    self.factor = numpy.array([self.factor,]).T

                    #self.mean = numpy.array([self.mean,]).T
                    #self.off = numpy.array([self.off,]).T

                self.real = self.amp*numpy.cos(self.phase)
                self.imag = self.amp*numpy.sin(self.phase)
                self.dataShape = dataShape

                #IPython.embed()

                self.delay = self.phase/(2.*numpy.pi*self.freq0)
                self.statusCalculated = True

            return True
               # self.amp = 
        else:
            return False


