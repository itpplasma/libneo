
import math

import scipy.signal
import numpy
import matplotlib.pylab as plt

import IPython


class correlator:
    def __init__(self, start = 1.5, end = 3.0, statusbar=None, statusprefix=''):
       #first index is always timebase
        self.time = []
        self.data = []
        self.refSignal = []
        
        self.statusParam = False
        self.statusBinning = False
        self.statusEqui = False
        self.statusData = False 
        self.statusRef = False 
        
        self.statusSingleFFT = False
        self.statusCC = False
        self.statusDual = False
        self.statusFFT = False

        self.setParams(start, end, 0.0, 0.0, 0.0)
        
        
    def Unload(self):

        self.statusParam = False
        self.statusBinning = False
        self.statusEqui = False
        self.statusData = False
        self.statusRef = False 

        self.statusSingleFFT = False
        self.statusCC = False
        self.statusDual = False
        self.statusFFT = False

    def setParams(self, start, end,  binning=0.0, interval=0.0, window=0.0):
	
	#in sec
        if (start > 0.0) & (start < end):
            self.start = start
            self.end = end
        else:
            print "No valid timerange"
            return

        if window <= 0.0 :
            self.window = end - start
        else:
            if window > end - start:
                self.window = end - start
            else:
                self.window = window


        if interval <= 0.0 :
            self.interval = end - start
        else:
            if interval > end -start:
                self.interval = end - start
            else:
                self.interval = interval

        # no binning takes places, if binning is 0
        if (binning == 0.0):
            print "No binning"
            self.binning = 0.0
            self.statusBinning = False
        else:
            self.binning = binning
            self.statusBinning = True
           # if (binning > 0.0) & ( binning < interval ) & ( binning < window ) :
             #   self.binning = binning
              #  self.statusBinning = True
 
        self.statusParam = True

    def checkTimebase(self, timebase):
        #IPython.embed()

        # calculate time difference
        distance = numpy.diff(timebase)
        stdDistance = numpy.std(distance)
        meanDistance = numpy.mean(distance)

        if (stdDistance/meanDistance < 0.25) :
            return True
        else:
            print 'Timebase not equidistant'
            return False


    def initializeData( self, timebase, data, time_ref = None, signal_ref = None , time_ext = None):

        #IPython.embed()
        if (timebase != None) | (data != None):
            if (time_ext == None):
                if self.checkTimebase(timebase):
                    self.time = timebase
                    self.data = data 
                else:
                    self.time = self.createEquiTimebase(timebase)       
                    self.data = self.transformToTimebase( timebase, data , self.time)
            else:
                if self.checkTimebase(time_ext):         
                    self.time = time_ext
                    self.data = self.transformToTimebase( timebase, data , self.time)
                else:
                    self.time = self.createEquiTimebase(time_ext)       
                    self.data = self.transformToTimebase( timebase, data , self.time)

            if self.statusBinning:    
                self.time,self.data = self.signalBinning(timebase, data, 1/self.binning)
                
            if (time_ref != None) & (signal_ref != None):
                self.refSignal = self.transformToTimebase( time_ref, signal_ref , self.time)
                self.statusRef = True
            
            self.dtime = numpy.mean(numpy.diff(self.time))
            self.dataDim = numpy.ndim(self.data)
            if self.dataDim == 1:
                self.nChannels = 1
            else:
                if self.dataDim == 2:
                    self.nChannels = numpy.size(self.data, axis=1)
            self.statusData = True

        else:
            print "no valid data"
                    
        
    def createEquiTimebase(self, timebase):
         
        distance = numpy.diff(timebase)
        stdDistance = numpy.std(distance)
        meanDistance = numpy.mean(distance)  
        idx = numpy.where( numpy.abs(distance-meanDistance) < stdDistance*2.0 )
        newMeanDiff = numpy.mean(distance[idx])
        return numpy.arange(timebase[0], timebase[-1], newMeanDiff)

    
    def transformToTimebase( self, time_in, signal_in, time_out):

        if (time_in != None) | (signal_in != None) | (time_out != None):
            if numpy.ndim(signal_in) == 1 :
                if numpy.size(signal_in,axis=0) == numpy.size(time_in) :
                    f=scipy.interpolate.interp1d( time_in, signal_in, bounds_error=False, fill_value=0 )
                    return f(time_out)
                else:
                    print "Timebase and data do not have the same size."
                    return 0
            #first index must be time
            if numpy.ndim(signal_in) == 2:
                if numpy.size(signal_in,axis=0) == numpy.size(time_in) :
                    channels = numpy.arange(numpy.size( signal_in, axis=1 ))
                    f=scipy.interpolate.interp2d( time_in,channels, numpy.transpose(signal_in), bounds_error=False, fill_value=0 )
                    return numpy.transpose(f(time_out,channels))
                else:
                    print "Timebase and data do not have the same size."
                    return 0
            print "Dimension of Data differ from 1 or 2."
            return 0
	else:
            print "one signal does not exist"
            return 0
	

    def signalBinning(self, time_in, signal_in, samplefreq = 1.0):

        if (time_in != None) | (signal_in != None):
            
            ntimes = numpy.size(time_in)
            samplingrate = 1.0/numpy.mean(numpy.diff(time_in))

            if  ntimes == numpy.size(signal_in,axis=0):

                print "Binning with ",samplefreq," Hz"
                bins = int(ntimes*(float(samplefreq)/samplingrate))
						
                signal_temp = signal_in
                slices= numpy.linspace(0, ntimes, bins+1, True).astype(numpy.int)
                counts = numpy.diff(slices)

                time_out = numpy.add.reduceat(time_in, slices[:-1]) / counts
                ntimes_new = numpy.size(time_out)

                if numpy.ndim(signal_in) == 1 :
                    signal_out =  numpy.zeros(ntimes)
                    return (time_out,numpy.add.reduceat(signal_in[:], slices[:-1]) / counts)

                if numpy.ndim(signal_in) == 2 :
                    if numpy.size(signal_in,axis=0) == numpy.size(time_in) :
                        signal_out=numpy.zeros((ntimes_new,numpy.size( signal_in, axis=1 )))
                        for i in range(numpy.size( signal_in, axis=1 )):
                            signal_out[:,i] = numpy.add.reduceat(signal_in[:,i], slices[:-1]) / counts
                        
                        return (time_out,signal_out)

                print "Dimension of Data differ from 1 or 2."    
                return (0,0)

# calculate fft of all channels 
    def singleFFT( self ):
        
        self.statusSingleFFT = False
        if self.statusData:

            nInt = int( (self.end - self.start) / self.interval )
            startIndex = numpy.argmin( numpy.abs( ( self.time - self.start ) ) )
            endIndex = numpy.argmin( numpy.abs( ( self.time - self.end ) ) )
            halfIntervalIndex = int(float(endIndex - startIndex) /2. /  (( self.end - self. start) / self.interval ))
            halfWindowIndex =  int(float(endIndex - startIndex) / 2. /  (( self.end - self. start) / self.window ))
        
            timebase = []
            ntimebase = 0
            dataFFT = []
            refSigFFT = []     
            freq = []
            meanChannel = []
            meanFFT = []

            for i in range(0, int(nInt)):
            
                mn = startIndex + (2*i+1)*halfIntervalIndex - halfWindowIndex 
                mx = startIndex + (2*i+1)*halfIntervalIndex + halfWindowIndex
               

                if (mn >= 0) & (mx <= numpy.size(self.time)) & ((mx-mn) > 0):
                                               
                    timebase.append(numpy.mean(self.time[mn:mx]))
                    ntimebase = ntimebase + 1
                    freq.append(numpy.fft.fftfreq(2*halfWindowIndex, self.dtime))
                
                #if ref signal avialable, take it, else first signal
                    if self.statusRef:
                        refSig = self.refSignal[mn:mx] - numpy.average( self.refSignal[mn:mx])
                    else:
                        refSig = self.data[mn:mx,0] - numpy.average( self.data[mn:mx,0] )
                        
                    #IPython.embed()
                    refSigFFT.append(numpy.fft.fft(refSig)/numpy.size(refSig))
                    

                    #if the signal has only one dimension
                    if self.dataDim == 1 :  
                        nData = numpy.size(self.data[mn:mx])
                        channel = self.data[mn:mx] - numpy.mean(self.data[mn:mx])
                        meanChannel.append(numpy.mean(self.data[mn:mx]))
                        dataFFT.append(numpy.fft.fft(channel)/nData)


                    else:
                        if self.dataDim == 2 :
                            nData = numpy.size(self.data[mn:mx,0])
                            allChannels = []
                            for ch in range(numpy.size(self.data,axis=1)):
                                channel = self.data[mn:mx,ch] - numpy.mean(self.data[mn:mx,ch])
                                meanChannel.append(numpy.mean(self.data[mn:mx,ch]))
                                allChannels.append(numpy.fft.fft(channel)/nData)
                            
                            dataFFT.append(allChannels)
                            meanFFT.append(meanChannel)               
                               # IPython.embed()
                               # allChannels.append(channelFFT*refFFT.conjugate())
                                  
                        
            self.statusSingleFFT = True
            self.freq = numpy.array(freq)
            self.dataFFT =numpy.array(dataFFT)
            self.meanFFT = numpy.array(meanFFT)
                    
                   
            if self.statusRef:
                self.refSigFFT = numpy.array(refSigFFT) 
            self.timebase = numpy.array(timebase)       
            self.ntimebase = int(ntimebase)
            return True
        else:
            print  "Data is not initiliazed"
            return False

        


    def dualLockIn( self ):  

        self.statusDual = False

        if self.statusRef == False:
            print "no reference"
            return False

        if self.statusData:

            nInt = int( (self.end - self.start) / self.interval )
            startIndex = numpy.argmin( numpy.abs( ( self.time - self.start ) ) )
            endIndex = numpy.argmin( numpy.abs( ( self.time - self.end ) ) )
            halfIntervalIndex = int(float(endIndex - startIndex) /2. /  (( self.end - self. start) / self.interval ))
            halfWindowIndex =  int(float(endIndex - startIndex) / 2. /  (( self.end - self. start) / self.window ))
        
            timebase = []
            ntimebase = 0
            dualPhase = []
            dualAmp = []
        
            for i in range(0, int(nInt)):
            
                mn = startIndex + (2*i+1)*halfIntervalIndex - halfWindowIndex 
                mx = startIndex + (2*i+1)*halfIntervalIndex + halfWindowIndex
               

                if (mn >= 0) & (mx <= numpy.size(self.time)) & ((mx-mn) > 0):
                                               
                    timebase.append(numpy.mean(self.time[mn:mx]))
                    ntimebase = ntimebase + 1
                              
                #if ref signal avialable, take it, else first signal
                    if self.statusRef:
                        refSig = (self.refSignal[mn:mx] - numpy.mean( self.refSignal[mn:mx]))/numpy.std(self.refSignal[mn:mx])
                    else:
                        refSig = (self.data[mn:mx,0] - numpy.mean( self.data[mn:mx,0] ))/numpy.std(self.data[mn:mx,0])
                        
                    nData = numpy.size(self.data[mn:mx,0])
                    orefSig = numpy.imag(scipy.signal.hilbert(refSig))
                    orefSig = (orefSig - numpy.mean( orefSig ))

                    #IPython.embed()

                    if self.dataDim == 1 :  
                        channel = ( self.data[mn:mx] - numpy.mean(self.data[mn:mx]) ) / numpy.std(self.data[mn:mx])
                        C1 =  numpy.sum( orefSig[:]*channel[:] )/ nData
                        C2 =  numpy.sum( refSig[:]*channel[:] )/nData
                        dualAmp.append(numpy.sqrt( C1 * C1 + C2 * C2 ) )
                        dualPhase.append(numpy.arctan2( C2,C1 ))

                    else:

                        if self.dataDim == 2 :

                            allChannelsPhase = []
                            allChannelsAmp = []

                            for ch in range(numpy.size(self.data,axis=1)):
                                channel = ( self.data[mn:mx,ch] - numpy.mean(self.data[mn:mx,ch]) ) / numpy.std(self.data[mn:mx,ch])
                                C1 =  numpy.sum( orefSig[:]*channel[:] ) / nData
                                C2 =  numpy.sum( refSig[:]*channel[:] ) / nData
                                allChannelsAmp.append(numpy.sqrt( C1 * C1 + C2 * C2 ) )
                                allChannelsPhase.append(numpy.array(numpy.arctan2( C1, C2 )))
                                #IPython.embed()

                            dualPhase.append( allChannelsPhase )
                            dualAmp.append( allChannelsAmp )


            self.phaseDual = numpy.array(dualPhase) 
            self.phaseDual = numpy.squeeze(numpy.unwrap( numpy.array(dualPhase) ))
            while (numpy.mean(self.phaseDual) < 0.0 ): 
                self.phaseDual += numpy.pi
            while (numpy.mean(self.phaseDual) > numpy.pi ):
                self.phaseDual -= numpy.pi
            self.ampDual = numpy.squeeze( numpy.absolute( numpy.array(dualAmp) ))
            self.timeDual = numpy.squeeze(self.phaseDual/2.0/math.pi/self.freq0 )

            self.timebase = numpy.array(timebase)       
            self.ntimebase = int(ntimebase)
            self.statusDual = True
            return True
        else:
            print  "Data is not initiliazed"
            return False



    def crossCorrelation(self):

        self.statusCC = False
        if self.statusData:

            nInt = int( (self.end - self.start) / self.interval )
            startIndex = numpy.argmin( numpy.abs( ( self.time - self.start ) ) )
            endIndex = numpy.argmin( numpy.abs( ( self.time - self.end ) ) )
            halfIntervalIndex = int(float(endIndex - startIndex) /2. /  (( self.end - self. start) / self.interval ))
            halfWindowIndex =  int(float(endIndex - startIndex) / 2. /  (( self.end - self. start) / self.window ))
        
            timebase = []
            ntimebase = 0
            ampCC = []
            timedelayCC = []
            refSigCC = []
        
            for i in range(0, int(nInt)):
            
                mn = startIndex + (2*i+1)*halfIntervalIndex - halfWindowIndex 
                mx = startIndex + (2*i+1)*halfIntervalIndex + halfWindowIndex
               

                if (mn >= 0) & (mx <= numpy.size(self.time)) & ((mx-mn) > 0):
                                               
                    timebase.append(numpy.mean(self.time[mn:mx]))
                    ntimebase = ntimebase + 1
                              
                #if ref signal avialable, take it, else first signal
                    if self.statusRef:
                        refSig = ( self.refSignal[mn:mx] - numpy.mean( self.refSignal[mn:mx]) ) / numpy.std(self.refSignal[mn:mx])
                    else:
                        refSig = ( self.data[mn:mx,0] - numpy.mean( self.data[mn:mx,0]) ) / numpy.std(self.data[mn:mx])
                     
                    nData = numpy.size(self.data[mn:mx,0])

                    refSigCC =  scipy.signal.correlate( refSig,refSig ) / nData               
                    

                    
                    #if the signal has only one dimension
                    if self.dataDim == 1 :  

                        channel = self.data[mn:mx] - numpy.mean(self.data[mn:mx])
                        corr = scipy.signal.correlate( channel, refSig,  mode='full')
                        ampCC.append( numpy.max( corr )/numpy.std(refSig) )
                        timedelayCC.append( ( numpy.argmax( corr[halfWindowIndex+1:-1] ) ) * self.dtime )

                    else:
                        if self.dataDim == 2 :

                            allChannelsTime = []   
                            allChannelsAmp = []

                            for ch in range(numpy.size(self.data,axis=1)):
                                channel = ( self.data[mn:mx,ch] - numpy.mean(self.data[mn:mx,ch]) ) / numpy.std(self.data[mn:mx,ch])
                                corr = scipy.signal.correlate(channel,refSig, mode='full')/ nData
                                allChannelsAmp.append( numpy.max( corr[nData:-1] )  )
                                allChannelsTime.append( ( numpy.argmax( corr[nData:-1] ) ) * self.dtime )                              

                            ampCC.append( allChannelsAmp )
                            timedelayCC.append( allChannelsTime )

                                          # allChannels.append(channelFFT*refFFT.conjugate())
                                  
                        
            self.timeCC = numpy.squeeze( numpy.array( timedelayCC ) )
            self.ampCC =numpy.squeeze( numpy.array( ampCC ) )
 
            self.timebase = numpy.array(timebase)       
            self.ntimebase = int(ntimebase)

            self.statusCC = True
            return True
        else:
            print  "Data is not initiliazed"
            return False


    def correlateFFT(self, minF=0.0, maxF=5000.0, takeRefFreq = False, normalise=True):

        if self.statusSingleFFT == False:
            if self.singleFFT() == False:
                print "Could not make FFT of single Signals"
                return False
                
        refResult = []
        singleResult = []
        absoluteResult = []
        meanResult = []

        singleAmp = []
        singlePhase =  []    
        absoluteAmp = []

        ampFFT = []
        phaseFFT = [] 
        timeFFT = [] 
        freq0 = []

        #selectedIndex = numpy.arange(numpy.abs(self.freq-minF).argmin(),numpy.abs(self.freq-maxF).argmin())

        #IPython.embed()
        #select index for correlation
        for t in range( self.ntimebase ):
            #dynamic range of freq
            if takeRefFreq:
                if self.statusRef:
                    selectedIndex = numpy.abs( self.refSigFFT[t,:int(numpy.size(self.refSigFFT[t])/2.0)] ).argmax()
                    freq = self.freq[t]
                    freq0.append(freq[selectedIndex])

                    print "Use frequency of ",freq[selectedIndex]
                else:
                    print "no Reference Signal"

            refBuffer=[]
            singleBuffer=[]
            absoluteBuffer=[]
            meanBuffer = []

            if self.dataDim == 1:
                refBuffer.append( self.dataFFT[t,selectedIndex]*self.refSigFFT[t,selectedIndex].conjugate() )
                singleResult.append( numpy.sqrt(self.dataFFT[t,selectedIndex].conjugate()*self.dataFFT[t,selectedIndex]) )

            if self.dataDim == 2:
                nData = numpy.size( self.dataFFT[t,0,:] )
                for ch in range( self.nChannels ):  
                    refBuffer.append(  self.dataFFT[t,ch,selectedIndex].conjugate()*self.refSigFFT[t,selectedIndex]/numpy.absolute(self.refSigFFT[t,selectedIndex])*2. )

#/numpy.sum(numpy.absolute(self.dataFFT[t,ch,:]))/numpy.sum(numpy.absolute(self.refSigFFT[t,selectedIndex])) ) 

                    singleBuffer.append(  self.dataFFT[t,ch,selectedIndex]*2.  ) 
                    meanBuffer.append(  numpy.mean(numpy.absolute(self.dataFFT[t,ch,:]) ) ) 


                refResult.append(refBuffer)
                singleResult.append(singleBuffer)
                meanResult.append(meanBuffer)

            
            singleAmp.append( numpy.absolute(singleResult)  )
            singlePhase.append( numpy.unwrap(numpy.angle(singleResult)) )

            if self.statusRef:
                ampFFT.append( numpy.absolute(refResult) )
                phaseFFT.append( numpy.angle(refResult) )
                # phaseFFT.append( numpy.unwrap(numpy.angle(refResult)) )
                timeFFT.append( numpy.unwrap(numpy.angle(refResult))/2.0/math.pi/freq0 ) 

        self.freq0 = numpy.array( freq0 )
        self.singleAmp = numpy.squeeze(numpy.array( singleAmp ))
        self.singlePhase =  numpy.squeeze(numpy.array( singlePhase ) )
        self.meanSig = numpy.squeeze(numpy.array( meanBuffer ) )
 
        if self.statusRef:
            self.ampFFT = numpy.squeeze( numpy.array( ampFFT ) )
            self.phaseFFT =  numpy.squeeze( numpy.unwrap(numpy.array( phaseFFT )) )
            while (numpy.mean(self.phaseFFT) < 0 ): 
                self.phaseFFT += numpy.pi
            while (numpy.mean(self.phaseFFT) > numpy.pi ):
                self.phaseFFT -= numpy.pi
            self.timeFFT = numpy.unwrap(self.phaseFFT) /2.0/math.pi/self.freq0  



       # def leastSquare(self):
