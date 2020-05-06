import sys
import os
import pylab
import math
import pdb
import ctypes
import time
import scipy.weave
import numpy
import copy
from scipy.weave import converters
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

        self.setParams(start, end, 0.0, 0.0, 0.0)
        
        
    def Unload(self):

        self.statusParam = False
        self.statusBinning = False
        self.statusEqui = False
        self.statusData = False
        self.statusRef = False 

        self.statusSingleFFT = False
    
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

        # calculate time difference
        distance = numpy.diff(timebase)
        stdDistance = numpy.std(distance)
        meanDistance = numpy.mean(distance)

        if (stdDistance/meanDistance < 0.1) :
            return True
        else:
            print 'Timebase not equidistant'
            return False


    def initializeData( self, timebase, data, time_ref = None, signal_ref = None, time_ext = None ):


        if (timebase != None) | (data != None):
            if self.checkTimebase(timebase):
                self.time = timebase
                self.data = data 

            else:
                self.time = self.createEquiTimebase(timebase)       
                self.data = self.transformToTimebase( timebase, data , self.time)
            
            if self.statusBinning:    
                self.time,self.data = self.signalBinning(timebase, data, 1/self.binning)
                
            if (time_ref != None) & (signal_ref != None):
                self.refSignal = self.transformToTimebase( time_ref, signal_ref , self.time)
                self.statusRef = True
            
            self.dtime = numpy.mean(numpy.diff(timebase))
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
    def singleFFT(self, normalise=True):
        
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
    
            #IPython.embed()
        
            for i in range(0, int(nInt)):
            
                mn = startIndex + (2*i+1)*halfIntervalIndex - halfWindowIndex 
                mx = startIndex + (2*i+1)*halfIntervalIndex + halfWindowIndex
               

                if (mn >= 0) & (mx <= numpy.size(self.time)) & ((mx-mn) > 0):
                                               
                    timebase.append(numpy.mean(self.time[mn:mx]))
                    ntimebase = ntimebase + 1
                    freq = numpy.fft.fftfreq(2*halfWindowIndex, self.dtime)[:halfWindowIndex]
                
                #if ref signal avialable, take it, else first signal
                    if self.statusRef:
                        refSig = self.refSignal[mn:mx] - numpy.average( self.refSignal[mn:mx])
                    else:
                        refSig = self.data[mn:mx,0] - numpy.average( self.data[mn:mx,0] )
                        
                    if normalise:
                        stdref = numpy.std(refSig)
                        if stdref > 0.0:
                            refFFT = numpy.fft.fft(refSig/stdref)[:halfWindowIndex]*2.0
                        else:
                            print "nothing to reference"
                            refFFT = numpy.fft.fft(refSig*0.0)[:halfWindowIndex]*2.0
                    else:
                        refFFT = numpy.fft.fft(refSig)[:halfWindowIndex]*2.0
                        refFFT = refFFT/numpy.sum(refFFT)*numpy.mean(self.refSignal[mn:mx])

                    refSigFFT.append(refFFT)
                    

                    #if the signal has only one dimension
                    if self.dataDim == 1 :  

                        channel = self.data[mn:mx] - numpy.mean(self.data[mn:mx])
                        if normalise:
                            stdchannel = numpy.std(channel)
                            if stdchannel > 0.0:
                                channelFFT = numpy.fft.fft(channel/stdchannel)[:halfWindowIndex]*2.0
                                channelFFT =  channelFFT/numpy.sum(channelFFT)
                            else:
                                channelFFT = numpy.fft.fft(channel*0.0)[:halfWindowIndex]*2.0
                        else:
                            channelFFT = numpy.fft.fft(channel)[:halfWindowIndex]*2.0
                            channelFFT = channelFFT/numpy.sum(channelFFT)*numpy.mean(self.data[mn:mx])

                        dataFFT.append(channelFFT)


                    else:
                        if self.dataDim == 2 :

                            allChannels = []
                            for ch in range(numpy.size(self.data,axis=1)):
                                channel = self.data[mn:mx,ch] - numpy.mean(self.data[mn:mx,ch])
                                if normalise:
                                    stdchannel = numpy.std(channel)
                                    if stdchannel > 0.0 :
                                        channelFFT = numpy.fft.fft(channel/stdchannel)[:halfWindowIndex]*2.0
                                        channelFFT =  channelFFT/numpy.sum(channelFFT)
                                    else:
                                        channelFFT = numpy.fft.fft(channel*0)[:halfWindowIndex]*2.0
                                else:
                                    channelFFT = numpy.fft.fft(channel)[:halfWindowIndex]*2.0
                                    channelFFT = channelFFT/numpy.sum(channelFFT)*numpy.mean(self.data[mn:mx,ch])

                                allChannels.append(channelFFT)
                            
                            dataFFT.append(allChannels)
                               # IPython.embed()
                               # allChannels.append(channelFFT*refFFT.conjugate())
                                  
                        
           # IPython.embed()
            self.statusSingleFFT = True
            self.freq = numpy.array(freq)
            self.dataFFT = numpy.array(dataFFT)
            if self.statusRef:
                self.refSigFFT = numpy.array(refSigFFT) 
            self.timebase = numpy.array(timebase)       
            self.ntimebase = int(ntimebase)
            return True
        else:
            print  "Data is not initiliazed"
            return False




    def correlateFFT(self, minF=0.0, maxF=5000.0, takeRefFreq = False, normalise=True):

        if self.statusSingleFFT == False:
            if self.singleFFT() == False:
                print "Could not make FFT of single Signals"
                return False
                
        result = []

        selectedIndex = numpy.arange(numpy.abs(self.freq-minF).argmin(),numpy.abs(self.freq-maxF).argmin())


        #IPython.embed()
        #select index for correlation
        for t in range( self.ntimebase ):

            #dynamic range of freq
            if takeRefFreq:
                #IPython.embed()
                selectedIndex = numpy.abs( self.refSigFFT[t,:] ).argmax()

            bufferRef=[]
            if self.dataDim == 1:
                bufferRef.append( self.dataFFT[t,selectedIndex]*self.refSigFFT[t,selectedIndex].conjugate() )

            if self.dataDim == 2:
                for ch in range( self.nChannels ):  
                    bufferRef.append( self.dataFFT[t,ch,selectedIndex]*self.refSigFFT[t,selectedIndex].conjugate() )
                   
            result.append(bufferRef)

        self.amp = numpy.abs(result)
        self.phase =  numpy.angle(result)

      #      freq = numpy.fft.fftfreq(windowindex, 1e-6)[:windowindex/2]
