#!/usr/bin/env python
# encoding: utf-8
"""
calc.py

Created by Alexander Bock on 2011-05-20.
Copyright (c) 2011 __MyCompanyName__. All rights reserved.
"""

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


class calculator:
	def __init__(self, statusbar=None, statusprefix=''):
		self.N1 = []
		self.N2 = []
		self.Te = []
		self.statusbar = statusbar
		self.statusprefix = statusprefix
		
		self.setParams(0, 9, 5, 5)

		#self.start = 0
		#self.end = 6000
		#self.interval = 10000
		#self.win = 1000
		#self.minF = 0
		#self.maxF = 100000
	
	
	def printStatus(self, status):
		s = self.statusprefix + status
		if self.statusbar != None:
			self.statusbar.SetStatusText(s)
		else:
			print s
	
	
	def setParams(self, start, end, interval, win, minF=0, maxF=500):
		self.start = start*1000 # convert to ms
		self.end = end*1000
		
		self.interval = interval # already in ms
		self.win = win
		
		self.minF = minF * 1000.0 # convert to Hz
		self.maxF = maxF * 1000.0
	

	def dualPhaseLockInAmplify(self, argTe, N1, N2, timebase, normalise=True, relative=False, debug=False):
		self.normalise = normalise
		self.relative = relative

		useN2 = True
		if numpy.average(N2) == 0 and numpy.std(N2) == 0:
			useN2 = False

		nInt = (self.end - self.start) / self.interval
		Te = copy.copy(argTe)
		startindex = int((self.start - timebase[0] * 1000) * 1000)
		endindex =     int((self.end - timebase[0] * 1000) * 1000)
		shiftindex = self.interval * 1000
		windowindex = self.win * 1000
		
		self.timebase = []
		for i in range(int(nInt)):
			self.timebase.append(self.start/1000.0 + i*self.interval/1000.0)
			
		chans = len(argTe)
		
		tmpN1 = N1
		if useN2: tmpN2 = N2
		
		self.corrAmp1 = []
		self.corrPhase1 = []
		if useN2:
			self.corrAmp2 = []
			self.corrPhase2 = []
		
		# weave code:
		# lock in
		code = """
		double toReturn = 0;
		double d_avgTe = 0;
		double varTe = 0, stdTe = 0;
		int i = 0;
		
		for (i = 0; i < windowindex; i++) {
			varTe += tmpTe[i] * tmpTe[i];
			d_avgTe += tmpTe[i];
		}
		d_avgTe = d_avgTe/windowindex;
		varTe = varTe/windowindex - d_avgTe*d_avgTe;
		stdTe = sqrt(varTe);

		for (i = 0; i < windowindex; i++) {
			toReturn += tmpTe[i]*ref[i];
		}

		if (rel == 1) {
			if (d_avgTe > 0)
				toReturn = toReturn / d_avgTe;
			else
				toReturn = 0;
		} else if (norm == 1) {
			if (stdTe > 0)
				toReturn = toReturn / stdTe;
			else
				toReturn = 0;
		}
		return_val = toReturn;
		"""
		
		self.printStatus('Calculating...')
		
		dTe = 1.0
		if relative and not normalise: # use fluctuation relative to temperature average
			rel = 1
		else:
			rel = 0
		
		if normalise: # normalise temperature to [-1, 1]
			norm = 1
		else:
			norm = 0

		progressbarLength = 50
		line = '|'
		dot = '\xc2\xb7'
		for i in range(1, int(nInt+1)):
			if i % 10 == 0:
				progress = int(float(i)/float(nInt+1.) * progressbarLength)
				sys.stdout.write("\b" * progressbarLength)
				sys.stdout.write(line*progress + dot*(progressbarLength-progress - 1))
				sys.stdout.flush()
			
			
			mn = startindex + i*shiftindex - windowindex
			mx = startindex + i*shiftindex
			
			# normalise to [-1, 1]			
			oN1 = N1[mn:mx]-pylab.average(N1[mn:mx])
			if useN2: oN2 = N2[mn:mx]-pylab.average(N2[mn:mx])
			oN1 /= (max(oN1) - min(oN1))/2.
		    	if useN2: oN2 /= (max(oN2) - min(oN2))/2.
			
			# shift by pi
			sN1 = pylab.fft(oN1); sN1[windowindex/2:] *= 1j; sN1[:windowindex/2] *= -1j; sN1 = pylab.ifft(sN1)
			if useN2: sN2 = pylab.fft(oN2); sN2[windowindex/2:] *= 1j; sN2[:windowindex/2] *= -1j; sN2 = pylab.ifft(sN2)
			
			# normalise to [-1, 1]			
			sN1 = sN1-pylab.average(sN1)
			sN1 /= numpy.std(sN1)
			if useN2: 
				sN2 = sN2-pylab.average(sN2)
				sN2 /= numpy.std(sN2)
			
			amp1, amp2 = [], []
			pha1, pha2 = [], []
			
			for ch in xrange(chans):
				tmpTe = Te[ch][mn:mx]
				avgTe = float(pylab.average(tmpTe))
				tmpTe -= avgTe # remove DC component
				
				tmpTe = tmpTe.astype(numpy.float64)
				
				ref = oN1.astype(numpy.float64)
				L1a = scipy.weave.inline(code, ['tmpTe', 'ref', 'windowindex', 'rel', 'norm'])
				ref = numpy.real(sN1).astype(numpy.float64)
				L1b = scipy.weave.inline(code, ['tmpTe', 'ref', 'windowindex', 'rel', 'norm'])
				a1 = math.sqrt(numpy.abs(L1a**2 + L1b**2))
				if L1a > 0:                 ph1 = math.atan(L1b/L1a)
				elif L1a < 0 and L1b >= 0:  ph1 = math.atan(L1b/L1a) + math.pi
				elif L1a < 0 and L1b  < 0:  ph1 = math.atan(L1b/L1a) - math.pi
				elif L1a == 0 and L1b > 0:  ph1 = - math.pi/2
				elif L1a == 0 and L1b < 0:  ph1 = math.pi/2
				elif L1a == 0 and L1b == 0: ph1 = 0
				
				if useN2: 
					ref = oN2.astype(numpy.float64)
					L2a = scipy.weave.inline(code, ['tmpTe', 'ref', 'windowindex', 'rel', 'norm'])
					ref = numpy.real(sN2).astype(numpy.float64)
					L2b = scipy.weave.inline(code, ['tmpTe', 'ref', 'windowindex', 'rel', 'norm'])
					a2 = math.sqrt(abs(L2a**2 + L2b**2))
					#print ch, 'weave:', L2a, L2b
					if L2a > 0:                 ph2 = math.atan(L2b/L2a)
					elif L2a < 0 and L2b >= 0:  ph2 = math.atan(L2b/L2a) + math.pi
					elif L2a < 0 and L2b  < 0:  ph2 = math.atan(L2b/L2a) - math.pi
					elif L2a == 0 and L2b > 0:  ph2 = - math.pi/2
					elif L2a == 0 and L2b < 0:  ph2 = math.pi/2
					elif L2a == 0 and L2b == 0: ph2 = 0
				
				amp1.append(a1); pha1.append(ph1)
				if useN2: amp2.append(a2); pha2.append(ph2)

				if debug:
					print ch, '1:', a1, ph1, L1a, L1b, avgTe
					print ch, '2:', a2, ph2, L2a, L2b, numpy.std(tmpTe)
				
			self.corrAmp1.append(amp1); self.corrPhase1.append(pha1)
			if useN2: self.corrAmp2.append(amp2); self.corrPhase2.append(pha2)
		
		self.corrAmp1 = numpy.transpose(self.corrAmp1);
		if useN2: self.corrAmp2 = numpy.transpose(self.corrAmp2); 
		self.corrPhase1 = numpy.transpose(self.corrPhase1);
		if useN2: self.corrPhase2 = numpy.transpose(self.corrPhase2);
		
		self.printStatus('')
		self.printStatus('Calculation finished.')
		return True
	
	
	def compareFFT(self, argTe, order, N1, N2, timebase, normalise=True, relative=False, stride = 1):

		self.normalise = normalise
		self.relative = relative

		chans = len(argTe)
		self.order = order
		self.stride = stride

		self.printStatus('Converting to suitable data type...')
		Te = numpy.array(argTe)
		self.printStatus('Converting to suitable data type finished.')

		nInt = (self.end - self.start) / self.interval
		startindex = int((self.start - timebase[0] * 1000) * 1000)
		endindex =     int((self.end - timebase[0] * 1000) * 1000)
		shiftindex = self.interval * 1000
		windowindex = int(self.win * 1000)

		self.timebase = []
		for i in range(int(nInt)):
			self.timebase.append(self.start/1000.0 + i*self.interval/1000.0)

		corr1, corr2 = [], []

		self.freqs1 = []
		self.freqs2 = []

		progressbarLength = 50
		line = '|'
		dot = '\xc2\xb7'
		for i in range(1, int(nInt+1)):
			if i % 10 == 0:
				progress = int(float(i)/float(nInt+1.) * progressbarLength)
				sys.stdout.write("\b" * progressbarLength)
				sys.stdout.write(line*progress + dot*(progressbarLength-progress - 1))
				sys.stdout.flush()			
			
			mn = startindex + i*shiftindex - windowindex
			mx = startindex + i*shiftindex

			freq = numpy.fft.fftfreq(windowindex, 1e-6)[:windowindex/2]

			ref = N1[mn:mx] - numpy.average(N1[mn:mx])
			IPython.embed()
			fft = numpy.fft.fft(ref)[:windowindex/2]
			cutoffindex = numpy.abs(freq-1000).argmin()
			fft[:cutoffindex] = 0
			maximum1 = numpy.abs(fft).argmax() # use this for comparing the channels with each other
			self.freqs1.append(freq[maximum1])

			ref = N2[mn:mx] - numpy.average(N2[mn:mx])
			fft = numpy.fft.fft(ref)[:windowindex/2]
			cutoffindex = numpy.abs(freq-1000).argmin()
			fft[:cutoffindex] = 0
			maximum2 = numpy.abs(fft).argmax() # use this for comparing the channels with each other
			self.freqs2.append(freq[maximum2])

			corrline1, corrline2 = [], []

			for ch in xrange(chans-stride):
				ch1 = order[ch]
				ch2 = order[ch+stride]
				if normalise:
					channel1 = Te[ch1, mn:mx] - numpy.average(Te[ch1, mn:mx])
					if numpy.std(channel1) > 0: channel1 = numpy.fft.fft(channel1/numpy.std(channel1))[:windowindex/2]
					else: channel1 = numpy.fft.fft(channel1*0)[:windowindex/2]

					channel2 = Te[ch2, mn:mx] - numpy.average(Te[ch2, mn:mx])				
					if numpy.std(channel2) > 0: channel2 = numpy.fft.fft(channel2/numpy.std(channel2))[:windowindex/2]
					else: channel2 = numpy.fft.fft(channel2*0)[:windowindex/2]
				else:
					channel1 = numpy.fft.fft(Te[ch1, mn:mx] - numpy.average(Te[ch1, mn:mx]))[:windowindex/2]
					channel2 = numpy.fft.fft(Te[ch2, mn:mx] - numpy.average(Te[ch2, mn:mx]))[:windowindex/2]
				corrline1.append(channel1[maximum1]*channel2[maximum1].conjugate())
				corrline2.append(channel1[maximum2]*channel2[maximum2].conjugate())
			
			corr1.append(corrline1)
			corr2.append(corrline2)

		self.correlation1 = numpy.array(corr1)
		self.correlation2 = numpy.array(corr2)
		self.printStatus('\nCorrelation finished.')
		return True


	def correlateNeighbors(self, argTe, order, timebase, normalise=True, relative=False, stride = 1):

		self.normalise = normalise
		self.relative = relative

		chans = len(argTe)
		self.order = order
		self.stride = stride

		self.printStatus('Converting to suitable data type...')
		Te = numpy.array(argTe)
		self.printStatus('Converting to suitable data type finished.')

		nInt = (self.end - self.start) / self.interval
		startindex = int((self.start - timebase[0] * 1000) * 1000)
		endindex =     int((self.end - timebase[0] * 1000) * 1000)
		shiftindex = self.interval * 1000
		windowindex = self.win * 1000

	#	print startindex,endindex

		self.timebase = []
		for i in range(int(nInt)):
			self.timebase.append(self.start/1000.0 + i*self.interval/1000.0)

		corr = []

		progressbarLength = 50
		line = '|'
		dot = '\xc2\xb7'
		for i in range(1, int(nInt+1)):
			if i % 10 == 0:
				progress = int(float(i)/float(nInt+1.) * progressbarLength)
				sys.stdout.write("\b" * progressbarLength)
				sys.stdout.write(line*progress + dot*(progressbarLength-progress - 1))
				sys.stdout.flush()
			
			
			mn = startindex + i*shiftindex - windowindex
			mx = startindex + i*shiftindex

			corrline = []

			for ch in xrange(chans-stride):
				ch1 = order[ch]
				ch2 = order[ch+stride]
				channel1 = Te[ch1, mn:mx] - numpy.average(Te[ch1, mn:mx])
				channel2 = Te[ch2, mn:mx] - numpy.average(Te[ch2, mn:mx])
				if normalise and numpy.std(channel1) > 0 and numpy.std(channel2) > 0:
					corrline.append(sum(channel1*channel2)/numpy.std(channel1)/numpy.std(channel2))
				elif normalise:
					corrline.append(0)
				else:
					corrline.append(sum(channel1*channel2))

			corr.append(corrline)

		self.correlation = numpy.array(corr)
		self.printStatus('\nCorrelation finished.')
		return self.correlation	


	def correlate(self, argTe, N1, N2, timebase, normalise=True, relative=False):

		self.normalise = normalise
		self.relative = relative

		chans = len(argTe)

		self.printStatus('Converting to suitable data type...')
		Te = numpy.array(argTe)
		self.printStatus('Converting to suitable data type finished.')

		nInt = (self.end - self.start) / self.interval
		startindex = int((self.start - timebase[0] * 1000) * 1000)
		endindex =     int((self.end - timebase[0] * 1000) * 1000)
		shiftindex = self.interval * 1000
		windowindex = self.win * 1000

		self.timebase = []
		for i in range(int(nInt)):
			self.timebase.append(self.start/1000.0 + i*self.interval/1000.0)

		corr1 = []
		corr2 = []

		progressbarLength = 50
		line = '|'
		dot = '\xc2\xb7'
		for i in range(1, int(nInt+1)):
			if i % 10 == 0:
				progress = int(float(i)/float(nInt+1.) * progressbarLength)
				sys.stdout.write("\b" * progressbarLength)
				sys.stdout.write(line*progress + dot*(progressbarLength-progress - 1))
				sys.stdout.flush()
			
			
			mn = startindex + i*shiftindex - windowindex
			mx = startindex + i*shiftindex

			corrline1, corrline2 = [], []

			ref1 = N1[mn:mx] - numpy.average(N1[mn:mx])
			ref1 /= numpy.std(ref1)
		
			ref2 = N2[mn:mx] - numpy.average(N2[mn:mx])
			ref2 /= numpy.std(ref2)

			for ch in xrange(chans):
				channel = Te[ch, mn:mx] - numpy.average(Te[ch, mn:mx])
				corrline1.append(sum(channel*ref1)/numpy.std(channel) if numpy.std(channel) > 0 else 0)
				corrline2.append(sum(channel*ref2)/numpy.std(channel) if numpy.std(channel) > 0 else 0)
			
			corr1.append(corrline1)
			corr2.append(corrline2)

		self.correlation1 = numpy.array(corr1)
		self.correlation2 = numpy.array(corr2)
		self.printStatus('\nCorrelation finished.')
		return [self.correlation1, self.correlation2]


