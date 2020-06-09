#############################################################################
# THIS FUNCTION EXTRACTS COIL DATA FROM SHOTFILES AND SAVES THEM            #
#############################################################################
# AUTHOR:  PHILIPP ULBL							    						#
# CREATED: 17.12.2019							   							#
#############################################################################

#!/usr/bin/env python
# line above sets the Python interpreter

import sys

sys.path.append('/afs/ipp/aug/ads-diags/common/python/lib')
import matplotlib.pylab as plt 			# load plotting library
import dd 								# load latest (!) dd library
# use "import dd_YYYYMMDD as dd" to safely load a specific
# version from /afs/ipp/aug/ads-diags/common/python/lib

import numpy as np
import os
import sys

def coil(shot, time, AVG = 0.01, path='', debug=True):
	"""
	function coil.py
	-------------------------------------------------------
	def coil(shot, time, AVG = 0.01, path='', debug=True):
	
	writes a coil file for a single shot and single time.
	currents are averaged over +-AVG/2 seconds.
	file will be located in a subdirectory for shot in path.
	if path is empty a directory coil will be created in home.
	if time has 2 elements, they are treated as tmin, tmax

	"""

	#overwrite path if it is empty
	if path == '':
		path = os.path.expanduser('~') + '/pyout/coil/'

	#time window for average coil current in seconds
	

#############################################################################
# IMPORT IDA: Te, ne
#############################################################################

	try:
		ncoil = 8

		maw = dd.shotfile('MAW', shot, experiment='AUGD') # instantiate and open shotfile
		IBl = []
		IBu = []
		#get signals from shotfile and put them into lists
		for i in range(ncoil):
		    IBl.append(maw('IBl' + str(i+1)))
		    IBu.append(maw('IBu' + str(i+1)))
	
		maw.close()					  # close shotfile (optional)
	except:
		print('\nINFO: no MAW file found, skipped coil data.')
		return False

#############################################################################
# EXTRACT AND WRITE 2 FILE
#############################################################################

	#create directory if not existent
	dirname = path + str(shot)
	if not os.path.exists(dirname):
	    os.mkdir(dirname, 0755)

	#get times
	times = IBl[0].time
	if len(np.atleast_1d(time)) == 1:
		#get time nearest to time
		diff = abs(times-time)
		ind = np.where(diff==min(diff))[0]
	elif len(time) == 2:
	    #get indices all times in time window specified by tmin and tmax
		ind = np.where(np.logical_and(times >= time[0], times <= time[1]))[0]#take 1st element of tuple
	else:
		print('\n INFO: len of time expected to be 1 or 2. skipped coil. ')
		return False

	#iterate over times
	for i in range(len(ind)):
		#time as string
		tname = '{:4.0f}'.format(times[ind[i]]*1e3)
		
		#print status message
		if debug == True:
			print('(' + str(i+1) + '/' + str(len(ind)) + ') coils: processing shot ' + str(shot) + ' at t=' + tname + 'ms')
	    #init data arrays
		datal = np.zeros((ncoil, 1))
		datau = np.zeros((ncoil, 1))

		#compute time windows for average
		#get time - AVG/2
		diff = abs(times-times[ind[i]]+(AVG/2))#minus changes to plus avg/2 because of bracket
		ind1 = np.where(diff == min(diff))[0][0]
		#get time + AVG/2
		diff = abs(times-times[ind[i]]-(AVG/2))
		ind2 = np.where(diff == min(diff))[0][0]

	    #fill data arrays with data from shot file
		for j in range(ncoil):
			datal[j] = np.mean(IBu[j].data[ind1:ind2]) #avg over times
			datau[j] = np.mean(IBl[j].data[ind1:ind2])

	    #concat arrays along columns
		prof = np.transpose(np.concatenate((datal, datau)))

	    #build filename
		fname = str(shot) + '.' + tname + '_coil'

	    #save file as txt with ending dat
		np.savetxt(dirname + '/' + fname + '.dat', prof)

	return True





