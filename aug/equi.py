#############################################################################
# THIS FUNCTION EXTRACTS G-EQDSK EQUILIBRIA FROM SHOTFILES AND SAVES THEM   #
#############################################################################
# AUTHOR:  PHILIPP ULBL													    #
# CREATED: 12.12.2019													    #
#############################################################################

#!/usr/bin/env python
# line above sets the Python interpreter

import sys
sys.path.append('/afs/ipp/aug/ads-diags/common/python/lib')
import matplotlib.pylab as plt 			# load plotting library
import dd 					# load latest (!) dd library
# use "import dd_YYYYMMDD as dd" to safely load a specific
# version from /afs/ipp/aug/ads-diags/common/python/lib

import eqtools_modified
import os

import numpy as np

def equi(shot, time, typ='EQH', path='', debug=True):
	"""
	function equi.py
	-------------------------------------------------------
	def equi(shot, time, typ='EQH', path='', debug=True):
	
	writes a gfile for single shot and single time
	equilibrium type default is EQH
	file will be located in a subdirectory for shot in path.
	if path is empty a directory will be created in home.
	if time has 2 elements, they are treated as tmin, tmax

	"""

#############################################################################
# LOAD EQUI
#############################################################################

	home = os.path.expanduser('~')

        #overwrite path if it is empty
	if path == '':
		path = home + '/pyout/equi/'
        
        path = path + str(shot)

	#create dir if not existent
	if not os.path.exists(path):
	    os.mkdir(path, 0755)

	#change dir to shot path
	os.chdir(path)

	try:
	    #load equi with eqtools
	    e = eqtools_modified.AUGData.AUGDDData(shot, shotfile=typ, tspline=True)
	except:
	    print('\nINFO: no ' + typ + ' file found, skipped equilibrium.')
	    return False

#############################################################################
# GET TIMES AND SAVE g-FILE
#############################################################################

	#get times in time window (between tmin and tmax)
	times = e.getTimeBase()
	if len(np.atleast_1d(time)) == 1:
		#get time nearest to time
		diff = abs(times-time)
		ind = np.where(diff==min(diff))[0]
	elif len(time) == 2:
	    #get indices all times in time window specified by tmin and tmax
		ind = np.where(np.logical_and(times >= time[0], times <= time[1]))[0]#take 1st element of tuple
	else:
		print('\n INFO: len of time expected to be 1 or 2. skipped equi. ')
		return False

	#for each time
	for i in range(ind.size):
		#time as string
		tname = '{:4.0f}'.format(int(times[ind[i]]*1e3))

		#print status message
		if debug == True:
			print('(' + str(i+1) + '/' + str(len(ind)) + ') equi: processing shot ' + str(shot) + ' at t=' + tname + 'ms')

		#build name of g-file
		gname = 'g' + str(shot) + '.' + tname + '_' + typ

		#write file with eqtools
		eqtools_modified.filewriter.gfile(e, times[ind[i]], nw=129, nh=129, name = gname, title=typ)
    
	return True
