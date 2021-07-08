#############################################################################
# THIS SCRIPT READS DATA FROM MIRNOV COILS AND WRITES ONE FILE PER COIL		#
# WITH TIME AND DATA.													    #
#############################################################################
# AUTHOR:  PHILIPP ULBL							    						#
# CREATED: 18.12.2019							    						#
#############################################################################

#!/usr/bin/env python
# line above sets the Python interpreter

import sys

sys.path.append('/afs/ipp/aug/ads-diags/common/python/lib')
import matplotlib.pylab as plt 			# load plotting library
import dd 					# load latest (!) dd library
# use "import dd_YYYYMMDD as dd" to safely load a specific
# version from /afs/ipp/aug/ads-diags/common/python/lib

#sys.path.append('/afs/ipp-garching.mpg.de/home/u/ulblp/.local/lib/python2.7/site-packages')
import map_equ 					# load map_equ library

import numpy as np
import os
import sys

import matplotlib.pyplot as plt

path = os.path.expanduser('~') + '/pyout/mirnov/'

#create dir if not existent
if not os.path.exists(path):
	os.mkdir(path, 0755)

#############################################################################
# PARAMETER: SHOT NUMBERS
#############################################################################

shots = [34548]
debug = True

#############################################################################
# READ DATA AND WRITE TO FILE
#############################################################################

for shot in shots:
	#create directory if not existent
	dirname = path + str(shot)
	if not os.path.exists(dirname):
		os.mkdir(dirname, 0755)

	try:
		mas = dd.shotfile('MAS', shot, experiment='AUGD') # instantiate and open shotfile
	except:
		print('\nINFO: no MAS file found, skipped shot ' + str(shot))
		continue

	names = mas.getObjectNames()
	signals = ()

	for n in names:
		#fetch poloidal fields
		if names[n][0:2] == 'Bp':
			signals = signals + (mas(names[n]),)

	for i in range(len(signals)):
		
		#print status message
		if debug == True:
			print('(' + str(i+1) + '/' + str(len(signals)) + ') mirnov: processing shot ' + str(shot) + ', signal: ' + signals[i].name)
	
		#concat arrays along columns
		data = np.transpose(np.stack((signals[i].time, signals[i].data)))

		#build filename
		fname = str(shot) + '_' + signals[i].name

		#save file as txt with ending dat
		np.savetxt(dirname + '/' + fname + '.dat', data)





