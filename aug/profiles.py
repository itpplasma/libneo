#############################################################################
# THIS FUNCTION EXTRACTS PROFILES FROM SHOTFILES AND SAVES THEM             #
#############################################################################
# AUTHOR:  PHILIPP ULBL							    						#
# CREATED: 12.12.2019							    						#
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

def profiles(shot, time, exper='AUGD', path='', area='rho_pol', debug=True):
	"""
	function profiles.py
	-------------------------------------------------------
	def profiles(shot, time, exper='AUGD', path='', area='rho_pol', debug=True):
	
	writes 4 profiles files (.dat) for Te, ne, Ti and vt
	file will be located in a subdirectory for shot in path.
	if path is empty a directory coil will be created in home.
	if time has 2 elements, they are treated as tmin, tmax
        area specifies the x-value of the profile: r, rho_pol or rho_tor

	"""

	#overwrite path if it is empty
	if path == '':
		path = os.path.expanduser('~') + '/pyout/profiles/'

#############################################################################
# IMPORT IDA: Te, ne
#############################################################################

	quant = ()
	typ = ()

        #try to find PED and then IDI
        try:
                ped = dd.shotfile('PED', shot, experiment=exper)
                Te = ped('TeFit')
                Te.name = 'Te'
                ne = ped('neFit')
                ne.name = 'ne'
                ped.close()
                #append to quant tuple
                quant = quant + (Te, ne)
                typ = typ + ('PED_' + exper, 'PED_' + exper)
        except:
                print('\nINFO: no PED file for Ti, vt found for experiment ' + exper)

	try:
		ida = dd.shotfile('IDA', shot, experiment='AUGD') # instantiate and open shotfile
		Te = ida('Te') 					  # grab electron temperature
		ne = ida('ne')
		ida.close()					  # close shotfile (optional)

		#put quantities in a tuple
		quant = quant + (Te, ne)
		typ = typ + ('IDA', 'IDA')
	except:
		print('\nINFO: no IDA file for Te, ne found.')
        
#############################################################################
# IMPORT PED: Ti, vt
#############################################################################

	#try to find PED and then IDI
	try:
		ped = dd.shotfile('PED', shot, experiment=exper)
		Ti = ped('TiFit')
		Ti.name = 'Ti'
		vt = ped('vTFit')
		vt.name = 'vt'
		ped.close()
		#append to quant tuple
		quant = quant + (Ti, vt)
		typ = typ + ('PED_' + exper, 'PED_' + exper)
	except:
	        print('\nINFO: no PED file for Ti, vt found for experiment ' + exper)
	try:
		idi = dd.shotfile('IDI', shot, experiment='AUGD')
		Ti = idi('TiFit')
		Ti.name = 'Ti'
		vt = idi('vTFit')
		vt.name = 'vt'
		idi.close()
		#append to quant tuple
		quant = quant + (Ti, vt)
		typ = typ + ('IDI', 'IDI')
	except:
	        print('\nINFO: no IDI file for  Ti, vt found.')

#############################################################################
# EXTRACT EQUILIBRIUM USED FOR MAPPING rho_pol to r or rho_tor
#############################################################################

	#create map equ instance
	eqm = map_equ.equ_map()
	#open eqm with cliste equilibrium
	eqm.Open(shot, 'EQH')

#############################################################################
# EXTRACT AND WRITE 2 FILE
#############################################################################

	#create directory if not existent
	dirname = path + str(shot)
	if not os.path.exists(dirname):
	    os.mkdir(dirname, 0755)

	#iterate quants
	for q in range(len(quant)):
		#get times
		times = quant[q].time
		
		if len(np.atleast_1d(time)) == 1:
			#get time nearest to time
			diff = abs(times-time)
			ind = np.where(diff==min(diff))[0]
		elif len(time) == 2:
	    	#get indices all times in time window specified by tmin and tmax
			ind = np.where(np.logical_and(times >= time[0], times <= time[1]))[0]#take 1st element of tuple
		else:
			print('\n INFO: len of time expected to be 1 or 2. skipped profiles. ')
			return False

	    #go through all times in time window
		for i in range(len(ind)):
			#time as string
			tname = '{:4.0f}'.format(times[ind[i]]*1e3)

	        #print status message
			if debug == True:
				print('(' + str(i+1) + '/' + str(len(ind)) + ') ' + quant[q].name + ': processing shot ' + str(shot) + ' at t=' + tname + 'ms') 
        
			#extract data and rho
			rho_pol = quant[q].area[ind[i],:]
			data = np.atleast_2d(quant[q].data)[ind[i],:] #atleast 2d because some PED data may be only 1D which breaks indexing
			#map to rho <=1 ... inside separatrix
			#data = data[rho_pol < 1]
			#rho_pol = rho_pol[rho_pol < 1]

                        #use r, rho_tor or rho_pol as x-axis
                        if area == 'r':
                        
		                #IF R AS AREA
			        #---------------------------------------------------------------------------
			        #calculate R out of rho_pol
			        rhoTheta2rz(self, rho, theta_in, t_in=None, coord_in='rho_pol', n_line=201)
			        x = eqm.rhoTheta2rz(rho_pol, 0, t_in=times[ind[i]], coord_in='rho_pol')
			        #get large and small radius
			        R = x[0][0][:][0]
			        Z = x[1][0][:][0]
			        r = R - R[0]
			        #---------------------------------------------------------------------------
                        
                        elif area == 'rho_tor':
                	
                                #IF RHO TOROIDAL AS AREA
			        #---------------------------------------------------------------------------
			        #calculate rho_tor from rho_pol
			        rho2rho(self, rho_in, t_in=None, coord_in='rho_pol', coord_out='rho_tor', extrapolate=False) 
			        r = eqm.rho2rho(rho_pol, t_in=times[ind[i]], coord_in='rho_pol', coord_out='rho_tor')
			
                        else:
                                r = rho_pol
                                r = np.tile(r, (1,1))

                        data = np.tile(data, (1,1))
			#build profile: 2 column structure of area+data
			prof = np.transpose(np.concatenate((r, data))) #use r, rho_pol or rho_tor from above

			#build filename
			fname = str(shot) + '.' + tname + '_' + quant[q].name + '_' + typ[q] + '_' + area

			#save file as txt with ending dat
			np.savetxt(dirname + '/' + fname + '.dat', prof)

	return True



