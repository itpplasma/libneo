#############################################################################
# THIS FUNCTION EXTRACTS POWER DATA FROM SHOTFILES AND SAVES THEM            #
#############################################################################
# AUTHOR:  PHILIPP ULBL							    						#
# CREATED: 30.03.2020						   							#
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

def power(shot, leng=1e4, path='', debug=True):
    """
    function power.py
    -------------------------------------------------------
    def power(shot, leng=1e4, path='', debug=True):

    writes 3 files in path which contain ascii data of
    radiated power in the main plasma and the total input power
    minus NBI shine through minus dw/dt.
    if path is empty a directory coil will be created in home.
    leng sets the length of the time traces (which are very
    huge by default)

    """

    #overwrite path if it is empty
    if path == '':
        path = os.path.expanduser('~') + '/pyout/power/'


#############################################################################
# IMPORT DATA
#############################################################################

    quant = list()
    name = list()

    try:
        bpd = dd.shotfile('BPD', shot, experiment='AUGD') # instantiate and open shotfile
        quant.append(bpd('Prad'))
        name.append('BPD_Prad')
        bpd.close()					  # close shotfile (optional)
    except:
        print('\nINFO: no BPD file found, skipped radiated power data.')

    try:
        tot = dd.shotfile('TOT', shot, experiment='AUGD') # instantiate and open shotfile
        quant.append(tot('P_TOT'))
        name.append('TOT_P_TOT')
        tot.close()					  # close shotfile (optional)
    except:
        print('\nINFO: no TOT file found, skipped total input power data.')

#############################################################################
# EXTRACT AND WRITE 2 FILE
#############################################################################

    #create directory if not existent
    dirname = path + str(shot)
    if not os.path.exists(dirname):
        os.makedirs(dirname, 0755)

    for i in range(0, len(quant)):
        
        print('('+str(i+1)+'/'+str(len(quant)) + ') writing ' + name[i])
        
        #get indices
        delta = np.floor(len(quant[i].time) / leng)
        
        if delta>0:
            ind = range(0,len(quant[i].time),int(delta))
        else:
            ind = range(0,len(quant[i].time))

        #concat arrays along columns
        prof = np.transpose(np.vstack((quant[i].time, quant[i].data)))

        #build filename
        fname = str(shot) + '_' + name[i]

        #save file as txt with ending dat
        np.savetxt(dirname + '/' + fname + '.dat', prof[ind,:])

    return True






