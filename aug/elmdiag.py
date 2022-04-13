#############################################################################
# THIS FUNCTION EXTRACTS ELMDIAG DATA FROM SHOTFILES AND SAVES THEM            #
#############################################################################
# AUTHOR:  PHILIPP ULBL							    						#
# CREATED: 20.03.2020						   							#
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

def elmdiag(shot, leng=1e4, path='', debug=True):
    """
    function elmdiag.py
    -------------------------------------------------------
    def elmdiag(shot, leng=1e4, path='', debug=True):

    writes 3 files in path which contain ascii data of the 
    divertor current (MAC, Ipolsoli), Dalpha radiation
    (POT, ELMi-Han) and the current of 1 RMP coil (Bl6).
    if path is empty a directory coil will be created in home.
    leng sets the length of the time traces (which are very
    huge by default)

    """

    #overwrite path if it is empty
    if path == '':
        path = os.path.expanduser('~') + '/pyout/elmdiag/'


#############################################################################
# IMPORT DATA
#############################################################################

    quant = list()
    name = list()

    try:
        mac = dd.shotfile('MAC', shot, experiment='AUGD') # instantiate and open shotfile
        quant.append(mac('Ipolsoli'))
        name.append('MAC_Ipolsoli')
        mac.close()					  # close shotfile (optional)
    except:
        print('\nINFO: no MAC file found, skipped divertor current data.')

    try:
        pot = dd.shotfile('POT', shot, experiment='AUGD') # instantiate and open shotfile
        quant.append(pot('ELMi-Han'))
        name.append('POT_ELMi-Han')
        pot.close()					  # close shotfile (optional)
    except:
        print('\nINFO: no POT file found, skipped Dalpha radiation data.')

    try:
        maw = dd.shotfile('MAW', shot, experiment='AUGD') # instantiate and open shotfile
        quant.append(maw('IBl6'))
        name.append('MAW_IBl6')
        maw.close()					  # close shotfile (optional)
    except:
        print('\nINFO: no MAW file found, skipped coil current data.')



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
        ind = range(0,len(quant[i].time),int(delta))

        #concat arrays along columns
        prof = np.transpose(np.vstack((quant[i].time, quant[i].data)))

        #build filename
        fname = str(shot) + '_' + name[i]

        #save file as txt with ending dat
        np.savetxt(dirname + '/' + fname + '.dat', prof[ind,:])

    return True






