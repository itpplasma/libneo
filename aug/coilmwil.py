#############################################################################
# THIS FUNCTION EXTRACTS COIL DATA FROM SHOTFILES USING PY LIB OF M WILLENSDORFER #
#############################################################################
# AUTHOR:  PHILIPP ULBL                                                                                                 #
# CREATED: 07.04.2019                                                                                                           #
#############################################################################

#!/usr/bin/env python
# line above sets the Python interpreter

import sys

sys.path.append('/afs/ipp/aug/ads-diags/common/python/lib')
sys.path.append('./mwilmodpy')

import MAWnew as MAW

import matplotlib.pylab as plt                  # load plotting library
import dd                                       # load latest (!) dd library
# use "import dd_YYYYMMDD as dd" to safely load a specific
# version from /afs/ipp/aug/ads-diags/common/python/lib

import numpy as np
import os
import sys

def coilmwil(shot, time, PSL=True, path='', debug=True):
    """
    function coilmwil(shot, time, PSL = true, path='', debug=True):
    -------------------------------------------------------
    
    writes a coil file for a single shot and single time using
    the py lib of M Willensdorfer to include effects of PSL.
    file will be located in a subdirectory for shot in path.
    if path is empty a directory coil will be created in home.
    if time has 2 elements, they are treated as tmin, tmax

    """

    home = os.path.expanduser('~')

    #overwrite path if it is empty
    if path == '':
        path = home + '/pyout/coilmwil/'

    #load instance of MAWnew class
    maw = MAW.MAW(shot)

    #get coil currents with or without PSL effects
    Icoil = maw.getIcoils(time, usePSL=PSL)
    #shape them to 1x16 array
    prof = np.array(np.concatenate((Icoil[0], Icoil[1])), ndmin=2)

    #create directory if not existent
    dirname = path + str(shot)
    if not os.path.exists(dirname):
        os.mkdir(dirname, 0755)
    
    tname = '{:4.0f}'.format(time*1e3)

    #build filename
    fname = str(shot) + '.' + tname + '_coil'
    if PSL == True:
        fname = fname + '_PSL'

    #save file as txt with ending dat
    np.savetxt(dirname + '/' + fname + '.dat', prof)
