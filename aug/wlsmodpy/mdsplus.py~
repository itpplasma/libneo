#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mdsplus.py
- MDSplus connection with some fallbacks

Created on Tue Mar 31 16:47:25 2020

@author: wls
"""

import os,sys

# find path to mdsplus (client) software
pcand = []
if 'MDSPLUS_DIR' in os.environ:     # specified
    pcand.append(os.environ['MDSPLUS_DIR'])
# any hint in python path
for p in sys.path:
    i = p.find('mdsplus')
    if i>=0: 
        pcand.append(p[:i+7])
# and popular locations
pcand.append('/usr/local/mdsplus')
pcand.append('/afs/ipp-garching.mpg.de/tok/soft/mdsplus/7.84.8/amd64_sles15')

have_mdsplus = False  
for p in pcand:
    # print('Checking MDSplus path %s' % p)
    sys.path.append(p+'/python')
    ld_lib_path = p+'/lib'
    if 'LD_LIBRARY_PATH' in os.environ:
        if ld_lib_path not in os.environ['LD_LIBRARY_PATH']:
            os.environ['LD_LIBRARY_PATH'] \
            = os.environ['LD_LIBRARY_PATH'] + ';' + ld_lib_path
    else:
        os.environ['LD_LIBRARY_PATH'] = ld_lib_path
    try:
        import MDSplus
    except:
        continue
    have_mdsplus = True
    print('Found MDSplus in: %s' % p)
    break

import numpy as np

# =========================================================================

class mdsplus(MDSplus.Connection):
    '''
      Replicate MDSplus.Connection class
      - use the mechanism above to find MDSplus module
      - fall back to alternative 'althostspec' if no connection can be made
      - allow all sorts of useful overloads
    '''

    def __init__(self, hostspec, althostspec='localhost'):
        try:
          super(mdsplus, self).__init__(hostspec)
          return
        except:
          print('No connection to host %s, trying %s' \
                % (hostspec, althostspec))
        self.socket = 0
        super(mdsplus, self).__init__(althostspec)


# =========================================================================

class mdsaug(mdsplus):
    '''
      Connection to AUG MDSplus server
      - can do everything an MDSplus.Connection can do
      - plus methods to read data from AUG archive
    '''
    
    def augsignal(self, shot, diag, signame,  \
                  experiment='AUGD', edition=0, t1=-9999.9, t2=9999.9, \
                  use_augsignal=False, raw=False):
        if use_augsignal:
            cmd = 'augsignal'
        else:
            cmd = 'augdiag'
        c = '_oshot=0, _oedition=0, _s = %s(%d, "%s", "%s", "%s", %d, %f, %f, _oshot, _oedition' \
            % (cmd, shot, diag, signame, experiment, edition, t1, t2)
        if raw:
            c += ', "raw"'
        c += ')'            
        # print("Command = >%s<" % c)
        timebase=None; areabase=None
        values = self.get(c)
        if values.ndim>0:
            timebase = self.get('dim_of(_s, 0)')
            if values.ndim>1:
                areabase = self.get('dim_of(_s, 1)')
                values = np.transpose(values)
                areabase = np.transpose(areabase)
        unit = self.get('units_of(_s)')
        oshot = self.get('_oshot')
        oedition= self.get('_oedition')
        return values, timebase, areabase, unit, oshot, oedition


# =================== test main ========================
if __name__ == '__main__':
    import matplotlib.pyplot as plt
    c = mdsaug('mdsplus')
    v,t,a,u,s,e = c.augsignal(37572, "CEZ", "vrot", t1=3.0, t2=3.2)
    print('Unit: %s  Actual shot: %d, Actual edition: %d' % (u,s,e))
    plt.plot(t, v)
    plt.show()


