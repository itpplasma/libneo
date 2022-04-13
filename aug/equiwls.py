import sys

sys.path.append('/afs/ipp/home/w/wls/.local/lib/python2.7/site-packages/pyPmds')
sys.path.append('./wlsmodpy/')

import os
import pmds as mds

from eqrzfd import eqrzfd, geqdsk
from eqmdsaug import eqmdsaug

def equiwls(shot, time, expr='AUGD', diag='EQH', ed=1, path='', debug=True):
    """
    function equiwls.py
    -------------------------------------------------------
    def equiwls(shot, time, expr='AUGD', diag='EQH', ed=1, path='', debug=True):
    
    using python library of W. Suttrop:
    writes a gfile for single shot and single time (in s)
    for experiment expr (default = AUGD) and diagnostics
    diag (default = EQH) and edition ed (default = 1)
    file will be located in a subdirectory for shot in path.
    if path is empty a directory will be created in home.
    
    """
    
    home = os.path.expanduser('~')
    #overwrite path if it is empty
    if path == '':
        path = os.path.expanduser('~') + '/pyout/equiwls/'
    
    #create path for shot
    path = path + str(shot) + '/'

    #create dir if not existent
    if not os.path.exists(path):
        os.mkdir(path, 0755)

    os.chdir(home + '/pyscripts/wlsmodpy/')

    #connect to MDSPlus
    mds.mdsconnect('mdsplus.aug.ipp.mpg.de')

    eq = eqmdsaug(shot, time, experiment=expr, diagnostic=diag, edition=ed)
    tname = '{:4.0f}'.format(int(time*1e3))
    fname = 'g' + str(shot) + '.' + tname + '_' + expr + '_' + diag + '_ed' + str(ed)

    eq.write_geqdsk(path + fname)
