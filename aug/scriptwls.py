#############################################################################
# THIS SCRIPT USES FUNCTIONS equiwls.py USING W. SUTTROPS PY LIBTO COLLEC   #
# DATA FOR GIVEN SHOTS AND TIMES                                            #
#############################################################################
# AUTHOR:  PHILIPP ULBL                                                     #
# CREATED: 07.04.2020                                                       #
#############################################################################

import sys
import os

from equiwls import equiwls

#put shotnumber and time pairs in list here
run = [[30835, 2.300],
       [30835, 3.000],
       [33120, 5.500],
       [33120, 5.635],
       [33133, 2.600],
       [33133, 3.000],
       [33353, 2.325],
       [33353, 2.670],
       [33353, 2.900]]

#possible types of equilibrium
conf = [['AUGD', 'IDE'],
        ['MICDU', 'EQB']]
       #['AUGD', 'EQH']

#editions
ed = range(4)

#go through all shots
for r in range(len(run)):
    
    shot = run[r][0]
    time = run[r][1]

    print('(' + str(r+1) + '/' + str(len(run)) + ') processing shot ' + str(shot) + ' @' + str(time) + 's')
    
    #go through all experiment, equilibria configs
    for d in range(len(conf)):

        #go through all versions
        for e in ed:

            #try to find equi files
            try:
                #switch off internal prints
                sys.stdout = open(os.devnull, 'w')
                equiwls(shot, time, expr = conf[d][0], diag = conf[d][1], ed = e)
                #switch on again prints
                sys.stdout = sys.__stdout__
                print('(x) equi found: ' + conf[d][0] + ', ' + conf[d][1] + ', ed=' + str(e))
            
            #runtime error occurs when nothing is found, skip this config
            #other errors will still get through which is nice to detect errors inside the functions
            except RuntimeError:
                #switch on again prints
                sys.stdout = sys.__stdout__
                print('(o) equi not found: '  + conf[d][0] + ', ' + conf[d][1] + ', ed=' + str(e))
                #pass



