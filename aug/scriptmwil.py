#############################################################################
# THIS SCRIPT USES FUNCTIONS equiwls.py USING W. SUTTROPS PY LIBTO COLLEC   #
# DATA FOR GIVEN SHOTS AND TIMES                                            #
#############################################################################
# AUTHOR:  PHILIPP ULBL                                                     #
# CREATED: 07.04.2020                                                       #
#############################################################################

import sys
import os

from coilmwil import coilmwil

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

#go through all shots
for r in range(len(run)):

    shot = run[r][0]
    time = run[r][1]

    print('(' + str(r+1) + '/' + str(len(run)) + ') processing shot ' + str(shot) + ' @' + str(time) + 's')

    coilmwil(shot, time, PSL=True)
    coilmwil(shot, time, PSL=False)


