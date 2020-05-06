#############################################################################
# THIS SCRIPT USES FUNCTIONS coil.py, profiles.py & equi.py TO COLLECT      #
# DATA FOR GIVEN SHOTS AND TIMES (SHOT+TIME PAIRS)	 					    				#
#############################################################################
# AUTHOR:  PHILIPP ULBL							    						#
# CREATED: 25.03.2020							    						#
#############################################################################

from coil import coil
from profiles import profiles
from equi import equi
from elmdiag import elmdiag
from power import power

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
#type of equilibrium
eqtype = 'EQH'
#experiment for PED
ex = 'ULBLP'

#go through all shots
for r in range(len(run)):

    shot = run[r][0]
    time = run[r][1]

    power(shot);
    elmdiag(shot);

    #save coilfile, profiles, shots
    coil(shot, time)
    profiles(shot, time, exper=ex)
    equi(shot, time, eqtype)
