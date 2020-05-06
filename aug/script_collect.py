#############################################################################
# THIS SCRIPT USES FUNCTIONS coil.py, profiles.py & equi.py TO COLLECT      #
# DATA FOR GIVEN SHOTS AND TIMES	 					    				#
#############################################################################
# AUTHOR:  PHILIPP ULBL							    						#
# CREATED: 17.12.2019							    						#
#############################################################################

from coil import coil
from profiles import profiles
from equi import equi
from elmdiag import elmdiag
from power import power

#put shotnumbers in list here
shot = [35712]
#put times in list here
time = [2.675, 4.55]
#flag if 2-element time should be interpreted as tmin and tmax. no effect if len(time) not 2
window = False
#type of equilibrium
eqtype = 'EQH'
#experiment for PED
ex = 'ULBLP'

#go through all shots
for s in range(len(shot)):

	#get times in window
	if window == True and len(time) == 2:
			#save coilfile, profiles, shots
			coil(shot[s], time)
			profiles(shot[s], time, exper=ex)
			equi(shot[s], time)
	#get specific times
	else:
		#go through all times
		for t in range(len(time)):
		
			#save coilfile, profiles, shots
			coil(shot[s], time[t])
			profiles(shot[s], time[t], exper=ex)
			equi(shot[s], time[t], eqtype)
                        #here one can include wls and mwil data (equi+coil)
        
        #get power and elmdiag data
        power(shot[s])
        elmdiag(shot[s])
