#############################################################################
# THIS SCRIPT USES FUNCTION elmdiag.py TO COLLECT DATA FOR THE PREVIEW OF   #
# GIVEN SHOTS                                                               
#############################################################################
# AUTHOR:  PHILIPP ULBL							    						#
# CREATED: 20.03.2020							    						#
#############################################################################

from elmdiag import elmdiag

#put shotnumbers in list here
shot=[30648, 30839, 31128, 31131, 32091, 33118, 33345, 33569, 34213, 34548, 34832]

for i in range(0, len(shot)):
    elmdiag(shot[i])
