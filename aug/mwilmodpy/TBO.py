
from IPython import embed
import numpy as np

import dd as dd

class TBOhelp:
    status = False


class TBO:
    def __init__( self , Shotnumber = None ):
        self.Status = False

        if Shotnumber != None :
            self.Load( Shotnumber )

    def __del__( self ):
        self.Unload( )
        del self.Status

#ETBO.Load(30839,Exoeriment='ECEI',Diagnostic='TBO')
    def Load( self ,  Shotnumber, Experiment='AUGD', Diagnostic='TBO', Edition = 0L, tBegin=0.0, tEnd=11.0):
        self.Unload()

        if Diagnostic == 'TBO':
            try:
                sf = dd.shotfile( Diagnostic, Shotnumber, Experiment, Edition)
                self.Shotnumber = Shotnumber
            except:
                print "Error reading shotfile" 
                return False

            try:
                self.tBegin = tBegin
                self.tEnd = tEnd

                $xr = sf.getSignalGroup( 'xr', tBegin=tBegin, tEnd=tEnd )

            except:
                print "Error reading ECI data"
                embed()
                sf.close()
                self.Status = False
                return False


    def Unload( self ):
        if self.Status:
            self.Status = False

#        zr	Signal-Group	z position of the absorption (m)#
#	t	Time-Base	Timebase for all signals (s)
#	yr	Signal-Group	toroidal coordinate of the absorption (m)
#	xr	Signal-Group	xr (R=sqrt(xr^2+yr^2)), x,y,z is a cartesian system (m)
#	xrmin	Signal-Group	x coordinate of the absorption path along the beam path (m)
#	yrmin	Signal-Group	toroidal coordinate of the absorption path (m)
#	zrmin	Signal-Group	z coordinate of the absorption path (m)
#	xrmax	Signal-Group	x coordinate of the absorption path (m)
#	yrmax	Signal-Group	toroidal coordinate of the absorption path (m)
#	zrmax	Signal-Group	z coordinate of the absorption path (m)
#	abspow	Signal-Group	absorbed power MW
#	rho_abs
