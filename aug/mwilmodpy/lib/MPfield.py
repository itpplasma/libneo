

import kk_mwillens as kk
import scipy.interpolate
from IPython import embed
#import eqi_map as fastkk
import numpy
import matplotlib.pylab as plt

class MPfield:
    status = False


    ## return coordinates coord, rzt is default
def getSurfPert( shot, time, dt = 0.01, q = None, rhop = 1.0, nTorAngle = 180., nPolAngle = 180., coord='ntgt', normalized = True ):

       
    if (rhop == None) & (q == None):
        return None
    else:
            #first read magnetic axis and separatrix
        output = kk.KK().kkrhorz( shot, time, [0.0,1.0] )
        RmagAxis_r = output.r[0]
        RmagAxis_z = output.z[0]
        Sep_r = output.r[1]
        Sep_z = output.z[1]

    if q != None:
        if (numpy.size(q) != 1):
            return None
                     
            #generate
        runningR = numpy.arange( RmagAxis_r, Sep_r,0.0005)
        runningz = numpy.repeat(RmagAxis_z,numpy.size(runningR))
            #embed()
        output = kk.KK().kkrzq( shot, time, runningR, runningz )
        f = scipy.interpolate.interp1d(-output.q,output.rhop, bounds_error=False, fill_value = float('NaN'))	
        print 'redefinition of rhop'
        rhop = f(q).tolist()

    if rhop != None:

        if (numpy.size(rhop) != 1) :
            return None

        if ( rhop > 1.0 ) | ( rhop < 0.0 ):
            return None

            # define running variables
        torAngle = numpy.arange(nTorAngle)*360./float(nTorAngle)
        polAngle = numpy.arange(nPolAngle)*360./float(nPolAngle)-180.

            #get rz of one surface
        print 'get rz of separatrix ', time
        output = kk.KK().kkrhorz( shot, time, [rhop],angle = polAngle)
        surf_r_sym = numpy.squeeze(output.r)
        surf_z_sym = numpy.squeeze(output.z)
            
            #reproduce array
        surf_r = numpy.tile(surf_r_sym,nTorAngle)
        surf_z = numpy.tile(surf_z_sym,nTorAngle)
        surf_torAngle = numpy.repeat(torAngle,nPolAngle)

        #embed()

            #get MP field
        print 'get MP field ', time
        coil = kk.KK().kkrztmprzt( shot, time, dt , surf_r, surf_z, surf_torAngle)

            ## reshape 1st dim torangle, 2nd dim polangle
        MP_r = numpy.reshape(coil.mpr,(nTorAngle,nPolAngle))
        MP_z = numpy.reshape(coil.mpz,(nTorAngle,nPolAngle))
        MP_t = numpy.reshape(coil.mpt,(nTorAngle,nPolAngle))

        #get magnetic field from eqilibrium
        field = kk.KK().kkrzBrzt( shot, time, surf_r_sym, surf_z_sym )
        B_r = numpy.reshape(numpy.tile(field.br,nTorAngle),(nTorAngle,nPolAngle))
        B_z = numpy.reshape(numpy.tile(field.bz,nTorAngle),(nTorAngle,nPolAngle))
        B_t = numpy.reshape(numpy.tile(field.bt,nTorAngle),(nTorAngle,nPolAngle))
            
        dB_r = MP_r+B_r
        dB_z = MP_z+B_z
        dB_t = MP_t+B_t

        IBI = numpy.sqrt( B_r*B_r + B_z*B_z + B_t*B_t  )
        dB = numpy.sqrt(dB_r*dB_r + dB_z*dB_z + dB_t*dB_t )

        x = numpy.reshape(surf_torAngle,(nTorAngle,nPolAngle))
        y = numpy.reshape(numpy.tile(polAngle,nTorAngle),(nTorAngle,nPolAngle))
        #embed()

        if coord == 'rzt':

            if normalized:
                return x,y, dB_r/IBI,dB_z/IBI,dB_t/IBI 
            else:  
                return x,y, dB_r,dB_z,dB_t
            
        if coord == 'abs':

            if normalized:
                return x,y,dB/IBI
            else:
                return x,y,dB

        if coord == 'ntgt':


            cosPhi =numpy.cos(y*numpy.pi/180.)
            sinPhi =numpy.sin(y*numpy.pi/180.)
            # normal to the flux surface
            MP_n = MP_r * cosPhi - MP_z * sinPhi
            #tangential to the flux surface
            MP_tg = MP_r * sinPhi + MP_z * cosPhi
            if normalized:
                return x,y,MP_n/IBI, MP_tg/IBI, MP_t/IBI
            else:
                return x,y,MP_n, MP_tg, MP_t

            #dB_tg = B_r * sinPhi + B_z * cosPhi
            #B_n = B_r * cosPhi - B_z * sinPhi





#poloidal cut at one or averaged
def getTorPert( shot, time, dt = 0.01, q = None, torAngle = 22.5, nR = 100., nz = 50., coord='ntgt',  confinedRegion=False, normalized = True ):

    nAngle = 60
    polAngle = numpy.arange(nAngle)*360./float(nAngle)
    
    out = kk.KK().kkrhorz( shot, time, [1.0],angle = 0.0)
    Rmag = numpy.squeeze(out.r)
    zmag = numpy.squeeze(out.z)

            #get rz of sep
    output = kk.KK().kkrhorz( shot, time, [1.0],angle = polAngle)
    surf_r_sym = numpy.squeeze(output.r)
    surf_z_sym = numpy.squeeze(output.z)
    

    Rcoord =  numpy.linspace( numpy.min(surf_r_sym)-0.05, numpy.max(surf_r_sym)+0.05, num=nR )
    zcoord =  numpy.linspace( numpy.min(surf_z_sym)-0.05, numpy.max(surf_z_sym)+0.05, num=nz )
    
                #reproduce array
    R = numpy.repeat( Rcoord, nz )
    z = numpy.tile( zcoord, nR )

        #get magnetic field from eqilibrium
    field = kk.KK().kkrzBrzt( shot, time, R, z )

    #embed()
    B_r = numpy.reshape(field.br,(nR,nz))
    B_z = numpy.reshape(field.bz,(nR,nz))
    B_t = numpy.reshape(field.bt,(nR,nz))

    #only one slice
    if (numpy.size(torAngle) == 1):

        torAngle =  numpy.repeat( torAngle, nR*nz )
        coil = kk.KK().kkrztmprzt( shot, time,dt, R, z,  torAngle)

            ## reshape 1st dim R, 2nd z and average the toroidal gnle
        MP_r = numpy.reshape(coil.mpr,(nR, nz))
        MP_z = numpy.reshape(coil.mpz,(nR, nz))
        MP_t = numpy.reshape(coil.mpt,(nR, nz))
        
    #for averaging
    else:
        
        torAngle = numpy.arange(360.)
        ntorAngle =  numpy.size(torAngle )
        Rtmp = numpy.tile( R, ntorAngle )
        ztmp = numpy.tile( z, ntorAngle )
        torAngletmp = numpy.repeat( torAngle, nR*nz )
        coil = kk.KK().kkrztmprzt( shot, time,dt, Rtmp, ztmp,  torAngletmp)
        
            ## reshape 1st dim R, 2nd z
        MP_r = numpy.average(numpy.reshape(coil.mpr,(ntorAngle,nR*nz)),axis=0)
        MP_z = numpy.average(numpy.reshape(coil.mpz,(ntorAngle,nR*nz)),axis=0)
        MP_t = numpy.average(numpy.reshape(coil.mpt,(ntorAngle,nR*nz)),axis=0)
     ## reshape 1st dim R, 2nd z
        MP_r = numpy.reshape(MP_r,(nR, nz))
        MP_z = numpy.reshape(MP_z,(nR, nz))
        MP_t = numpy.reshape(MP_t,(nR, nz))

    #embed()

    rhop_out = numpy.reshape(kk.KK().kkrzptfn( shot, time, R, z,).rho_p,(nR, nz))
    R_out = numpy.reshape(R,(nR, nz))
    z_out = numpy.reshape(z,(nR, nz))

    dB_r = MP_r+B_r
    dB_z = MP_z+B_z
    dB_t = MP_t+B_t

    if confinedRegion:
        idx = numpy.where(rhop_out > 1.0 )
        B_r[idx] = 0.0
        B_z[idx] = 0.0
        B_t[idx] = 0.0

        MP_r[idx] = 0.0
        MP_z[idx] = 0.0
        MP_t[idx] = 0.0

    IBI = numpy.sqrt( B_r*B_r + B_z*B_z + B_t*B_t  )
    IdBI = numpy.sqrt(dB_r*dB_r + dB_z*dB_z + dB_t*dB_t )
    IMPI = numpy.sqrt( MP_r*MP_r + MP_z*MP_z + MP_t*MP_t  )

    
## output are the different values
    if coord == 'rzt':

        if normalized:
            return R_out,z_out, dB_r/IBI,dB_z/IBI,dB_t/IBI 
        else:  
            return R_out,z_out, dB_r,dB_z,dB_t
            
    if coord == 'abs':

        if normalized:
            return R_out,z_out,IdBI/IBI
        else:
            return R_out,z_out,IdBI

    if coord == 'ntgt':


        phi = numpy.arctan((zmag-z_out) / (Rmag-R_out))
        cosPhi =numpy.cos(phi)
        sinPhi =numpy.sin(phi)
        # normal to the flux surface
        MP_n = MP_r * cosPhi - MP_z * sinPhi
            #tangential to the flux surface
        MP_tg = MP_r * sinPhi + MP_z * cosPhi

        embed()
        if normalized:
            return R_out,z_out,MP_n/IBI, MP_tg/IBI, MP_t/IBI
        else:
            return R_out,z_out,MP_n, MP_tg, MP_t


    
    
