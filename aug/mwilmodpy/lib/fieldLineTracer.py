import dd
import Fieldline
import numpy
import kk
from IPython import embed
import EQU

class fieldLineTracer(object):
   
    ## load equilibrium
    def __init__(self, pulseNumber, diag='EQH', exp='AUGD'):
        object.__init__(self)
	self.pulseNumber = pulseNumber
        self.EQ = EQU.EQU()
        self.EQ.Load(pulseNumber,exp, diag)
        sf = dd.shotfile(diag, pulseNumber,exp)
	R = self.EQ.R
	z =  self.EQ.z
        self.time =  self.EQ.time 
        self.NR = self.EQ.NR
	self.Nz = self.EQ.Nz
        self.zmax = float(z.max())
	self.zmin = float(z.min())
	self.Rmax = float(R.max())
	self.Rmin = float(R[R!=0.0].min())
	self.R0 = 1.65
	
        self.Ntime = self.EQ.Ntime
        self.psi = self.EQ.Psi
	del sf
        try:
            self.Btor = dd.shotfile('MAI', pulseNumber)('BTFABB')
        except:
            self.Btor = dd.shotfile('MAI', pulseNumber)('BTF')*1.01

    def getPsi(self, time):

        psi = self.EQ.getPsi(time)
	return Fieldline.axiSymmetric.magneticField(float(numpy.interp(time, self.Btor.time, self.Btor.data)), self.R0, self.Rmin, self.Rmax, self.zmin, self.zmax, self.NR, self.Nz, psi.flatten().tolist())


    # start coordinate R,z,phi tor in degree, kind 'rhopanglephi'
    def getFieldline(self, time ,start=[2.0,0.0,0.0], coord_in='Rzphi'):
        
        Psi= self.getPsi(time)

        ## calculate rz from the data
        if  coord_in == 'rhopanglephi':
            output = kk.KK().kkrhorz( self.pulseNumber, time, [start[0]],angle=start[1] )
            R = numpy.squeeze(output.r)
            z = numpy.squeeze(output.z)
            phi = start[2]/180.*float(numpy.pi)
        else:
            R = start[0]
            z = start[1]
            phi = start[2]/180.*float(numpy.pi)

        stopCriterion = Fieldline.trace.stopCriterionSteps(1300)
        #embed()
        print 'Start values: ',start, R,z,phi
   #     if z < 0.0:
    #        direc=1
    #    else:
    #        direc=2
       ### embed()
        
        fieldline1 = Fieldline.trace.fieldline(Psi, float(R),float(z),float(phi),0.01,stopCriterion, 2)
        fieldline2 = Fieldline.trace.fieldline(Psi, float(R),float(z),float(phi),0.01,stopCriterion, 1)

        fieldline = Fieldline.trace.fieldline(Psi, float(R),float(z),float(phi),0.01,stopCriterion, 3)

        x = numpy.hstack((numpy.array(fieldline1.x)[::-1],numpy.array(fieldline2.x)))
        y = numpy.hstack((numpy.array(fieldline1.y)[::-1],numpy.array(fieldline2.y)))
        z = numpy.hstack((numpy.array(fieldline1.z)[::-1],numpy.array(fieldline2.z)))
        R = numpy.hstack((numpy.array(fieldline1.R)[::-1],numpy.array(fieldline2.R)))
        phi = numpy.hstack((numpy.array(fieldline1.phi)[::-1],numpy.array(fieldline2.phi)))

        #return x, y, z, R, phi
        return fieldline2.x,fieldline2.y,fieldline2.z,fieldline2.R,fieldline2.phi
        #return phi and the toroidal angle
 
