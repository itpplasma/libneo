import numpy

import matplotlib.pylab as plt
import IPython
import kk_mwillens as kk



class plasmaParameters:

    def __init__( self, shotnumber=None ):

        self.status = False

        if shotnumber != None:
            self.Load(shotnumber)
            


    def Load( self , shotnumber ):
        
        self.Unload()

        if shotnumber != None:
            self.shotnumber = shotnumber
            self.status = True

    def Unload( self ):

        if self.status:
            del self.shotnumber

        self.status = False

    def veperp(self, time, ne, Te, v_tor, v_pol):
    
        if ((self.status) & (ne!=None) & (Te != None) & (rhop != None)& (time != None)):
            if numpy.size(time) == 1:
                if (numpy.shape(ne) == numpy.shape(Te)) & (numpy.size(numpy.shape(Te)) == 1): 

                    outmaj=kk.KK().kkrhorz(self.shotnumber, time,  rhop, angle=0.0, exp='AUGD', diag='EQH', ed=0)
                    Rin = outmaj.r[0]
                    zin = outmaj.z[0]
                    out=kk.KK().kkrzbrzt(self.shotnumber, time, Rin, zin)
                    B = numpy.reshape(numpy.sqrt( out.br*out.br + out.bt*out.bt + out.bz*out.bz),rshape)
                    

                    return veperp

    def nuStar( self, time, ne, Te, rhop, Zeff = 1.5, eps = 1.65/0.5 ):

        if (self.status) & (ne!=None) & (Te != None) & (rhop != None)& (time != None):
#first case only one time
            if numpy.size(time) == 1:
                if (numpy.shape(ne) == numpy.shape(Te)) & (numpy.size(numpy.shape(Te)) == 1):                                                 
                 
                  
                    outmaj=kk.KK().kkrhorz(self.shotnumber, time, [0.0,1.0], angle=0.0, exp='AUGD', diag='EQH', ed=0)
                    R = outmaj.r[0]
                    a = outmaj.r[1]-R
                    eps = a/R
                    
                    q = numpy.abs(kk.KK().kkrhopfq(self.shotnumber, time, rhop).q)
                    lnLambda = 31.3 - numpy.log( numpy.sqrt(ne) / Te  )
                    
                    coll_e = 6.921e-18*(q*R*Zeff*ne*lnLambda)/(Te*Te*eps**(3./2.))

                    return coll_e
                   # IPython.embed()
                       
                else:
                    print 'ne and Te should have the same shape'

            else:
                if (numpy.shape(ne) == numpy.shape(Te)) & (numpy.shape(Te)[0] == numpy.size(time)): 
                    coll_e = numpy.zeros_like(ne)  
                    for i in range(numpy.size(time)):
                        outmaj=kk.KK().kkrhorz(self.shotnumber, time[i], [0.0,1.0], angle=0.0, exp='AUGD', diag='EQH', ed=0)   
                        R = outmaj.r[0]
                        a = outmaj.r[1]-R
                        eps = a/R           
                        q = numpy.abs(kk.KK().kkrhopfq(self.shotnumber, time[i], rhop).q)
                        lnLambda = 31.3 - numpy.log( numpy.sqrt(ne[i,:]) / Te[i,:]  )
                        coll_e[i,:] = 6.921e-18*(q*R*Zeff*ne[i,:]*lnLambda)/(Te[i,:]*Te[i,:]*eps**(3./2.))

                    return coll_e

                else:

                    print 'shapes do not fit'

        else:

            print 'ne and Te not available'

            

            

    def beta( self, time, ne, Te, rhop ):

        if (self.status) & (ne!=None) & (Te != None) & (rhop != None) & (time != None):
#first case only one time
            if numpy.size(time) == 1:
                if (numpy.shape(ne) == numpy.shape(Te)) & (numpy.size(numpy.shape(Te)) == 1):                                                 
                    angle = numpy.arange(1,360,1.)
                    r = numpy.zeros( (numpy.size(angle),numpy.size(rhop)) )
                    z = numpy.zeros_like(r)

                    for i in range(numpy.size(angle)):
                        out = kk.KK().kkrhorz(self.shotnumber, time, rhop, angle=angle[i], exp='AUGD', diag='EQH', ed=0)
                        r[i,:] = out.r
                        z[i,:] = out.z

               #     IPython.embed()
                    Rin =  numpy.ndarray.flatten(r)
                    zin =  numpy.ndarray.flatten(z)
                    rshape = numpy.shape(r)
                    out=kk.KK().kkrzbrzt(self.shotnumber, time, Rin, zin)
                    B = numpy.reshape(numpy.sqrt( out.br*out.br + out.bt*out.bt + out.bz*out.bz),rshape)

                    B_av = numpy.average(B,axis=0)
                   
                    beta = 0.00251*ne*Te/(B_av*B_av)*1.6022*10**-3.*10**-19

                    return beta 
                   # IPython.embed()
                       
                else:

                    print 'ne and Te should have the same shape'

            else:

                print 'only one timepoint is allowed'

        else:

            print 'ne and Te not available'
        

    def rhoStar( self, time, T, rhop ):

        if (self.status) &  (T != None) & (rhop != None) & (time != None):
#first case only one time
            if numpy.size(time) == 1:
                if numpy.size(numpy.shape(T)) == 1:                                                 
                    angle = numpy.arange(1,360,1.)
                    r = numpy.zeros( (numpy.size(angle),numpy.size(rhop)) )
                    z = numpy.zeros_like(r)

                    for i in range(numpy.size(angle)):
                        out = kk.KK().kkrhorz(self.shotnumber, time, rhop, angle=angle[i], exp='AUGD', diag='EQH', ed=0)
                        r[i,:] = out.r
                        z[i,:] = out.z

                    #IPython.embed()
                    Rin =  numpy.ndarray.flatten(r)
                    zin =  numpy.ndarray.flatten(z)
                    rshape = numpy.shape(r)
                    out=kk.KK().kkrzbrzt(self.shotnumber, time, Rin, zin)
                    B = numpy.reshape(numpy.sqrt( out.br*out.br + out.bt*out.bt + out.bz*out.bz),rshape)
                    B_av = numpy.average(B,axis=0)

                    outmaj=kk.KK().kkrhorz(self.shotnumber, time, [0.0,1.0], angle=0.0, exp='AUGD', diag='EQH', ed=0)

                    R = outmaj.r[0]
                    a = outmaj.r[1]-R
                   
                    rhostar = 0.00646*numpy.sqrt(T*10.**-3)/(a*B_av)
                    return rhostar

                else:

                    print 'ne and Te should have the same shape'

            else:

                print 'only one timepoint is allowed'

        else:

            print 'ne and Te not available'
        
def eta( Te, ne, Zeff = 1.5 ):
    if ((numpy.all(Te) != None) &  (numpy.all(ne) != None)):  
        lnLambda = 31.3 - numpy.log( numpy.sqrt(ne) / Te  )
        eta=5.e-5 * lnLambda/Te**(3./2.)
        return eta
