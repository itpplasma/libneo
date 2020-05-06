import scipy.interpolate
import math
import scipy.signal
import scipy
import numpy
import matplotlib.pylab as plt
import IPython
#import kk_abock as kk2
import kk_mwillens as kk
            


class particleTransport:


    def __init__( self, shot = None, time = None ):
      
        self.status = False
        self.status_gs = False
        self.status_Vder = False

        if (shot != None) & (time != None) :
            self.Load( shot, time)


    def Load( self, shot, time):
        
        if (shot != None) & (time != None) :
            self.shot = shot
            self.time = time
            self.status = True
        else:
            print 'no shotnumber and time'
    
        
    def Unload( self ):

        if self.status:

            del self.shot
            del self.time

        if self.status_gs:

            del self.g1
            del self.g2

        self.status = False
        self.status_gs = False
        self.status_Vder = False



# calculate grad of rho in R
    def calcGs(self, rhop, exp = 'AUGD', diag = 'EQH'):

        if (rhop == None):
            return False
# for the average
        angle = numpy.arange(1,360,1.)
        r = numpy.zeros( (numpy.size(angle),numpy.size(rhop)) )
        z = numpy.zeros_like(r)
        R = numpy.zeros_like(r)
        #diffR = numpy.zeros((numpy.size(angle),numpy.size(rho)-1) )
        gradR = numpy.zeros((numpy.size(angle),numpy.size(rhop)) )
            
#calculate the difference in rho    
        #diffRho = numpy.diff(rho)
        gradRho = numpy.gradient(rhop)
        try:
            out = kk.KK().kkrhorz(self.shot, self.time, [0.0], angle=0.0, exp='AUGD', diag='EQH', ed=0)
            r0 = out.r
            z0 = out.z
#loop for every angle because kk eo
            #IPython.embed()

            for i in range(numpy.size(angle)):
                #       if self.status_rhot:
                out = kk.KK().kkrhorz(self.shot, self.time, rhop, angle=angle[i], exp='AUGD', diag='EQH', ed=0)
                r[i,:] = r0 - out.r
                z[i,:] = z0 - out.z
                R[i,:] = numpy.sqrt(r[i,:]*r[i,:]+ z[i,:]* z[i,:])
                gradR[i,:] = numpy.gradient( R[i,:] )
            
            gradRho = gradR/gradRho
            self.g1 = numpy.average(gradRho,axis=0)
            self.g2 = numpy.average(gradRho**2.,axis=0)

            self.status_gs = True
            return True

        except:

            return False



    def calcVder(self, rhop, exp = 'AUGD', diag = 'EQH'):

        if (rhop == None):
            return False

        try:
            out1 = kk.KK().kkrhorz(self.shot, self.time, rhop, angle=0.0, exp='AUGD', diag='EQH', ed=0)
            out2 = kk.KK().kkrzptfn(self.shot, self.time, out1.r, out1.z, exp='AUGD', diag='EQH', ed=0)
            out3 = kk.KK().kkpfv(self.shot, self.time,  out2.pf, exp='AUGD', diag='EQH', ed=0)
            
            self.dVdRho = numpy.gradient(out3.v)/numpy.gradient(rhop)

            self.status_Vder = True
            print 
            return True

        except:

            return False


    def nonHarmonicAnalysis(self, time, rhop, ne):

        if (rhop == None) | (ne == None):
            return False

        meanRhop = numpy.mean( rhop,axis=0 )

        if not self.status_Vder:
            if not self.calcVder( meanRhop ):
                print 'Something went wrong in calculating the Volume derivative'

        if not self.status_gs:
            if not self.calcGs( meanRhop ):
                print 'Something went wrong in calculating the Volume derivative'
           
     
        meanne = numpy.mean( ne,axis=0 )
        nChannel = numpy.size(meanne)
        ntime = numpy.size(time)
        

        g1 =  numpy.reshape(numpy.tile( self.g1,(ntime) ) ,(ntime,nChannel) )
        g2 =  numpy.reshape(numpy.tile( self.g2,(ntime) ) ,(ntime,nChannel) )
        dVdRho = numpy.reshape(numpy.tile( self.dVdRho,(ntime) ) ,(ntime,nChannel) )
        t = numpy.transpose( numpy.reshape(numpy.tile( time ,(nChannel) ) ,(nChannel,ntime)) )
        

        x0 = numpy.gradient( meanne ) / numpy.gradient( numpy.mean( rhop,axis=0 ) ) / meanne
        x = numpy.gradient( ne )[1] / numpy.gradient( rhop )[1] / ne
        y = 1.0/dVdRho/g2/ne * scipy.integrate.cumtrapz( rhop, numpy.gradient(dVdRho*ne)[0]/numpy.gradient(t)[0], initial = 0.0, axis = 1 )


        D1 = numpy.zeros_like(meanne)
        D2 = numpy.zeros_like(meanne)
        v1 = numpy.zeros_like(meanne)
        v2 = numpy.zeros_like(meanne)
        
        
        for i in range(nChannel):

            A1 = numpy.vstack([x[:,i]-x0[i], numpy.ones(ntime)]).transpose()
            A2 = numpy.vstack([x[:,i]-x0[i], numpy.zeros(ntime)]).transpose()
            D1[i], v1[i] = numpy.linalg.lstsq(A1, y[:,i])[0]
            D2[i], v2[i] = numpy.linalg.lstsq(A2, y[:,i])[0]
            v2[i] = x0[i]*D2[i]*g2[0,i]/g1[0,i]

        IPython.embed()


    #    y = 1.0/self.dVdRho/self.g2/ne * scipy.integrate.cumtrapz( rhop, numpy.gradient(self.dVdRho*ne)[0]/, initial = 0.0, axis = 1 )
        
        #x = 

"""    def calculateHarmonicDv(self, amp, phase, rhop = None, exp = 'AUGD', diag = 'EQH'):


        if (numpy.size(amp) == numpy.size(phase) & numpy.size(phase) == numpy.size(rhop)):
            
#limit the rho to 1.0  
            if (rhop != None):
                idx = numpy.where( rhop <= 1.0 ) 
                self.rho = rhop[idx]
            if (rhot != None):
                idx = numpy.where( rhot <= 1.0 ) 
                    
            self.rho = rhot[idx] 
            self.status_rhot = True

            if (shot != None) & (freq != None) & (time != None) :
                self.shot = shot
                self.freq = freq
                self.time = time
            else:
                print "shotnumber or/and frequency is missing."
                return False
                        
            self.amp = amp[idx]
            self.phase = phase[idx]          
                        
        else:
        
            print "Data for amplitude or phase or rhop have not the same size"
            return False


        self.D = numpy.zeros_like(self.amp)
        self.v = numpy.zeros_like(self.amp)
           
#### Calculate the flux surface average rho gradient
    
        if self.rho[0] == 0.0 :
            amp = self.amp
            phase =self.phase
            rho = self.rho
            newRho = numpy.diff(self.rho)/2.+self.rho[:-1]
        else:
            amp = numpy.append(self.amp[0],self.amp)
            phase = numpy.append(self.phase[0]+0.1,self.phase)
            rho = numpy.append([0.0],self.rho )
            newRho = numpy.diff(numpy.append([0.0],self.rho))/2.+numpy.append([0.0],self.rho[:-1]) 

            ##for amplitudes and phases inbetween
        f = scipy.interpolate.interp1d( rho, amp )
        newAmp = f( newRho )			
        g = scipy.interpolate.interp1d( rho, phase )
        newPhase = g( newRho ) 
        


            out1 = kk.KK().kkrhorz(self.shot, self.time, rho, angle=0.0, exp='AUGD', diag='EQH', ed=0)
            out2 = kk.KK().kkrzptfn(self.shot, self.time, out1.r, out1.z, exp='AUGD', diag='EQH', ed=0)

  
            out3 = kk.KK().kkpfv(self.shot, self.time,  out2.pf, exp='AUGD', diag='EQH', ed=0)
            
           # IPython.embed()
            
            diffV = numpy.diff(out3.v)
            
            dVdRho = diffV / diffRho
            dVdRho2 = numpy.gradient(out3.v)/numpy.gradient(rho)
            #dVdRho2[0] = 0.0 
            
            gradAmp = numpy.diff(amp)  /  numpy.diff(rho)
            gradPhase = numpy.diff(phase) /  numpy.diff(rho) 
           
            w = 2.*numpy.pi*self.freq
          #  dVRho
  

            X = scipy.integrate.cumtrapz(rho,dVdRho2*amp*numpy.cos( phase ))
            Y = scipy.integrate.cumtrapz(rho,dVdRho2*amp*numpy.sin( phase ))


            J=g2*dVdRho

            D1 = -w*( Y * numpy.sin(newPhase)  + X * numpy.cos(newPhase) ) / (J*gradPhase*newAmp)
            v1 = -w*( numpy.sin(newPhase) * (gradAmp*Y - gradPhase* newAmp*X) + numpy.cos(newPhase) * (gradAmp*X + gradPhase*newAmp*Y) ) / (J*gradPhase*newAmp**2.)*g1

            self.D1 = D1
            self.v1 = v1
            self.newRho = newRho
            self.g2 = g2
            self.g1 = g1




# use formula from Emiliano,PhD page 57, sign of N is changed
            gradAmp2 = numpy.gradient(amp)  /  numpy.gradient(rho)
            gradAmp2[0] = 0.0
            gradPhase2 = numpy.gradient(phase) /  numpy.gradient(rho)
            gradPhase2[0] = 0.0

            gradRho = gradR/numpy.gradient(rho)
            
    #        g1 = numpy.average(gradRho,axis=0)       
    #        g2 = numpy.average(gradRho**2.,axis=0)
            
            J=g2*dVdRho
            X = scipy.integrate.cumtrapz(rho,dVdRho2*amp*numpy.cos( phase ) )
            #IPython.embed()
            M11 = gradAmp*numpy.cos(newPhase) - newAmp*gradPhase*numpy.sin(newPhase)
            #M11 = gradAmp2*numpy.cos(phase) - amp*gradPhase2*numpy.sin(phase)

            M21 = -gradAmp*numpy.sin(newPhase) - newAmp*gradPhase*numpy.cos(newPhase)
            #M21 = -gradAmp2*numpy.sin(phase) - amp*gradPhase2*numpy.cos(phase)

            M12 = -g1/g2*newAmp*numpy.cos(newPhase)
            #M12 = -g1/g2*amp*numpy.cos(phase)
            M22 = g1/g2*newAmp*numpy.sin(newPhase)
            #M22 = g1/g2*amp*numpy.sin(phase)

            M=numpy.array([[M11,M12],[M21,M22]])

            N1 = scipy.integrate.cumtrapz(rho,dVdRho2*amp*numpy.sin( phase ) )
            N2 = scipy.integrate.cumtrapz(rho,dVdRho2*amp*numpy.cos( phase ) )
            
            N=numpy.array([-w*N1/dVdRho/g2,-w*N2/dVdRho/g2])
            
            D2 = numpy.zeros_like(M11)
            v2 = numpy.zeros_like(M11)
            D3 = numpy.zeros_like(M11)
            v3 = numpy.zeros_like(M11)
            D4 = numpy.zeros_like(M11)
            v4 = numpy.zeros_like(M11)
            D5 = numpy.zeros_like(M11)
            v5 = numpy.zeros_like(M11)
            for i in range(numpy.size(M11)):
                try:
                    D2[i],v2[i] = numpy.linalg.solve(M[:,:,i], N[:,i])
                    D3[i],v3[i] = numpy.linalg.solve(M[:,:,i].transpose(), N[:,i])
                    D4[i],v4[i] = numpy.linalg.solve(M[:,:,i], -N[:,i])
                    D5[i],v5[i] = numpy.linalg.solve(M[:,:,i].transpose(), -N[:,i])
                except:
                    D2[i] = 0.0
                    v2[i] = 0.0

            self.D2 = D2
            self.v2 = v2
            IPython.embed()
            #v3 =  ( w*X + D2 * amp * J * gradPhase2 * numpy.cos( phase ) + D2 * J * gradAmp2 * numpy.sin(phase) ) / ( amp * J *  numpy.sin(phase) )
        #    ( w*X + D.*AMP.*J.*dPHA.*cos(PHA) + D.*J.*dAMP.*sin(PHA) ) ./ ...
#  (AMP.*J.*sin(PHA));

           # IPython.embed()
                                                        

#  cos(PHA).* (dAMP.*X+dPHA.*AMP.*Y) ) ./ (J.*dPHA.*AMP.^2);


      #  % different does non-equidistant gradient
#dAMP=different(AMP,d.rho);
#dPHA=different(PHA,d.rho);


#%% TAKENAGA FORMULA with GENERALISED GEOMETRY (see A. Salmi EPS'14)
#w=2*pi*d.freq;

#X=cumtrapz(d.rho,d.dVdrho.*AMP.*cos(PHA));
#Y=cumtrapz(d.rho,d.dVdrho.*AMP.*sin(PHA));

#J=d.avegradrho2.*d.dVdrho;

#D=-w*( Y.*sin(PHA)  + X.*cos(PHA) ) ./ (J.*dPHA.*AMP);
#V=-w*( ...
#  sin(PHA).* (dAMP.*Y-dPHA.*AMP.*X) + ...
#  cos(PHA).* (dAMP.*X+dPHA.*AMP.*Y) ) ./ (J.*dPHA.*AMP.^2);

#% another way of getting V when D is known
#V2= ( w*X + D.*AMP.*J.*dPHA.*cos(PHA) + D.*J.*dAMP.*sin(PHA) ) ./ ...
#  (AMP.*J.*sin(PHA));


   
"""
