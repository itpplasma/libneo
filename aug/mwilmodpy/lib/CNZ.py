import numpy as np

from scipy.interpolate import LSQUnivariateSpline,interp1d,UnivariateSpline


from IPython import embed
import matplotlib.pyplot as plt

import sys
sys.path.append("/afs/ipp-garching.mpg.de/home/m/mwillens/Programme/python/lib/")
import dd
import kk_mwillens as kk
import eqi_map as fastkk
import CXRS



def shiftFuncResi( shift, x,  f, ferror,coeff ):
    return  (f - splev(x-shift,coeff))/ferror



class CNZhelp:
    status = False

class CNZ:
    def __init__( self , Shotnumber = None , Exp_tor='CXRS',Exp_pol='ULP',Ed_tor=1,Ed_pol=1,tBegin=0.0,tEnd=10.0 ,useTor4Er=False):
        self.Status = False
        if Shotnumber != None:
	    self.Load( Shotnumber,Exp_tor=Exp_tor,Exp_pol=Exp_pol,Ed_tor=Ed_tor,Ed_pol=Ed_pol,tBegin=tBegin,tEnd=tEnd,useTor4Er=useTor4Er )

    def Load( self ,  Shotnumber, Exp_tor='CXRS', Exp_pol='ULP', Diagnostic='CNZ', Edition = 0L, eqExp = 'AUGD', eqDiag = 'EQH',Rshift_pol=0.00,Rshift_tor=-0.0025,Ed_tor=1,Ed_pol=1,plotInte=False,tBegin=0.0,tEnd=10.0,useTor4Er=False):


        Pol=CXRS.CXRS()
        Pol.Load(Shotnumber,Experiment=Exp_pol,Diagnostic='CNZ',loadAllRhop=True,Rshift=Rshift_pol, Edition=Ed_pol,tBegin=tBegin,tEnd=tEnd)
        Tor=CXRS.CXRS()
        Tor.Load(Shotnumber,Experiment=Exp_tor,Diagnostic='CNZ',loadAllRhop=True,Rshift=Rshift_tor, Edition=Ed_tor,tBegin=tBegin,tEnd=tEnd)
        idxTorTime = np.ones_like(Tor.time,dtype='bool')
        idxPolTime = np.ones_like(Pol.time,dtype='bool')

        if np.all(Tor.time==Pol.time) == False:
            print('something wrong, timebase is not the same, using pol timebase')
            idxTorTime = np.zeros_like(Tor.time,dtype='bool')
            idxPolTime = np.zeros_like(Pol.time,dtype='bool')
            for i,t in enumerate(Pol.time):
                if (t in Tor.time):
                    idxPolTime[i]=True

            for i,t in enumerate(Tor.time):
                if (t in Pol.time[idxPolTime]):
                    idxTorTime[i]=True
                    
            #embed()
            #return
        
        useIdxTor=Tor.inte.mean(axis=0)>1.e14
        vtor_time=Tor.time[idxTorTime]
        vtor=Tor.vrot[idxTorTime][:,useIdxTor]
        vtor_err=Tor.vrot_unc[idxTorTime][:,useIdxTor]
        vtor_Rall = Tor.Rall[idxTorTime][:,useIdxTor]
        vtor_R = Tor.R[useIdxTor]
        vtor_zall = Tor.zall[idxTorTime][:,useIdxTor]
        vtor_z = Tor.z[useIdxTor]       
        vtor_inte = Tor.inte[idxTorTime][:,useIdxTor]
        vtor_inte_err = Tor.err_inte[idxTorTime][:,useIdxTor]
        vtor_Ti = Tor.Ti[idxTorTime][:,useIdxTor]
        vtor_Ti_err = Tor.Ti_unc[idxTorTime][:,useIdxTor]
        vtor_inte_mean = np.mean(vtor_inte)
        #get the signal with the largest intensity
        vtor_inte_maxCh = np.mean(vtor_inte[:,vtor_inte.argmax(axis=1)])
        vtor_rhop = Tor.rhop[:,useIdxTor]
        vtor_rmaj = Tor.rmaj[idxTorTime][:,useIdxTor]
        print "vtor ",vtor_R
        
        useIdxPol=Pol.inte.mean(axis=0)>1.e14
        vpol=Pol.vrot[idxPolTime][:,useIdxPol]
        vpol_time=Pol.time[idxPolTime]
        vpol_err=Pol.vrot_unc[idxPolTime][:,useIdxPol]
        vpol_Rall = Pol.Rall[idxPolTime][:,useIdxPol]
        vpol_R = Pol.R[useIdxPol]
        vpol_zall = Pol.zall[idxPolTime][:,useIdxPol]
        vpol_z = Pol.z[useIdxPol]
        vpol_Ti = Pol.Ti[idxPolTime][:,useIdxPol]
        vpol_Ti_err = Pol.Ti_unc[idxPolTime][:,useIdxPol]
        vpol_inte = Pol.inte[idxPolTime][:,useIdxPol]
        vpol_inte_err = Pol.err_inte[idxPolTime][:,useIdxPol]
        vpol_inte_mean = np.mean(vpol_inte)
        #get the signal with the largest intensity
        vpol_inte_maxCh = np.mean(vpol_inte[:,vpol_inte.argmax(axis=1)])
        vpol_rhop = Pol.rhop[idxPolTime][:,useIdxPol]
        vpol_rmaj = Pol.rmaj[idxPolTime][:,useIdxPol]
        print "vpol ",vpol_R
        self.time = vpol_time
        #check shift

        #xShift[t],xTraceCov = leastsq(shiftFuncResi,[0.0] ,args=(self.winData[i]['x'][t,useChan],self.winData[i]['data'][t,useChan]/dataMean,np.ones_like( self.winData[i]['x'][t,useChan] ),splCoeffData),full_output=1)[0:2]
        #nots = np.array([1.65,1.75,1.82,2.05])
#IPython.embed()
#LSSpline = LSQUnivariateSpline(fitdata_R[sort_index],fitdata_Te[sort_index],knots )
        #Fit

        #embed()
        outB = kk.KK().kkrzBrzt(Shotnumber,Pol.time,vpol_R,vpol_z)
        pol_Bt = outB.bt
        pol_Bp = np.sqrt(outB.br**2.+outB.bz**2.)
        pol_IBI =np.sqrt(outB.br**2.+outB.bz**2.+ outB.bt**2.)      
        outB = kk.KK().kkrzBrzt(Shotnumber,Tor.time,vtor_R,vtor_z)
        tor_Bt = outB.bt
        tor_Bp = np.sqrt(outB.br**2.+outB.bz**2.)
        tor_IBI =np.sqrt(outB.br**2.+outB.bz**2.+ outB.bt**2.)      

        GQH=dd.shotfile('GQH',Shotnumber)
        Raus = GQH('Raus',tBegin=vtor_time.min(),tEnd=vtor_time.max())
        self.Rsep = np.interp(vtor_time,Raus.time,Raus.data )
        
        inte_Fit=[]
        inte_derFit=[]
        inte_pFit=[]
        inte_pderFit=[]
        Ti_Fit=[]
        Ti_pFit=[]
        Ti_derFit=[]
        Ti_pderFit=[]
        vtor_Fit=[]
        xinte=vpol_rmaj #np.hstack([vpol_rmaj,vtor_rmaj])
        inte=vpol_inte/vpol_inte_maxCh#np.hstack([vpol_inte/vpol_inte_maxCh,vtor_inte/vtor_inte_maxCh])
        inte_err=vpol_inte_err/vpol_inte_maxCh#np.hstack([vpol_inte_err/vpol_inte_maxCh,vtor_inte_err/vtor_inte_maxCh])
        self.rmaj=vpol_rmaj
        self.rhop=vpol_rhop

        self.diaTerm, self.vpolBtor,  self.vpolBtor_err, self.vtorBpol, self.vtorBpol_err = [],[],[],[],[]
            
        self.Er, self.Er_err , self.vExB, self.vExB_err = [],[],[],[]

        self.vdia,  self.vxBoB, self.vpolBtoroB,  self.vtorBpoloB = [],[],[],[]
      

        
        for t in np.arange(vtor_time.size):
            idxInteSort=np.argsort(xinte[t])
            knots4inte=np.linspace(xinte[t].min()+0.003,xinte[t].max()-0.003,4,endpoint=True)
            #inte_Fit.append( UnivariateSpline(xinte[t][idxInteSort] ,inte[t][idxInteSort], np.abs(1./inte_err[t][idxInteSort]),s=30 ))
            
            rmaj_run = np.linspace(vpol_rmaj[t].min()-0.001,vpol_rmaj[t].max()+0.001,500)
            pInte=np.polyfit(xinte[t][idxInteSort],np.log(inte[t][idxInteSort]),3,w= np.abs(1./inte_err[t][idxInteSort]))

            inte_pFit.append( interp1d(rmaj_run,np.exp(np.poly1d(pInte)(rmaj_run)) ))
            #embed()
            inte_pderFit.append( interp1d(rmaj_run,np.gradient(np.exp(np.poly1d(pInte)(rmaj_run)))/np.gradient(rmaj_run) ) )
                              
                              #rmaj_run np.gradient(np.exp(inte_pFit))
            #inte_pderFit.append( np.poly1d(pInte[:-1]*np.arange(pInte.size-1,0,-1)))
            #inte_derFit.append(inte_Fit[-1].derivative())
            #print np.abs(1/inte_err[t][idxSort]) 

            idxVpolSort=np.argsort( vpol_rmaj[t] )
            xvpol =  vpol_rmaj[t,idxVpolSort]
            knots4vpol=np.linspace(xvpol.min()+0.007,xvpol.max()-0.007,3,endpoint=True)
            
            #Ti_Fit.append( LSQUnivariateSpline(xvpol ,vpol_Ti[t][idxVpolSort], knots4vpol ,np.abs(1./vpol_Ti_err[t][idxVpolSort]),s=10. ))
            Ti_Fit.append(UnivariateSpline(xvpol ,vpol_Ti[t][idxVpolSort], np.abs(1./vpol_Ti_err[t][idxVpolSort]),s=30. ))

            pTi=np.polyfit(xvpol ,vpol_Ti[t][idxVpolSort],2,w= np.abs(1./vpol_Ti_err[t][idxVpolSort]) )
            Ti_pFit.append(np.poly1d(pTi))
            Ti_pderFit.append(np.poly1d(pTi[:-1]*np.arange(pTi.size-1,0,-1)))
            Ti_derFit.append(Ti_Fit[-1].derivative())

            idxVtorSort=np.argsort( vtor_rhop[t] )
            xvtor =  vtor_rmaj[t,idxVtorSort]

            #knots4vtor
            knots4vtor=np.linspace(xvtor.min()+0.001,xvtor.max()-0.005,3,endpoint=True)
            #embed()
            vtor_Fit.append(np.poly1d(np.polyfit(xvtor ,vtor[t][idxVtorSort], 1,w=np.abs(1./vtor_err[t][idxVtorSort]) )))
            #knots4vtor =
            #p1 = np.poly1d(np.polyfit(xvtor ,vtor[t][idxVtorSort], 1,w=np.abs(1./vtor_err[t][idxVtorSort]) ))
            #p2 = np.poly1d(np.polyfit(xvtor ,vtor[t][idxVtorSort], 2,w=np.abs(1./vtor_err[t][idxVtorSort]) ))
            #p0 = np.poly1d(np.polyfit(xvtor ,vtor[t][idxVtorSort], 0,w=np.abs(1./vtor_err[t][idxVtorSort]) ))
            if plotInte:
            
                plt.errorbar(vtor_rmaj[t],vtor_inte[t]/vtor_inte_maxCh,vtor_inte_err[t]/vtor_inte_maxCh,color='r',fmt='o')
                plt.errorbar(vpol_rmaj[t],vpol_inte[t]/vpol_inte_maxCh,vpol_inte_err[t]/vpol_inte_maxCh,color='r',fmt='x')
                plt.plot(rmaj_run,inte_pFit[-1](rmaj_run),'r-')
                plt.show()
            
            #plt.errorbar(vpol_rmaj[t],vpol_Ti[t],vpol_Ti_err[t],color='r',fmt='o')
            #plt.plot(xvpol,Ti_Fit[-1](xvpol),'r-')
            #plt.plot(xvpol,Ti_pFit[-1](xvpol),'b-')
            #plt.show()
            #p_arr = np.array([p0(vpol_rmaj[t].min()),p2(vpol_rmaj[t].min()),p1(vpol_rmaj[t].min()),vtor[t][idxVtorSort][0]])
            
            #vtor_Fit.append( UnivariateSpline(np.append(vpol_rmaj[t].min(),xvtor) , np.append(p_arr.mean(),vtor[t][idxVtorSort]), np.append(1./(p_arr.std()),np.abs(1./vtor_err[t][idxVtorSort])),s=20.  ))

            
     #       plt.errorbar(vtor_rmaj[t],vtor[t],vtor_err[t],color='b',fmt='o')
     #       plt.plot(vpol_rmaj[t],p0(vpol_rmaj[t]),'b--')
     #       plt.plot(vpol_rmaj[t],p2(vpol_rmaj[t]),'b-')
     #       plt.plot(vpol_rmaj[t],p1(vpol_rmaj[t]),'b:')
     #       plt.plot(vpol_rmaj[t],vtor_Fit[-1](vpol_rmaj[t]),'r-')


    #plt.xlim([vpol_rmaj[t].min(),vpol_rmaj[t].max()])
    #plt.show()
             
             #Ti_Fit.append( LSQUnivariateSpline(xinte[t][idxInteSort] ,vtor_Ti[t][idxInteSort], knots4inte ,np.abs(1./vtor_Ti_err[t][idxInteSort]) ))

   
            #plt.plot(vpol_rhop[t],inte_Fit[-1](inte_Fit[-1](vpol_rhop[t]) )
            #plt.show()
            
            self.diaTerm.append( inte_pderFit[-1](vpol_rmaj[t])/inte_pFit[-1](vpol_rmaj[t])*Ti_pFit[-1](vpol_rmaj[t])+Ti_pderFit[-1](vpol_rmaj[t]))
            self.vpolBtor.append(-pol_Bt[t]*vpol[t])
            self.vpolBtor_err.append( -pol_Bt[t]*vpol_err[t])
            self.vtorBpol.append(+pol_Bp[t]*vtor_Fit[-1](vpol_rmaj[t]))
            
            if useTor4Er:
                self.Er.append(self.diaTerm[-1]+self.vpolBtor[-1]+self.vtorBpol[-1])
            else:
                self.Er.append(self.diaTerm[-1]+self.vpolBtor[-1])
                
            self.Er_err.append(np.sqrt( self.vpolBtor_err[-1]**2. + self.vtorBpol[-1]**2.) )

            self.vExB.append(self.Er[-1] / pol_IBI[t] )
            self.vExB_err.append(self.Er_err[-1] / pol_IBI[t] )

            self.vdia.append( self.diaTerm[-1] / pol_IBI[t] )
            self.vxBoB.append( (self.vpolBtor[-1] +self.vtorBpol[-1] )/ pol_IBI[t] )
            self.vpolBtoroB.append(  (self.vpolBtor[-1] )/ pol_IBI[t] )
            self.vtorBpoloB.append( (self.vtorBpol[-1] )/ pol_IBI[t] )
            
            #embed()
#            plt.title('%.3fs'%Pol.time[t])
#            plt.plot(vpol_rmaj[t],self.vdia,'r-',label='dia')
#            plt.errorbar(vpol_rmaj[t],self.vxBoB,color='b',fmt='o')
      
#            plt.plot(vpol_rmaj[t], self.vtorBpoloB,'g-',label='vtorBpol')
#            plt.errorbar(vtor_rmaj[t],tor_Bp[t]*vtor[t]/tor_IBI[t],tor_Bp[t]*vtor_err[t]/tor_IBI[t],color='g',fmt='o')

#            plt.errorbar(vpol_rmaj[t],self.vExB,self.vExB_err,color='k',fmt='o')
            
#            plt.xlim([vpol_rmaj[t].min(),vpol_rmaj[t].max()])
#            plt.legend()
#            plt.show()

        self.diaTerm, self.vpolBtor,  self.vpolBtor_err, self.vtorBpol =  np.array(self.diaTerm), np.array(self.vpolBtor),  np.array(self.vpolBtor_err), np.array(self.vtorBpol)
            
        self.Er, self.Er_err , self.vExB, self.vExB_err = np.array(self.Er), np.array(self.Er_err) , np.array(self.vExB), np.array(self.vExB_err)

        self.vdia,  self.vxBoB, self.vpolBtoroB,  self.vtorBpoloB = np.array(self.vdia),  np.array(self.vxBoB), np.array(self.vpolBtoroB),  np.array(self.vtorBpoloB)
      


        #embed()
        #uperpSpline =UnivariateSpline(rhop.data[i][idxSort] ,uperp.data[i][idxSort]/uperpmean ,np.abs(uperpmean/uperp_er.data[i][idxSort]),s=9.,ext=3 )
        #uperpFit=uperpSpline(rhoFit)
