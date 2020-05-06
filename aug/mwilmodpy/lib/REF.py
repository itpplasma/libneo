import numpy
import dd
import scipy.interpolate
import kk_mwillens as kk
from IPython import embed
import matplotlib.pylab as plt
import tools
import PHI

class REFhelp:
    status = False

class REF:

    def __init__( self , Shotnumber = None, Experiment = 'AUGD' ):
        
        self.Status = False
        self.Status_LFS = False
        self.Status_HFS = False
        self.Status_allrhop = False
        if Shotnumber != None:
            self.Load( Shotnumber )

    def __del__( self ):
        self.Unload( )
        self.Status = False
        self.Status_LFS = False
        self.Status_HFS = False
        self.Status_allrhop = False		

    def Unload( self ):
        if self.Status:
            del self.time	
            del self.ne
            del self.R
            del self.z
            self.Status = False
		

    def Load( self ,  Shotnumber, Experiment='AUGD', Diagnostic='RPS', Edition = 0L, tBegin=-1.0, tEnd=12.0, LFSorHFS='LFS', loadAllRhop=False, eqExp = 'AUGD', eqDiag = 'EQH',RICAntenna=1 ):
        
        self.Unload()
        self.Shotnumber = Shotnumber
        if Diagnostic == 'RPS':

            self.Shotnumber = Shotnumber

            if LFSorHFS == 'LFS':
                
                #embed()
                try:
                    sf_LFS = dd.shotfile( Diagnostic, Shotnumber, Experiment, Edition)
                except:
                    print "Error reading shotfile" 
                    return False

                output_LFS = sf_LFS( 'neb_LFS' )
                timeTMP_LFS = output_LFS.time
                index = numpy.where( ( timeTMP_LFS > tBegin ) & ( timeTMP_LFS < tEnd ) )[0]
                self.Diagnostic = Diagnostic
             
                self.R = numpy.copy(output_LFS.area[index])              
                z  = sf_LFS.getParameter('AuxInfo', 'z_lfs',dtype=numpy.float64).data
                self.z = numpy.zeros_like(self.R)
                self.z[:] = z
   
                self.ne = numpy.copy(output_LFS.data[index])
                self.time = numpy.copy(timeTMP_LFS[index])
                self.ntime = numpy.size(self.time)

                self.nChannels = numpy.size(self.R[0,:])
                print 'LFS data read'
                #embed()
                #del output
                self.Status = True
                self.Status_LFS = True
                sf_LFS.close()

            elif LFSorHFS == 'HFS':

                try:
                    sf_HFS = dd.shotfile( Diagnostic, Shotnumber, Experiment, Edition)
                except:
                    print "Error reading shotfile" 
                    return False

                output_HFS = sf_HFS( 'neb_HFS' )
                timeTMP_HFS = output_HFS.time
                index = numpy.where( ( timeTMP_HFS > tBegin ) & ( timeTMP_HFS < tEnd ) )[0]
                
                self.Diagnostic = Diagnostic
                self.R = numpy.copy(output_HFS.area[index])
                z  = sf_HFS.getParameter('AuxInfo', 'z_hfs',dtype=numpy.float64).data
                self.z = numpy.zeros_like(self.R)
                self.z[:] = z

                self.ne = numpy.copy(output_HFS.data[index])
                self.time = numpy.copy(timeTMP_HFS[index])
                self.ntime = numpy.size(self.time)
                
                self.nChannels = numpy.size(self.R[0,:])
                print 'HFS data read'
                #embed()
                self.Status = True
                self.Status_HFS = True
                sf_HFS.close()
            else:
                print 'somethin went wrong'

            if loadAllRhop:
					#check if the entire dataset is read
                if ( tBegin == -1.0 ) & ( tEnd == 12.0 ):

                    self.eqExp = eqExp
                    self.eqDiag = eqDiag
                    self.rhop = kk.KK().kkrzptfn(Shotnumber, self.time, self.R, self.z, exp=eqExp, diag=eqDiag )
                    self.Status_allrhop = True
                else:
                    
                    print 'It is not possible to read Tomas kk with limited timerange, but time is not eqidistant'
                    self.Status_allrhop = False
                    return False
            
#some RPS failed
            #embed()
            self.RPS_exceptions()

            return True
                #self.Status = True	
                

        elif Diagnostic == 'FRS':
                    
            try:
                sf = dd.shotfile( Diagnostic, Shotnumber, Experiment, Edition)
            except:
                print "Error reading FRS  shotfile" 
                return False

            try:
                #embed()
                output = sf( 'ne' )
                
                timeTMP = output.time
                
                index = numpy.where( ( timeTMP > tBegin ) & ( timeTMP < tEnd ) )[0]
                
                self.ne = numpy.copy(output.data[index])
                self.time = numpy.copy(output.time[index])
                self.ntime = numpy.size(self.time)
                self.R = numpy.copy(output.area[index])
                self.Diagnostic = Diagnostic
            #    self.rhop = numpy.copy(sf( 'rhop' ).data)[index]
            #    self.Status_allrhop = True
                #del output
                self.Shotnumber=Shotnumber
                self.z = numpy.zeros_like(self.R)
                self.nChannels = numpy.size(self.R[0,:])
                self.Status = True
                sf.close()
            except:
                print "Error reading FRS data" 
                embed()
                return False
                        
        elif (Diagnostic == 'RIL') | (Diagnostic == 'RIC'):
   

            if Shotnumber==34852:
                import scipy.io as sc
                data=sc.loadmat("/afs/ipp/u/augd/rawfiles/RIF/3485/2/profiles/density_profile_valid_34852_ant%d.mat"%RICAntenna)
                self.time = numpy.squeeze(data['time'])
                self.ntime = numpy.size(self.time)
                self.ne = data['dp_ne']
                self.z = data['dp_z']
                self.R = data['dp_r']
                self.rhop = data['dp_rhop']
                self.Status_allrhop = True
                self.Shotnumber=Shotnumber
                self.Diagnostic = Diagnostic
                self.Status = True
                return True  

            try:
                sf = dd.shotfile( Diagnostic, Shotnumber, Experiment, Edition)
            except:
                print "Error reading %s shotfile"%Diagnostic
                return False
         
            if (RICAntenna == 1) | (RICAntenna == 4) | (RICAntenna == 8):
                antenna = RICAntenna
            else:
                print 'Only Antenna 1, 4 or 8 are available'
                print 'use Antenna 1'
                antenna = 1

            self.Diagnostic = Diagnostic
            if (Diagnostic == 'RIL'):
                output_R = sf( 'R_Ant%d'%antenna )
                timeTMP_HFS = output_R.time
                index = numpy.where( ( timeTMP_HFS > tBegin ) & ( timeTMP_HFS < tEnd ) )[0]
                self.time = output_R.time[index]
                self.ntime = numpy.size(self.time)
                self.R = output_R.data[index]
                self.z = sf( 'Z_Ant%d'%antenna ).data[index]
                self.rhop = sf( 'RhoAnt%d'%antenna ).data[index]
                self.Shotnumber=Shotnumber
                self.Status_allrhop = True
                ne = sf.getParameter('CFG_NE', 'NeVal',dtype=numpy.float64).data
                self.ne = numpy.array([ne,]*self.ntime)
            else:
                output_ne = sf( 'Ne_Ant%d'%antenna )
                timeTMP_HFS = output_ne.time
                index = numpy.where( ( timeTMP_HFS > tBegin ) & ( timeTMP_HFS < tEnd ) )[0]
                self.time = output_ne.time[index]
                self.ntime = numpy.size(self.time)
                self.ne = output_ne.data[index]
                self.z = sf( 'Z_Ant%d'%antenna ).data[index]
                self.R = sf( 'R_Ant%d'%antenna ).data[index]
                self.rhop = sf( 'RhoAnt%d'%antenna ).data[index]
                self.Status_allrhop = True
            
            self.Status = True
            return True      

            try:
                sf = dd.shotfile( Diagnostic, Shotnumber, Experiment, Edition)
            except:
                print "Error reading %s shotfile"%Diagnostic
                return False
         
            if (RICAntenna == 1) | (RICAntenna == 4) | (RICAntenna == 8):
                antenna = RICAntenna
            else:
                print 'Only Antenna 1, 4 or 8 are available'
                print 'use Antenna 1'
                antenna = 1

            self.Diagnostic = Diagnostic
            if (Diagnostic == 'RIL'):
                output_R = sf( 'R_Ant%d'%antenna )
                self.time = output_R.time
                self.ntime = numpy.size(self.time)
                self.Shotnumber=Shotnumber
                self.R = output_R.data
                self.z = sf( 'Z_Ant%d'%antenna ).data
                self.rhop = sf( 'RhoAnt%d'%antenna ).data
                self.Status_allrhop = True
                ne = sf.getParameter('CFG_NE', 'NeVal',dtype=numpy.float64).data
                self.ne = numpy.array([ne,]*self.ntime)
            else:
                output_ne = sf( 'Ne_Ant%d'%antenna )
                self.time = output_ne.time
                self.ntime = numpy.size(self.time)
                self.ne = output_ne.data
                self.z = sf( 'Z_Ant%d'%antenna ).data
                self.R = sf( 'R_Ant%d'%antenna ).data
                self.rhop = sf( 'RhoAnt%d'%antenna ).data
                self.Status_allrhop = True
            

            self.Shotnumber = Shotnumber
            self.Status = True
            return True      

    def RPS_exceptions(self):
        if self.Status:
            if self.Shotnumber == 33568:
                tskip = [3.9,4.3]
            elif self.Shotnumber == 33570:
                tskip = [3.08,5.43]
            elif self.Shotnumber ==33346:
                tskip = [3.9,12.0]        
            else:
                return

            #embed()
            idxOut = numpy.where((self.time >= tskip[0])&(self.time <= tskip[1]))[0]
            self.time = numpy.delete(self.time,idxOut)
            self.ne = numpy.delete(self.ne,idxOut,axis=0)
            self.R = numpy.delete(self.R,idxOut,axis=0)
            self.z = numpy.delete(self.z,idxOut,axis=0)
            if self.Status_allrhop:
                self.rhop = numpy.delete(self.rhop,idxOut,axis=0)
            self.ntime = numpy.size(self.time)
            
            return
                        

    def getphi( self ):
        if self.Diagnostic == 'RPS':
            return PHI.RPS()
        else:
            return PHI.FRS()
            
    def getRsep( self,  Experiment = 'AUGD', Diagnostic='EQH'):
        
        if self.Status:

            Rmin=1.8
            Rmax=2.2

	    ntimes = self.ntime
            R = self.getData( self.time,value='R')
            z = self.getData( self.time,value='z')
            p = numpy.polyfit(numpy.ravel(R), numpy.ravel(z),  1)
            lin = numpy.poly1d(p)

           

            Rscan =numpy.arange(Rmin,Rmax,0.001)
            zscan=lin(Rscan)
            out_temp = numpy.zeros((numpy.size(self.time),numpy.size(Rscan)))
            	
            output = kk.KK().kkrzptfn( self.Shotnumber, self.time, Rscan, zscan, exp= Experiment, diag=Diagnostic, ed=0)
            Rsep = numpy.zeros((numpy.size(self.time)))
            zsep = numpy.zeros((numpy.size(self.time)))

           

            for i in numpy.arange(ntimes):
                Rsep[i] = numpy.squeeze(numpy.interp([1.0],output.rho_p[i,:],Rscan))
                zsep[i] = numpy.squeeze(numpy.interp([1.0],output.rho_p[i,:],zscan))

            self.Rsep,self.zsep=Rsep,zsep

            return Rsep
            #out_temp[i,:] = output.rho_p	
            
            #rhop = self.getData( self.time,value='rhop')
            #embed()
            #Rsep = self.getTimetrace( selectData_in = 'rhop', data_in = [1.0], selectData_out = 'R' )	
            #self.R_sep = self.Rall - numpy.reshape(numpy.repeat(LIB_Rsep,self.nChannels),(self.ntime,self.nChannels))
            #return LIB_Rsep


 #   def getRsep( self ):
 #       if self.Status_allrhop:
 #           ntimes = self.ntime
 #           data = numpy.zeros((ntimes))
 #           data_in = [1.0]
 #           x = self.rhop
 #           y = self.R#

#            if x[0,0] < x[0,-1]:
#                data = tools.interpolSimpleArray(x,y,data_in)
#            else:
#                data = tools.interpolSimpleArray(x[:,::-1],y[:,::-1],data_in)
                
#            self.R_sep = self.R - numpy.reshape(numpy.repeat(data,self.nChannels),(ntimes,self.nChannels))
            
#            return data



	#Default Te, indexbased
    def getData( self , timepoints = None, arange = None, value='ne', Experiment = 'AUGD', Diagnostic='EQH' ):
        ntimes = 0

        if self.Status:
            
            if not ( (value == 'ne') | (value == 'rhop')| (value == 'time')| (value == 'R')| (value == 'z') ):
                print 'value must be Te or ne or rhop'
                return None

            if (value == 'rhop'):
                if not self.Status_allrhop:
                    print 'not all rhop values loaded, every point has to be loaded: very slow'

            if (timepoints == None) & (range == None) :
                print 'no timepoints or range are given, return None'
                return None

            if timepoints != None:
                            
                timepoints = numpy.array(timepoints)
                            
                if arange != None:
                    print 'either range or timepoints must be given'
                    return False

                ntimes = numpy.size(timepoints)
                if ntimes > 1:
					#get min from the time matrix of the timepoints vector, 1. self.time is copied ntimes, then timepoints is subtracted row wise.
                
                    search_array = numpy.sum([numpy.array([self.time,]*ntimes).T,-numpy.array([timepoints,]*1)],axis=0).T
					# the indices of the closest points are calculated
                    idx = numpy.nanargmin(numpy.abs(search_array),axis=1)
				
					#if (value == 'Te') | (value == 'ne') | (value == 'rhop'):
					#ne/Te/rhop should have the same area base
						#out_temp = numpy.zeros((ntimes,numpy.size(self.Te[0,:])))
                                
                elif ntimes == 1:
                    idx = numpy.argmin( numpy.abs( self.time - timepoints ) )
                    if (value == 'time'): 
                        return self.time[idx]
 
			

            if arange != None:
                arange = numpy.array(arange)
                
                if not numpy.size(arange) == 2:
                    print 'range must have two elements'
                    return None
                            
                if ( (arange[0]>arange[1]) & (arange[0] >= self.time.min()) & (arange[1] <= self.time.max())):
                    print 'second Element must be bigger than first and range must be in range of tBegin or tEnd'
                    return None				
                            
                if timepoints != None:
                    print 'either range or timepoints must be given'
                    return None
				
                idx = numpy.squeeze(numpy.where( (self.time[:] >= arange[0]) & (self.time[:] <= arange[1] ) ))
				
            if (value == 'time'): 
                return self.time[idx]

	
            if (value == 'ne') | (value == 'rhop') | (value == 'R') | (value == 'z'):
					#ne/Te/rhop should have the same area base
								
                if value == 'ne':
                    return  self.ne[idx,:]
                elif value == 'R':
                    return  self.R[idx,:]
                elif value == 'z':
                    return  self.z[idx,:]

                if (value == 'rhop'):
                    if (self.Status_allrhop):
                        if  (Experiment == self.eqExp) & (Diagnostic == self.eqDiag):
                            return self.rhop[idx,:]
                        else:
                            print 'Equlibrium Experiment differ, calculate new'
                            rhop_out = fastkk.eqi_map().kkrzptfn(self.Shotnumber, self.time, self.R, self.z, exp=Experiment, diag=Diagnostic )
                            self.eqDiag = Diagnostic
                            self.eqExp = Experiment
                            return rhop_out[idx,:]

                    else:
                        nindex = numpy.size(idx)
                        if nindex > 1: 
                            out_temp = numpy.zeros((nindex,numpy.size(self.ne[0,:])))
                            for i in numpy.arange(nindex):	
                                output = kk.KK().kkrzptfn( self.Shotnumber, self.time[idx[i]], self.R[idx[i],:], self.z[idx[i],:], exp= Experiment, diag=Diagnostic, ed=0)
                                out_temp[i,:] = output.rho_p	
                                return out_temp
                        else:
                            output = kk.KK().kkrzptfn( self.Shotnumber, self.time[idx], self.R[idx,:], self.z[idx,:], exp= Experiment, diag=Diagnostic, ed=0)
                            return output.rho_p	


    def get(self, time, eqExp = 'AUGD', eqDiag='EQH'):

        if self.Status:		
            if (numpy.size(time)== 1) | ( numpy.size(time) > 2):
                output = RPShelp()
                output.ne = self.getData( timepoints = time, value='ne')
                output.rhop = self.getData( timepoints = time, value='rhop')
                output.time = self.getData( timepoints = time, value='time')
                return output

            elif numpy.size(time)== 2:
                output = RPShelp()
                output.ne = self.getData( arange = time, value='ne')
                output.rhop = self.getData( arange = time, value='rhop')
                output.time = self.getData( arange = time, value='time')
                return output
            else:
                print 'no Data loaded'
                return None





    def getne(self, timepoint):

        if self.Status:		
            if (numpy.size(timepoint)== 1) | ( numpy.size(timepoint) > 2):
                return self.getData( timepoints = timepoint, value='ne')
            
            elif numpy.size(timepoint)== 2:
                return  self.getData( arange = timepoint, value='ne')

            return None
        else:
            print 'no Data loaded'
            return None



    def getRhop(self, timepoint, eqExp = 'AUGD', eqDiag='EQH'):

        if self.Status:		
            if (numpy.size(timepoint)== 1) | ( numpy.size(timepoint) > 2):
                return self.getData( timepoints = timepoint, value='rhop')
            
            elif numpy.size(timepoint)== 2:
                return  self.getData( timepoints  = timepoint, value='rhop')

            return None
        else:
            print 'no Data loaded'
            return None


    def getR(self, timepoint):

        if self.Status:		
            if (numpy.size(timepoint)== 1) | ( numpy.size(timepoint) > 2):
                return self.getData( timepoints = timepoint, value='R')
            
            elif numpy.size(timepoint)== 2:
                return  self.getData( arange = timepoint, value='R')

            return None
        else:
            print 'no Data loaded'
            return None

    def getz(self, timepoint):

        if self.Status:		
            if (numpy.size(timepoint)== 1) | ( numpy.size(timepoint) > 2):
                return self.getData( timepoints = timepoint, value='z')
            
            elif numpy.size(timepoint)== 2:
                return  self.getData( arange = timepoint, value='z')

            return None
        else:
            print 'no Data loaded'
            return None


    def getTime(self, timepoint):

        if self.Status:		
            if (numpy.size(timepoint)== 1) | ( numpy.size(timepoint) > 2):
                return self.getData( timepoints = timepoint, value='time')

            elif numpy.size(timepoint)== 2:
                return  self.getData( arange = timepoint, value='time')

            return None
        else:
            print 'no Data loaded'
            return None


## only one value for ne, get R and Z
    def getTimetrace( self , time_in = None, selectData_in = 'ne', data_in = [1.5e19], selectData_out = 'R' ,kind = '1D'):

        if self.Status:

			#if no time is given, every point will be used to calc timetrace
            if (numpy.all(time_in) == None) :
                time_in = self.time
			
            if selectData_in == 'ne':
                x = self.getne(time_in)
            if selectData_in == 'Te':
                x = self.getTe(time_in)
            if (selectData_in == 'Rhop')| (selectData_in == 'rhop'):
                x = self.getRhop(time_in)
			

            if selectData_out == 'R':
                y = self.getR(time_in)
            if selectData_out == 'z':
                y = self.getz(time_in)
            if (selectData_out == 'Rhop') | (selectData_out == 'rhop'):
                y = self.getRhop(time_in)
            if selectData_out == 'ne':
                y = self.getne(time_in)
            if selectData_out == 'Te':
                y = self.getTe(time_in)

            ntimes = numpy.size(time_in)
            data = numpy.zeros((ntimes,numpy.size(data_in)))
			#embed()
            if kind == '1D':
                for i in numpy.arange(ntimes):
                    try:
                        idxSort=numpy.argsort(x[i])
                        data[i,:] =numpy.interp(data_in,x[i][idxSort],y[i][idxSort])	
                    except:
                        print 'upsi daisy'
                        #embed()

                return numpy.squeeze(data)
				
            elif kind == '2D':
                time = numpy.reshape( numpy.repeat( self.time, ( self.nChannels ) ) ,( self.ntime, self.nChannels) )
                f =  scipy.interpolate.interp2d(time, x, y, kind='linear', fill_value = float('NaN'), bounds_error=False)
                return f(time_in,data_in)
            


    def __call__( self , time , rhop = None, eqExp = 'AUGD', eqDiag='EQH' ):

        if self.Status:
            if rhop == None:
                return self.getne( time )
            else:
                ne_out = self.getne( time )
                rhop_out = self.getRhop( time )
                ntimes = numpy.size(time)
        if ntimes > 1:
            ne_temp = numpy.zeros((ntimes,numpy.size(rhop)))
            for i in numpy.arange(ntimes) :						
                if (rhop_out[i,1]<= rhop_out[i,0]):
                    ne = scipy.interpolate.interp1d(rhop_out[i,::-1] , ne_out[i,::-1],  bounds_error=False, fill_value = float('inf'))
                else:
                    ne = scipy.interpolate.interp1d(rhop_out[i,:] , ne_out[i,:1],  bounds_error=False, fill_value = float('inf'))
                    ne_temp[i,:] = ne( rhop )
                    return ne_temp
            else:
                if rhop_out[i,1] <= rhop_out[i,0]:
                    ne = scipy.interpolate.interp1d(rhop_out[i,::-1] , ne_out[i,::-1],  bounds_error=False, fill_value = float('inf'))
                else:
                    ne = scipy.interpolate.interp1d(rhop_out[i,:] , ne_out[i,:1],  bounds_error=False, fill_value = float('inf'))
                    return ne( rhop )					

        else:
            print 'no Data loaded'
            return None


