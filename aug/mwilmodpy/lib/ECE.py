import numpy
import dd as dd
import kk_mwillens as kk
import eqi_map as fastkk
from scipy.interpolate import interp1d
import IPython

class ECEhelp:
    status = False


class ECE:
    def __init__( self ,  Experiment = 'AUGD', Diagnostic = 'CEC', Shotnumber = None ):
        self.Status = False
        self.StatusRz = False
        self.exists_Rmaj = False
        self.exists_rhop = False
        self.Status_allrhop = False

        if Shotnumber != None :
            self.Load( Shotnumber )
		

    def __del__( self ):
        self.Unload( )
        del self.Status

    def Load( self ,  Shotnumber, Experiment='AUGD', Diagnostic='CEC', Edition = 0L, tBegin=-1.0, tEnd=12.0, loadAllRhop=False, eqExp = 'AUGD', eqDiag = 'EQH', Rshift=None,Binning = 20.0,debug=False,onlyHeader=False):
        self.Unload()
        if (Diagnostic == 'CEC') or (Diagnostic == 'RMD') or (Diagnostic == 'RMC'):
            try:
                if  Diagnostic == 'RMC':
                    DiagnosticIn = 'RMD'
                else:
                    DiagnosticIn = Diagnostic

                sf = dd.shotfile( DiagnosticIn, Shotnumber, Experiment, Edition)
                self.Shotnumber = Shotnumber
            except:
                print "Error reading shotfile" 
                return False

            

#####################################################################################################################################
##################### READING calibration of CEC or RMD
#####################################################################################################################################
            try: 
                self.tBegin = tBegin
                self.tEnd = tEnd

                if Diagnostic != 'CEC':
                    self.Calfact = sf.getParameter('parms-A', 'calfact',dtype=numpy.float64).data
                    self.Shift = sf.getParameter('parms-A', 'shift',dtype=numpy.float64).data
                    self.Multi00 = numpy.concatenate( (sf.getParameter('eCAL-A1', 'MULTIA00',dtype=numpy.float64).data,sf.getParameter('eCAL-A2', 'MULTIA00',dtype=numpy.float64 ).data),axis=0)
                    self.Shift00 = numpy.concatenate( (sf.getParameter('eCAL-A1', 'SHIFTB00',dtype=numpy.float64).data,sf.getParameter('eCAL-A2', 'SHIFTB00',dtype=numpy.float64).data),axis=0)
            except:
                print "Error reading Calibration shotfile" 
                if debug:
                    IPython.embed()
                return False

#####################################################################################################################################
##################### READING AreaBase of CEC or RMD
#####################################################################################################################################
            try:
                print "get Areabase z"

                print "Reading z"
                self.z = sf( 'z-A', tBegin=tBegin, tEnd=tEnd ).data
                print "get Areabase R"
                print "Reading R"
                self.R = sf( 'R-A', tBegin=tBegin, tEnd=tEnd ).data 
                if (numpy.all(Rshift) != None):
                    self.R += Rshift
                    self.Rshift = Rshift

    #				print "get timebase of Areabase"
                self.rztime = sf.getTimeBase( 'R-A', tBegin=tBegin, tEnd=tEnd )

                self.nChannels = numpy.size(self.z[0,:])
            except:
                print "Error reading area base" 
                if debug:
                    IPython.embed()
                return False
#####################################################################################################################################
##################### READING Parameters of CEC or RMD
#####################################################################################################################################
            try:
                self.Freq = sf.getParameter('parms-A', 'f',dtype=numpy.float64).data
                self.Btot = sf.getParameter('parms-A', 'Btot',dtype=numpy.float64).data
                self.Ifgroup = sf.getParameter('parms-A', 'IFGROUP',dtype=numpy.float64).data
                self.Calfact = sf.getParameter('parms-A', 'calfact',dtype=numpy.float64).data
                self.avail  = sf.getParameter('parms-A', 'AVAILABL',dtype=numpy.float64).data
            except:
                print "Error reading parameters" 
                if debug:
                    IPython.embed()
                return False


            try:
                self.Sideband = sf.getParameter('METHODS', 'SIDEBAND').data
                self.Freqlo = sf.getParameter('METHODS', 'FREQLO',dtype=numpy.float64).data
                self.nIfgroups = int(sf.getParameter('METHODS', 'IFGROUPS',dtype=numpy.float64) .data)
            except:
                print "Error reading extra parameters" 

            try:
                self.snr  = sf.getParameter('parms-A', 'SNR_cali',dtype=numpy.float64).data 
            except:
                print "Error reading snr parameters" 

            if onlyHeader:
                return


            if (Diagnostic == 'RMD')| (Diagnostic == 'CEC'):

#####################################################################################################################################
##################### READING Te of CEC or RMD
#####################################################################################################################################
                try:
                    print "Reading Time"
                    self.time = sf.getTimeBase( 'Trad-A', tBegin=tBegin, tEnd=tEnd )
                    self.ntime = numpy.size(self.time)
                    print "get Trad-A"
                    print "Reading Te"
                    self.Te = sf( 'Trad-A', tBegin=tBegin, tEnd=tEnd ).data
                    sf.close()
                except:
                    print "Error reading Te" 
                    if debug:
                        IPython.embed()
                    return False     
            else:
		#reading RMC data
                sf.close()

#####################################################################################################################################
##################### Reading RMC
#####################################################################################################################################

                try:
                    sf_RMC = dd.shotfile( 'RMC' , Shotnumber, 'AUGD', Edition)
                    
                except:
                    print "Error opening RMC shotfile" 
                    if debug:
                        IPython.embed()
                    return

                try:	
                    self.time = sf_RMC.getTimeBase( 'Trad-A1', tBegin=tBegin, tEnd=tEnd )
                    self.ntimes = numpy.size(self.time)
                        
                    # read raw rmc data Trad-A1 and Trad-A2 and sort it			
                    self.Te = numpy.concatenate( (sf_RMC('Trad-A1', tBegin=tBegin, tEnd=tEnd,dtype=numpy.float64).data,sf_RMC('Trad-A2',tBegin=tBegin,tEnd=tEnd,dtype=numpy.float64).data),axis=1 )[:,:]
                    sf_RMC.close()
                    print "RMC read"

                    self.Te = numpy.multiply(self.Multi00,self.Te,dtype=numpy.float64)
                    #		self.Te = numpy.add(self.Shift00,self.Te,dtype=numpy.float64)
                    if (tBegin < -0.1):
                        SelectTime = numpy.where( self.time <= 0.0 ) 
                        CalcShift = numpy.reshape(numpy.mean(self.Te[SelectTime,:],axis=1,dtype=numpy.float64),(self.nChannels))
                        self.Te = numpy.add(-CalcShift,self.Te,dtype=numpy.float64)
                    else:		
                        self.Te = numpy.add(self.Shift00,self.Te,dtype=numpy.float64)
                except:
                    print "Error reading Te RMC" 
                    if debug:
                        IPython.embed()
                    return

            self.samplingrate = 1.0/numpy.mean(numpy.diff(self.time))
            self.phi = 191.25


            self.Status = True

#####################################################################################################################################
##################### Binning of Te 
#####################################################################################################################################
            
            if (numpy.all(Binning) != None) & (Binning*1.e3 < self.samplingrate):
                self.Binning(Binning)
                self.samplingrate = 1.0/numpy.mean(numpy.diff(self.time))

#####################################################################################################################################
##################### Get R and z for every timepoint
#####################################################################################################################################

            self.Rall = numpy.zeros_like(self.Te)
            self.zall = numpy.zeros_like(self.Te)
                               
            for i in numpy.arange(self.nChannels):
                self.Rall[:,i] = numpy.interp(self.time,self.rztime,self.R[:,i])
                self.zall[:,i] = numpy.interp(self.time,self.rztime,self.z[:,i]) 

            self.max_Te = numpy.max(self.Te)
                      
					#calcule all rhop values
                        
            if loadAllRhop:
                print "Calculate rhop values"
                self.funcLoadAllRhop(eqExp = eqExp, eqDiag = eqDiag )
			
            self.ChannelNr =  numpy.arange(1,numpy.size(self.Te[0,:])+1)
            self.ChannelsSorted = False
            self.Shotnumber = Shotnumber			
            self.SortECEChannels()
            print "close"
            return True
        else:
            print 'Diagnotic should be CEC or RMD or RMC, but is:' ,Diagnostic


    def Unload( self ):
        if self.Status:
            if self.StatusRz:
                self.UnloadRz()
            self.Status = False
            del self.tBegin
            del self.tEnd
            del self.time	
            del self.Te
            del self.nChannels
            del self.R
            del self.z
            del self.phi
            del self.Rall
            del self.zall
            del self.rztime
            del self.Shotnumber
            del self.max_Te
                    
                    
            del self.Freq
            del self.Btot
            del self.Ifgroup
            del self.Calfact
            del self.avail
            try:
                del self.nIfgroups
            except:
                print 'self.nIfgroups not deleting'

            del self.Sideband
            del self.Freqlo
			
            if self.Status_allrhop:
                del self.rhop
                del self.eqExp
                del self.eqDiag
                self.Status_allrhop = False

            if self.ChannelsSorted:
                del self.SortedIndex
                try:
                    del self.IfGroupIndex
                except:
                    print 'no IfGroupIndex'

                self.ChannelsSorted = False


    def UnloadRz( self ):
        if (self.StatusRz) & (self.Status == False):
            self.StatusRz = False
            del self.R
            del self.z 
            del self.rztime
            del self.phi

    def LoadRz( self ,  Shotnumber, Experiment='AUGD', Diagnostic='CEC', Edition = 0L, tBegin=-1.0, tEnd=12.0,  eqExp = 'AUGD', eqDiag = 'EQH' ):
        self.UnloadRz()
        if Diagnostic == 'CEC' or Diagnostic == 'RMD':
            try:
                sf = dd.shotfile( Diagnostic, Shotnumber, Experiment, Edition)
                self.Shotnumber = Shotnumber
            except:
                print "Error reading shotfile" 
                return False

            print "Reading z"
            self.z = sf( 'z-A', tBegin=tBegin, tEnd=tEnd ).data
			#				print "get Areabase R"
            print "Reading R"
            self.R = sf( 'R-A', tBegin=tBegin, tEnd=tEnd ).data
#				print "get timebase of Areabase"
            self.rztime = sf.getTimeBase( 'R-A', tBegin=tBegin, tEnd=tEnd )

            sf.close()
 
            self.phi = 191.25
                
            self.StatusRz = True


    def getRsep( self ):

        if self.Status_allrhop:
            ntimes = self.ntime
            data = numpy.zeros((ntimes))
            data_in = [1.0]
                
            if self.ChannelsSorted:
                rhop = self.rhop[:,self.SortedIndex]
                Rall = self.Rall[:,self.SortedIndex]
            else:
                rhop = self.rhop
                Rall = self.Rall
                    
            print 'tools is neccessar'
            import sys
            sys.path.append("/afs/ipp-garching.mpg.de/home/m/mwillens/Programme/python/lib/")
            import tools
            data = tools.interpolSimpleArray(rhop[:,::-1],Rall[:,::-1],data_in)

            self.R_sep = self.Rall - numpy.reshape(numpy.repeat(data,self.nChannels),(ntimes,self.nChannels))
            
            return data


        #single time point
    def getRzphi(self, time ):
        if self.StatusRz | self.Status:

            idx =  numpy.nanargmin(numpy.abs(self.rztime - time))
            R_out= numpy.squeeze(self.R[idx,:]) 
            z_out= numpy.squeeze(self.z[idx,:]) 
            phi = numpy.repeat([self.phi],numpy.size(R_out))
            return R_out, z_out, phi


    def funcLoadAllRhop( self, eqExp = 'AUGD', eqDiag = 'EQH', ed=0, eqi_map=False ):
		
        if self.Status:

            if (((self.tBegin==-1.0) & (self.tEnd == 12.0)) & eqi_map):
                self.eqExp = eqExp
                self.eqDiag = eqDiag
				
                R =self.Rall  # numpy.zeros_like(self.Te)
                z =self.zall  # numpy.zeros_like(self.Te)
                #IPython.embed()
                self.rhop = fastkk.eqi_map().kkrzptfn(self.Shotnumber, self.time, R, z, exp=eqExp, diag=eqDiag )
                self.q = fastkk.eqi_map().kkrhotpq(self.Shotnumber, self.time, self.rhop, exp=eqExp, diag=eqDiag )
                self.Status_allrhop = True
            else:
                self.eqExp = eqExp
                self.eqDiag = eqDiag       
                out = kk.KK().kkrzptfn(self.Shotnumber, self.time, self.Rall, self.zall, exp=eqExp, diag=eqDiag, ed=ed)
                self.rhop = out.rho_p
                self.rhot = out.rho_t
                #self.q = out.q
                del out
                        
                self.Status_allrhop = True
			
## only one value for ne, get R and Z
    def getTimetrace( self , time_in = None, selectData_in = 'Te', data_in = [1.5e19], selectData_out = 'R' ,kind = '1D'):

        if self.Status:

			#if no time is given, every point will be used to calc timetrace
            if (numpy.all(time_in) == None) :
                time_in = self.time
			
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
            if selectData_out == 'Te':
                y = self.getTe(time_in)


            R = self.getR(time_in)
            Te= self.getTe(time_in)
            ntimes = numpy.size(time_in)
            data = numpy.zeros((ntimes,numpy.size(data_in)))
			#embed()
            if kind == '1D':
                #IPython.embed()
                for i in numpy.arange(ntimes):
                    try:
                        idxRSort=numpy.argsort(R[i])
                        maxR=numpy.argmax(numpy.gradient(numpy.gradient(Te[i][idxRSort])/numpy.gradient(R[i][idxRSort]))/numpy.gradient(R[i][idxRSort]))
                        minR=numpy.argmin(numpy.gradient(Te[i][idxRSort[:maxR]])/numpy.gradient(R[i][idxRSort[:maxR]]))
                        useIdx=idxRSort[minR:maxR]
                        idxSort=useIdx[numpy.argsort(x[i][useIdx])]
                        data[i,:] =numpy.interp(data_in,x[i][idxSort],y[i][idxSort])	
                    except:
                        print 'uppsala'

                return numpy.squeeze(data)
				
            elif kind == '2D':
                time = numpy.reshape( numpy.repeat( self.time, ( self.nChannels ) ) ,( self.ntime, self.nChannels) )
                f =  scipy.interpolate.interp2d(time, x, y, kind='linear', fill_value = float('NaN'), bounds_error=False)
                return f(time_in,data_in)





    def get(self, time=None, sorted = True, eqExp = 'AUGD', eqDiag='EQH',doInterpolation =  False):

        if self.Status:	
            if (self.ChannelsSorted) & (sorted == True):
                idx = self.SortedIndex
            else:
                idx = numpy.arange(self.nChannels)
	
            if numpy.all(time) == None:
                doInterpolation = False
                time=self.time
                
                
            if (numpy.size(time)== 1):
                output = ECEhelp()
                output.Te = (self.getData( timepoints = time, value='Te',doInterpolation = doInterpolation))[idx]
                output.rhop = (self.getData( timepoints = time, value='rhop', eqExp = eqExp, eqDiag=eqDiag,doInterpolation = doInterpolation))[idx]
                output.time = self.getData( timepoints = time, value='time',doInterpolation = doInterpolation)
                return output

            elif numpy.size(time)== 2:
                output = ECEhelp()
                output.Te = (self.getData( range = time, value='Te',doInterpolation = doInterpolation))[:,idx]
                output.rhop = (self.getData( range = time, value='rhop', eqExp = eqExp, eqDiag=eqDiag,doInterpolation = doInterpolation))[:,idx]
                output.time = self.getData( range = time, value='time',doInterpolation = doInterpolation)
                return output
            elif numpy.size(time) > 2:
                output = ECEhelp()
                output.Te = (self.getData( timepoints = time, value='Te',doInterpolation = doInterpolation))[:,idx]
                output.rhop = (self.getData( timepoints = time, value='rhop', eqExp = eqExp, eqDiag=eqDiag,doInterpolation = doInterpolation))[:,idx]
                output.time = self.getData( timepoints = time, value='time',doInterpolation = doInterpolation)
                
            else:
                print 'no Data loaded'
                return None


    def getTe(self, timepoint=None, sorted = True,doInterpolation =  False):

        if self.Status:	

            if (self.ChannelsSorted) & (sorted == True):
                idx = self.SortedIndex
            else:
                idx = numpy.arange(self.nChannels)
                
            if numpy.all(timepoint) == None:
                doInterpolation = False
                timepoint=self.time

            if (numpy.size(timepoint) == 1) :
                return (self.getData( timepoints = timepoint, value='Te',doInterpolation = doInterpolation))[idx]     
            elif numpy.size(timepoint) == 2:
                return  (self.getData( range = timepoint, value='Te',doInterpolation = doInterpolation))[:,idx]
            elif (numpy.size(timepoint) > 2):
                return  (self.getData( timepoints = timepoint, value='Te',doInterpolation = doInterpolation))[:,idx]
            return None
        else:
            print 'no Data loaded'
            return None


    def getR(self, timepoint=None, sorted = True,doInterpolation =  False):

        if self.Status:	

            if numpy.all(timepoint) == None:
                doInterpolation = False
                timepoint=self.time

            if (self.ChannelsSorted) & (sorted == True):
                idx = self.SortedIndex
            else:
                idx = numpy.arange(self.nChannels)
                
            if (numpy.size(timepoint) == 1) :
                return (self.getData( timepoints = timepoint, value='R',doInterpolation = doInterpolation))[idx]     
            elif numpy.size(timepoint) == 2:
                return  (self.getData( range = timepoint, value='R',doInterpolation = doInterpolation))[:,idx]
            elif (numpy.size(timepoint) > 2):
                return  (self.getData( timepoints = timepoint, value='R',doInterpolation = doInterpolation))[:,idx]
            return None
        else:
            print 'no Data loaded'
            return None

    def getRhop(self, timepoint=None, sorted = True, eqExp = 'AUGD', eqDiag='EQH',doInterpolation =  False):

        if self.Status:		

            if numpy.all(timepoint) == None:
                doInterpolation = False
                timepoint=self.time

            if (self.ChannelsSorted) & (sorted == True):
                idx = self.SortedIndex
            else:
                idx = numpy.arange(self.nChannels)
                    
            if (numpy.size(timepoint)==1) :
                return (self.getData( timepoints = timepoint, value='rhop', eqExp = eqExp, eqDiag=eqDiag,doInterpolation = doInterpolation))[idx] 
            elif (numpy.size(timepoint) == 2):
                return  (self.getData( range = timepoint, value='rhop', eqExp = eqExp, eqDiag=eqDiag,doInterpolation = doInterpolation))[:,idx]
            elif (numpy.size(timepoint) > 2) :
                return (self.getData( timepoints = timepoint, value='rhop', eqExp = eqExp, eqDiag=eqDiag,doInterpolation = doInterpolation))[:,idx] 
            return None
        else:
            print 'no Data loaded'
            return None

    def getq(self, timepoint=None, sorted = True, eqExp = 'AUGD', eqDiag='EQH',doInterpolation =  False):

        if self.Status:		

            if numpy.all(timepoint) == None:
                doInterpolation = False
                timepoint=self.time

            if (self.ChannelsSorted) & (sorted == True):
                idx = self.SortedIndex
            else:
                idx = numpy.arange(self.nChannels)
                    
            if (numpy.size(timepoint)==1) :
                return (self.getData( timepoints = timepoint, value='q', eqExp = eqExp, eqDiag=eqDiag,doInterpolation = doInterpolation))[idx] 
            elif (numpy.size(timepoint) == 2):
                return  (self.getData( range = timepoint, value='q', eqExp = eqExp, eqDiag=eqDiag,doInterpolation = doInterpolation))[:,idx]
            elif (numpy.size(timepoint) > 2) :
                return (self.getData( timepoints = timepoint, value='q', eqExp = eqExp, eqDiag=eqDiag,doInterpolation = doInterpolation))[:,idx] 
            return None
        else:
            print 'no Data loaded'
            return None



    def MapToRhop(self, timepoint=None, eqExp = 'AUGD', eqDiag='EQH',doInterpolation = False ):
        if self.Status:	
            return self.getRhop( timepoint, eqExp = eqExp, eqDiag=eqDiag,doInterpolation = doInterpolation)


	#Default Te, indexbased
    def getData( self , timepoints = None, range = None, value='Te', eqExp = 'AUGD', eqDiag='EQH',doInterpolation =  False):
        ntimes = 0

        if self.Status:

            if not ( (value == 'Te') | (value == 'rhop')| (value == 'time')| (value == 'R')| (value == 'z') |  (value == 'q')):
                print 'value must be Te or ne or rhop'
                return None

            if (value == 'rhop'):
                if not self.Status_allrhop:
                    print 'not all rhop values loaded, every point has to be loaded: very slow'

            if (numpy.all(timepoints) == None) & (numpy.all(range) == None) :
                print 'no timepoints or range are given, return None'
                return None

            
            if doInterpolation == False:
                if numpy.all(timepoints) != None:
				
                    timepoints = numpy.array(timepoints)
                    if range != None:
                        print 'either range or timepoints must be given'
                        return False

                    ntimes = numpy.size(timepoints)
                    if ntimes > 1:
                        try:
                            if numpy.all(timepoints==self.time):
                                idx=numpy.arange(self.ntime)
                            else:
					#get min from the time matrix of the timepoints vector, 1. self.time is copied ntimes, then timepoints is subtracted row wise.
                                search_array = numpy.sum([numpy.array([self.time,]*ntimes).T,numpy.multiply(timepoints,-1.)],axis=0).T
                            # the indices of the closest points are calculated
                                idx = numpy.nanargmin(numpy.abs( numpy.reshape( numpy.concatenate(search_array),(ntimes,self.ntime) )),axis=1)
				
					#if (value == 'Te') | (value == 'ne') | (value == 'rhop'):
					#ne/Te/rhop should have the same area base
						#out_temp = numpy.zeros((ntimes,numpy.size(self.Te[0,:])))
                        except:
                            print 'failure in timepoints'
                            IPython.embed()
                    elif ntimes == 1:
                        idx = numpy.argmin( numpy.abs( self.time - timepoints ) )
                        if (value == 'time'): 
                            return self.time[idx]
 
            else:
                if (value == 'time'): 
                    return numpy.array(timepoints)
                
            if numpy.all(range) != None:
                range = numpy.array(range)
                        
                if not numpy.size(range) == 2:
                    print 'range must have two elements'
                    return None
                if ( (range[0]>range[1]) & (range[0] >= self.tBegin) & (range[1] <= self.tEnd)):
                    print 'second Element must be bigger than first and range must be in range of tBegin or tEnd'
                    return None				
                if timepoints != None:
                    print 'either range or timepoints must be given'
                    return None
                
                idx = numpy.squeeze(numpy.where( (self.time[:] >= range[0]) & (self.time[:] <= range[1] ) ))
		 	
                if (value == 'time'): 
                    return self.time[idx]

                doInterpolation = False
           

            
            if doInterpolation == True:
                if (value == 'Te') |  (value == 'rhop') |  (value == 'R') |  (value == 'z')|  (value == 'q'):
                    timepoints = numpy.array(timepoints)
                    if range != None:
                        print 'either range or timepoints must be given'
                        return False
                    ntimes = numpy.size(timepoints)

                    if value == 'Te':
 #                       TeOut = numpy.zeros((ntimes,self.nChannels))
                        func = interp1d(self.time,self.Te,axis=0,bounds_error=False  )
                        return func(timepoints)
                    elif value == 'R':
                        func = interp1d(self.time,self.Rall,axis=0,bounds_error=False )
                        return func(timepoints)
                    elif value == 'z':
                        func = interp1d(self.time,self.zall,axis=0,bounds_error=False )
                        return func(timepoints)

                    if (value == 'rhop') | (value == 'q'):
                        if (self.Status_allrhop):
                            if  (eqExp == self.eqExp) & (eqDiag == self.eqDiag):
                                if (value == 'rhop'):
                                    func = interp1d(self.time,self.rhop,axis=0,bounds_error=False )
                                else:
                                    func = interp1d(self.time,self.q,axis=0,bounds_error=False )
                                return func(timepoints)
                            else:
                                print 'Equlibrium Experiment differ, calculate new'
                                self.funcLoadAllRhop( eqExp=eqExp, eqDiag=eqDiag)
                                if (value == 'rhop'):
                                    func = interp1d(self.time,self.rhop,axis=0,bounds_error=False )
                                else:
                                    func = interp1d(self.time,self.q,axis=0,bounds_error=False )
                                return  func(timepoints)
                        else:
                            Rfunc = interp1d(self.time,self.Rall,axis=0,bounds_error=False )
                            zfunc = interp1d(self.time,self.zall,axis=0,bounds_error=False)
                            
                            output = kk.KK().kkrzptfn( self.Shotnumber, timepoints, Rfunc(timepoints), zfunc(timepoints), exp= eqExp, diag=eqDiag, ed=0)
                          
                            if (value == 'rhop'):
                                return output.rho_p
                            else:
                                return output.q





                print 'hallo'
            else:
                if (value == 'Te') |  (value == 'rhop') |  (value == 'R') |  (value == 'z')|  (value == 'q'):
					#ne/Te/rhop should have the same area base
				
                    if value == 'Te':
                        return  self.Te[idx,:]

                    if value == 'R':
                        return  self.Rall[idx,:]

                    if value == 'z':
                        return  self.zall[idx,:]
				
                    if (value == 'rhop')|  (value == 'q'):
                        if (self.Status_allrhop):
                            if  (eqExp  == self.eqExp) & (eqDiag == self.eqDiag):
                                if (value == 'rhop'):
                                    return self.rhop[idx,:]
                                else:
                                    return self.q[idx,:]
                            else:
                                print 'Equlibrium Experiment differ, calculate new'
                                self.funcLoadAllRhop( eqExp=eqExp, eqDiag=eqDiag)
                                if (value == 'rhop'):
                                    return self.rhop[idx,:]
                                else:
                                    return self.q[idx,:]
                        else:
                            output = kk.KK().kkrzptfn( self.Shotnumber, self.time[idx], self.Rall[idx,:], self.zall[idx,:], exp= eqExp, diag=eqDiag, ed=0)
                            if (value == 'rhop'):
                                return output.rho_p
                            else:
                                return output.q
  
                            
			
                return None
		




    def SortECEChannels( self ):

        if self.Status:
			
            sort_ECE = numpy.zeros(numpy.size(self.Freq), dtype=[('freq',float),('idx',int),('avail',float),('ifgroup',float)])
            sort_ECE['freq'] = self.Freq
            sort_ECE['idx']= numpy.arange(int(numpy.size(self.Freq)))
            sort_ECE['avail']= self.avail 
            sort_ECE['ifgroup']= self.Ifgroup

            a = numpy.sort( sort_ECE, axis=0, order='freq' ) 
            
            self.SortedIndex =  a['idx'][numpy.where(a['avail'] == 1) ]
            #self.ChannelNr =  numpy.arange(1,numpy.size(self.avail)+1)[self.SortedIndex]
            

            self.IfGroupIndex = []
		
            try:
                for i in range(1,self.nIfgroups+1):
                    self.IfGroupIndex.append(a['idx'][numpy.where( (a['ifgroup'] == i) & (a['avail'] == 1) )])
            except:
                print 'no If GroupIndex sorting'
            
            self.ChannelsSorted = True
            self.nChannelsSorted = self.SortedIndex.size




    def __call__( self , timepoints=None, sorted = True ):
        if self.Status:	
            return self.getTe( timepoints)
        
 



    def Binning( self, samplefreq = 10.0 ):				   
        if self.Status: 
            newtime, self.Te = dataBinning( self.time, self.Te, samplefreq = samplefreq )
            try:
                newtime, self.Rall = dataBinning( self.time, self.Rall, samplefreq = samplefreq ) 
                newtime, self.zall = dataBinning( self.time, self.zall, samplefreq = samplefreq )
            except:
                pass

            if self.Status_allrhop:
                self.time,self.rhop = dataBinning(  self.time, self.rhop, samplefreq = samplefreq )
            else:
                self.time = newtime
                
            self.ntime = numpy.size( self.time) 


# return new timebase and Data
def dataBinning( time, data, samplefreq = 1.0 ):			
            						      
    print "binning with ",samplefreq," kHz"
    ntimes= numpy.size(time)
    samplingrate = 1.0/numpy.mean(numpy.diff(time))
    dataShape =numpy.array(numpy.shape(data))  
    #get the time index
    idxOfTime = numpy.squeeze(numpy.where(dataShape == ntimes))
    # if more index with the number of times exists, take the first one
    if numpy.size(idxOfTime) > 1:
        idxOfTime = idxOfTime[0]

    bins = int(ntimes*(float(samplefreq)*1.0e3/samplingrate))

    slices= numpy.linspace(0, ntimes, bins+1, True).astype(int)
    counts = numpy.diff(slices)

    #calculate new timebase
    newTime = numpy.add.reduceat(time, slices[:-1]) / counts
    newNtimes = numpy.size(newTime)

    #create new shape
    newDataShape = dataShape
    #replace old shape
    numpy.put(newDataShape, idxOfTime, newNtimes)
    #create new Data array
    newData = numpy.zeros( (newDataShape)  )

    #simplify array such as the first index is always the timebase
    newData = numpy.swapaxes(newData,0,idxOfTime)
    data = numpy.swapaxes( data,0,idxOfTime )

    storeShape = numpy.shape( newData )

    # rehape data to two dimensions
    data = numpy.reshape(data,(ntimes,numpy.size(data)/ntimes))
    newData = numpy.reshape(newData,(newNtimes,numpy.size(newData)/newNtimes))

    for i in numpy.arange(numpy.shape(data)[1]):
        newData[:,i] = numpy.add.reduceat(data[:,i], slices[:-1]) / counts

#shape back
    newData = numpy.reshape(newData,(storeShape))
    #swap back to original shape
    newData = numpy.swapaxes(newData,0,idxOfTime)

    return newTime,newData

