import numpy
import dd
import scipy.interpolate
import kk_mwillens as kk
import eqi_map as fastkk
import IPython
import tools

class CXRShelp:
    status = False

class CXRS:
	def __init__( self ,  Experiment = 'AUGD', Shotnumber = None ):
		self.Status = False
		self.StatusRz = False
		self.Status_allrhop = False
		if Shotnumber != None:
			self.Load( Shotnumber )

	def __del__( self ):
		self.Unload( )
		self.Status = False
		self.StatusRz = False
		self.Status_allrhop = False		

	def Unload( self ):

		if self.Status:
			del self.Shotnumber
			del self.tBegin
			del self.tEnd
			del self.time
			del self.Ti
                        del self.vrot
                        del self.inte
			del self.R
			del self.z
			del self.phi	
			self.Status = False
                        self.StatusRz = False
		

	def Load( self ,  Shotnumber, Experiment='AUGD', Diagnostic='CEZ', Edition = 0L, tBegin=-1.0, tEnd=12.0, loadAllRhop=False ,eqExp = 'AUGD', eqDiag = 'EQH',Rshift=0.00):
		self.Unload()
		
		if (Diagnostic == 'CEZ') | (Diagnostic == 'CMZ')| (Diagnostic == 'CPZ')| (Diagnostic == 'CNZ')| (Diagnostic == 'CUZ'):
			try:
				sf = dd.shotfile( Diagnostic, Shotnumber, Experiment, Edition)
                                self.Shotnumber = Shotnumber
			except:
				print "Error reading shotfile" 
				return False

                        try:
                            self.time = sf.getTimeBase( 'Ti_c', tBegin=tBegin, tEnd=tEnd )
                            self.ntime = numpy.size(self.time)
                            self.Ti = sf( 'Ti_c', tBegin=tBegin, tEnd=tEnd ).data
                            self.Ti_unc = sf( 'err_Ti_c', tBegin=tBegin, tEnd=tEnd ).data
                            self.inte = sf( 'inte', tBegin=tBegin, tEnd=tEnd ).data
                            self.err_inte = sf( 'err_inte', tBegin=tBegin, tEnd=tEnd ).data
                            self.vrot_unc = sf( 'err_vrot', tBegin=tBegin, tEnd=tEnd ).data
                            self.vrot = sf( 'vrot', tBegin=tBegin, tEnd=tEnd ).data
                            if (Diagnostic == 'CEZ') | (Diagnostic == 'CMZ')| (Diagnostic == 'CPZ')| (Diagnostic == 'CNZ'):
                                self.R = numpy.squeeze(sf( 'R', tBegin=tBegin, tEnd=tEnd ).data)+Rshift
                                self.nChannels = numpy.size(self.R)

                                self.Rall =  numpy.reshape(numpy.tile( self.R,(self.ntime) ) ,(self.ntime,self.nChannels) )
                                self.z = numpy.squeeze(sf( 'z', tBegin=tBegin, tEnd=tEnd ).data)
                                self.zall =  numpy.reshape(numpy.tile( self.z,(self.ntime) ) ,(self.ntime,self.nChannels) )
                                self.phi = sf( 'phi', tBegin=tBegin, tEnd=tEnd ).data-67.5+360.
                                self.phiall =  numpy.reshape(numpy.tile( self.phi,(self.ntime) ) ,(self.ntime,self.nChannels) )
                                self.StatusRz = True

                            self.Status = True
                        except: 
                            print "Error reading data" 
                            sf.close()
                            return False			
					
					
                        self.tBegin = tBegin
			self.tEnd = tEnd

					#calcule all rhop values
			if loadAllRhop:
                            self.eqExp = eqExp
                            self.eqDiag = eqDiag
					#check if the entire dataset is read
                            if (tBegin==-1.0) & (tEnd==12.0):
                                self.eqExp = eqExp
                                self.eqDiag = eqDiag
                                self.rhop = fastkk.eqi_map().kkrzptfn(self.Shotnumber, self.time, self.R, self.z, exp=eqExp, diag=eqDiag )
                                self.Status_allrhop = True
                            else:
                                out = kk.KK().kkrzptfn(self.Shotnumber, self.time, self.Rall, self.zall, exp=eqExp, diag=eqDiag )
                                self.rhop = out.rho_p
                                del out
                                out = kk.KK().kkrhorz(self.Shotnumber, self.time, self.rhop, angle=[0.0], exp=eqExp, diag=eqDiag )
                                self.rmaj = out.r
                                self.Status_allrhop = True
                                    
					
			sf.close()
			return True


        def UnloadRz( self ):
            if (self.StatusRz) & (self.Status == False):
                self.StatusRz = False
                del self.nChannels
                del self.R
                del self.Rall
                del self.z 
                del self.zall
                del self.phi
                del self.phiall 



        def LoadRz( self ,  Shotnumber, Experiment='AUGD', Diagnostic='CEZ', Edition = 0L, tBegin=-1.0, tEnd=12.0,  eqExp = 'AUGD', eqDiag = 'EQH' , Rshift = 0.0 ):
            self.UnloadRz()
            if (Diagnostic == 'CEZ') or (Diagnostic == 'CMZ'):
                try:
                    sf = dd.shotfile( Diagnostic, Shotnumber, Experiment, Edition)
                    self.Shotnumber = Shotnumber
                except:
                    print "Error reading shotfile" 
                    return False

                self.R = sf( 'R', tBegin=tBegin, tEnd=tEnd ).data+Rshift
                self.z = sf( 'z', tBegin=tBegin, tEnd=tEnd ).data
                self.phi = sf( 'phi', tBegin=tBegin, tEnd=tEnd ).data#-67.5+360.

                idx = numpy.where(self.phi != 0.)
 
                self.R = self.R[idx]
                self.z = self.z[idx]
                self.phi = self.phi[idx] - 67.5 + 360.

                sf.close()
                self.StatusRz = True


        def getRsep( self ):
            if self.Status_allrhop:
                ntimes = self.ntime
                data = numpy.zeros((ntimes))
                data_in = [1.0]
                x = self.rhop
                y = self.Rall
               
                data = tools.interpolSimpleArray(x,y,data_in)
                self.R_sep = self.Rall - numpy.reshape(numpy.repeat(data,self.nChannels),(ntimes,self.nChannels))
                return data
                
			

	def get(self, time):

		if self.Status:		
			if (numpy.size(time)== 1) | ( numpy.size(time) > 2):
				output = CXRShelp()
				output.Ti = self.getData( timepoints = time, value='Ti')
				output.rhop = self.getData( timepoints = time, value='rhop')
				output.time = self.getData( timepoints = time, value='time')
				return output

			elif numpy.size(time)== 2:
				output = CXRShelp()
				output.Ti = self.getData( range = time, value='Ti')
				output.rhop = self.getData( range = time, value='rhop')
				output.time = self.getData( range = time, value='time')
				return output
		else:
			print 'no Data loaded'
			return None


	def getTi(self, timepoint):

		if self.Status:		
			if (numpy.size(timepoint)== 1) | ( numpy.size(timepoint) > 2):
				return self.getData( timepoints = timepoint, value='Ti')

			elif numpy.size(timepoint)== 2:
				return self.getData( range = timepoint, value='Ti')

			return None
		else:
			print 'no Data loaded'
			return None

	def getVrot(self, timepoint):

		if self.Status:		
			if (numpy.size(timepoint)== 1) | ( numpy.size(timepoint) > 2):
				return self.getData( timepoints = timepoint, value='vrot')

			elif numpy.size(timepoint)== 2:
				return self.getData( range = timepoint, value='vrot')

			return None
		else:
			print 'no Data loaded'
			return None



	def getRhop(self, timepoint):

		if self.Status:		
			if (numpy.size(timepoint)== 1) | ( numpy.size(timepoint) > 2):
				return self.getData( timepoints = timepoint, value='rhop')

			elif numpy.size(timepoint)== 2:
				return  self.getData( range = timepoint, value='rhop')

			return None
		else:
			print 'no Data loaded'
			return None

	def getTime(self, timepoint):

		if self.Status:		
			if (numpy.size(timepoint)== 1) | ( numpy.size(timepoint) > 2):
				return self.getData( timepoints = timepoint, value='time')

			elif numpy.size(timepoint)== 2:
				return  self.getData( range = timepoint, value='time')

			return None
		else:
			print 'no Data loaded'
			return None



	#Default Te, indexbased
	def getData( self , timepoints = None, range = None, value='Ti', Experiment = 'AUGD', Diagnostic='EQH' ):
		ntimes = 0

		if self.Status:

			if not  (value == 'Ti') | (value == 'rhop')| (value == 'time')| (value == 'vrot')| (value == 'inte') :
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
                            
                            if range != None:
                                print 'either range or timepoints must be given'
                                return False

                            ntimes = numpy.size(timepoints)
                            if ntimes > 1:
					#get min from the time matrix of the timepoints vector, 1. self.time is copied ntimes, then timepoints is subtracted row wise.
                                search_array = numpy.sum([numpy.array([self.time,]*ntimes),numpy.multiply(timepoints,-1.)],axis=0)
					# the indices of the closest points are calculated
                                idx = numpy.nanargmin(numpy.abs( numpy.reshape( numpy.concatenate(search_array),(ntimes,self.ntime) )),axis=1)
				
					#if (value == 'Te') | (value == 'ne') | (value == 'rhop'):
					#ne/Te/rhop should have the same area base
						#out_temp = numpy.zeros((ntimes,numpy.size(self.Te[0,:])))
                                
                            elif ntimes == 1:
                                idx = numpy.argmin( numpy.abs( self.time - timepoints ) )
                                if (value == 'time'): 
                                    return self.time[idx]

                            else:
                                print 'Something went wrong with time points'
                                IPython.embed()
                                return None
				

                        if range != None:
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

	
                        if ((value == 'Ti') | (value == 'ne') | (value == 'rhop') | (value == 'vrot')):
					#ne/Te/rhop should have the same area base
				
                            if value == 'Ti':
                                return  self.Ti[idx,:]

                            if value == 'vrot':
                                return  self.vrot[idx,:]
								

                            if value == 'inte':
                                return  self.inte[idx,:]
				
                            if (value == 'rhop'):
                                if (self.Status_allrhop):
                                    if  (Experiment == self.eqExp) & (Diagnostic == self.eqDiag):
                                        return self.rhop[idx,:]
                                    else:
                                        print 'Equlibrium Experiment differ, calculate new'
                                        #rhop_out = kk.KK().kkrzptfn(self.Shotnumber, self.time, self.R, self.z, exp=Experiment, diag=Diagnostic )
                                        output = kk.KK().kkrzptfn( self.Shotnumber, self.time, self.R[:], self.z[:], exp= Experiment, diag=Diagnostic, ed=0)
                                        self.eqDiag = Diagnostic
                                        self.eqExp = Experiment
         
                                        return output.rho_p 

                                else:
                                    nindex = numpy.size(idx)
                                    if nindex > 1: 
                                        out_temp = numpy.zeros((nindex,numpy.size(self.Ti[0,:])))
                                        output = kk.KK().kkrzptfn( self.Shotnumber, self.time, self.R, self.z, exp= Experiment, diag=Diagnostic, ed=0)
                                        out_temp = output.rho_p	
                                        return out_temp
                                    else:
                                        
                                        output = kk.KK().kkrzptfn( self.Shotnumber, self.time[idx], self.R[:], self.z[:], exp= Experiment, diag=Diagnostic, ed=0)
                                        return output.rho_p	

                                        
                return None

