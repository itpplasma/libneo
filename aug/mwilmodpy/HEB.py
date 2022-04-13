import numpy
import dd
import scipy.interpolate
import kk
import tools

#import matplotlib.pylab as plt
from IPython import embed

class HEB:
    def __init__( self ,  Experiment = 'AUGD', Shotnumber = None ):
        self.Status = False
        if Shotnumber != None:
            self.Load( Shotnumber ,Experiment = Experiment)

    def __del__( self ):
        self.Unload( )
        del self.Status

    def Load( self ,   Shotnumber, Experiment = 'MGRIEN',Diagnostic = 'HEP', Edition = 0L, tBegin=-1.0, tEnd=12.0 , Rshift=0.00,select=[0,-1]):	
        self.Unload()
        if Diagnostic == 'HEP':
            try:
                sf = dd.shotfile( Diagnostic, Shotnumber, Experiment, Edition)
            except:
                print "Error reading shotfile" 
                #embed()
                return False

            try:        
                time = sf.getTimeBase( 'ne' )
                index = numpy.where( ( time > tBegin ) & ( time < tEnd ) )
                output = sf( 'ne' )
                self.Rshift = Rshift
                self.ne =  numpy.squeeze(output.data[index,select[0]:select[1]])
                self.ne_unc = numpy.squeeze(sf( 'ne_err' ).data[index,select[0]:select[1]])
                self.time = output.time[index]
                self.ntime = numpy.size(self.time)
                self.rhop = output.area[index,select[0]:select[1]]
                del output
                self.Te = numpy.squeeze(sf( 'Te' ).data[index,select[0]:select[1]])
                self.Te_unc = numpy.squeeze(sf( 'Te_err' ).data[index,select[0]:select[1]])
                self.rhop = numpy.squeeze(sf( 'RHO_POL' ).data[index,select[0]:select[1]])
               # self.Te_green = numpy.squeeze(sf( 'Te_green' ).data)[index]
               # self.Te_comb = numpy.squeeze(sf( 'Te_comb' ).data)[index]
               # self.ne_green = numpy.squeeze(sf( 'ne_green' ).data)[index]
               # self.ne_comb = numpy.squeeze(sf( 'ne_comb' ).data)[index]
                
                self.R = numpy.squeeze(sf( 'R' ).data)[select[0]:select[1]]+Rshift
                self.nChannels = numpy.size(self.R)
                self.Rall =  numpy.reshape(numpy.tile( self.R,(self.ntime) ) ,(self.ntime,self.nChannels) )
                self.z = numpy.squeeze(sf( 'z' ).data)[select[0]:select[1]]
                self.zall = numpy.reshape(numpy.tile( self.z,(self.ntime) ) ,(self.ntime,self.nChannels) )
                self.phiCh = numpy.squeeze(sf( 'phi' ).data)[select[0]:select[1]]
                self.phi = numpy.squeeze(sf( 'phi' ).data).mean()
            except:
                print "Loading data failed" 
                embed()             
            
 #260.64
		#	self.phiRad = self.phi*numpy.pi/180.
			#self.lib_dat=numpy.squeeze(sf( 'lib_dat' ).data)[index]
            self.Shotnumber = Shotnumber
            self.tBegin = tBegin
            self.tEnd = tEnd
            self.Status = True
            sf.close()
            return True

    def Unload( self ):
        if self.Status:
            self.Status = False
            del self.time
            del self.ntime
            del self.rhop
            del self.ne
            del self.ne_unc
            del self.Te
            del self.Te_unc
 #           del self.Te_green 
 #           del self.Te_comb 
 #           del self.ne_green
 #           del self.ne_comb 


			#del self.lib_dat
            del self.nChannels
            del self.R
            del self.Rall
            del self.z
            del self.Shotnumber
                        


    def getne( self , timepoints ):
        if self.Status:    			
            return  self.getData(timepoints, selectData='ne' )

    def getTe( self , timepoints ):
        if self.Status:    			
            return  self.getData(timepoints, selectData='Te' )

    def getnewiUnc( self , timepoints ):
        if self.Status:    			
            return  self.getData(timepoints, selectData='ne' ),self.getData(timepoints, selectData='ne_unc' )

    def getTewiUnc( self , timepoints ):
        if self.Status:    			
            return  self.getData(timepoints, selectData='Te' ),self.getData(timepoints, selectData='Te_unc' )

    def getRhop( self , timepoints ):
        if self.Status:    			
            return  self.getData(timepoints, selectData='rhop' )

    def getRsep( self , timepoints ):
        if self.Status:    			
            return  self.getTimetrace(timepoints, selectData_in = 'rhop', data_in = [1.0], selectData_out = 'R' )


    def getR( self , timepoints ):
        if self.Status:    			
            return  self.getData(timepoints, selectData='R' )

    def getz( self , timepoints ):
        if self.Status:    			
            return  self.getData(timepoints, selectData='z' )

    def getData( self , timepoints, selectData='ne' ):
        ntimes = 0
        if self.Status:
            ntimes = numpy.size(timepoints)
            data = numpy.zeros( (ntimes,self.nChannels)) 
            if ntimes > 1:
                for i in range(ntimes):
                    
                    idx = numpy.argmin( numpy.abs( self.time - timepoints[i] ) )

                    if selectData == 'ne':
                        data[i,:] = self.ne[idx,:]
                    if selectData == 'Te':
                        data[i,:] = self.Te[idx,:]
                    if selectData == 'neunc':
                        data[i,:] = self.ne_unc[idx,:]
                    if selectData == 'Te_unc':
                        data[i,:] = self.Te_unc[idx,:]
                    if selectData == 'R':
                        data[i,:] = self.Rall[idx,:]
                    if selectData == 'z':
                        data[i,:] = self.zall[idx,:]
                    if selectData == 'rhop':
                        data[i,:] = self.rhop[idx,:]				
            else:
                idx = numpy.argmin(numpy.abs(self.time-timepoints))
                if selectData == 'ne':
                    data = self.ne[idx,:]
                if selectData == 'Te':
                    data = self.Te[idx,:]
                if selectData == 'ne_unc':
                    data = self.ne[idx,:]
                if selectData == 'Te_unc':
                    data = self.Te[idx,:]
                if selectData == 'R':
                    data = self.Rall[idx,:]
                if selectData == 'z':
                    data = self.zall[idx,:]
                if selectData == 'rhop':
                    data = self.rhop[idx,:]	
			     			
            return data 


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
                ntimes = numpy.size(time)
				#get the time points from the get routines
                ne_out = self.getne( time )
                rhop_out = self.getRhop( time )
                ntimes = numpy.size(time)
                if ntimes > 1:
                    ne_temp = numpy.zeros((ntimes,numpy.size(rhop)))
                    for i in numpy.arange(ntimes) :						
                        if rhop_out[i,1]<= rhop_out[i,0]:
                            ne = scipy.interpolate.interp1d(rhop_out[i,::-1] , ne_out[i,::-1],  bounds_error=False, fill_value = float('nan'))
                        else:
                            ne = scipy.interpolate.interp1d(rhop_out[i,:] , ne_out[i,:1],  bounds_error=False, fill_value = float('nan'))
                            ne_temp[i,:] = ne( rhop )
                            return ne_temp
                else:
                    if rhop_out[i,1] <= rhop_out[i,0]:
                        ne = scipy.interpolate.interp1d(rhop_out[i,::-1] , ne_out[i,::-1],  bounds_error=False, fill_value = float('nan'))
                    else:
                        ne = scipy.interpolate.interp1d(rhop_out[i,:] , ne_out[i,:1],  bounds_error=False, fill_value = float('nan'))
                        return ne( rhop )					

        else:
            print 'no Data loaded'
            return None

