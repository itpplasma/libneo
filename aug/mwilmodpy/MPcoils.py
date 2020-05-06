import numpy
import kk_mwillens as kk
import scipy.interpolate
from IPython import embed
import eqi_map as fastkk
from StringIO import StringIO
# structure nRows,nCoils,...
#import matplotlib.pylab as plt
#from mpl_toolkits.mplot3d import Axes3D
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')


class MPcoilhelp:
    status = False

class MPcoils:

#reading 
    def __init__( self , nRows = 2, nCoils = 8, labelRows = ['l','u']):
        self.Status = False
        
        ## generate Stringarray
        if nRows != numpy.size(labelRows):
            print 'number of Rows do not have the same'
            return False

        self.nRows =  nRows
        self.nCoils =  nCoils
        # generate Stringarray to read coil stuff#
        coilStr = numpy.chararray((nRows, nCoils),itemsize=2)
        for i in range(nRows):
            for j in range(nCoils):
                coilStr[i,j] = labelRows[i]+str(j+1) 
            
	self.coilStr = coilStr

    def __del__( self ):
        del self.Status

        
    def LoadCoord( self, path='/afs/ipp-garching.mpg.de/home/m/mwillens/Programme/python/data/Bcoils/', preFile='B',postFile='.asc'):
        
#generate file to read
        fileCoilStr = path + preFile +  self.coilStr + postFile
        
#X=hallo[:,0]*numpy.cos(hallo[:,2]/2./numpy.pi)
#Y=hallo[:,0]*numpy.sin(hallo[:,2]/2./numpy.pi)
        #coordPol = [] #numpy.array(numpy.shape(fileCoilStr))
        #coordCar = [] #numpy.array(numpy.shape(fileCoilStr))

        tmpRows=[]
        for i in range(self.nRows):
            tmpCoils=[]
            for j in range(self.nCoils):
                #R[m] z[m] phi[degree]
                tmpCoils.append(numpy.genfromtxt(fileCoilStr[i,j]))
                #ax.plot3D(coordCar[i,j,:,0],coordCar[i,j,:,1],coordCar[i,j,:,2])
            tmpRows.append(tmpCoils)

        coordPol =  numpy.squeeze(tmpRows)
         #X[m] Y[m] Z[m]
        coordCar = numpy.zeros_like(coordPol)
        coordCar[:,:,:,0] = coordPol[:,:,:,0] * numpy.cos(coordPol[:,:,:,2] )
        coordCar[:,:,:,1] = coordPol[:,:,:,0] * numpy.sin(coordPol[:,:,:,2] )
        coordCar[:,:,:,2] = coordPol[:,:,:,1] 

        self.coordPol = coordPol
        self.coordCar = coordCar

  #  ax.plot3D(coordCar[i,j,:,0],coordCar[i,j,:,1],coordCar[i,j,:,2])


