"""

    Routine data binning. Data must be eqidistant


"""

__author__='Matthias Willensdorfer (mwillens@ipp.mpg.de)'
__date__='20 Nov 2014'
__version__='1.0'


#from IPython import embed
#import matplotlib.pylab as plt
from numpy import *
import numpy as np
from IPython import embed

# return new timebase and Data
def dataBinning( time, data, samplefreq = 1.0 ):			
            						      
    print "binning with ",samplefreq," kHz"
    ntimes= size(time)
    samplingrate = 1.0/mean(diff(time))
    dataShape =array(shape(data))  
    #get the time index
    idxOfTime = squeeze(where(dataShape == ntimes))
    # if more index with the number of times exists, take the first one
    if size(idxOfTime) > 1:
        idxOfTime = idxOfTime[0]

    bins = int(ntimes*(float(samplefreq)*1.0e3/samplingrate))

    slices= linspace(0, ntimes, bins+1, True).astype(int)
    counts = diff(slices)

    #calculate new timebase
    newTime = add.reduceat(time, slices[:-1]) / counts
    newNtimes = size(newTime)

    #create new shape
    newDataShape = dataShape
    #replace old shape
    put(newDataShape, idxOfTime, newNtimes)
    #create new Data array
    newData = zeros( (newDataShape)  )

    #simplify array such as the first index is always the timebase
    newData = swapaxes(newData,0,idxOfTime)
    data = swapaxes( data,0,idxOfTime )

    storeShape = shape( newData )

    # rehape data to two dimensions
    data = reshape(data,(ntimes,size(data)/ntimes))
    newData = reshape(newData,(newNtimes,size(newData)/newNtimes))

    for i in range(shape(data)[1]):
        newData[:,i] = add.reduceat(data[:,i], slices[:-1]) / counts

#shape back
    newData = reshape(newData,(storeShape))
    #swap back to original shape
    newData = swapaxes(newData,0,idxOfTime)

    return newTime,newData
       
def rotatexy(xy, phiIn = 67.5, phiUnit = 'degree'):
 
    if phiUnit == 'degree':
        fac = pi/180.
    else:
        fac = 1.0

    return squeeze([xy[0]*sin(phiIn*fac)+xy[1]*cos(phiIn*fac),-xy[0]*cos(phiIn*fac)+xy[1]*sin(phiIn*fac)])   


def rzphi2xyz(rzphi, phiUnit = 'rad'):
    
    if phiUnit == 'degree':
        fac = pi/180.
    else:
        fac = 1.0

    dataShape =array(shape(rzphi)) 
    idxOfSpace = squeeze(where(dataShape == 3))

    if size(idxOfSpace) > 1:
        idxOfSpace = idxOfSpace[0]

    if size(idxOfSpace) < 1:
        return False

    if size(dataShape)>1:
         #simplify array such as the first index is always the timebase
        rzphi = swapaxes(rzphi,0,idxOfSpace)
        # store shape
        storeShape = shape( rzphi )
        # rehape data to two dimensions
        data = reshape(rzphi,(3,size(rzphi)/3.))
        output = zeros_like(data)
        output[0,:] = data[0,:]*cos(data[2,:]*fac)
        output[1,:] = data[0,:]*sin(data[2,:]*fac)
        output[2,:] = data[1,:]
        output = reshape(output,(storeShape) )
        output = swapaxes(output,0,idxOfSpace)
        
        return output

        
    elif (size(dataShape) == 1):
        return squeeze([rzphi[0]*cos(rzphi[2]*fac),
                        rzphi[0]*sin(rzphi[2]*fac),
                        rzphi[1]
                        ])
def Rotate3D(alpha, beta, gamma):

    import numpy as np
    ones = np.ones_like(alpha)
    zeros = np.zeros_like(alpha)

    # Trig.
    sin_a = np.sin(alpha)
    sin_b = np.sin(beta)
    sin_g = np.sin(gamma)

    cos_a = np.cos(alpha)
    cos_b = np.cos(beta)
    cos_g = np.cos(gamma)

    ## xachse
    ##  1   0    0
    ##  0   co   -si
    ##  0   si   co
    Ralpha=np.zeros((3,3,alpha.size))
    Ralpha[0,:]=[ ones, zeros,zeros]
    Ralpha[1,:]= [zeros,cos_a,-sin_a]
    Ralpha[2,:]= [zeros,sin_a,cos_a]
    ## yachse
    ##  co   0   si
    ##  0    1   0
    ##  -si  0   co
    Rbeta=np.zeros((3,3,alpha.size))
    Rbeta[0,:]= [cos_b,zeros,sin_b]
    Rbeta[1,:]= [zeros,ones,zeros]
    Rbeta[2,:]=  [-sin_b,zeros,cos_b]  
    ## zachse
    ##  co   -si  0
    ##  si   co   0
    ##  0    0    1
    Rgamma=np.zeros((3,3,alpha.size))
    Rgamma[0,:]= [cos_g,-sin_g,zeros]
    Rgamma[1,:]= [sin_g,cos_g,zeros]
    Rgamma[2,:]=  [zeros,zeros,ones]         

    return np.einsum('ijf,jkf,klf->ilf',Ralpha,Rbeta,Rgamma)



def makewire(xIn,yIn,zIn,r=0.01,nPhi=100):
    xIn = np.ravel(xIn)
    yIn = np.ravel(yIn)
    zIn = np.ravel(zIn)
    nSize = xIn.size
    
    xOut = np.zeros((nSize,nPhi))
    yOut = np.zeros_like(xOut)
    zOut = np.zeros_like(xOut)
    phiRun=np.linspace(0.,2*np.pi,nPhi)
    xDiff=np.gradient(xIn)
    yDiff=np.gradient(yIn)
    zDiff=np.gradient(zIn)
    I=np.sqrt(xDiff*xDiff+yDiff*yDiff+zDiff*zDiff)
    vecIn=np.array([xDiff/I,yDiff/I,zDiff/I])
    xvec=[np.ones_like(xDiff),np.zeros_like(xDiff),np.zeros_like(zDiff)]
    yvec=[np.zeros_like(xDiff),np.ones_like(xDiff),np.zeros_like(zDiff)]
    zvec=[np.zeros_like(xDiff),np.zeros_like(xDiff),np.ones_like(zDiff)]

    
    alpha = np.arccos(np.sum(vecIn*yvec,axis=0))
    beta = np.arccos(np.sum(vecIn*zvec,axis=0))
    gamma = np.arccos(np.sum(vecIn*xvec,axis=0)) 

    xrad = np.array([r*np.cos(phiRun),]*nSize)
    yrad = np.array([r*np.sin(phiRun),]*nSize)
    zrad = np.zeros_like(xrad)
    rad = np.array([xrad,yrad,zrad])
    rotMatrix=Rotate3D(alpha, beta, gamma)
    outRot = np.einsum('ijk,jkf->ikf',rotMatrix,rad)

    xOut[:],yOut[:],zOut[:] = (outRot[0].T+xIn).T,(outRot[1].T+yIn).T,(outRot[2].T+zIn).T

    return xOut,yOut,zOut
    



def interpolSimpleArray( x, y, x_in ):

    x =  squeeze(x)
    y = squeeze(y)
    x_in = squeeze(x_in)
    nx = shape(x)[1] 
    ny = shape(y)[1]
    nt = shape(y)[0]
    y_out = zeros_like(repeat(x_in,nt))
    y_out[:] = y[:,0]

    if size( x_in ) != 1 :
        print 'input must have index 1'
        return float('NaN')
        

    if ( nx == ny ) :

        
        idx = nanargmin( abs( x - x_in ), axis=1 )      
        idxOut = where( (idx > 1 ) & (idx < nx - 2)) 
        
        runIdx = arange(nt)
        #embed()
        if size(idxOut) > 0: 
            idxOut = squeeze( idxOut )
            runIdx = runIdx[idxOut]
        
 
        idxEqual = where( (x[runIdx,idx[runIdx]] == x[runIdx,idx[runIdx]+1])  ) 
 
        if size(idxEqual) > 0: 
            idxEqual = squeeze( idxEqual )
            idx[idxEqual] = idx[idxEqual] + 1
            
        idx1 = zeros_like(idx) 
        idx2 = zeros_like(idx) 
        idxPos = where( (x[ runIdx, idx[runIdx] ] < x_in ) & ( x_in < x[ runIdx, idx[runIdx]+1 ] ) ) 

        if size(idxPos) > 0: 
            idxPos = squeeze( idxPos )
            idx2[idxPos] = idx[idxPos] + 1
            idx1[idxPos] = idx[idxPos]


        
        idxNeg =  where( (x[ runIdx, idx[runIdx]-1 ] < x_in ) & ( x_in < x[ runIdx, idx[runIdx] ] ) ) 
        if size(idxNeg) > 0:
            idxNeg = squeeze( idxNeg )
            idx2[idxNeg] = idx[idxNeg]
            idx1[idxNeg] = idx[idxNeg] - 1

 
        runIdx = runIdx[ where(idx1[runIdx] != 0) ]
        y_out[runIdx] = ( y[runIdx,idx2[runIdx]] - y[runIdx,idx1[runIdx]] ) / ( x[runIdx,idx2[runIdx]] - x[runIdx,idx1[runIdx]]  ) * ( x_in -  x[runIdx,idx1[runIdx]] ) + y[runIdx,idx1[runIdx]] 
        
        return y_out

    else:

        print 'x and y input must equally sized'
        return float('NaN')









def interpolSimple( x, y, x_in ):

    x =  squeeze(x)
    y = squeeze(y)
    x_in = squeeze(x_in)
    nx = size(x) 
    ny = size(y)
 
    if size( x_in ) != 1 :
        print 'input must have index 1'
        return float('NaN')
        

    if ( nx == ny ) :
        idx1 = nanargmin( abs( x - x_in ) )

        
        if ( idx1 < 1 ) | (idx1 >= nx - 2) :
         #   print 'input must have index'
          #  embed()   
            return float('NaN')

        if x[ idx1 ] == x[ idx1+1 ]:
            idx1 = idx1 + 1

#is x in lower or higher
        if ( x[ idx1 ] < x_in ) & ( x_in < x[ idx1+1 ] ) :
            idx2 = idx1 + 1
            idx1 = idx1 
        elif ( x[idx1 -1 ] < x_in ) & ( x_in < x[ idx1 ] ) :
            idx2 = idx1 
            idx1 = idx1 - 1
        else:
            return y[idx1]
           ### embed()
         #   return float('NaN')
            
        y_out = ( y[idx2] - y[idx1] ) / ( x[idx2] - x[idx1]  ) * ( x_in -  x[idx1] ) + y[idx1] 
        if y_out==0.0 :
            print 'y_out'
            embed()
        return y_out

    #take the next inde
    else:

        print 'x and y input must equally sized'
        return float('NaN')
