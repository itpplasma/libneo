import matplotlib.pylab as plt   # For plotting graphs.
import os                         # For issuing commands to the OS.
import numpy as np


class MKMOVIE:

  def __init__(self,x,y,psym='b-',xlab='X',ylab='Y',title='Plot'):

# Plot

    plt.ion()
    movfig = plt.figure(1,figsize=(9,7))
    plt.clf()
    movsub = movfig.add_subplot(111)
    movfig.subplots_adjust(left=0.15,bottom=0.12,right=0.95,top=0.9)
    mov, = movsub.plot(x[0],y[0],psym)
    movsub.set_xlabel(xlab)
    movsub.set_ylabel(ylab)
    movsub.set_title(title)
    for i in range(1,len(y)) :
      mov.set_xdata(x[i])
      mov.set_ydata(y[i])
      plt.draw()
      filename = str('%03d' % i) + '.png'
      plt.savefig(filename, dpi=100)
      print 'Wrote file', filename

    command = ('mencoder',
               'mf://*.png',
               '-mf',
               'type=png:w=800:h=600:fps=25',
               '-ovc',
               'lavc',
               '-lavcopts',
               'vcodec=mpeg4',
               '-oac',
               'copy',
               '-o',
               'output.avi')

    encoder='/afs/@cell/common/soft/visualization/mencoder/svn-2011-05-17/@sys/bin/mencoder'
    os.spawnvp(os.P_WAIT, encoder, command)

def __main__():

  print 'Executing on', os.uname()
  print 'Initializing data set...'   # Let the user know what's happening.

# Initialize variables needed to create and store the example data set.
  nt = 100   # Number of frames we want in the movie.
  x = np.arange(-10,10,0.01)   # Values to be plotted on the x-axis.
  mean = -6                 # Initial mean of the Gaussian.
  stddev = 0.2              # Initial standard deviation.
  meaninc = 0.1             # Mean increment.
  stddevinc = 0.1           # Standard deviation increment.

# Create an array of zeros and fill it with the example data.
  y = np.zeros((nt,len(x)))  
  x = np.repeat([x],nt,axis=0)
  print x.shape, y.shape
  for i in range(nt) :
    y[i] = (1/np.sqrt(2*np.pi*stddev))*np.exp(-((x[i]-mean)**2)/(2*stddev))
    mean = mean + meaninc
    stddev = stddev + stddevinc

  MKMOVIE(x,y)

if __name__ == "__main__":
  __main__()
