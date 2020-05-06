#!/usr/bin/env python
#
# --------
# rprof.py
# --------
# radial profiles vs. normalised polodal flux (1 time point each)
# with optional error bars

import numpy as np
import matplotlib.pyplot as plt
import pickle
from netCDF4 import Dataset

# -------- generic profile class --------------------

class RadialProfile:

# manually preset a profile
  def __init__ (self, yunit, value, errorbar, xunit, psinorm):
    self.psinorm = psinorm
    self.value = value
    self.xunit = xunit
    self.yunit = yunit
    self.errorbar = errorbar

# plot a profile, solid line by default
  def plot(self, scale, nonzero=True, label=None, nearestpsi=None,
           color='black', linestyle='solid', linewidth=2, marker='None', markersize=4):
    x = self.psinorm
    y = self.value
    if nearestpsi is not None:     # pick only data nearest to given positions
      ix = np.digitize(nearestpsi, x)
      ix = np.unique(ix)
      x = x[ix]
      y = y[ix]
    plt.plot(x, y/scale, label=label,
             color=color, linestyle=linestyle, linewidth=linewidth,
             marker=marker, markersize=markersize)


# plot a profile with error bars, markers by default
  def plot_errors (self, scale, label=None, nearestpsi=None,
           color='black', linestyle = 'None', linewidth=2, marker='x', markersize=4):

    if self.errorbar==[]:
      return
    x = self.psinorm
    y = self.value
    errbar = self.errorbar
# sort data in ascending order
    ix = np.argsort(x)
    x = x[ix]
    y = y[ix]
    errbar = errbar[ix]
# option: pick only data nearest to given positions
    if nearestpsi is not None:
      ix = np.digitize(nearestpsi, x)
      ix = np.unique(ix)
      wh = (ix<len(x)-1)
      ix = ix[wh]
      x = x[ix]
      y = y[ix]
      errbar = errbar[ix]
    plt.errorbar(x, y/scale, yerr=errbar/scale, 
         label=label, linestyle=linestyle, color=color, marker=marker, markersize=markersize)
      

# dump instance to (open) pickle file 'output'
  def dump_pickle(self, output):
    pickle.dump(self, output, -1)


# ------- One profile, by concatenation of two others ----------

class Concatenate (RadialProfile):

  def __init__ (self, rprof1, rprof2):
    psinorm  = np.append(rprof1.psinorm, rprof2.psinorm)
    value    = np.append(rprof1.value, rprof2.value)
    errorbar = np.append(rprof1.errorbar, rprof2.errorbar)
    ix = np.argsort(psinorm)
    self.psinorm = psinorm[ix]
    self.value = value[ix]
    try:
      self.errorbar = errorbar[ix]
    except:
      pass
#      self.errorbar = None
    self.xunit = rprof1.xunit
    self.yunit = rprof1.yunit


# ------- One profile out of a peqdsk file ----------

class peqdsk (RadialProfile):

# read from p-eqdsk file, handle f must be open
  def __init__ (self, f):
      s = f.readline()
      sn = s.split()
      npsi = int(sn[0])
      self.psinorm = np.zeros((npsi))
      self.value = np.zeros((npsi))
      self.errorbar = []
      self.xunit = sn[1]
      self.yunit = sn[2]
      for i in range(npsi):
          s = f.readline()
          sn = s.split()
          self.psinorm[i]=float(sn[0])
          self.value[i]=float(sn[1])
          # derivative ignored for the moment
          # self.deriv[i]=float(sn[2])    


# -----------------------------------------------------------
# read entire peqdsk file and return a dictionary of profiles

def read_peqdsk(fn):
  profiles = {}    # dictionary
  try:
    f = open(fn, 'r')
    while True:
      try:
        p = peqdsk(f)
        quan =  (p.yunit.split('('))[0]
        profiles[quan] = p
      except:
        break
    f.close()
  except:
    pass
  return(profiles)


# -----------------------------------------------------------
# read entire pickle file (assuming it contains only profiles)

def read_pickle(fn):
  profiles = {}    # dictionary
  try:
    f = open(fn, 'rb')
    while True:
      try:
        p = pickle.load(f)
        quan =  (p.yunit.split('('))[0]
        print('Reading quantity: %s' % quan)
        profiles[quan] = p
      except:
        break
    f.close()
  except:
    return None
  return(profiles)


# -----------------------------------------------------------
# write a dictionary of profiles 'pd' to netcdf file 'fn'

def write_netcdf(fn, pd, comment='created by rprof.py', title='no title', version=1):

  try:
    file = Dataset(fn, 'w', format='NETCDF4')
    file.description = comment
#    file = NetCDFFile(fn, 'w', comment)
  except:
    print("Could not open netcdf file: %s" % fn)
    return None

  file.title = title
  file.version = version

  for key, p in pd.iteritems():
    file.createDimension('dim_'+key, p.psinorm.size)

    value = file.createVariable(key, 'f8', ('dim_'+key,))
    value[:] = p.value
    value.units = p.yunit

    psinorm = file.createVariable(key+'_psinorm', 'f8', ('dim_'+key,))
    psinorm[:] = p.psinorm
    psinorm.units = p.xunit

    errorbar = file.createVariable(key+'_error', 'f8', ('dim_'+key,))
    errorbar[:] = p.errorbar
    errorbar.units = p.yunit

  file.close()



# ------- test main program -------------------------

# if True:
if __name__ == '__main__':  

#  d3d_profs = read_peqdsk('../data/d3d/eqdsk/p164362.03400')
#  p = d3d_profs['ne']
#  p.plot(0.1, color='red', label='DIII_D 164277 2500 ms')

  aug_profs = read_pickle('/afs/ipp-garching.mpg.de/u/wls/public/Modeling_kit/AUG/prof/pickle/Profiles_33353_2.9.pkl')
  p = aug_profs['neData']
  p.plot_errors (1e19, color='blue')

  p = aug_profs['neFit']
  p.plot (1e19, color='blue', label='AUG 33353 t=2.325s')

  plt.xlabel(p.xunit)
  plt.ylabel(p.yunit)
  plt.legend(loc='lower left')
  plt.show()

