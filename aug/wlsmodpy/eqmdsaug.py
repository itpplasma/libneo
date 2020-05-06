#!/usr/bin/env python
#
# -----------
# eqmdsaug.py
# -----------
# read AUG equilibrium using MDSplus

import sys

sys.path.append('/afs/ipp/home/w/wls/.local/lib/python2.7/site-packages/pyPmds')

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colrs
import pmds as mds
#/afs/ipp/home/w/wls/.local/lib/python2.7/site-packages/pyPmds
from scipy.interpolate import interp1d


from eqrzfd import eqrzfd, geqdsk


# ...................................................
class eqmdsaug(eqrzfd):
  """
     Read ASDEX Upgrade equilibrium using MDSplus.
     If a remote connection is to be used it must be open already.
  """
  
  def __init__ (self, shot, time, experiment='AUGD', diagnostic='EQI', edition=0):
    s1 = 'augdiag('+str(shot)+',"'+diagnostic+'","'
    s2 = '","'+experiment+'",'+str(edition)+','+str(time)+','+str(time)+')'
    
    # read flux matrix
    sh = mds.mdsvalue('_s=shape(_psi='+s1+'PFM'+s2+')')
    nr = sh[1]
    nz = sh[2]
    dummy = mds.mdsvalue('_p0=zero(_s[1:2])')  # trick for 2D result
    R = mds.mdsvalue('_r=dim_of(_psi,1)')
    z = mds.mdsvalue('_z=dim_of(_psi,2)')

    # number of internal flux labels
    Lpf = mds.mdsvalue('_lpf='+s1+'Lpf'+s2)
    npsi = 1 + (int(Lpf[0][0]) % 10000)
    npso = 1 + (int(Lpf[0][0]) / 10000)

    # set up arrays and geometry coefficients
    # super(eqmdsaug, self).__init__(nr, nz, npsi+npso, R[0, 0], R[-1, 0], z[0, 0], z[-1, 0])
    eqrzfd.__init__(self, nr, nz, npsi+npso, R[0, 0], R[-1, 0], z[0, 0], z[-1, 0],
                    descriptor = experiment+':'+diagnostic+'('+str(edition)+')',
                    shot = shot, time = time)

    # flux matrix, convert to Vs/radians
    self.PsiRz = mds.mdsvalue('_p0+_psi').T / (2*np.pi)
    
    # index arrays to sort out 1D flux profiles
    inner = np.arange(npsi-1,-1,-1)
    even_inner = np.arange(npsi+npsi-2,-1,-2)
    odd_inner  = np.arange(npsi+npsi-1,-1,-2)

    outer = np.arange(npsi+2,npsi+npso+1,1)
    even_outer = np.arange(npsi+npsi+4,2*(npsi+npso+1),2)
    odd_outer = np.arange(npsi+npsi+5,2*(npsi+npso+1)+1,2)
    
    # flux quantities (1D profiles)
    Qpsi = mds.mdsvalue('_qpsi='+s1+'Qpsi'+s2)
    Psi  = mds.mdsvalue('dim_of(_qpsi,1)')  / (2*np.pi)

    # put 1D profiles into order:  mag.axis..separatrix..SOL
    self.Psi = np.append(Psi[inner,0],Psi[outer,0])

    # do not expect Qpsi on open flux surfaces
    self.Qpsi = np.append(Qpsi[inner,0], np.zeros(len(outer)))

    # p/p' and J/J' are stored in a funny way in AUG shot files:
    # even elements contain value, odd elements contain derivative wrt psi
    Pres = mds.mdsvalue(s1+'Pres'+s2)
    self.Pres   = np.append(Pres[even_inner,0],Pres[even_outer,0])
    self.pprime = np.append(Pres[odd_inner,0], Pres[odd_outer,0]) * 2*np.pi

    Fpol = 2.0e-7 * mds.mdsvalue(s1+'Jpol'+s2)
    self.Fpol = np.append(Fpol[even_inner,0], Fpol[even_outer,0]) 
    self.FFprime = np.append(Fpol[odd_inner,0],Fpol[odd_outer,0]) * self.Fpol * 2*np.pi

    # Psimag, Psibdry  -- use libkk
    cc = 'augconv('+str(shot)+',"'+diagnostic+'","'+experiment+'",'+str(edition)+','+str(time)
    print 'cc = ', cc+',28,[0,1],,_pflux)'
    dummy = mds.mdsvalue(cc+',28,[0,1],,_pflux)')
    pflux = mds.mdsvalue('_pflux') / (2*np.pi)
    self.Psimag = pflux[0][0]
    self.Psibdry = pflux[1][0]

    # Rmaxis, zmaxis
    raxsepo = mds.mdsvalue(cc+',30,[0,1],0.0,_zaxsepo)')
    zaxsepo = mds.mdsvalue('_zaxsepo')
    self.Rmaxis = raxsepo[0][0]
    self.zmaxis = zaxsepo[0][0]

    # Limiter contour
    Rlim, zlim = read_aug_pfc('aug/augiii_vessel.asc')
    self.set_limiter (Rlim, zlim)

    # Plasma current
    jplasma = self.ToroidalCurrentDensity()
    self.Iplasma = np.sum(jplasma * self.dR * self.dz)
    
    # plasma boundary (try how close to the separatrix we can sit and still get a closed contour)
    bfluxes = [0.99999, 0.99993, 0.9999, 0.9993, 0.999, 0.993, 0.99, 0.93, 0.9]
    for bflux in bfluxes:
      c = self.fluxcontours(self.Psimag + np.array([bflux])*(self.Psibdry-self.Psimag))
      if c<>[]:
        print "Found closed boundary at ", bflux*100.0, "% flux"
        self.Rbdry = c[0][:,0]
        self.zbdry = c[0][:,1]
        break

# end of eqmdsaug.__init__()
    


def read_aug_pfc(fn):
  """
  read ASDEX Upgrade pfc coordinates from an ascii file
  """
  try:
    f = open(fn, 'r')
  except:
    print "Cannot open file ", fn
    return None, None
  Rv = []
  zv = []
  for line in f:
    columns = line.split()
    if len(columns)>0:
      Rv.append(float(columns[0]))
      zv.append(float(columns[1]))
  f.close()
  return Rv, zv



# ---------------------------------------------------------------------
# demo


if __name__ == '__main__':
# if True:

  mds.mdsconnect('mdsplus.aug.ipp.mpg.de')
#  mds.mdsconnect('sxaug31:8001')
#  mds.mdsconnect('gate2:8001')

  
  if False:
    eq = geqdsk("testdata/eqdsk/d3d/g164362.04545_828")
    ofn = 'g164362.04545_test'

  if False:
    eq = eqmdsaug(33353, 4.25, experiment='AUGD', diagnostic='EQI', edition=0)
    ofn = 'testdata/eqdsk/aug/g33353_4250_eqi'
    
  if False:
    eq = eqmdsaug(34214, 3.2, experiment='AUGE', diagnostic='IDE', edition=2)
    ofn = 'g34214_3200_IDE2_Sauter'

  if False:
    eq = eqmdsaug(34214, 3.355, experiment='MICDU', diagnostic='EQI', edition=1)
    ofn = 'g34214_3355_EQI'

  if False:
    eq = eqmdsaug(34548, 1.0, experiment='MICDU', diagnostic='EQB', edition=0)
    ofn = 'testdata/eqdsk/aug/g34548.1000_EQB'

  if False:
    eq = eqmdsaug(31529, 4.700, experiment='MICDU', diagnostic='EQI', edition=3)
    ofn = 'g31529_4700_MICDU_EQI_Ed3'

  if False:
    eq = eqmdsaug(33120, 5.5, experiment='MICDU', diagnostic='EQB', edition=1)
    ofn = 'g33120_5500_MICDU_EQI_Ed1'

  if False:
    eq = eqmdsaug(33120, 5.635, experiment='MICDU', diagnostic='EQB', edition=1)
    ofn = 'g33120_5635_MICDU_EQI_Ed1'

  if False:
    eq = eqmdsaug(34398, 3.55, experiment='AUGD', diagnostic='IDE', edition=0)
    ofn = 'g34398_3550_IDE'

  if False:
    eq = eqmdsaug(34398, 4.10, experiment='AUGD', diagnostic='IDE', edition=1)
    ofn = 'g34398_4100_IDE'

  if False:
    eq = eqmdsaug(37438, 2.5, experiment='AUGD', diagnostic='EQH', edition=0)
    ofn = 'g37438_2500_EQH'

  if True:
    eq = eqmdsaug(37438, 3.8, experiment='AUGD', diagnostic='EQH', edition=0)
    ofn = 'g37438_3800_EQH'

  if False:
    eq = eqmdsaug(36285, 3.0, experiment='AUGD', diagnostic='EQH', edition=0)
    ofn = 'g36285_3000_EQH'

    
    
#  mds.mdsdisconnect()
    
    
  if False:
    plt.axis("equal")
    eq.plot_Rz(color='blue')

    psin = np.array([0.9999])
    c = eq.fluxcontours(eq.Psimag + psin*(eq.Psibdry-eq.Psimag))
    for c1 in c:
      plt.plot(c1[:,0], c1[:,1], color='red')  

    plt.show()
    
  print "Toroidal plasma current: ", eq.Iplasma, " A"

  # plot current density
  if False:
    jplasma = eq.ToroidalCurrentDensity()
    
    plt.axis("equal")
    levels = np.linspace(np.amin(jplasma), np.amax(jplasma), 300)
    plt.contour(eq.R[1:-1], eq.z[1:-1], jplasma.T, levels=levels)
    plt.plot(eq.Rlim, eq.zlim)
    plt.show()
  
    
# write geqdsk file
  eq.write_geqdsk(ofn)
