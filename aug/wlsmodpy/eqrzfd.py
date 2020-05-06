#!/usr/bin/env python
#
# ---------
# eqrzfd.py
# ---------
# 2-D ideal MHD force equilibrium (vs. R,z), finite differences with rectangular grid

import sys
import re   # regular expressions for g-file descriptor line
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as colrs
from scipy.interpolate import interp1d

from contourc import contourc
from matrixutil import bresenham_line, floodfill


# ---------------------------------------------------
# aux routines
# ---------------------------------------------------

def read_array_n5 (f, n, d):
  """Read n floating point values into 1D array d from file f.
     Expect 5 floats per line.
  """
  col=6
  for i in range(0,n):
    if col>5:
      s = f.readline()
      col=1
    o = (col-1)*16
    d[i] = float(s[o+0:o+16])
    col+=1
  return

# ...................................................

def write_array_n5 (f, d, n=5):
  """Write 1D array d to file f.
     Write n floats per line (default: 5).
  """
  col=0
  for v in d:
    f.write('{: 16.9e}'.format(v))
    col+=1
    if col>=n:
      s = f.write('\n')
      col=0
  if col>0:
    f.write('\n')

# ...................................................

def write_interpol_profile(f, xn, x, y):
  fint = interp1d(x, y, kind='cubic')
  y = fint(xn)
  write_array_n5 (f, y)

    
# ---------------------------------------------------
# EQRZFD
# ---------------------------------------------------

class eqrzfd:
  """Generic class for 2D ideal MHD equilibrium on a finite difference rectangular (R,z) grid."""

# ...................................................
  def __init__ (self, nr, nz, npsi, R1, Rm, z1, zn,
                descriptor=None, date=None, shot=None, time=None ):
    
    """
    Create empty 2D finite difference equilibrium with given grid sizes:
    nr    number of radial grid points (nr>2)
    nz    number of vertical grid points  (nz>2)
    npsi  number of grid points for 1D flux functions
    R1    radius of innermost / first [0] radial grid point
    Rm    radius of outermost / last [nr-1] radial grid point
    z1    vertical position of lowest / first [0] vertical grid point
    zn    vertical position of highest / last [nz-1] vertical grid point
    optional annotation:
      descriptor, date, shot, time (in s)
    """
    self.R = np.linspace(R1, Rm, nr)
    self.z = np.linspace(z1, zn, nz)
    self.PsiRz = np.zeros((nr,nz))
    self.Psi = np.zeros((npsi))
    self.Fpol = np.zeros((npsi))
    self.Pres = np.zeros((npsi))
    self.pprime = np.zeros((npsi))
    self.FFprime = np.zeros((npsi))
    self.Qpsi = np.zeros((npsi))

    self.descriptor = descriptor
    self.date = date
    self.shot = shot
    self.time = time
    
    # geometrical coefficients
    self.dR = (Rm-R1) / (nr-1)
    self.dz = (zn-z1) / (nz-1)    
    s = self.dz / self.dR
    s2 = np.square(s)
    R_1_2 = 0.5*(self.R[1:] + self.R[0:-1])
    
    # 1D coefficients vs. R
    sp = s2 / R_1_2[1:]
    sm = s2 / R_1_2[0:-1]
    delta = -sp -sm -2.0/ self.R[1:-1]

    # 2D versions (i.e. same for all z)
    nz1 = np.ones(nz-2)
    self.sp2  = np.outer(sp, nz1)
    self.sm2  = np.outer(sm, nz1)
    self.Rm2  = np.outer(1.0 / self.R[1:-1], nz1)
    self.delta2 = np.outer(delta, nz1)

    # no limiter and plasma boundary yet
    self.Rlim = None
    self.zlim = None
    self.Rbdry = None
    self.zbdry = None

    # entire grid is plasma domain (for now, until limiter boundary loaded)
    self.PlasmaDomain = 2 * np.ones([nr-2, nz-2], dtype='int8')


# ...................................................


  def write_geqdsk(self, fn):
    """
      Save equilibrium as geqdsk file with path/name fn.
    """
    f = open(fn, 'w')

    nw = self.PsiRz.shape[0]
    nh = self.PsiRz.shape[1]
    f.write('{:12}{:12}{:8}{:8}{:8}{:4d}{:4d}{:4d}\n'.format(self.descriptor, self.date, ' #'+str(self.shot),
                                                          '  '+str(int(1000*self.time))+'ms', '        ', 3, nw, nh))
    
    # line 2
    rdim = self.R[-1] - self.R[0]         # "Horizontal dimension in meter of computational box"
    zdim = self.z[-1] - self.z[0]         # "Vertical dimension in meter of computational box"  
    rcentr = self.Rmaxis                # "R in meter of vacuum toroidal magnetic field BCENTR"
    rleft = self.R[0]                   # "Minimum R in meter of rectangular computational box"
    zmid =  0.5*(self.z[-1] + self.z[0])  # "Z of center of computational box in meter"
    write_array_n5 (f, [rdim, zdim, rcentr, rleft, zmid])

    # line 3
    bcentr = self.Fpol[0] / rcentr      # "Vacuum toroidal magnetic field in Tesla at RCENTR"
    write_array_n5 (f, [self.Rmaxis, self.zmaxis, self.Psimag, self.Psibdry, bcentr])

    # line 4
    write_array_n5 (f, [self.Iplasma, self.Psimag, 0.0, self.Rmaxis, 0.0])

    # line 5
    write_array_n5 (f, [self.zmaxis, 0.0, self.Psibdry, 0.0, 0.0])

    # flux functions must be interpolated onto an equidistant grid
    # with nw elements (a g-file format quirk)
    Psimag = self.Psimag
    Psibdry = self.Psibdry
    if Psimag>Psibdry:
      if Psimag > np.amax(self.Psi):
        print "Ajusting Psimag from ", Psimag, " to ", np.amax(self.Psi)
        Psimag = np.amax(self.Psi)
      if Psibdry < np.amin(self.Psi):
        print "Ajusting Psibdry from ", Psibdry, " to ", np.amin(self.Psi)
        Psibdry = np.amin(self.Psi)
    else:      
      if Psimag < np.amin(self.Psi):
        print "Ajusting Psimag from ", Psimag, " to ", np.amin(self.Psi)
        Psimag = np.amin(self.Psi)
      if Psibdry > np.amax(self.Psi):
        print "Ajusting Psibdry from ", Psibdry, " to ", np.amax(self.Psi)
        Psibdry = np.amax(self.Psi)    
    xn = np.linspace(Psimag, Psibdry, num=nw, endpoint=True)
    write_interpol_profile (f, xn, self.Psi, self.Fpol)
    write_interpol_profile (f, xn, self.Psi, self.Pres)
    write_interpol_profile (f, xn, self.Psi, self.FFprime)
    write_interpol_profile (f, xn, self.Psi, self.pprime)
    write_array_n5 (f, self.PsiRz.T.reshape(nh*nw))
    write_interpol_profile (f, xn, self.Psi, self.Qpsi)

    # size of boundary and limiter contours
    if (not self.Rbdry is None) and (not self.zbdry is None):
      nbbs = self.Rbdry.size
    else:
      nbbs = 0
    if (not self.Rlim is None) and (not self.zlim is None):
      limitr = self.Rlim.size
    else:
      limitr = 0
    f.write('{:4d}{:4d}\n'.format(nbbs, limitr))

    # plasma boundary
    if nbbs>0:
      boundary = np.zeros(2*nbbs)
      boundary[0::2] = self.Rbdry
      boundary[1::2] = self.zbdry
      write_array_n5 (f, boundary)

    # limiter contour
    if limitr>0:
      Rzlim = np.zeros(2*limitr)
      Rzlim[0::2] = self.Rlim
      Rzlim[1::2] = self.zlim
      write_array_n5 (f, Rzlim)

    # done.
    f.close()

    
# ...................................................
  def plot_Rz (self, color='blue', nint=10):
    """Plot 2D equilibrium (R,z) poloidal cross section using Matplotlib.
       color='black'    plot color 
       nint = 10        number of internal surfaces to plot
    """
    cc = colrs.ColorConverter()
    co = cc.to_rgb(color)    # original RGB tupel
    cm = max(co)
    cb = (co[0]/cm, co[1]/cm, co[2]/cm)    # bright version 
    cd = (co[0]/1.5, co[1]/1.5, co[2]/1.5)    # dark version

    # flux matrix is compulsory
    if self.Psibdry > self.Psimag:   # contour wants increasing levels
      levels = self.Psimag + (self.Psibdry - self.Psimag) * np.linspace(0,1,nint)
    else:
      levels = self.Psibdry + (self.Psimag - self.Psibdry) * np.linspace(0,1,nint)
    if nint>0:
      plt.contour(self.R, self.z, self.PsiRz.T, levels=levels, 
                  linestyles='solid', colors=color, linewidths=0.3)
    # boundary and vessel contour are optional
    if (not self.Rbdry is None) and (not self.zbdry is None):
      plt.plot(self.Rbdry, self.zbdry, color=cb)
    if (not self.Rlim is None) and (not self.zlim is None):
      plt.plot(self.Rlim, self.zlim, color=cd)


# ...................................................
  def set_limiter(self, Rlim, zlim):
    """
      Set limiter contour, defined as line segments (Rlim[i], zlim[i]).
      Calculate dependent quantities:
      - flux matrix mask for plasma interior
    """
    self.Rlim = np.array(Rlim)
    self.zlim = np.array(zlim)

    # For later use, identify grid elements that are bounded by the limiter.
    # We consider a grid of the size of the current density matrix
    self.PlasmaDomain = np.zeros([self.R.size-2, self.z.size-2], dtype='int8')

    # limiter contour, expressed as indices to current density matrix
    R1 = self.R[1:-1]
    z1 = self.z[1:-1]
    R2 = np.outer(R1, np.ones(self.Rlim.size))
    z2 = np.outer(z1, np.ones(self.zlim.size))
    Rlim2 = np.outer(np.ones(R1.size), self.Rlim)
    zlim2 = np.outer(np.ones(z1.size), self.zlim)
    ir = np.argmin(abs(Rlim2-R2),axis=0)
    iz = np.argmin(abs(zlim2-z2),axis=0)

    for i in range(ir.size-1):
      line = bresenham_line(ir[i],iz[i],ir[i+1],iz[i+1])
      for p in line:
        self.PlasmaDomain[p[0],p[1]] = 1   # limiter contour
        # print "Bresenham: ", p[0],p[1]

    # limiter interior = plasma domain
    # assume grid center is inside plasma domain
    floodfill (self.PlasmaDomain, (R1.size/2, z1.size/2), 2)

# ...................................................
  def fluxcontours(self, psi):
    """
      Find (R,z) contours with flux stored in the 'psi' array.
      Returns a list of 
      For each flux, the first closed contour is returned
    """
    Rv, zv = np.meshgrid(self.R, self.z, sparse=False, indexing='ij')
    r = contourc(Rv, zv, self.PsiRz, psi)
    c = []    
    for level in r:
      for pk in level:
          Rc = pk[:, 0]
          zc = pk[:, 1]
          # check if a closed contour (first and last points identical)
          if pk[0, 0]==pk[-1, 0] and pk[0, 1]==pk[-1, 1]:
            c.append(pk)
    return c

# ...................................................
  def DeltaStar(self):
    """ 
      calculate elliptic differential operator on flux matrix.
      Prealculated coefficients (sp2, sm2, Rm2, delta2) are passed by reference
    """
    rho = self.sp2 * self.PsiRz[2:,1:-1] + self.sm2 * self.PsiRz[:-2,1:-1] \
        + self.Rm2 * (self.PsiRz[1:-1,2:] + self.PsiRz[1:-1,:-2]) \
        + self.delta2 * self.PsiRz[1:-1,1:-1]
    return rho

# ...................................................
  def ToroidalCurrentDensity(self):
    """ 
      calculate toroidal current density by differentiation of the poloidal flux matrix 
    """    
    rho = self.DeltaStar()
    # toroidal current density, on entire (R,z) grid
    jtor = - rho  / (self.dz * self.dz) / (4e-7*np.pi)
    # toroidal current density in plasma domain only
    jplasma = (jtor * 0.5* self.PlasmaDomain)
    return jplasma
  
# ..................................................
  def BRz(self):
    """
      poloidal field components B_R, Bz, on intermediate grid points
    """
    nz1 = np.ones(self.z.size-1)
    R1 = 0.5* (self.R[0:-1] + self.R[1:])
    Rm2  = np.outer(1.0 / R1, nz1)
    BR = -0.5* (self.PsiRz[0:-1,1:] + self.PsiRz[1:,1:] \
                - self.PsiRz[0:-1,0:-1] - self.PsiRz[1:,0:-1] ) * Rm2 / self.dz 
    Bz =  0.5* (self.PsiRz[1:,0:-1] + self.PsiRz[1:,1:] \
                - self.PsiRz[0:-1,0:-1] - self.PsiRz[0:-1,1:] ) * Rm2 / self.dR
    return BR, Bz



# -------------------------------------------------
# GEQDSK 
# -------------------------------------------------

class geqdsk (eqrzfd):
  """Read equilibrium from file fn in geqdsk (EFIT) format."""
  
  def __init__ (self, fn):
    f = open(fn, 'r')

    s1 = f.readline()  # line 1: descriptors and (R,z) matrix size
    descriptor = s1[0:11]
    date = s1[11:22]
    try:
      shot = int(re.findall(r'\d+', s1[22:32])[0])
    except:
      shot = 0
    try:
      time = float(re.findall(r'\d+', s1[32:48])[0]) / 1000.0
      print "time=", time
    except:
      time = 0.0
    dum = s1[49:52]
    nw = int(s1[52:56]);
    nh = int(s1[56:]);

    s = f.readline()  # line 2
    rdim = float(s[0:16])     # "Horizontal dimension in meter of computational box"
    zdim = float(s[16:32])     # "Vertical dimension in meter of computational box"
    rcentr = float(s[32:48])   # "R in meter of vacuum toroidal magnetic field BCENTR"
    rleft = float(s[48:64])    # "Minimum R in meter of rectangular computational box"
    zmid = float(s[64:80])     # "Z of center of computational box in meter"

    # allocate space now that we know how much. geqdsk: npsi = nw :-(
    eqrzfd.__init__ (self, nw, nh, nw, rleft, rleft+rdim, zmid-0.5*zdim, zmid+0.5*zdim,
                     descriptor = descriptor, date = date, shot = shot, time = time)

    s = f.readline()  # line 3
    self.Rmaxis = float(s[0:16])    # "R of magnetic axis in meter"
    self.zmaxis = float(s[16:32])    # "Z of magnetic axis in meter"
    self.Psimag = float(s[32:48])    # "poloidal flux at magnetic axis in Weber /rad"
    self.Psibdry = float(s[48:64])    # "poloidal flux at the plasma boundary in Weber /rad"
    bcentr = float(s[64:80])       # "Vacuum toroidal magnetic field in Tesla at RCENTR"

    s = f.readline()  # line 4
    self.Iplasma = float(s[0:16])   # "Plasma current in Ampere"

    s = f.readline()  # line 5     # nothing needed from this line

    self.Psi = np.linspace(self.Psimag, self.Psibdry, nw)
    read_array_n5 (f, nw, self.Fpol)
    read_array_n5 (f, nw, self.Pres)
    read_array_n5 (f, nw, self.FFprime)
    read_array_n5 (f, nw, self.pprime)
    PsiRz1 = np.zeros((nw*nh))
    read_array_n5 (f, nw*nh, PsiRz1)
    self.PsiRz = PsiRz1.reshape(nh,nw).T

# allow q profile, boundary and vessel contour as optional
    try:
      read_array_n5 (f, nw, self.Qpsi)
      try:
        s = f.readline() 
        nbbs = int(s[0:5]);
        limitr = int(s[5:10]);
        boundary = np.zeros((2*nbbs))
        read_array_n5 (f, 2*nbbs, boundary)
        vessel = np.zeros((2*limitr))
        read_array_n5 (f, 2*limitr, vessel)
        self.Rbdry = boundary[0::2]
        self.zbdry = boundary[1::2]
        self.set_limiter (vessel[0::2], vessel[1::2])
      except:
        pass
    except:
      pass
    f.close()

# end of geqdsk.__init__()
    



# ---------------------------------------------------------------------
# demo

# if True:
if __name__ == '__main__':

  # eq = geqdsk('testdata/eqdsk/d3d/g164362.04545_828')
  eq = geqdsk('testdata/eqdsk/aug/g31128.2300_mkit')
  # eq = geqdsk('g34214_3200_IDE2_Hager')

  # toroidal plasma current density    
  jplasma = eq.ToroidalCurrentDensity()

  # total toroidal plasma current 
  Iplasma = np.sum(jplasma * eq.dR * eq.dz)
  print "Toroidal plasma current: ", eq.Iplasma, " A (equil), ", Iplasma, " A (calc)"
    
  if True:    
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111, aspect='equal')
  
  # plot equilibrium and some contours
  if False:
    eq.plot_Rz(color='blue')
    psin = np.array([0.5, 0.7, 0.9])
    c = eq.fluxcontours(eq.Psimag + psin*(eq.Psibdry-eq.Psimag))
    for c1 in c:
      plt.plot(c1[:,0], c1[:,1], color='red')  
    plt.axis("equal")
    plt.show()

  # plot current density
  if False:
    levels = np.linspace(np.amin(jtor), np.amax(jtor), 300)
    plt.contour(eq.R[1:-1], eq.z[1:-1], jtor.T, levels=levels)
    plt.show()

