#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 10:17:41 2016

@author: Christopher Albert
"""

__all__ = ['write_boozer_head', 'append_boozer_block_head',
           'append_boozer_block', 'convert_to_boozer', 'BoozerFile']


def write_boozer_head(filename, version, shot: int, m0b, n0b, nsurf_fmt, nfp, psi_tor_a, aminor, Rmajor):
  import getpass

  global_variables_head_line=" m0b   n0b  nsurf  nper    flux [Tm^2]        a [m]          R [m]"

  with open(filename, 'w') as outfile:
    outfile.write("CC Boozer-coordinate data file\n")
    outfile.write("CC Version: " + version + '\n')
    outfile.write("CC Author: " + getpass.getuser() + '\n')
    outfile.write("CC shot: {:4d}".format(shot) + '\n')
    outfile.write(global_variables_head_line + '\n')
    outfile.write("{:4d}  {:4d}  {:3d} {:4d}  ".format(m0b, n0b, nsurf_fmt, nfp))
    outfile.write("{:12.6e}   {:12.6e}   {:12.6e}".format(psi_tor_a, aminor, Rmajor))
    outfile.write("\n")

def append_boozer_block_head(filename, s, iota, bsubvB, bsubuB, pprime, vp, enfp):
  import numpy as np
  import scipy.constants

  mu0 = scipy.constants.mu_0

  _append_boozer_block_head(filename, s, iota,
    -2.0*np.pi/mu0*bsubvB/enfp, -2.0*np.pi/mu0*bsubuB, pprime,
    -4.0*np.pi**2*vp/enfp)

def _append_boozer_block_head(filename, s, iota, Jpol_divided_by_nper, Itor, pprime, sqrt_g_00):
  with open(filename, 'a') as outfile:
    outfile.write('        s               iota           Jpol/nper          '+
    'Itor            pprime         sqrt g(0,0)\n')
    outfile.write('                                          [A]           '+
    '[A]             [Pa]         (dV/ds)/nper\n')
    outfile.write(' {:16.8e}'.format(s))
    outfile.write(' {:16.8e}'.format(iota))
    outfile.write(' {:16.8e}'.format(Jpol_divided_by_nper))
    outfile.write(' {:16.8e}'.format(Itor))
    outfile.write(' {:16.8e}'.format(pprime))
    outfile.write(' {:16.8e}'.format(sqrt_g_00))
    outfile.write('\n')

def append_boozer_block(filename, mb, nb, rmnb, zmnb, vmnb, bmnb, enfp):
  from numpy import array, int32
  _append_boozer_block(filename, mb, array(nb/enfp, dtype=int32),
      rmnb.real, -rmnb.imag,
      zmnb.real, -zmnb.imag,
      vmnb.real, -vmnb.imag,
      bmnb.real, -bmnb.imag)

def _append_boozer_block(filename, mb, nb, rmnc, rmns, zmnc, zmns, vmnc, vmns, bmnc, bmns):
  with open(filename, 'a') as f:
    f.write('    m    n      rmnc [m]         rmns [m]         zmnc [m]  '+
            '       zmns [m]         vmnc [ ]         vmns [ ]         '+
            'bmnc [T]         bmns [T]\n')
    for k in range(len(mb)):
      f.write(' {:4d} {:4d}'.format(mb[k], nb[k]))
      f.write(' {:16.8e} {:16.8e}'.format(rmnc[k], rmns[k]))
      f.write(' {:16.8e} {:16.8e}'.format(zmnc[k], zmns[k]))
      f.write(' {:16.8e} {:16.8e}'.format(vmnc[k], vmns[k]))
      f.write(' {:16.8e} {:16.8e}'.format(bmnc[k], bmns[k]))
      f.write('\n')


def convert_to_boozer(infile, ks, outfile, uv_grid_multiplicator: int = 6):
  """
  input:
  ------
  infile: string, name(+path) of the input file (including extension).
  ks: integer, number of the flux surface to handle. Must be in the range
    [1, n-2] where n is the number of flux surfaces in the file. This is
    because of the interpolation that is done.
  outfile: string, name of the output file.
  uv_grid_multiplicator: integer, the function creates a grid in vmec u
    and v coordinates for calculation of the boozer variables. This
    factor determines how many more points should be used compared to
    the maximum m/n-mode number. For example, a value of 6 would mean,
    60 points in u and 24 in v direction will be used if 10 and 4 are
    the maximum modes numbers for m and n respectively.
    The default value should work for some cases(?), but might need to
    be increased for others. [Default: 6]

  output:
  -------
  none

  sideeffects:
  ------------
  creates outfile.
  """
  # Import modules.

  # Look for a fourier transformation module.
  try:
    from fourier import fourierseries as fourierseries1
    print('Using f2py fourierseries .so')
  except:
    try:
      from fourier_win import fourierseries as fourierseries1
      print('Using Win32 f2py fourierseries DLL')
    except:
      def fourierseries1(fmn, u, v, m, n):
        ycpl = np.sum(fmn*np.exp(1j*(m*u + n*v)))
        return np.real(ycpl)
      print('Using Python fourierseries (SLOW!)')
  # Standard stuff.
  import math
  import numpy as np
  import scipy.interpolate as ip
  import string
  import sys
  import time

  def nextpow2(i):
    """Return next largest power of 2.

    Return the smallest n > i, with n of the form n = 2^k.

    input:
    ------
    i: number (if integer or float should not matter), for which to
      calculate the next power of 2.

    output:
    -------
    integer, as defined above.

    sideeffects:
    ------------
    None
    """
    n = 1
    while n < i: n *= 2
    return n

  def fourierseries(fmn, u, v, m, n):
    if type(u) == float or u.size == 1:
      ycpl = np.sum(fmn*np.exp(1j*(m*u + n*v)))
      return np.real(ycpl)
    else:
      y = np.zeros(u.shape)
      for k in range(u.size):
        y[k] = fourierseries1(fmn, u[k], v[k], m, n)
      return y

  pi = np.pi

  nl     = 3          # number of points on each side for Lagrange interpolation
  plot   = False

  t = time.time()

  low = True # use lower mode number for NetCDF (not Nyquist one)
  if infile.endswith('nc'):
    print('Reading NetCDF file {} at flux surface {}'.format(infile,ks))
    import scipy.io.netcdf as ncdf
    f = ncdf.netcdf_file(infile)
    data = f.variables

    enrho = np.copy(data['ns'].data)
    mlow = np.array(np.copy(data['xm'].data),int)
    nlow = np.array(-np.copy(data['xn'].data),int)
    m = np.array(np.copy(data['xm_nyq'].data),int)
    n = np.array(-np.copy(data['xn_nyq'].data),int)
    empol = int(np.max(np.abs(m)) + 1)
    entor = int(np.max(np.abs(n)))
    empoll = int(np.max(np.abs(mlow)) + 1)
    entorl = int(np.max(np.abs(nlow)))
    phip = -np.copy(data['phipf'].data)/(2.0*np.pi)
    pres = np.copy(data['presf'].data)
    empmnt= np.copy(data['mnmax_nyq'].data)
    empmntl= np.copy(data['mnmax'].data)
    buco =  np.copy(data['buco'].data)
    bvco =  np.copy(data['bvco'].data)
    iota =  np.copy(data['iotas'].data)
    enfp =  np.copy(data['nfp'].data)
    vp =    np.copy(data['vp'].data)

    try:
      rmnl = np.copy(data['rmnc'].data - 1.0j*data['rmns'].data)
      zmnl = np.copy(data['zmnc'].data - 1.0j*data['zmns'].data)
      lmnl = np.copy(data['lmnc'].data - 1.0j*data['lmns'].data)
      bsubumn = np.copy(data['bsubumnc'].data - 1.0j*data['bsubumns'].data)
      bsubvmn = np.copy(data['bsubvmnc'].data - 1.0j*data['bsubvmns'].data)
      bsubsmn = np.copy(data['bsubsmnc'].data - 1.0j*data['bsubsmns'].data)
      bsupumn = np.copy(data['bsupumnc'].data - 1.0j*data['bsupumns'].data)
      bsupvmn = np.copy(data['bsupvmnc'].data - 1.0j*data['bsupvmns'].data)
    except:  # Stellarator-symmetric case
      print('Stellarator-symmetric case')
      rmnl = np.copy(data['rmnc'].data)
      zmnl = np.copy(- 1.0j*data['zmns'].data)
      lmnl = np.copy(- 1.0j*data['lmns'].data)
      bsubumn = np.copy(data['bsubumnc'].data)
      bsubvmn = np.copy(data['bsubvmnc'].data)
      bsubsmn = np.copy(- 1.0j*data['bsubsmns'].data)
      bsupumn = np.copy(data['bsupumnc'].data)
      bsupvmn = np.copy(data['bsupvmnc'].data)

    del data
    f.close()
    # use only modes where all quantities are defined
    condi = (np.abs(m)<empoll) & (np.abs(n)<=entorl)

    if low:
      m = mlow
      n = nlow
      empol = empoll
      entor = entorl
      empmnt = empmntl
      rmn = rmnl
      zmn = zmnl
      lmn = lmnl
      bsubumn = bsubumn[:,condi]
      bsubvmn = bsubvmn[:,condi]
      bsubsmn = bsubsmn[:,condi]
      bsupumn = bsupumn[:,condi]
      bsupvmn = bsupvmn[:,condi]
    else:
      rmn = np.zeros(bsubumn.shape, complex)
      zmn = np.zeros(bsubumn.shape, complex)
      lmn = np.zeros(bsubumn.shape, complex)
      rmn[:,condi] = rmnl
      zmn[:,condi] = zmnl
      lmn[:,condi] = lmnl

  else:
    print('Reading text file {} at flux surface {}'.format(infile,ks))
    with open(infile) as f:
      lines = f.readlines()

    block = np.fromstring(lines[0], sep=' ')
    gamma = block[0]
    enfp  = int(block[1])
    enrho = int(block[2])
    block = np.fromstring(lines[1], sep=' ')
    empol = int(block[0])
    entor = int(block[1])
    empmnt= int(block[2])
    block = np.fromstring(lines[2], sep=' ')
    eiasym = int(block[0])
    phiedge  = block[1]

    lfourier = lines[3:enrho*empmnt+3]
    lprofile = lines[enrho*empmnt+3:enrho*empmnt+int(math.ceil((enrho-1)*12/5))+4]
    fourier = np.fromstring(string.join(lfourier),sep=' ').reshape(enrho,empmnt,16)
    profile = np.fromstring(string.join(lprofile), sep=' ')[0:12*(enrho-1)]\
        .reshape(enrho-1,12)

    # Fourier quantities
    rmn = fourier[:,:,0] - 1.0j*fourier[:,:,2]
    zmn = fourier[:,:,3] - 1.0j*fourier[:,:,1]
    bsupumn = fourier[:,:,4] - 1.0j*fourier[:,:,6]
    bsupvmn = fourier[:,:,5] - 1.0j*fourier[:,:,7]
    lmn = fourier[:,:,9] - 1j*fourier[:,:,8]
    bsubsmn = fourier[:,:,15] - 1.0j*fourier[:,:,12]
    bsubumn = fourier[:,:,10] - 1.0j*fourier[:,:,13]
    bsubvmn = fourier[:,:,11] - 1.0j*fourier[:,:,14]

    # Profile quantities
    iota = profile[:,0] # TODO: this seems to be full mesh, need Lagrange
    mass = profile[:,1]
    pres = profile[:,2]
    phip = profile[:,3]
    buco = profile[:,4]
    bvco = profile[:,5]
    phi  = profile[:,6]
    vp   = profile[:,7]
    overr = profile[:,8]
    jcuru = profile[:,9]
    jcurv = profile[:,10]
    specw = profile[:,11]

    # Modenumbers
    k = 0
    m = np.zeros(empmnt,dtype=int)
    n = np.zeros(empmnt,dtype=int)
    for mk in range(empol):
      nmin0 = -entor
      if mk < 1:
        nmin0 = 0
      for nk in range(nmin0,entor+1):
        m[k] = mk
        n[k] = -nk*enfp
        k = k+1

    # Adjustments
    iota = np.insert(iota, 0, 0.0)
    phip = np.insert(phip, 0, 0.0)
    buco = np.insert(buco, 0, 0.0)
    bvco = np.insert(bvco, 0, 0.0)
    vp = np.insert(vp, 0, 0.0)

  ns = enrho - 1
  ds = 1.0/ns
  s   = (np.arange(0,ns)+0.5)*ds
  sf  = (np.arange(0,ns+1))*ds

  # Reduce number of points for radial interpolation if near axis or
  # near outer border.
  if(ks < 3 or ks > ns-3):
    nl = 2

  if(ks < 2 or ks > ns-2):
    nl = 1

  print('Defining parameters and functions')

  cond1   = (m != 0)
  cond2   = (n != 0)
  m1 = m[cond1]; n1 = n[cond1]
  m2 = m[cond2]; n2 = n[cond2]

  s           = np.insert(s, 0, 0.0)
  dpsitords   = -phip

  # Full mesh quantities
  ppoly  = ip.lagrange(sf[ks-nl:ks+nl],pres[ks-nl:ks+nl])
  pspoly = np.polyder(ppoly)
  psval  = np.polyval(ppoly,s[ks])

  # Radial interpolation of the quantities. Also get derivative from interpolation.
  # Half to full grid?
  rmnval = []; rsmnval = []
  zmnval = []; zsmnval = []
  lmnval = []; lsmnval = []
  for km in range(empmnt):
    rpoly = ip.lagrange(sf[ks-nl:ks+nl],rmn[ks-nl:ks+nl,km])
    rspoly = np.polyder(rpoly)
    rmnval.append(np.polyval(rpoly, s[ks]))
    rsmnval.append(np.polyval(rspoly, s[ks]))

    zpoly = ip.lagrange(sf[ks-nl:ks+nl],zmn[ks-nl:ks+nl,km])
    zspoly = np.polyder(zpoly)
    zmnval.append(np.polyval(zpoly, s[ks]))
    zsmnval.append(np.polyval(zspoly, s[ks]))

    lpoly = ip.lagrange(sf[ks-nl:ks+nl],lmn[ks-nl:ks+nl,km])
    lspoly = np.polyder(lpoly)
    lmnval.append(np.polyval(lpoly, s[ks]))
    lsmnval.append(np.polyval(lspoly, s[ks]))

  rmnval = np.array(rmnval); zmnval = np.array(zmnval); lmnval = np.array(lmnval)
  rsmnval= np.array(rsmnval);zsmnsval=np.array(zsmnval);lsmnval=np.array(lsmnval)

  # Cylindrical coordinates
  def r(u,v): return fourierseries(rmnval,u,v,m,n)
  def drdu(u,v): return fourierseries(1j*m*rmnval,u,v,m,n)
  def drds(u,v): return fourierseries(rsmnval,u,v,m,n)

  def z(u,v): return fourierseries(zmnval,u,v,m,n)
  def dzdu(u,v): return fourierseries(1j*m*zmnval,u,v,m,n)
  def dzds(u,v): return fourierseries(zsmnsval,u,v,m,n)

  # Stream function
  def lam(u,v): return fourierseries(lmnval,u,v,m,n)
  def dlamdu(u,v): return fourierseries(1j*m*lmnval,u,v,m,n)
  def dlamdv(u,v): return fourierseries(1j*n*lmnval,u,v,m,n)
  def dlamds(u,v): return fourierseries(lsmnval,u,v,m,n)

  def bsupu(u,v,fsi): return fourierseries(bsupumn[fsi,:],u,v,m,n)
  def bsupv(u,v,fsi): return fourierseries(bsupvmn[fsi,:],u,v,m,n)

  def bsubs(u,v,fsi): return fourierseries(bsubsmn[fsi,:],u,v,m,n)
  def bsubu(u,v,fsi): return fourierseries(bsubumn[fsi,:],u,v,m,n)
  def bsubv(u,v,fsi): return fourierseries(bsubvmn[fsi,:],u,v,m,n)

  def bmod2(u,v,fsi): return bsupu(u,v,fsi)*bsubu(u,v,fsi) +\
      bsupv(u,v,fsi)*bsubv(u,v,fsi)
  def bmod(u,v): return np.sqrt(bmod2(u,v,ks))

  # Metric tensor
  def G(u,v): return drdu(u,v)*dzds(u,v)-drds(u,v)*dzdu(u,v)
  def sqrtg(u,v): return np.abs(r(u,v)*G(u,v))

  # Alternative definition of stream function
  def dlamdu0(u,v): return sqrtg(u,v)/dpsitords[ks]*bsupv(u,v,ks)-1.0
  def dlamdv0(u,v): return (-sqrtg(u,v)/(iota[ks]*dpsitords[ks])*bsupu(u,v,ks)+1.0)*iota[ks]

  # VMEC magnetic coordinates
  def uf(u,v,fsi): return u + lam(u,v)
  def sqrtgf(u,v): return sqrtg(u,v)/np.abs(1+dlamdu(u,v))
  def bsupuf(u,v,fsi): return (1+dlamdu(u,v,fsi))*bsupu(u,v,fsi)+\
      dlamdv(u,v,fsi)*bsubv(u,v,fsi)

  # Boozer coordinates
  bsubuB   = np.real(bsubumn[:,0])
  bsubvB   = np.real(bsubvmn[:,0])
  bsubuB2  = buco
  bsubvB2  = bvco
  pprime   = 0.0

  # output flux surface quantities
  # TODO: check sign convention. This one matches Strumberger output
  # TODO: check dV/ds and bvco enfp factor. Documentation says it's already included
  #       but it seems it needs to be added to get the correct output
  # TODO: pprime
  append_boozer_block_head(outfile, s[ks], iota[ks], bsubvB[ks], bsubuB[ks], pprime, vp[ks], enfp)

  def sqrtgB(u,v): return np.abs(dpsitords[ks]*(iota[ks]*bsubuB[ks]+
                                  bsubvB[ks])/bmod2(u,v,ks))


  def hcheck(u,v): return sqrtgf(u,v)/sqrtgB(u,v) - 1.0

  # Boozer conversion
  print('Boozer conversion after Nuehrenberg/Zille')

  hmn1 = (bsubumn[ks,cond1]-1j*m1*bsubuB[ks]*lmnval[cond1])/\
          (1j*m1*(bsubuB[ks]+bsubvB[ks]/iota[ks]))
  hmn2 = (bsubvmn[ks,cond2]-1j*n2*bsubuB[ks]*lmnval[cond2])/\
          (1j*n2*(bsubuB[ks]+bsubvB[ks]/iota[ks]))

  hmn = np.zeros(m.shape, dtype='complex')
  hmn[cond2] = hmn2
  hmn[cond1] = hmn1

  def H(u,v): return fourierseries(hmn,u,v,m,n)
  def dHdu(u,v): return fourierseries(1j*m*hmn,u,v,m,n)
  def dHdv(u,v): return fourierseries(1j*n*hmn,u,v,m,n)

  def dth(u,v): return fourierseries(lmnval+hmn, u, v, m, n)
  def dph(u,v): return fourierseries(hmn/iota[ks], u, v, m, n)

  # Calculate Boozer modes
  m0b = 2*empol-1
  n0b = 4*entor+1
  mb = np.zeros(int((m0b-1)*n0b+(n0b-1)/2+1),dtype=int)
  nb = np.zeros(int((m0b-1)*n0b+(n0b-1)/2+1),dtype=int)
  k=0
  for mk in range(m0b):
    nmin0 = -2*entor
    if mk < 1:
      nmin0 = 0
    for nk in range(nmin0,2*entor+1):
      mb[k] = mk
      nb[k] = nk*enfp
      k = k+1

  nu = uv_grid_multiplicator*np.max(np.abs(m))-1
  nv = uv_grid_multiplicator*np.max(np.abs(n))+1
  du = 2.0*pi/nu
  dv = 2.0*pi/nv
  up = np.arange(0,2*pi,du)
  vp = np.arange(0,2*pi,dv)
  [U,V] = np.meshgrid(up,vp)
  U = U.flatten().T; V = V.flatten().T
  THB = U + dth(U,V)
  PHB = V + dph(U,V)
  R = r(U,V)
  Z = z(U,V)
  B = bmod(U,V)
  H1 = H(U,V)
  HCHECK = hcheck(U,V)
  Jbinv = np.abs((1.0+dlamdu(U,V)+dHdu(U,V))*(1.0+dHdv(U,V)/iota[ks])\
          -dHdu(U,V)/iota[ks]*(dlamdv(U,V)+dHdv(U,V)))

  print("Computing Boozer modes")
  dthmnb = np.zeros(mb.shape, complex)
  dphmnb = np.zeros(mb.shape, complex)
  rmnb = np.zeros(mb.shape, complex)
  zmnb = np.zeros(mb.shape, complex)
  bmnb = np.zeros(mb.shape, complex)
  hmnb = np.zeros(mb.shape, complex)
  hcheckmnb = np.zeros(mb.shape, complex)
  for km in range(len(mb)):
    efun = 2.0/(2.0*np.pi)**2*np.exp(-1.0j*(mb[km]*THB + nb[km]*PHB))*du*dv
    if (mb[km]==0) and (nb[km]==0):
      efun = efun/2.0
    dthmnb[km] = np.sum(efun*Jbinv*(THB-U))
    dphmnb[km] = np.sum(efun*Jbinv*(PHB-V))
    rmnb[km] = np.sum(efun*Jbinv*R)
    zmnb[km] = np.sum(efun*Jbinv*Z)
    bmnb[km] = np.sum(efun*Jbinv*B)
    hmnb[km] = np.sum(efun*Jbinv*H1)
    hcheckmnb[km] = np.sum(efun*Jbinv*HCHECK)

  vmnb = -enfp*dphmnb/(2*np.pi)

  append_boozer_block(outfile, mb, nb, rmnb, zmnb, vmnb, bmnb, enfp)

  elapsed = time.time() - t
  print('Elapsed time: {} s'.format(elapsed))


  #~ plot = True
  if plot:
    import matplotlib.pyplot as plt
    # Check contravariant B_Boozer
    u = np.linspace(-np.pi,np.pi,100)
    v = 0

    plt.figure(1)
    plt.clf()
    plt.subplot(1,2,1)
    plt.title('dlamdu')
    plt.plot(u, [dlamdu(ui,v) for ui in u])
    plt.plot(u, [dlamdu0(ui,v) for ui in u],'--')
    plt.legend(['dlamdu','dlamdu0'])

    plt.subplot(1,2,2)
    plt.title('dlamdv')
    plt.plot(u, [dlamdv(ui,v) for ui in u])
    plt.plot(u, [dlamdv0(ui,v) for ui in u],'--')
    plt.legend(['dlamdv','dlamdv0'])

    plt.figure(3)
    plt.clf()
    plt.title('sqrtg')
    plt.plot(u, [sqrtg(ui,v) for ui in u])
    plt.plot(u, [sqrtgf(ui,v) for ui in u],'--')
    plt.plot(u, [sqrtgB(ui,v) for ui in u],'-.')

    plt.figure(4)
    plt.clf()
    plt.subplot(1,2,1)
    plt.title('bsubuB')
    plt.plot(s, bsubuB)
    plt.plot(s, bsubuB2,'--')
    plt.subplot(1,2,2)
    plt.title('bsubvB')
    plt.plot(s, bsubvB)
    plt.plot(s, bsubvB2,'--')

  #  plt.figure(4)
  #  plt.subplot(1,2,1)
  #  plt.title('bsupu')
  #  plt.plot(s, bsupumn)
  #  plt.subplot(1,2,2)
  #  plt.title('bsupv')
  #  plt.plot(s, bsupvmn)
  #  plt.show()


  #  plt.figure(5)
  #  plt.subplot(1,2,1)
  #  plt.title('rmn')
  #  plt.plot(rmnval[(n==0)&(m<5)])
  #  plt.subplot(1,2,2)
  #  plt.title('zmn')
  #  plt.plot(zmnval[(n==0)&(m<5)])
  #  plt.show()

class BoozerFile:
  """Storing the information in a boozer file.

  This class is designed to store the information of a boozer file.

  Note that there are different versions of boozer files.
  So far this class can read two of them. The first one is those of
  Strumberger, the second is called 'other' until a better name is
  found.
  They differ in the number of modes present in the file, the number of
  columns ('other' has only the real values) and in the order of the
  modes.
  Note that the functions were written with the order of 'Strumberger'
  in mind, thus they might not work with 'other' (except for
  'read_boozer').
  """
  # Note: order due to history, only Strumberger was handled in the beginning.
  file_versions_to_read = {'Strumberger': 0, 'other': 1}

  file_version = file_versions_to_read['other']

  def read_boozer(self, filename: str):
    """Reads information from a file, whose name is given as a string.

    This routine will read the content of file 'filename', assuming it
    is a boozer file (thus expecting a specific format).
    The comments are available in a list, each line an entry.
    Global fields are available as simple elements.
    Fields that depend on radius are lists, those that depend on radius
    and mode number lists of lists.
    """

    with open(filename) as f:
      lines = f.readlines()

    self.comments = []
    self.m0b = 0.0
    self.n0b = 0.0
    self.nsurf = 0
    self.nper = 0
    self.flux = 0.0 #[Tm^2]
    self.a = 0.0
    self.R = 0.0 # [m]

    self.s = []
    self.iota = []
    self.Jpol_divided_by_nper = []
    self.Itor = []
    self.pprime = []
    self.sqrt_g_00 = []

    self.m = []
    self.n = []
    self.rmnc = [] # [m]
    self.rmns = [] # [m]
    self.zmnc = [] # [m]
    self.zmns = [] # [m]
    self.vmnc = []
    self.vmns = []
    self.bmnc = [] # [T]
    self.bmns = [] # [T]

    # Get header information, e.g. comments and sizes.
    for lineindex in range(len(lines)):
      line = lines[lineindex]
      if (len(line) > 2):
        # Comments start with 'CC' followed by a whitespace.
        if (line.split()[0] == 'CC'):
          self.comments.append(line[2+1:]) # 2+1 to ignore 'CC' and the following whitespace.

        if (lineindex > 1):
          if ((lines[lineindex-1].split())[0] == 'm0b'):
            self.m0b = int((lines[lineindex].split())[0])
            self.n0b = int((lines[lineindex].split())[1])
            self.nsurf = int((lines[lineindex].split())[2])
            self.nper = int((lines[lineindex].split())[3])
            self.flux = float((lines[lineindex].split())[4])
            self.a = float((lines[lineindex].split())[5])
            self.R = float((lines[lineindex].split())[6])
            break

    self.check_which_file_version()

    blocklines = []

    blockindex = -1
    for lineindex in range(len(lines)):
      line = lines[lineindex]
      if (line.split()[0] == 's'):
        blockindex += 1
        blocklines.append([])

      if (blockindex >= 0 and len(line) > 0):
        blocklines[blockindex].append(line)

    if (len(blocklines) != self.nsurf):
      print('m0b = ' + str(self.m0b))
      print('n0b = ' + str(self.n0b))
      print('nsurf = ' + str(self.nsurf))
      print('nper = ' + str(self.nper))
      print('flux = ' + str(self.flux))
      print('a = ' + str(self.a))
      print('R = ' + str(self.R))
      print(str(len(blocklines)) + ' != ' + str(self.nsurf))
      raise Exception

    head_number_of_lines = 4
    if (self.file_version == self.file_versions_to_read['Strumberger']):
      # Strumberger format +1 for zero mode, and n modes go from -n to n,
      # so also +1 for the zero mode.
      expected_block_length = (self.m0b + 1)*(2*self.n0b + 1)
    elif (self.file_version == self.file_versions_to_read['other']):
      # other format: m = 0 has only n+1 positive modes, while other m
      # modes have positivie and negative n, so 2n+1 for each
      expected_block_length = 1*(self.n0b + 1) + self.m0b*(2*self.n0b + 1)

    for i in range(self.nsurf):
      if (len(blocklines[i]) != expected_block_length + head_number_of_lines):
        print("ERROR: number of lines in block does not match expectation (modes+header)")
        print(str(len(blocklines[i])) + ' != ' + str(expected_block_length) + ' + ' + str(head_number_of_lines))
        raise Exception

      self.m.append([])
      self.n.append([])
      self.rmnc.append([])
      self.rmns.append([])
      self.zmnc.append([])
      self.zmns.append([])
      self.vmnc.append([])
      self.vmns.append([])
      self.bmnc.append([])
      self.bmns.append([])

      for j in range(0, expected_block_length + head_number_of_lines):
        if (j == 2):
          line_split = blocklines[i][j].split();
          self.s.append(float(line_split[0]))
          self.iota.append(float(line_split[1]))
          self.Jpol_divided_by_nper.append(float(line_split[2]))
          self.Itor.append(float(line_split[3]))
          self.pprime.append(float(line_split[4]))
          self.sqrt_g_00.append(float(line_split[5]))

        if (j > 3):
          line_split = blocklines[i][j].split();
          # Consider different formats(?)
          if (len(line_split) == 10):
            self.m[i].append(int(line_split[0]))
            self.n[i].append(int(line_split[1]))
            self.rmnc[i].append(float(line_split[2]))
            self.rmns[i].append(float(line_split[3]))
            self.zmnc[i].append(float(line_split[4]))
            self.zmns[i].append(float(line_split[5]))
            self.vmnc[i].append(float(line_split[6]))
            self.vmns[i].append(float(line_split[7]))
            self.bmnc[i].append(float(line_split[8]))
            self.bmns[i].append(float(line_split[9]))
          elif (len(line_split) == 6):
            self.m[i].append(int(line_split[0]))
            self.n[i].append(int(line_split[1]))
            self.rmnc[i].append(float(line_split[2]))
            self.rmns[i].append(0.0)
            self.zmnc[i].append(float(line_split[3]))
            self.zmns[i].append(0.0)
            self.vmnc[i].append(float(line_split[4]))
            self.vmns[i].append(0.0)
            self.bmnc[i].append(float(line_split[5]))
            self.bmns[i].append(0.0)
          else:
            print('Unknown format of the boozer file. Number of elements: ' + str(len(line_split)))
            print('Can handle 6 and 10.')
            raise Exception

    self._read_from = filename

  def check_which_file_version(self):
    """
    Check comments to determine which type of boozer file this is.

    So far two versions of boozer files can be read:
    - Strumberger
    - an unnamed version
    """
    for c in self.comments:
      if 'Strumberger' in c:
        self.file_version = self.file_versions_to_read['Strumberger']
        return

  def get_rbeg(self):
    rbeg = []
    for radial_position in self.rmnc:
      rbeg.append(sum(radial_position))

    return rbeg

  def get_zbeg(self):
    zbeg = []
    for radial_position in self.zmnc:
      zbeg.append(sum(radial_position))

    return zbeg

  def get_bbeg(self):
    bbeg = []
    for radial_position in self.bmnc:
      bbeg.append(sum(radial_position))

    return bbeg

  def get_iota(self):

    return self.iota

  def __init__(self, filename: str):
    """Init routine which takes a string, representing the file to read.
    """
    if filename.endswith('nc'):
      self.convert_vmec_to_boozer(filename)
    else:
      self.read_boozer(filename)

  def write(self, filename: str):
    write_boozer_head(filename, '', 0, self.m0b, self.n0b, self.nsurf,
      self.nper, self.flux, self.a, self.R)

    for i in range(self.nsurf):
      _append_boozer_block_head(filename, self.s[i], self.iota[i],
        self.Jpol_divided_by_nper[i], self.Itor[i], self.pprime[i],
        self.sqrt_g_00[i])
      _append_boozer_block(filename, self.m[i], self.n[i],
        self.rmnc[i], self.rmns[i],
        self.zmnc[i], self.zmns[i],
        self.vmnc[i], self.vmns[i],
        self.bmnc[i], self.bmns[i])

  def contours_in_r_z_plane(self, phi: float, nplotsurf: int, outfile: str, add_last_flux_surface:bool = False):
    """Write outfile with contours based on data, at specified toroidal angle.

    Output data is first column for R, second for Z, with contours
    seperated by empty lines. This is intended for plotting with
    gnuplot, where the empty lines will cause different contours not to
    be connected. File can also be load'ed with octave/matlab, but then
    seperation into contours is lost

    Example:
    --------
      contours_in_r_z_plane(phi=0.0, nplotsurf=20, outfile="rz_contours.dat")

    Input:
    ------
    phi: toroidal angle (in units of pi) at which to compute the contours.
    nplotsurf: number of flux surfaces to calculate (output will contain
      one contour less).
    outfile: string, name of the file in which to store the data.

    output:
    -------
    none

    sideeffects:
    ------------
    Creates a file to which the data is written.

    limitations:
    ------------
    With current implementation first contour is innermost flux surface.
    Subsequent contours, may be further innwards or on top of first
    contour.
    """

    import math

    hrho = 1.0/float(nplotsurf)
    phi = phi*math.pi

    nmodes = (self.m0b+1)*(2*self.n0b+1)
    modfactor = 30
    nt = self.m0b*modfactor

    rho_tor = hrho
    s_plot = rho_tor**2
    s = 0.0
    with open(outfile, 'w') as outfile:

      for k in range(self.nsurf):
        s_old = s
        s = self.s[k]

        if (s > s_plot):
          [Rnew, Znew] = self.get_R_Z(np = nt, phi = phi, ind = k)
          # Add first point at end, to close the lines.
          Rnew.append(Rnew[0])
          Znew.append(Znew[0])

          if k > 0:
            [Rold, Zold] = self.get_R_Z(np = nt, phi = phi, ind = k-1)
            Rold.append(Rold[0])
            Zold.append(Zold[0])
          else:
            Rold = Rnew
            Zold = Znew

          w=(s_plot-s_old)/(s-s_old)
          for i in range(nt+1):
            outfile.write(' {:16.8e}'.format(Rnew[i]*w + Rold[i]*(1.0-w)))
            outfile.write(' {:16.8e}'.format(Znew[i]*w + Zold[i]*(1.0-w)))
            outfile.write('\n')

          outfile.write('\n')

          rho_tor = rho_tor + hrho
          s_plot = rho_tor**2

      if (add_last_flux_surface):
        [Rnew, Znew] = self.get_R_Z(np = nt, phi = phi, ind = -1)
        Rnew.append(Rnew[0])
        Znew.append(Znew[0])
        for i in range(nt+1):
          outfile.write(' {:16.8e}'.format(Rnew[i]))
          outfile.write(' {:16.8e}'.format(Znew[i]))
          outfile.write('\n')

  def get_radial_cut_of_minor_radius(self, poloidal_angle: float):
    """Get radial cut of the minor radius at a specified angle.
    """

    import math

    # \todo correct values for the axis.
    #raxis = self.R
    raxis = sum(self.rmnc[0]) - self.a*math.sqrt(self.s[0])
    zaxis = sum(self.zmnc[0])

    phi = 0

    radii = []

    for k in range(self.nsurf):

      R = sum(rr * math.cos(m*poloidal_angle - self.nper*n*phi) + rc * math.sin(m*poloidal_angle - self.nper*n*phi) for rr, rc, m, n in zip(self.rmnc[k], self.rmns[k], self.m[k], self.n[k]))
      Z = sum(zr * math.cos(m*poloidal_angle - self.nper*n*phi) + zc * math.sin(m*poloidal_angle - self.nper*n*phi) for zr, zc, m, n in zip(self.zmnc[k], self.zmns[k], self.m[k], self.n[k]))

      radii.append( math.sqrt((R - raxis)**2 + (Z - zaxis)**2) )

    return radii

  def get_x_point_range(self, inner_limitf: float, outer_limitf: float):
    """
    Determines the range of flux surfaces 'involved' in forming the
    'x-point'.
    At the moment this is defined as the flux surfaces where the
    derivative of the radius (regarding to index) is negative.
    Own function, so it might be replaced by somethin more elaborate.
    """
    import math
    from numpy import diff

    # Two step apporoach:
    # - first determine the extrema of the radius over index courve, the
    #   index of the o-point is then assumed to be the mean of the two
    #   indices.
    # - Determine the two radii that have (approximately) the same
    #   radius, these are thought to be the field lines that form the
    #   x-point, and thus form the outer part of the magnetic island.

    radiip = self.get_radial_cut_of_minor_radius(math.pi)
    radiip_diff = diff(radiip)

    inner_limit = max(int(inner_limitf*len(radiip)), 1)
    outer_limit = min(int(outer_limitf*len(radiip)), len(radiip)-1)

    x_point_range = [-1, -1];

    for k in range(inner_limit, outer_limit):
      if radiip_diff[k] < 0 and radiip_diff[k-1] > 0:
        x_point_range[0] = k
      if radiip_diff[k] > 0 and radiip_diff[k-1] < 0:
        x_point_range[1] = k

    k_0 = int(abs(x_point_range[1] + x_point_range[0])/2)

    delta_inner = +1.0e9
    delta_outer = +1.0e9
    k_inner = -5
    k_out = +5
    for k in range(inner_limit, outer_limit):
      # Restrict the test to regions outside the extrema, to avoid that
      # a point around k_0 is actually the closest one in radial
      # distance.
      if k < x_point_range[0] and abs(radiip[k] - radiip[k_0]) < delta_inner:
        delta_inner = abs(radiip[k] - radiip[k_0])
        k_inner = k
      if k > x_point_range[1] and abs(radiip[k] - radiip[k_0]) < delta_outer:
        delta_outer = abs(radiip[k] - radiip[k_0])
        k_outer = k

    x_point_range[0] = k_inner
    x_point_range[1] = k_outer

    return x_point_range

  def calculate_island_width(self, x_point_range: list):
    """
    Own function, so it might be replaced by somethin more elaborate.

    input:
    ------
    x_point_range: list with two elementes, giving the index of the
      inner and out flux surface that bound the island.
      If the list should have more elementes, they are ignored.

    output:
    -------
    float, calculated width of the island.
    """
    radii0 = self.get_radial_cut_of_minor_radius(0.0)
    return radii0[x_point_range[1]] - radii0[x_point_range[0]]

  def get_island_width(self):
    """Calculate and return the island width for the magnetic field.
    """

    # Limits in which to search for island/x-point.
    inner_limit = 0.1
    outer_limit = 0.7

    x_point_range = self.get_x_point_range(inner_limit, outer_limit)

    island_width = self.calculate_island_width(x_point_range)

    return island_width

  def get_rho_poloidal(self):
    """Calculate and return rho poloidal, by integrating q(s).

    \todo Find a better way to set the point closest to the axis?

    This function will determine rho_poloidal from s grid and
    iota = 1/q = d\psi_{pol} / d\psi_{tor}.
    For this it will integrate iota over s. The integral over the whole
    region serves as normalization constant, so that the values will be
    <= 1.
    \int_{0}^{s} ds' iota(s') /(\int_{0}^{a} ds' iota(s')) = \psi_{\pol} (s) \psi_{\pol}(a) = \rho_{pol}^2
    As s = \rho_{tor}^2, thus a connection between \rho_{tor} and \rho_{pol}
    has been established.

    output:
    -------
    list with values for rho_poloidal as a function.
    """

    from math import sqrt
    from numpy import array
    from scipy.integrate import simps

    iota = array(self.iota)

    psi_pol_a = simps(iota, self.s)# Order of arguments is y, x.

    rho_poloidal = []

    rho_poloidal.append(0.0)

    for k in range(1+1, len(self.iota)+1):
      rho_poloidal.append(sqrt(simps(iota[0:k], self.s[0:k])/ psi_pol_a))

    # Interpolate first point
    rho_poloidal[0] = sqrt(self.s[0]/self.s[1])*rho_poloidal[1]

    return rho_poloidal

  def write_rho_toroidal_vs_rho_poloidal(self, filename: str):
    """Write rho_tor as function of  rho_pol to a file.

    input:
    ------
    filename: string, that contains name (and path) of the file to write.

    output:
    -------
    none

    sideeffects:
    ------------
    creates file, overwrites file if it already exists.
    """

    from math import sqrt

    rho_poloidal = self.get_rho_poloidal()
    rho_toroidal = [sqrt(f) for f in self.s]

    with open(filename, 'w') as f:
      for [p, t] in zip(rho_poloidal, rho_toroidal):
        f.write("{:13.6e} {:13.6e}\n".format(p, t))

      f.write("\n")

  def get_R_Z(self, np:int = 100, phi:float = 0.0, ind:int = -1):
    """
    Get list of R and Z coordinates for given theta values.

    Returns a list with two lists for R and Z, repectively, that
    contains the values for flux surface with given index for equaly
    spaced theta values.

    input:
    ------
    np: integer, number of points to use for theta grid.
    phi: float, the phi value to use for the calculations. So far this
      can only be a single value. Defaults to 0.0.
    ind: integer, index of the flux surface to use for the calculation.
      Defaults to -1, i.e. the outermost flux surface.

    output:
    -------
    List with two elements R and Z, which are also lists (of floats).
    """
    from math import cos, sin, pi

    R = [0.0 for i in range(np)]
    Z = [0.0 for i in range(np)]

    htheta = 2.0*pi/float(np)
    theta = [htheta*x for x in range(0,np)]
    for i in range(np):
      R[i] = sum(rr * cos(m*theta[i] - self.nper*n*phi) + rc * sin(m*theta[i] - self.nper*n*phi) for rr, rc, m, n in zip(self.rmnc[ind], self.rmns[ind], self.m[ind], self.n[ind]))
      Z[i] = sum(zr * cos(m*theta[i] - self.nper*n*phi) + zc * sin(m*theta[i] - self.nper*n*phi) for zr, zc, m, n in zip(self.zmnc[ind], self.zmns[ind], self.m[ind], self.n[ind]))

    return [R, Z]

  def get_dR_dl(self, ind:int = -1, np:int = 100):
    """
    Get derivative of major radius R as a function of arc length l.

    For a given radial surface index calculate the derivative of the
    major radius R as a function of arc length l in a poloidal plane.

    Note: 'usual order' of arc length is outboard side, top, inboard
    side, bottom. Not sure how this is affected by sign conventions.

    input:
    ------
    ind: integer, index of the surface to use. Defaults to -1, i.e. the
        outermost surface.
    np: integer, number of poloidal points to use for calculatiing the
        derivative. Defaults to 100.
    """

    from math import sqrt

    phi = 0.0
    R = [0.0 for i in range(np)]
    Z = [0.0 for i in range(np)]
    l = [0.0 for i in range(np)]
    dR_dl = [0.0 for i in range(np)]

    [R, Z] = self.get_R_Z(np, phi, ind)

    # Calculate l
    l[0] = 0
    for i in range(1, np):
      l[i] = l[i-1] + sqrt((R[i]-R[i-1])**2 + (Z[i]-Z[i-1])**2)

    # Calculate dR/dl with forward difference
    for i in range(1, np-1):
      dR_dl[i] = (R[i+1] - R[i]) / (l[i+1] - l[i])
    # Last point calculated explicitly, note order in denumerator, this
    # is necessary to get the correct sign due to the jump in l.
    dR_dl[np-1] = (R[0] - R[np-1]) / (l[np-1] - l[0])

    return [dR_dl, l]


  def convert_vmec_to_boozer(self, filename: str, uv_grid_multiplicator: int = 6):
    """Intended to create a boozer file object from a vmec file.

    This function is thougt for creating a boozer file object from a
    vmec file.
    At the same time, it was thought that some speed improvement might
    be obtained, from not having the vmec file to be read again for each
    flux surface (as is the case for convert_to_boozer).

    input:
    ------
    filename: string, name (and maybe path), of the vmec file to convert.
    uv_grid_multiplicator: integer, serves as multiplicator for the m/n
      grid of vmec. Higher values might be necessary for having a
      reasanoble boundary. [6]
    """

    # Look for a fourier transformation module.
    try:
      from fourier import fourierseries as fourierseries1
      print('Using f2py fourierseries .so')
    except:
      try:
        from fourier_win import fourierseries as fourierseries1
        print('Using Win32 f2py fourierseries DLL')
      except:
        def fourierseries1(fmn, u, v, m, n):
          ycpl = np.sum(fmn*np.exp(1j*(m*u + n*v)))
          return np.real(ycpl)
        print('Using Python fourierseries (SLOW!)')

    # Standard stuff.
    import math
    import numpy as np
    import scipy.interpolate as ip
    import scipy.io.netcdf as ncdf
    import string
    import sys
    import time

    import getHeaderDataVMEC

    print("Input is netcdf file, converting to boozer.")

    t_start = time.time()

    [self.nper, self.flux, self.a, self.R, self.m0b, self.n0b] = getHeaderDataVMEC.getHeadDataVmecNc(filename)

    n = ncdf.netcdf_file(filename)
    data = n.variables

    self.nsurf = data['bmns'].shape[0]


    def nextpow2(i):
      """Return next largest power of 2.

      Return the smallest n > i, with n of the form n = 2^k.

      input:
      ------
      i: number (if integer or float should not matter), for which to
        calculate the next power of 2.

      output:
      -------
      integer, as defined above.

      sideeffects:
      ------------
      None
      """
      n = 1
      while n < i: n *= 2
      return n

    def fourierseries(fmn, u, v, m, n):
      if type(u) == float or u.size == 1:
        ycpl = np.sum(fmn*np.exp(1j*(m*u + n*v)))
        return np.real(ycpl)
      else:
        y = np.zeros(u.shape)
        for k in range(u.size):
          y[k] = fourierseries1(fmn, u[k], v[k], m, n)
        return y

    pi = np.pi

    low = True # use lower mode number for NetCDF (not Nyquist one)

    self.s = []
    self.iota = []
    self.Jpol_divided_by_nper = []
    self.Itor = []
    self.pprime = []
    self.sqrt_g_00 = []

    self.m = []
    self.n = []
    self.rmnc = [] # [m]
    self.rmns = [] # [m]
    self.zmnc = [] # [m]
    self.zmns = [] # [m]
    self.vmnc = []
    self.vmns = []
    self.bmnc = [] # [T]
    self.bmns = [] # [T]

    enrho = np.copy(data['ns'].data)
    mlow = np.array(np.copy(data['xm'].data),int)
    nlow = np.array(-np.copy(data['xn'].data),int)
    m = np.array(np.copy(data['xm_nyq'].data),int)
    n = np.array(-np.copy(data['xn_nyq'].data),int)
    empol = int(np.max(np.abs(m)) + 1)
    entor = int(np.max(np.abs(n)))
    empoll = int(np.max(np.abs(mlow)) + 1)
    entorl = int(np.max(np.abs(nlow)))
    phip = -np.copy(data['phipf'].data)/(2.0*np.pi)
    pres = np.copy(data['presf'].data)
    empmnt= np.copy(data['mnmax_nyq'].data)
    empmntl= np.copy(data['mnmax'].data)
    buco =  np.copy(data['buco'].data)
    bvco =  np.copy(data['bvco'].data)
    iota =  np.copy(data['iotas'].data)
    enfp =  np.copy(data['nfp'].data)
    vp =    np.copy(data['vp'].data)

    try:
      rmnl = np.copy(data['rmnc'].data - 1.0j*data['rmns'].data)
      zmnl = np.copy(data['zmnc'].data - 1.0j*data['zmns'].data)
      lmnl = np.copy(data['lmnc'].data - 1.0j*data['lmns'].data)
      bsubumn = np.copy(data['bsubumnc'].data - 1.0j*data['bsubumns'].data)
      bsubvmn = np.copy(data['bsubvmnc'].data - 1.0j*data['bsubvmns'].data)
      bsubsmn = np.copy(data['bsubsmnc'].data - 1.0j*data['bsubsmns'].data)
      bsupumn = np.copy(data['bsupumnc'].data - 1.0j*data['bsupumns'].data)
      bsupvmn = np.copy(data['bsupvmnc'].data - 1.0j*data['bsupvmns'].data)
    except:  # Stellarator-symmetric case
      print('Stellarator-symmetric case')
      rmnl = np.copy(data['rmnc'].data)
      zmnl = np.copy(- 1.0j*data['zmns'].data)
      lmnl = np.copy(- 1.0j*data['lmns'].data)
      bsubumn = np.copy(data['bsubumnc'].data)
      bsubvmn = np.copy(data['bsubvmnc'].data)
      bsubsmn = np.copy(- 1.0j*data['bsubsmns'].data)
      bsupumn = np.copy(data['bsupumnc'].data)
      bsupvmn = np.copy(data['bsupvmnc'].data)

    # use only modes where all quantities are defined
    condi = (np.abs(m)<empoll) & (np.abs(n)<=entorl)

    if low:
      m = mlow
      n = nlow
      empol = empoll
      entor = entorl
      empmnt = empmntl
      rmn = rmnl
      zmn = zmnl
      lmn = lmnl
      bsubumn = bsubumn[:,condi]
      bsubvmn = bsubvmn[:,condi]
      bsubsmn = bsubsmn[:,condi]
      bsupumn = bsupumn[:,condi]
      bsupvmn = bsupvmn[:,condi]
    else:
      rmn = np.zeros(bsubumn.shape, complex)
      zmn = np.zeros(bsubumn.shape, complex)
      lmn = np.zeros(bsubumn.shape, complex)
      rmn[:,condi] = rmnl
      zmn[:,condi] = zmnl
      lmn[:,condi] = lmnl

    ns = enrho - 1
    ds = 1.0/ns
    s   = (np.arange(0,ns)+0.5)*ds
    sf  = (np.arange(0,ns+1))*ds

    cond1   = (m != 0)
    cond2   = (n != 0)
    m1 = m[cond1]; n1 = n[cond1]
    m2 = m[cond2]; n2 = n[cond2]

    s           = np.insert(s, 0, 0.0)
    dpsitords   = -phip

    # Boozer coordinates
    bsubuB   = np.real(bsubumn[:,0])
    bsubvB   = np.real(bsubvmn[:,0])
    bsubuB2  = buco
    bsubvB2  = bvco
    pprime   = 0.0

    for ind in range(1, self.nsurf):
      print('Processing flux surface {}/{}'.format(ind, self.nsurf-1))
      t = time.time()

      # number of points on each side for Lagrange interpolation
      nl = 3
      # Reduce number of points for radial interpolation if near axis or
      # near outer border.
      if(ind < 3 or ind > ns-3):
        nl = 2

      if(ind < 2 or ind > ns-2):
        nl = 1

      # Full mesh quantities
      ppoly  = ip.lagrange(sf[ind-nl:ind+nl],pres[ind-nl:ind+nl])
      pspoly = np.polyder(ppoly)
      psval  = np.polyval(ppoly,s[ind])

      # Radial interpolation of the quantities. Also get derivative from interpolation.
      # Half to full grid?
      rmnval = []; rsmnval = []
      zmnval = []; zsmnval = []
      lmnval = []; lsmnval = []
      for km in range(empmnt):
        rpoly = ip.lagrange(sf[ind-nl:ind+nl],rmn[ind-nl:ind+nl,km])
        rspoly = np.polyder(rpoly)
        rmnval.append(np.polyval(rpoly, s[ind]))
        rsmnval.append(np.polyval(rspoly, s[ind]))

        zpoly = ip.lagrange(sf[ind-nl:ind+nl],zmn[ind-nl:ind+nl,km])
        zspoly = np.polyder(zpoly)
        zmnval.append(np.polyval(zpoly, s[ind]))
        zsmnval.append(np.polyval(zspoly, s[ind]))

        lpoly = ip.lagrange(sf[ind-nl:ind+nl],lmn[ind-nl:ind+nl,km])
        lspoly = np.polyder(lpoly)
        lmnval.append(np.polyval(lpoly, s[ind]))
        lsmnval.append(np.polyval(lspoly, s[ind]))

      rmnval = np.array(rmnval); zmnval = np.array(zmnval); lmnval = np.array(lmnval)
      rsmnval= np.array(rsmnval);zsmnsval=np.array(zsmnval);lsmnval=np.array(lsmnval)

      # Cylindrical coordinates
      def r(u,v): return fourierseries(rmnval,u,v,m,n)
      def drdu(u,v): return fourierseries(1j*m*rmnval,u,v,m,n)
      def drds(u,v): return fourierseries(rsmnval,u,v,m,n)

      def z(u,v): return fourierseries(zmnval,u,v,m,n)
      def dzdu(u,v): return fourierseries(1j*m*zmnval,u,v,m,n)
      def dzds(u,v): return fourierseries(zsmnsval,u,v,m,n)

      # Stream function
      def lam(u,v): return fourierseries(lmnval,u,v,m,n)
      def dlamdu(u,v): return fourierseries(1j*m*lmnval,u,v,m,n)
      def dlamdv(u,v): return fourierseries(1j*n*lmnval,u,v,m,n)
      def dlamds(u,v): return fourierseries(lsmnval,u,v,m,n)

      def bsupu(u,v,fsi): return fourierseries(bsupumn[fsi,:],u,v,m,n)
      def bsupv(u,v,fsi): return fourierseries(bsupvmn[fsi,:],u,v,m,n)

      def bsubs(u,v,fsi): return fourierseries(bsubsmn[fsi,:],u,v,m,n)
      def bsubu(u,v,fsi): return fourierseries(bsubumn[fsi,:],u,v,m,n)
      def bsubv(u,v,fsi): return fourierseries(bsubvmn[fsi,:],u,v,m,n)

      def bmod2(u,v,fsi): return bsupu(u,v,fsi)*bsubu(u,v,fsi) +\
          bsupv(u,v,fsi)*bsubv(u,v,fsi)
      def bmod(u,v): return np.sqrt(bmod2(u,v,ind))

      # Metric tensor
      def G(u,v): return drdu(u,v)*dzds(u,v)-drds(u,v)*dzdu(u,v)
      def sqrtg(u,v): return np.abs(r(u,v)*G(u,v))

      # Alternative definition of stream function
      def dlamdu0(u,v): return sqrtg(u,v)/dpsitords[ind]*bsupv(u,v,ind)-1.0
      def dlamdv0(u,v): return (-sqrtg(u,v)/(iota[ind]*dpsitords[ind])*bsupu(u,v,ind)+1.0)*iota[ind]

      # VMEC magnetic coordinates
      def uf(u,v,fsi): return u + lam(u,v)
      def sqrtgf(u,v): return sqrtg(u,v)/np.abs(1+dlamdu(u,v))
      def bsupuf(u,v,fsi): return (1+dlamdu(u,v,fsi))*bsupu(u,v,fsi)+\
          dlamdv(u,v,fsi)*bsubv(u,v,fsi)

      # store flux surface quantities
      # TODO: check sign convention. This one matches Strumberger output
      # TODO: check dV/ds and bvco enfp factor. Documentation says it's already included
      #       but it seems it needs to be added to get the correct output
      # TODO: pprime
      self.s.append(s[ind])
      self.iota.append(iota[ind])
      self.Jpol_divided_by_nper.append(-2.0*pi/mu0*bsubvB[ind]/enfp)
      self.Itor.append(-2.0*pi/mu0*bsubuB[ind])
      self.pprime.append(pprime)
      self.sqrt_g_00.append(-4.0*pi**2*vp[ind]/enfp)

      def sqrtgB(u,v): return np.abs(dpsitords[ind]*(iota[ind]*bsubuB[ind]+
                                      bsubvB[ind])/bmod2(u,v,ind))


      def hcheck(u,v): return sqrtgf(u,v)/sqrtgB(u,v) - 1.0

      # Boozer conversion
      print('Boozer conversion after Nuehrenberg/Zille')

      hmn1 = (bsubumn[ind,cond1]-1j*m1*bsubuB[ind]*lmnval[cond1])/\
              (1j*m1*(bsubuB[ind]+bsubvB[ind]/iota[ind]))
      hmn2 = (bsubvmn[ind,cond2]-1j*n2*bsubuB[ind]*lmnval[cond2])/\
              (1j*n2*(bsubuB[ind]+bsubvB[ind]/iota[ind]))

      hmn = np.zeros(m.shape, dtype='complex')
      hmn[cond2] = hmn2
      hmn[cond1] = hmn1

      def H(u,v): return fourierseries(hmn,u,v,m,n)
      def dHdu(u,v): return fourierseries(1j*m*hmn,u,v,m,n)
      def dHdv(u,v): return fourierseries(1j*n*hmn,u,v,m,n)

      def dth(u,v): return fourierseries(lmnval+hmn, u, v, m, n)
      def dph(u,v): return fourierseries(hmn/iota[ind], u, v, m, n)

      # Calculate Boozer modes
      m0b = 2*empol-1
      n0b = 4*entor+1
      mb = np.zeros(int((m0b-1)*n0b+(n0b-1)/2+1),dtype=int)
      nb = np.zeros(int((m0b-1)*n0b+(n0b-1)/2+1),dtype=int)
      k=0
      for mk in range(m0b):
        nmin0 = -2*entor
        if mk < 1:
          nmin0 = 0
        for nk in range(nmin0,2*entor+1):
          mb[k] = mk
          nb[k] = nk*enfp
          k = k+1

      nu = uv_grid_multiplicator*np.max(np.abs(m))-1
      nv = uv_grid_multiplicator*np.max(np.abs(n))+1
      du = 2.0*pi/nu
      dv = 2.0*pi/nv
      up = np.arange(0,2*pi,du)
      vp = np.arange(0,2*pi,dv)
      [U,V] = np.meshgrid(up,vp)
      U = U.flatten().T; V = V.flatten().T
      THB = U + dth(U,V)
      PHB = V + dph(U,V)
      R = r(U,V)
      Z = z(U,V)
      B = bmod(U,V)
      H1 = H(U,V)
      HCHECK = hcheck(U,V)
      Jbinv = np.abs((1.0+dlamdu(U,V)+dHdu(U,V))*(1.0+dHdv(U,V)/iota[ind])\
              -dHdu(U,V)/iota[ind]*(dlamdv(U,V)+dHdv(U,V)))

      print("Computing Boozer modes")
      dthmnb = np.zeros(mb.shape, complex)
      dphmnb = np.zeros(mb.shape, complex)
      rmnb = np.zeros(mb.shape, complex)
      zmnb = np.zeros(mb.shape, complex)
      bmnb = np.zeros(mb.shape, complex)
      hmnb = np.zeros(mb.shape, complex)
      hcheckmnb = np.zeros(mb.shape, complex)
      for km in range(len(mb)):
        efun = 2.0/(2.0*np.pi)**2*np.exp(-1.0j*(mb[km]*THB + nb[km]*PHB))*du*dv
        if (mb[km]==0) and (nb[km]==0):
          efun = efun/2.0
        dthmnb[km] = np.sum(efun*Jbinv*(THB-U))
        dphmnb[km] = np.sum(efun*Jbinv*(PHB-V))
        rmnb[km] = np.sum(efun*Jbinv*R)
        zmnb[km] = np.sum(efun*Jbinv*Z)
        bmnb[km] = np.sum(efun*Jbinv*B)
        hmnb[km] = np.sum(efun*Jbinv*H1)
        hcheckmnb[km] = np.sum(efun*Jbinv*HCHECK)

      vmnb = -enfp*dphmnb/(2*np.pi)

      self.m.append(mb)
      self.n.append(array(nb/enfp, dtype=int32))
      self.rmnc.append(+rmnb.real)
      self.rmns.append(-rmnb.imag)
      self.zmnc.append(+zmnb.real)
      self.zmns.append(-zmnb.imga)
      self.vmnc.append(+vmnb.real)
      self.vmns.append(-vmnb.imag)
      self.bmnc.append(+bmnb.real)
      self.bmns.append(-bmnb.imag)

      elapsed = time.time() - t
      print('Elapsed time: {} s'.format(elapsed))

    elapsed = time.time() - t_start
    print('Total time: {} s'.format(elapsed))


def main():
  import sys

  from . import getHeaderDataVMEC

  if (len(sys.argv) < 3):
    print("Usage:")
    print("./boozer.py infilename numberfluxsurfacesminustwo")
  else:
    infile = sys.argv[1]
    nsurf = int(sys.argv[2])
    wout_name = infile + '.bc'

    if (len(sys.argv) >= 4):
      uv_grid_multiplier = int(sys.argv[3])
    else:
      uv_grid_multiplier = 6

    shot = -1

    [nfp, psi_tor_a, aminor, Rmajor, m0b, n0b] = getHeaderDataVMEC.getHeadDataVmecNc(infile)

    write_boozer_head(wout_name, '01', shot, m0b, n0b, nsurf, nfp, psi_tor_a, aminor, Rmajor)

    for ind in range(1, nsurf+1):
      convert_to_boozer(infile, ind, wout_name, uv_grid_multiplier)


if __name__ == "__main__":
  main()
