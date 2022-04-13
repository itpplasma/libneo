#!/usr/bin/env python
#
# --------------------
# aug2equilf.py
# --------------------
#
# read AUG equilibrium and PF coil data 
# to produce input files for "equilf" (Lackner forward [predictive] equilibrium code):
# .in       - input control
# .coils    - PF coil locations and currents
# .limiter  - limiter contours
# .fluxprof - Pres, Fpol flux functions


import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colrs
import pmds as mds
from scipy.interpolate import interp1d

from eqrzfd import eqrzfd, geqdsk   # generic 2D free boundary equilibrium
from eqmdsaug import eqmdsaug       # ditto, read from AUG MDSplus server
from augpfcoils import pfcircuits, pfcoils_geo, pfcoils_nturns


DTL = 0.02   # integration length for self inductance
Dataprefix = os.path.expanduser('~')+'/public/Modeling_kit/AUG/'
maidiag = 'MAY'


# ------------------ read signal --------------------------

def read_value_1(shot, tbegin, tend, diag, signame, experiment='AUGD', edition=0):
  cmd = 'augdiag('+str(shot)+',"'+diag+'","'+signame+'","'+experiment+'",'+str(edition)+',' \
         +str(tbegin)+','+str(tend)+')'
  trace = mds.mdsvalue(cmd)
  ntime = trace.shape[0]
  v     = np.sum(trace, axis=0) / ntime
  return(v)
  

def read_value(shot, tbegin, tend, diag, signame, \
                experiment='AUGD', edition=0, tzero_beg=None, tzero_end=None):
  value = read_value_1(shot, tbegin, tend, diag, signame, experiment, edition)
  if (tzero_beg is not None) and (tzero_end is not None):
    offset = read_value_1(shot, tzero_beg, tzero_end, diag, signame, experiment, edition)
  else:
    offset = 0.0
  return(value - offset)




# ----------------- write case data --------------------

def write_case(f, casename, eq, Iplasma, betpol, prescribe_fluxprof=True):

  f.write('{0}\n'.format(casename))
  if prescribe_fluxprof:   # use prescribed flux functions
    IW = 2
  else:                    # Gruber (parabolic) profile parametrisation
    IW = 0
  C2 = Iplasma/1e6 * (1-betpol)
  C1 = Iplasma/1e6 - C2
  print "Case: ", casename, " C1=", C1, " C2=", C2
  # FVAC,IW,C1,C2,C3,C4,LL    - For Gruber's formula, need C3=C4=0 so that CURR=C1+C2
  f.write('{:6.4f} {:1d} {:8.6f} {:8.6f}  0.0  0.0 1.00\n'.format(eq.Fpol[-1], IW, C1, C2))
  # RAXIS0, ZAXIS0, ILDES, RDES0, ZDES0
  f.write('{:8.6f} {:8.6f} 3  0.0 -0.4\n'.format(eq.Rmaxis, eq.zmaxis))

  # write .fluxprof file
  if prescribe_fluxprof:
    fn2 = Dataprefix+'equilf/'+casename+'.fluxprof'
    f2 = open(fn2, 'w')
    npsi = len(eq.Psi)
    Psin = (eq.Psi - eq.Psimag) / (eq.Psibdry - eq.Psimag) 
    f2.write('{:3d}\n'.format(npsi))
    for i in range(npsi):
      f2.write('{:12.9f} {:12.5f} {:12.9f}\n'.format(Psin[i], eq.Pres[i], eq.Fpol[i]))    
    f2.close()
    print "Flux profile file "+fn2+" written"

# end of write_case()


def process_equil(f, shot, tbegin, tend, casename, eqexp='AUGD', eqdiag='IDE', egdiag='IDG', \
                  eqedition = 0, prescribe_fluxprof=True):
  time = 0.5*(tbegin+tend)
  eq = eqmdsaug(shot, time, experiment=eqexp, diagnostic=eqdiag, edition=eqedition)
  Iplasma = read_value(shot, tbegin, tend, "FPC", "IpiFP")

  if egdiag is not None:
    betpol = read_value(shot, tbegin, tend, egdiag, "betpol", experiment=eqexp)
    print "betpol (from "+egdiag+" shotfile) = ", betpol
  else:
    PSIN = (eq.Psi-eq.Psi[0]) / (eq.Psi[-1]-eq.Psi[0])
    PsiRzN = (eq.PsiRz[1:-1,1:-1]-eq.Psi[0]) / (eq.Psi[-1]-eq.Psi[0])
#  wh=np.logical_or(PsiRzN<0.0, PsiRzN>1.0)
    wh=(PsiRzN>1.0)
    PsiRzN[wh] = 0.0   # outside plasma

    Pres_CS  = interp1d(PSIN, eq.Pres, kind='linear', fill_value='extrapolate')
    nz1 = np.ones(eq.z.size-2)
    R2  = np.outer(eq.R[1:-1], nz1)
    dWmhd = 1.5 *2*np.pi*R2 *0.5*eq.PlasmaDomain * Pres_CS(PsiRzN)
    dWmhd[wh]= 0.0
    Wmhd = np.sum(dWmhd*eq.dR*eq.dz)
    print "Wmhd: ", 1e-6*Wmhd, " MJ"

    dArea= 0.5*eq.PlasmaDomain
    dArea[wh] = 0.0
    Area = np.sum(dArea*eq.dR*eq.dz)
    betpol = 1.35*Wmhd/Iplasma**2/(4*np.pi*1e-7)*2/Area
    print "betpol (calc from pressure profile) = ", betpol

  print eqexp+':'+eqdiag+' '+str(shot)+'('+str(eqedition)+')'
  print "  Rmag = ", eq.Rmaxis, "  zmag = ", eq.zmaxis, " m"
  print "  betpol = ", betpol
  print "  Plasma current (MA)   FPC: ", Iplasma/1e6, " Equilibrium: ", eq.Iplasma/1e6, " MA"
  write_case(f, casename, eq, Iplasma, betpol, prescribe_fluxprof)



# --------- make all files for an equilf run, i.e. one shot/time -----------

def prepare_equilf_run(shot, tbegin, tend, \
                       eqexp = None, eqdiag='EQH',  eqedition=0, \
                       eqbexp= None, eqbdiag='EQB', eqbedition=0, \
                       ideexp= None, idediag='IDE', ideedition=0  ):

  time = 0.5*(tbegin+tend)
  shottime = str(shot)+'.'+str(int(1000*time))
  Equilprefix = Dataprefix + 'equilf/'+shottime

# ---------------------- read PF coil currents ----------------------
  Icircuits = {}
  Icoils = {}

  for c in pfcircuits.keys():
    Icircuits[c] = read_value(shot, tbegin, tend, maidiag, c, tzero_beg=-10, tzero_end=-9)

  for c in pfcircuits.keys():
    coils = pfcircuits[c]
    print c, ' I=', Icircuits[c]/1e3, ' kA'  
    for d in coils:
      try:
        Icoils[d] += Icircuits[c]
      except:
        Icoils[d] = Icircuits[c]
  #    print '  ', d, Icoils[d]/1e3


  # ----------------- write .coils file --------------------

  f = open(Equilprefix+'.coils', 'w')

  # passive coils (PSL)
  for c in pfcoils_geo.keys():
    if c=='PSLo' or c=='PSLu':
      Rzt = np.array(pfcoils_geo[c])
      for i in range(Rzt.shape[0]):
        f.write('{: 7.4f} {: 7.4f} {: 7.4f}\n'.format(Rzt[i,0], Rzt[i,1], DTL))
  f.write('{: 7.1f} {: 7.1f} {: 7.1f}\n'.format(0.0, 0.0, 0.0))

  # "symmetric" conductors (vertical field = radial position control coils) => V2o,u
  for c in pfcoils_geo.keys():
    if c=='V2o' or c=='V2u':
      signturns = 1.0    # / pfcoils_nturns[d]
      Rzt = np.array(pfcoils_geo[c])
      for i in range(Rzt.shape[0]):
        f.write('{: 7.4f} {: 7.4f} {: 7.4f}\n'.format(Rzt[i,0], Rzt[i,1], signturns*Rzt[i,2]))
  f.write('{: 7.1f} {: 7.1f} {: 7.1f}\n'.format(0.0, 0.0, 0.0))

  # "anti-symmetric" conductors (radial field = vertical position control coils) => V2o,u
  for c in pfcoils_geo.keys():
    if c=='V2o' or c=='V2u':
      signturns = 1.0    # / pfcoils_nturns[c]
      if c=='V2u': signturns = -signturns
      Rzt = np.array(pfcoils_geo[c])
      for i in range(Rzt.shape[0]):
        f.write('{: 7.4f} {: 7.4f} {: 7.4f}\n'.format(Rzt[i,0], Rzt[i,1], signturns*Rzt[i,2]))
  f.write('{: 7.1f} {: 7.1f} {: 7.1f}\n'.format(0.0, 0.0, 0.0))

  # "unsymmetric" conductors (all coils outside FD grid with given currents) 
  for c in pfcoils_geo.keys():
    if c<>'PSLo' and c<>'PSLu' and c<>'CoIo' and c<>'CoIu':
      Icond = Icoils[c] / 1e6
      Rzt = np.array(pfcoils_geo[c])
      for i in range(Rzt.shape[0]):
        f.write('{: 7.4f} {: 7.4f} {: 7.4f}\n'.format(Rzt[i,0], Rzt[i,1], Icond*Rzt[i,2]))
  f.write('{: 7.1f} {: 7.1f} {: 7.1f}\n'.format(0.0, 0.0, 0.0))

  f.close()
  print "Coils file written"


  # ----------------- write .limiter file --------------------

  eq = eqmdsaug(shot, time, diagnostic='EQI')
  f = open(Equilprefix+'.limiter', 'w')
  nlim = len(eq.Rlim)
  f.write('{:3d}\n'.format(nlim))
  for i in range(nlim):
    f.write('{:12.9f} {:12.9f}\n'.format(eq.Rlim[i], eq.zlim[i]))
  f.close()
  print "Limiter file written"


  # ----------------- write .in file --------------------

  f = open(Equilprefix+'.in', 'w')
  f.write('{0}\n'.format(shottime))
  f.write('{0}\n'.format(shottime))
  f.write('1  2  0  20  10  5.4\n')                     # ILOG,IFLUXMAT,IPLOT,NLI,NLO,SCALEP
  f.write('6  7  0.80000  2.3500  -1.5000  1.5000\n')   # IJ,NJ,RMIN,RMAX,ZMIN,ZMAX
  f.write('30 0.0001\n')                                # MAXIT,ERROR

  if eqdiag is not None:
    process_equil(f, shot, tbegin, tend, 'g'+shottime+'_EQH_fwd',\
                  eqexp=eqexp, eqdiag=eqdiag, egdiag='GQH', eqedition=eqedition)
  if eqbdiag is not None:
    process_equil(f, shot, tbegin, tend, 'g'+shottime+'_EQB_fwd',\
                  eqexp=eqbexp, eqdiag=eqbdiag, egdiag=None, eqedition=eqbedition)
  if idediag is not None:
    process_equil(f, shot, tbegin, tend, 'g'+shottime+'_IDE_fwd',\
                  eqexp=ideexp, eqdiag=idediag, egdiag='IDG', eqedition=ideedition)

  print "Input file written"
  f.close()



# ========================= main =====================================

#prepare_equilf_run(34213, 2.86, 2.87, \
#                   eqexp = 'AUGD',  eqdiag='EQH',  eqedition=0, \
#                   eqbexp= 'MICDU', eqbdiag='EQB', eqbedition=1, \
#                   ideexp= None,  idediag= None, ideedition=0  )

prepare_equilf_run(34398, 3.50, 3.60, \
                   eqexp = 'AUGD',  eqdiag='EQH',  eqedition=0, \
                   eqbexp= None,  eqbdiag=None, eqbedition=0, \
                   ideexp= None,  idediag= None, ideedition=0  )

prepare_equilf_run(34398, 4.05, 4.15, \
                   eqexp = 'AUGD',  eqdiag='EQH',  eqedition=0, \
                   eqbexp= None,  eqbdiag=None, eqbedition=0, \
                   ideexp= None,  idediag= None, ideedition=0  )

#prepare_equilf_run(34832, 5.45, 5.50, \
#                   eqexp = 'AUGD',  eqdiag='EQH',  eqedition=0, \
#                   eqbexp= 'MICDU', eqbdiag='EQB', eqbedition=1, \
#                   ideexp= 'AUGD',  idediag='IDE', ideedition=2  )

#prepare_equilf_run(34834, 3.8, 3.85, \
#                   eqexp = 'AUGD',  eqdiag='EQH',  eqedition=0, \
#                   eqbexp= 'MICDU', eqbdiag='EQB', eqbedition=2, \
#                   ideexp= 'AUGD',  idediag='IDE', ideedition=2  )

#prepare_equilf_run(34834, 6.35, 6.4, \
#                   eqexp = 'AUGD',  eqdiag='EQH',  eqedition=0, \
#                   eqbexp= 'MICDU', eqbdiag='EQB', eqbedition=1, \
#                   ideexp= 'AUGD',  idediag='IDE', ideedition=2  )

#prepare_equilf_run(34835, 3.0, 3.05, \
#                   eqexp = 'AUGD',  eqdiag='EQH',  eqedition=0, \
#                   eqbexp= 'MICDU', eqbdiag='EQB', eqbedition=1, \
#                   ideexp= 'AUGD',  idediag='IDE', ideedition=2  )

#prepare_equilf_run(34835, 4.95, 5.0, \
#                   eqexp = 'AUGD',  eqdiag='EQH',  eqedition=0, \
#                   eqbexp= 'MICDU', eqbdiag='EQB', eqbedition=1, \
#                   ideexp= 'AUGD',  idediag='IDE', ideedition=2  )
