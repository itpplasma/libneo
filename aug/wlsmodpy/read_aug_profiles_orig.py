#!/usr/bin/env python
#
# --------------------
# read_aug_profiles.py
# --------------------
#
# - make "modeling kit" from AUG shotfiles:
#    2D equilibrium (if exist) from
#         MICDU:EQB
#         AUGD:IDE
#      -> write to geqdsk file
#    1D profiles (if exist) from
#         WLS:PED
#         MICDU:PED
#         AUGD:IDA
# - plot "review_sheet" (to file)


import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colrs
import pmds as mds
from scipy.interpolate import interp1d

from eqrzfd import eqrzfd, geqdsk   # generic 2D free boundary equilibrium
from eqmdsaug import eqmdsaug       # ditto, read from AUG MDSplus server

import rprof          # generic 1D profile
import rprof_mds      # ditto, read from AUG MDSplus server


matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['text.usetex'] = False
matplotlib.rcParams['figure.figsize'] = [18,12]
matplotlib.rcParams['figure.subplot.left'] = 0.08
matplotlib.rcParams['figure.subplot.right'] = 0.98
matplotlib.rcParams['figure.subplot.bottom'] = 0.08
matplotlib.rcParams['figure.subplot.top'] = 0.98
matplotlib.rcParams['figure.subplot.wspace'] = 0.25


# ----------------- read PED experimental data (with errors) ------------------------

def read_ped_data (shot, quant, experiment='AUGD', edition='0'):
  # measured data
  dcmd = 'augdiag('+str(shot)+',"PED","'+quant+'Data","'+experiment+'",'+str(edition)+')'
  ecmd = 'augdiag('+str(shot)+',"PED","d'+quant+'Data","'+experiment+'",'+str(edition)+')'
  Data = rprof_mds.rprof_mds(aug, shot, quant+"Data", "(_data="+dcmd+")[,0]", 
                   "(_errbar="+ecmd+")[,0]", "psinorm", "dim_of(_data,0)^2")
  return(Data)

# ---------------- read PED fit profile --------------------------------------------


def read_ped_fit (shot, tbegin, tend, quant, experiment='AUGD', edition='0'):
  dcmd = 'augdiag('+str(shot)+',"PED","'+quant+'Fit","'+experiment+'",'+str(edition)+',' \
         +str(tbegin)+','+str(tend)+')'
  ecmd = 'augdiag('+str(shot)+',"PED","d'+quant+'Fit","'+experiment+'",'+str(edition)+',' \
         +str(tbegin)+','+str(tend)+')'

  print dcmd
  
  Data =   mds.mdsvalue("(_fit="+dcmd+")")
  Errbar = mds.mdsvalue("(_fiterr="+ecmd+")")
  PsiNorm  = mds.mdsvalue("dim_of(_fit,0)^2")

  ntime   = Data.shape[1]
  Data1   = np.sum(Data,   axis=0) / ntime
  Errbar1 = np.sum(Errbar, axis=0) / ntime
  
  Prof = rprof.RadialProfile(quant, Data1, Errbar1, 'psinorm', PsiNorm)
  return(Prof)


# ---------------- read IDA fit profile --------------------------------------------

def read_ida_fit (shot, tbegin, tend, quant, experiment='AUGD', edition='0'):
  dcmd = 'augdiag('+str(shot)+',"IDA","'+quant+'","'+experiment+'",'+str(edition)+',' \
         +str(tbegin)+','+str(tend)+')'
  ecmd = 'augdiag('+str(shot)+',"IDA","'+quant+'_unc","'+experiment+'",'+str(edition)+',' \
         +str(tbegin)+','+str(tend)+')'

  Data =   mds.mdsvalue("(_fit="+dcmd+")")
  Errbar = mds.mdsvalue("(_fiterr="+ecmd+")")
  PsiNorm  = mds.mdsvalue("dim_of(_fit,1)^2")

  ntime   = Data.shape[1]
  Data1   = np.sum(Data,   axis=1) / ntime
  Errbar1 = np.sum(Errbar, axis=1) / ntime
  
  Prof = rprof.RadialProfile(quant, Data1, Errbar1, 'psinorm', PsiNorm[:,0])
  return(Prof)


# ---------------- save a profile as asci file ------------------------------------

def save_profile_as_ascii(fn, p, scalefac, nonnegative=False):
  f = open(fn, 'w')
  l = len(p.value)
  for i in range(l):
    f.write(' {: 8.6f} '.format( np.sqrt(p.psinorm[i]) ))
    v = scalefac*p.value[i]
    if nonnegative and v<0:
      v=0
    f.write(' {: 14.7e}'.format(v))
    v = scalefac*p.errorbar[i]
    f.write(' {: 14.7e}'.format(v))
    f.write("\n")
  f.close()


def save_q_profile(fn, eq):
    psinorm = (eq.Psi - eq.Psimag) / (eq.Psibdry - eq.Psimag)
    psinorm[0] = 0.0
    Qprof = rprof.RadialProfile('1', eq.Qpsi, eq.Qpsi*0.0, '1', psinorm)
    save_profile_as_ascii(fn, Qprof, 1.0, nonnegative=False)


# ---------------------------- main program ---------------------------------------

# experimentcode 'None' (default) means don't read this shot file


def process_slice(shot, tbegin, tend, \
                  pedexp=None, pededition=0, \
                  idaexp=None, idaedition=0, \
                  eqbexp=None, eqbdiag='EQB', eqbedition=0, \
                  ideexp=None, ideedition=0, \
                  fnsuffix = ''):

    time = 0.5*(tbegin+tend)

    Dataprefix = os.path.expanduser('~')+'/public/Modeling_kit/AUG/'
    shottime = str(shot)+'.'+str(int(1000*time))
    profofn = Dataprefix + 'prof/ascii/'+shottime+'_'

    # IDE equilibrium
    if ideexp is not None:
      eq_ide = eqmdsaug(shot, time, experiment=ideexp, diagnostic='IDE', edition=ideedition)
      eq_ide.write_geqdsk (Dataprefix+'eqdsk/g'+shottime+'_IDE'+fnsuffix)
      save_q_profile(profofn+'q_EQB', eq_ide)
    else:
      eq_ide = None

    # EQB equilibrium
    if eqbexp is not None:
      eq_eqb = eqmdsaug(shot, time, experiment=eqbexp, diagnostic=eqbdiag, edition=eqbedition)
      eq_eqb.write_geqdsk (Dataprefix+'eqdsk/g'+shottime+'_EQB'+fnsuffix)
      save_q_profile(profofn+'q_EQB', eq_eqb)
    else:
      eq_eqb = None
      
    # IDA profiles
    if idaexp is not None:
      Te_ida_fit = read_ida_fit (shot, tbegin, tend, 'Te', \
                                   experiment=idaexp, edition=idaedition)
      if Te_ida_fit is not None:
        save_profile_as_ascii(profofn+'Te_IDA'+fnsuffix, Te_ida_fit, 1.0, nonnegative=True)

      ne_ida_fit = read_ida_fit (shot, tbegin, tend, 'ne', \
                                   experiment=idaexp, edition=idaedition)
      if ne_ida_fit is not None:
        save_profile_as_ascii(profofn+'ne_IDA'+fnsuffix, ne_ida_fit, 1.0, nonnegative=True)
    else:
      Te_ida_fit = None    
      ne_ida_fit = None    

    # PED profiles
    if pedexp is not None:
      Te_ped_fit = read_ped_fit (shot, tbegin, tend, 'Te', \
                                 experiment=pedexp, edition=pededition)
      if Te_ped_fit is not None:
        save_profile_as_ascii(profofn+'Te_PED'+fnsuffix, Te_ped_fit, 1.0e3, nonnegative=True)

      ne_ped_fit = read_ped_fit (shot, tbegin, tend, 'ne', \
                                 experiment=pedexp, edition=pededition)

      if ne_ped_fit is not None:
        save_profile_as_ascii(profofn+'ne_PED'+fnsuffix, ne_ped_fit, 1.0e3, nonnegative=True)

      Ti_ped_fit = read_ped_fit (shot, tbegin, tend, 'Ti', \
                                 experiment=pedexp, edition=pededition)

      if Ti_ped_fit is not None:
        save_profile_as_ascii(profofn+'Ti_PED'+fnsuffix, Ti_ped_fit, 1.0e3, nonnegative=True)


      vt_ped_fit = read_ped_fit (shot, tbegin, tend, 'vT', \
                                 experiment=pedexp, edition=pededition)

      if vt_ped_fit is not None:
        save_profile_as_ascii(profofn+'vt_PED'+fnsuffix, vt_ped_fit, 1.0e3, nonnegative=False)
    else:
      Te_ped_fit = None    
      ne_ped_fit = None    
      Ti_ped_fit = None
      vt_ped_fit = None    


    # plot results

    if True:
    # multi-panel plot
      rmin = 0.0
      rmax = 1.1
      fig = plt.figure()

    #electron temperature
      plt.subplot(231)
      if Te_ped_fit is not None:
        Te_ped_fit.plot_errors(1.0, label='PED', color='red')
      if Te_ida_fit is not None:
        Te_ida_fit.plot_errors(1e3, label='IDA', color='blue')
      plt.xlabel('$\Psi_n$')
      plt.ylabel('electron temperature [keV]')
      plt.legend(loc='best')
      plt.axis([rmin, rmax, 0.01, 4.0])

    #ion temperature
      plt.subplot(232)
      if Ti_ped_fit is not None:
        Ti_ped_fit.plot_errors(1.0, label='PED', color='red')
      plt.xlabel('$\Psi_n$')
      plt.ylabel('ion temperature [keV]')
      plt.legend(loc='best')
      plt.axis([rmin, rmax, 0.01, 4.0])

    #electron density
      plt.subplot(234)
      if ne_ped_fit is not None:
        ne_ped_fit.plot_errors(1e16, label='PED', color='red')
      if ne_ida_fit is not None:
        ne_ida_fit.plot_errors(1e19, label='IDA', color='blue')
      plt.xlabel('$\Psi_n$')
      plt.ylabel('electron density [keV]')
      plt.legend(loc='best')
      plt.axis([rmin, rmax, 0.01, 7.0])

    #ion toroidal rotation
      plt.subplot(235)
      if vt_ped_fit is not None:
        vt_ped_fit.plot_errors(1.0, label='PED', color='red')
      plt.xlabel('$\Psi_n$')
      plt.ylabel('toroidal rotation frequency [krad/s]')
#      plt.legend(loc='best')
      plt.axis([rmin, rmax, -20.0, 100.0])

    #equilibrium pressure
      plt.subplot(233)
      if eq_ide is not None:
        Psinorm = (eq_ide.Psi-eq_ide.Psimag) / (eq_ide.Psibdry-eq_ide.Psimag)
        plt.plot(Psinorm, eq_ide.Pres*1e-3, label='IDE', color='black')
      if eq_eqb is not None:
        Psinorm = (eq_eqb.Psi-eq_eqb.Psimag) / (eq_eqb.Psibdry-eq_eqb.Psimag)
        plt.plot(Psinorm, eq_eqb.Pres*1e-3, label='EQB', color='magenta')
      if ne_ida_fit is not None:
        plt.plot(ne_ida_fit.psinorm, 2*1.602e-19*1e-3*ne_ida_fit.value*Te_ida_fit.value, \
                 label='IDA, 2*Te', color='blue')
      if ne_ped_fit is not None:
        plt.plot(ne_ped_fit.psinorm, 1.602e-19*ne_ped_fit.value*1e3*(Te_ped_fit.value+Ti_ped_fit.value), \
                 label='PED, Ti+Te', color='red')
      plt.xlabel('$\Psi_n$')
      plt.ylabel('plasma pressure [kPa]')
      plt.legend(loc='best')
      plt.axis([rmin, rmax, -20.0, 100.0])

    # toroidal current profile 
      plt.subplot(236)
      if eq_ide is not None:
        jtor = eq_ide.ToroidalCurrentDensity()
        nz = int(jtor.shape[1]/2)  # midplane section
        plt.plot(eq_ide.R[1:-1], 1e-6*jtor[:,nz], label='IDE', color='black')
      if eq_eqb is not None:
        jtor = eq_eqb.ToroidalCurrentDensity()
        nz = int(jtor.shape[1]/2)  # midplane section
        plt.plot(eq_eqb.R[1:-1], 1e-6*jtor[:,nz], label='EQB', color='magenta')
      plt.xlabel('R [m]')
      plt.ylabel('jtor [MA/m^2]')
      plt.legend(loc='best')

    # plt.show()
    fig.savefig(Dataprefix+'plots/prof_'+str(shot)+'_'+str(int(1000*time))+fnsuffix+'.eps', bbox_inches='tight')


# ---------------------------- main ----------------------------------


# ........... dithers between ELMy and suppression ............

# ELM suppression
#process_slice(33120, 5.46, 5.54, \
#              pedexp='MICDU', pededition=1, \
#              idaexp='AUGD', idaedition=3, \
#              eqbexp='MICDU', eqbedition=1, \
#              ideexp='AUGD', ideedition=1)

# small ELMs
#process_slice(33120, 5.625, 5.645, \
#              pedexp='MICDU', pededition=2, \
#              idaexp='AUGD', idaedition=3, \
#              eqbexp='MICDU', eqbedition=1, \
#              ideexp='AUGD', ideedition=1)


# ........... "outlier" with small Ti-ped ........
#process_slice(34213, 2.86, 2.87, \
#              pedexp='MICDU', pededition=1, \
#              idaexp='AUGD', idaedition=1, \
#              eqbexp='MICDU', eqbedition=1, \
#              ideexp=None, ideedition=1)


# ........... IB exit threshold ............
# Delta phi = -45 degrees

# ........... rigid MP rotation ............
#process_slice(34214, 2.66, 2.75, \
#              pedexp='WLS', pededition=1, \
#              idaexp='AUGD', idaedition=1, \
#              eqbexp=None, eqbedition=1, \
#              ideexp='AUGD', ideedition=1)

# ........... rigid MP rotation ............
process_slice(34548, 5.595, 5.645, \
              pedexp='MICDU', pededition=1, \
              idaexp='AUGD', idaedition=1, \
              eqbexp='MICDU', eqbedition=1, \
              ideexp='AUGD', ideedition=1)


#process_slice(34832, 5.45, 5.5, \
#              pedexp='WLS', pededition=1, \
#              idaexp='AUGD', idaedition=1, \
#              eqbexp='MICDU', eqbedition=1, \
#              ideexp='AUGD', ideedition=2)

# Delta Phi = +45 degrees
#process_slice(34834, 6.35, 6.4, \
#              pedexp='WLS', pededition=1, \
#              idaexp='AUGD', idaedition=1, \
#              eqbexp='MICDU', eqbedition=1, \
#              ideexp='AUGD', ideedition=2)

# Delta Phi = +135 degrees
#process_slice(34834, 3.8, 3.85, \
#              pedexp='WLS', pededition=1, \
#              idaexp='AUGD', idaedition=1, \
#              eqbexp='MICDU', eqbedition=2, \
#              ideexp='AUGD', ideedition=2)

# ........... triangularity scan ............
# delta_u ~ 0.23
#process_slice(34835, 3.0, 3.05, \
#              pedexp='WLS', pededition=1, \
#              idaexp='AUGD', idaedition=1, \
#              eqbexp='MICDU', eqbedition=1, \
#              ideexp='AUGD', ideedition=2)

# delta_u ~ 0.07
#process_slice(34835, 4.95, 5.0, \
#              pedexp='WLS', pededition=2, \
#              idaexp='AUGD', idaedition=1, \
#              eqbexp='MICDU', eqbedition=1, \
#              ideexp='AUGD', ideedition=2)

# .... q95 ramp for pump-out ..........

# EQB equilibria with SOL currents (ed=1)
if False:
  process_slice(34381, 2.35, 2.45, \
                pedexp='MICDU', pededition=4, \
                eqbexp='MICDU', eqbedition=1, \
                idaexp='AUGD', idaedition=1)

  process_slice(34381, 3.85, 3.95, \
                pedexp='MICDU', pededition=3, \
                eqbexp='MICDU', eqbedition=1, \
                idaexp='AUGD', idaedition=1)

  process_slice(34381, 5.20, 5.30, \
                pedexp='MICDU', pededition=1, \
                eqbexp='MICDU', eqbedition=1, \
                idaexp='AUGD', idaedition=1)

# EQB equilibria without SOL currents (ed=3)
if False:
  process_slice(34381, 2.35, 2.45, \
                eqbexp='MICDU', eqbedition=3, \
                fnsuffix = '_noSOLcur')

  process_slice(34381, 3.85, 3.95, \
                eqbexp='MICDU', eqbedition=3, \
                fnsuffix = '_noSOLcur')

  process_slice(34381, 5.20, 5.30, \
                eqbexp='MICDU', eqbedition=3, \
                fnsuffix = '_noSOLcur')

if False:
  process_slice(34381, 2.35, 2.45, \
                pedexp='NLEUTHOL', pededition=9, \
                eqbexp='MICDU', eqbedition=1, \
                fnsuffix = '_nleuthol')

  process_slice(34381, 3.85, 3.95, \
                pedexp='NLEUTHOL', pededition=10, \
                eqbexp='MICDU', eqbedition=1, \
                fnsuffix = '_nleuthol')

  process_slice(34381, 5.20, 5.30, \
                pedexp='NLEUTHOL', pededition=11, \
                eqbexp='MICDU', eqbedition=1, \
                fnsuffix = '_nleuthol')

if False:
  process_slice(35712, 2.66, 2.69, \
                pedexp='MWILLENS', pededition=16, \
                eqbexp='MICDU', eqbedition=2, \
                fnsuffix = '_ped16')
  process_slice(35712, 4.5, 4.6, \
                pedexp='MWILLENS', pededition=15, \
                fnsuffix = '_ped15')
