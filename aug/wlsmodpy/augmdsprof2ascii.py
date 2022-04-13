#!/usr/bin/env python

'''
  Read ASDEX Upgrade (fitted) radial profiles through MDSplus
  and save them in various (ASCII-, netCDF-) formats
'''

import numpy as np
import matplotlib.pyplot as plt
import mdsplus
import rprof


class rprof_augmds(rprof.RadialProfile):
 
  def __init__(self, conn, shot, diag, signame, tbegin, tend, \
               experiment='AUGD', edition=0, verbose=False):
    '''
        AUG fitted profiles, read through MDSplus
        (RadialProfile subclass)
        - read data and error bars from diag='PED' or diag='IDA'
        - average over time from tbegin to tend
        - make rprof instance
        - return None if unsuccessful
    '''

    if diag=='IDA':
        esigname = signame + '_unc'
    elif diag=='PED':
        esigname = 'd'+signame
    else:
        print('Diagnostic is %s, but must be either "IDA" or "PED"' % diag)
        return None

    try:    # profile data
        data2, tbase, rhop2, unit, oshot, oedition \
            = conn.augsignal(shot, diag, signame,  \
            experiment=experiment, edition=edition, t1=tbegin, t2=tend)
    except:
        print("Cannot read %d %s:%s (%d)" % (shot,experiment,diag,edition))
        return None

    try:    # error bars
        errbar2, tbase_, rhop_, unit_, oshot_, oedition_ \
             = conn.augsignal(shot, diag, esigname, \
             experiment=experiment, edition=edition, t1=tbegin, t2=tend)
    except:
        print("Cannot read %d %s:%s (%d)" % (shot,experiment,diag,edition))
        return None

    if np.any(tbase_.shape != tbase.shape) and np.any(tbase_ != tbase):
        print("Signal (%s) and Error bars (%s) do not have the same time base." % (signame, esigname))
        return None
    if np.any(rhop_.shape != rhop2.shape) and np.any(rhop_ != rhop2):
        print("Signal (%s) and Error bars (%s) do not have the same radial base." % (signame, esigname))
        return None

    ntime = data2.shape[0]
    data1 = np.sum(data2, axis=0) / ntime    
    if rhop2.ndim>1:   # 2-dimensional area base as in IDA
        psin1 = np.square( np.sum(rhop2, axis=0) / ntime )
    else:   # 1-dimensional are base as in PED
        psin1 = np.square(rhop2)
    errbar1 = np.sum(errbar2, axis=0) / ntime
    
    if verbose:
        print("Read %d %s:%s (%d) %s" % (shot,experiment,diag,edition,signame))
        print("  Time range requested: %5.3f - %5.3f s, found: %5.3f - %5.3f s" \
              % (tbegin,tend,np.min(tbase), np.max(tbase)))
        print("  Psi(normalised) range found: %5.3f - %5.3f s" \
              % (np.min(psin1), np.max(psin1)))    

    super(rprof_augmds, self).__init__(signame, data1, errbar1, 'psinorm', psin1)  
    


# ------------------- main -----------------------------------------------------

if __name__ == '__main__':
  import argparse
  import getpass
  import sys

  parser = argparse.ArgumentParser(description='Read AUG profiles through MDSplus and write as ASCII.')

  parser.add_argument('shot', metavar='shot', type=int, nargs='?',
                      help='AUG shot number')  
  parser.add_argument('-s', '--server', dest='server', action='store', \
                      default='mdsplus.aug.ipp.mpg.de', help='MDSplus server')
  parser.add_argument('-t1', '--tbegin', dest='tbegin', action='store', type=float, \
                      default='2.0', help='start of time interval of interest [s]')
  parser.add_argument('-t2', '--tend', dest='tend', action='store', type=float, \
                      default='3.0', help='end of time interval of interest [s]')

  parser.add_argument('-d', '--diag', dest='diag', action='store', \
                      default='IDA', help='shotfile diagnostic code (PED, IDA)')
  parser.add_argument('-exp', '--experiment', dest='experiment', action='store', \
                      default='AUGD', help='shotfile experiment code (AUGD, self, etc.)')
  parser.add_argument('-ed', '--edition', dest='edition', type=int, action='store', \
                      default=0, help='shotfile edition number')

  parser.add_argument('-p', '--plot', dest='plot_profiles', action="store_true",
                      help="plot profiles on screen")
  parser.add_argument('-v', '--verbose', dest='verbose', action="store_true",
                      help="print progress messages (not only errors)")

  args = parser.parse_args()

# inject parameters for testing in Spyder etc.
  if True:
      args.shot = 34548
      args.experiment = 'AUGD'
      args.diag = 'IDA'
      args.edition = 1
      args.tbegin = 5.6
      args.tend = 5.7
      args.verbose=True
      args.plot_profiles=True

  if args.shot is None:
    print('No shot number given! Abort.')
    sys.exit()

  if args.experiment=='self' or args.experiment=='SELF':
    args.experiment=getpass.getuser()

  print('Shot: %d  Diagnostic: %s  Experiment: %s  Edition: %d' \
          % (args.shot, args.diag, args.experiment, args.edition))
  print('Requested time interval: %6.3f - %6.3f s' % (args.tbegin, args.tend))
    
  print('Connecting to MDSplus server: %s' % args.server)
  conn = mdsplus.mdsaug(args.server)

  if args.diag=='PED':
    TeFit  = rprof_augmds(conn, args.shot, args.diag, 'TeFit',  args.tbegin, args.tend, \
                          experiment=args.experiment, edition=args.edition, \
                          verbose=args.verbose)
    TeData = rprof_augmds(conn, args.shot, args.diag, 'TeData', args.tbegin, args.tend, \
                          experiment=args.experiment, edition=args.edition, \
                          verbose=args.verbose)
    neFit  = rprof_augmds(conn, args.shot, args.diag, 'neFit',  args.tbegin, args.tend, \
                          experiment=args.experiment, edition=args.edition, \
                          verbose=args.verbose)
    neData = rprof_augmds(conn, args.shot, args.diag, 'neData', args.tbegin, args.tend, \
                          experiment=args.experiment, edition=args.edition, \
                          verbose=args.verbose)
    TiFit  = rprof_augmds(conn, args.shot, args.diag, 'TiFit',  args.tbegin, args.tend, \
                          experiment=args.experiment, edition=args.edition, \
                          verbose=args.verbose)
    TiData = rprof_augmds(conn, args.shot, args.diag, 'TiData', args.tbegin, args.tend, \
                          experiment=args.experiment, edition=args.edition, \
                          verbose=args.verbose)
    vtFit  = rprof_augmds(conn, args.shot, args.diag, 'vTFit',  args.tbegin, args.tend, \
                          experiment=args.experiment, edition=args.edition, \
                          verbose=args.verbose)
    vtData = rprof_augmds(conn, args.shot, args.diag, 'vTData', args.tbegin, args.tend, \
                          experiment=args.experiment, edition=args.edition, \
                          verbose=args.verbose)
    vpFit  = rprof_augmds(conn, args.shot, args.diag, 'vPFit',  args.tbegin, args.tend, \
                          experiment=args.experiment, edition=args.edition, \
                          verbose=args.verbose)
    vpData = rprof_augmds(conn, args.shot, args.diag, 'vPData', args.tbegin, args.tend, \
                          experiment=args.experiment, edition=args.edition, \
                          verbose=args.verbose)

    if args.plot_profiles:
        plt.subplot(231)
        if TeData is not None:
            TeData.plot_errors(1e3, color='grey', linewidth=0)
        if TeFit is not None:
            TeFit.plot(1e3, color='red')
        plt.xlabel('psi-norm')
        plt.ylabel('Te [keV]')
        plt.subplot(232)
        if TiData is not None:
            TiData.plot_errors(1e3, color='grey', linewidth=0)
        if TiFit is not None:
            TiFit.plot(1e3, color='red')
        plt.xlabel('psi-norm')
        plt.ylabel('Ti [keV]')
        plt.subplot(233)
        if vtData is not None:
            vtData.plot_errors(1e3, color='grey', linewidth=0)
        if vtFit is not None:
            vtFit.plot(1e3, color='red')
        plt.xlabel('psi-norm')
        plt.ylabel('omega_tor [krad/s]')
        plt.subplot(234)
        if neData is not None:
            neData.plot_errors(1e19, color='grey', linewidth=0)   
        if neFit is not None:
            neFit.plot(1e19, color='red')
        plt.xlabel('psi-norm')
        plt.ylabel('ne [10^19 m^-3]')
        plt.subplot(236)
        if vpData is not None:
            vpData.plot_errors(1e3, color='grey', linewidth=0)
        if vpFit is not None:
            vpFit.plot(1e3, color='red')
        plt.xlabel('psi-norm')
        plt.ylabel('omega_pol [krad/s]')
        plt.show()


  elif args.diag=='IDA':
    TeFit  = rprof_augmds(conn, args.shot, args.diag, 'Te', args.tbegin, args.tend, \
                experiment=args.experiment, edition=args.edition, \
                verbose=args.verbose)
    TeData = None
    neFit  = rprof_augmds(conn, args.shot, args.diag, 'ne', args.tbegin, args.tend, \
                experiment=args.experiment, edition=args.edition, \
                verbose=args.verbose)
    neData = None
    if args.plot_profiles:
        plt.subplot(121)
        # TeData.plot_errors(1e3, color='grey', linewidth=0)
        if TeFit is not None:
          TeFit.plot(1e3, color='red')
        plt.xlabel('psi-norm')
        plt.ylabel('Te [keV]')
        plt.subplot(122)
        # neData.plot_errors(1e19, color='grey', linewidth=0)   
        if neFit is not None:
            neFit.plot(1e19, color='red')
        plt.xlabel('psi-norm')
        plt.ylabel('ne [10^19 m^-3]')
        plt.show()
    
  else:
    print('Unsupported diagnostic: %s -- Abort.' % args.diag)
    sys.exit()


