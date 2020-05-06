#!/usr/bin/env python3

'''
  augmdsprof2netcdf.py

  Read ASDEX Upgrade (fitted) radial profiles through MDSplus
  and save them in netCDF formats
'''

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import mdsplus
import rprof


class rprof_augmds(rprof.RadialProfile):

  @classmethod
  def from_shotfile(cls, conn, shot, diag, signame, tbegin, tend, \
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
        errbar2, tbase_, rhop_, unit_, oshot, oedition \
             = conn.augsignal(shot, diag, esigname, \
             experiment=experiment, edition=edition, t1=tbegin, t2=tend)
    except:
        print("Cannot read %d %s:%s (%d)" % (shot,experiment,diag,edition))
        return None

    if np.any(tbase_.shape != tbase.shape) and np.any(tbase_ != tbase):
        print("Signal (%s) and Error bars (%s) do not have the same time base."\
              % (signame, esigname))
        return
    if np.any(rhop_.shape != rhop2.shape) and np.any(rhop_ != rhop2):
        print("Signal (%s) and Error bars (%s) do not have the same radial base."\
              % (signame, esigname))
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
        print("  Actual edition: %d  Unit: %s" % (oedition, unit))
        print("  Time range requested: %5.3f - %5.3f s, found: %5.3f - %5.3f s"\
              % (tbegin,tend,np.min(tbase), np.max(tbase)))
        print("  Psi(normalised) range found: %5.3f - %5.3f s" \
              % (np.min(psin1), np.max(psin1)))

    c = cls(unit, data1, errbar1, 'psinorm', psin1)
    c.shot = oshot
    c.edition = oedition
    c.time = 0.5*(tbegin + tend)
    return c


# -----------------------------------------------------------------------------

def plot_profiles(pd):
    # plot parameter
    matplotlib.rcParams['figure.figsize'] = [14,10]
    matplotlib.rcParams['figure.subplot.left'] = 0.08
    matplotlib.rcParams['figure.subplot.right'] = 0.96
    matplotlib.rcParams['figure.subplot.bottom'] = 0.08
    matplotlib.rcParams['figure.subplot.top'] = 0.96
    matplotlib.rcParams['figure.subplot.wspace'] = 0.15
    matplotlib.rcParams['figure.subplot.hspace'] = 0.2

    # determine plot layout
    k = pd.keys()
    sk = set([x[:2] for x in k])   # unique set of "Te", "Ti", "ne", ...
    npanels = len(sk)
    maxny = 2    # maximum number of traces stacked
    nxplot = int((npanels+maxny-1) / maxny)
    nyplot = int((npanels+nxplot-1) / nxplot)
    
    # iterate plot quantities
    i = 0
    for s in sk:
        i+=1
        plt.subplot("%1d%1d%1d"% (nyplot, nxplot, i))
        plt.title(s)
        if s[0]=='n':
            scale = 1e19
            yunit = '1e19 m^-3'
        elif s[0]=='T':
            scale = 1e3
            yunit = 'keV'
        elif s[0]=='v':
            scale = 1e3
            yunit = 'krad/s'
        else:
            scale = 1
        if s+'Data' in k:
            p = pd[s+'Data']
            p.plot_errors(scale, label='Data')
        if s+'Fit' in k:
            p = pd[s+'Fit']
            p.plot(scale, color='red', label='Fit')
        plt.xlabel('psi-norm')
        plt.ylabel(yunit)
        plt.legend(loc='best')
    plt.show()

#    plt.subplot(231)
#    if TeData is not None:
#        TeData.plot_errors(1e3, color='grey', linewidth=0)
#    if TeFit is not None:
#        TeFit.plot(1e3, color='red')
#    plt.xlabel('psi-norm')
#    plt.ylabel('Te [keV]')
#    plt.show()


# ------------------- main ----------------------------------------------------

if __name__ == '__main__':
  import os        # .path.expanduser(),
  import argparse
  import getpass   # .getuser()
  import sys

  parser = argparse.ArgumentParser( \
        description='Read AUG profiles through MDSplus and write to file(s).')

  parser.add_argument('shot', metavar='shot', type=int, nargs='?',
                      help='AUG shot number')

# ................ MDSplus options ............................................
  parser.add_argument('-s', '--server', dest='server', action='store', \
                      default='localhost', help='MDSplus server')
  parser.add_argument('-as', '--alt-server', dest='altserver', \
                      action='store', \
                      default='mdsplus.aug.ipp.mpg.de', \
                      help='alternative MDSplus server')

  # ............. shot file options ...........................................
  parser.add_argument('-t1', '--tbegin', dest='tbegin', action='store', \
                      type=float, default='2.0', \
                      help='start of time interval of interest [s]')
  parser.add_argument('-t2', '--tend', dest='tend', action='store', \
                      type=float, default='3.0', \
                      help='end of time interval of interest [s]')

  parser.add_argument('-d', '--diag', dest='diag', \
                      action='store', default='IDA', \
                      help='shotfile diagnostic code (PED, IDA)')
  parser.add_argument('-exp', '--experiment', dest='experiment', \
                      action='store', default='AUGD', \
                      help='shotfile experiment code (AUGD, self, etc.)')
  parser.add_argument('-ed', '--edition', dest='edition', \
                      type=int, action='store', default=0, \
                      help='shotfile edition number')

# .............. output options ...............................................
  parser.add_argument('-p', '--plot', dest='plot_profiles', action="store_true",
                      help="plot profiles on screen")
  parser.add_argument('-np', '--netcdf_path', dest='ncpath', action='store', \
                      default='~/public/Modeling_kit/AUG/prof', \
                      help='netCDF output file path, None for no netcdf output')

  parser.add_argument('-v', '--verbose', dest='verbose', action="store_true",
                      help="print progress messages (not only errors)")

  args = parser.parse_args()

# inject parameters for testing in Spyder etc.
  if True:
      args.shot = 34548
      args.experiment = 'WLS'
      args.diag = 'PED'
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

  # open MDSplus connection
  conn = mdsplus.mdsaug(args.server, althostspec=args.altserver)

  pd = {}  # start with empty profile dictionary

# ................. PED .......................................................
  if args.diag=='PED':
    signames = {'TeFit': 'TeFit', 'TeData': 'TeData', \
                'neFit': 'neFit', 'neData': 'neData', \
                'TiFit': 'TiFit', 'TiData': 'TiData', \
                'vtFit': 'vTFit', 'vtData': 'vTData', \
                'vpFit': 'vPFit', 'vpData': 'vPData' }
    for k in signames.keys():
        p = rprof_augmds.from_shotfile(conn, args.shot, args.diag, signames[k],  \
                         args.tbegin, args.tend, \
                         experiment=args.experiment, edition=args.edition, \
                         verbose=args.verbose)
        if p is not None:
            oshot = p.shot; oedition = p.edition
            pd[k] = p

# ................. IDA .......................................................
  elif args.diag=='IDA':
    signames = {'TeFit': 'Te', \
                'neFit': 'ne' }
    for k in signames.keys():
        p = rprof_augmds.from_shotfile(conn, args.shot, args.diag, signames[k], \
                         args.tbegin, args.tend, \
                         experiment=args.experiment, edition=args.edition, \
                         verbose=args.verbose)
        if p is not None:
            oshot = p.shot; oedition = p.edition
            pd[k] = p

# ............... unkown diagnostic ...........................................
  else:
    print('Unsupported diagnostic: "%s", Abort.' % args.diag)
    sys.exit()

  if pd=={}:
    print('No profiles collected. Nothing to plot or write. Abort.')
    sys.exit()

# ................ plotting ...................................................
  if args.plot_profiles:
    plot_profiles(pd)

# ................ netCDF file ................................................
  if not args.ncpath in ['None', 'none', 'NONE']:
    ncpath = os.path.expanduser(args.ncpath)
    time = 0.5*(args.tbegin+args.tend)
    ofn = '%s/%05d.%05d_%s.nc' \
          % (ncpath, args.shot, 1000*time, args.diag)
    username = getpass.getuser()
    if pd != {}:
      rprof.write_netcdf(ofn, pd, \
          comment='created by augmdsprof2netcdf ('+username+')', \
          title='AUG %s t=%6.3f s %5s:%3s (%d)' \
          % (oshot, time, args.experiment, args.diag, oedition), \
          version=1)

# -----------------------------------------------------------------------------