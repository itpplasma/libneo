#!/usr/bin/env python

'''
  Read ASDEX Upgrade radial profiles through MDSplus
  and save them in various (ASCII-) formats
'''

import numpy as np
import matplotlib.pyplot as plt
import MDSplus
import rprof

tree = 'augdp'

shot = 34835
tbegin = 3.0
tend = 3.05



# ---------------- read profile --------------------------------------------

def mdsvalue(conn, cmd, verbose=0):
  if verbose>1:
    print('mdsvalue: %s' % cmd)  
  try:
    data = np.asarray( conn.get(cmd) )
    return data
  except:
    return None


def read_2dsignal_augmds(conn, diag, signame, tbegin, tend, \
                         experiment='AUGD', edition=0, verbose=0):
  cmd = 'augdiag($SHOT,"'+diag+'","'+signame+'","'+experiment+'",'+str(edition)+',' \
         +str(tbegin)+','+str(tend)+')'
  data = mdsvalue(conn, "(_data="+cmd+")", verbose=verbose)
  if data is None:
      print("Signal %s not found." % signame)
      return None, None, None
  if verbose>0:
    print("Signal: %s" % signame)
  tbase = mdsvalue(conn, "dim_of(_data,0)", verbose=verbose)
  if tbase is None:
      print("No timebase for signal %s." % signame)
      return None, None, None
  if verbose>0:
    print("  time: %d sample(s) (t = %6.3f - %6.3f s)" \
          % (len(tbase), np.min(tbase), np.max(tbase)))
  rbase = mdsvalue(conn, "dim_of(_data,1)", verbose=verbose)
  if rbase is None:
      print("No radial base for signal %s." % signame)
      return None, None, None
  if verbose>0:
    print("  rhop: %d profile points (rho-pol = %6.3f - %6.3f)" \
          % (len(rbase), np.min(rbase), np.max(rbase)))
  return data, tbase, rbase


class rprof_augmds(rprof.RadialProfile):

  @classmethod
  def from_ped(cls, conn, signame, tbegin, tend, experiment='AUGD', diag='PED', edition=0, verbose=0):
    data2, tbase, rhop = read_2dsignal_augmds(conn, diag, signame, tbegin, tend, \
                                experiment=experiment, edition=edition, verbose=verbose)
    if data2 is None or tbase is None or rhop is None:
      return None

    esigname = "d"+signame
    errbar2, tbase_, rhop_ = read_2dsignal_augmds(conn, diag, esigname, tbegin, tend, \
                                experiment=experiment, edition=edition, verbose=verbose)
    if errbar2 is None:
      return None

    if np.any(tbase_.shape != tbase.shape) and np.any(tbase_ != tbase):
      print("Signal (%s) and Error bars (%s) do not have the same time base." % (signame, esigname))
      return None

    if np.any(rhop_.shape != rhop.shape) and np.any(rhop_ != rhop):
      print("Signal (%s) and Error bars (%s) do not have the same radial base." % (signame, esigname))
      return None    

    ntime = data2.shape[1]
    data1 = np.sum(data2, axis=1) / ntime
    psin = np.square(rhop)
    if errbar2 is None:
      errbar1 = None
    else:
      errbar1 = np.sum(errbar2, axis=1) / ntime
    return cls(signame, data1, errbar1, 'psinorm', psin)
  

  @classmethod
  def from_ida(cls, conn, signame, tbegin, tend, experiment='AUGD', diag='IDA', edition=0):
    data = read_signal_augmds(conn, diag, signame, tbegin, tend, experiment=experiment, edition=edition)
    if data is None:
      return None
    psin = mdsvalue(conn, "dim_of(_data,1)^2") 
    if psin is None:
      return None
    errbar = read_signal_augmds(conn, diag, signame+'_unc', tbegin, tend, experiment=experiment, edition=edition)
    ntime = data.shape[1]
    print("Shapes: ", data.shape, psin.shape)
    data1 = np.sum(data, axis=1) / ntime
    psin1 = np.sum(psin, axis=1) / ntime
    if errbar is None:
      errbar1 = None
    else:
      errbar1 = np.sum(errbar, axis=1) / ntime
    return cls(signame, data1, errbar1, 'psinorm', psin1)


# ------------------- plot -----------------------------------------------------



# ------------------- main -----------------------------------------------------

if True:
  conn = MDSplus.Connection('mdsplus')
  conn.openTree(tree, 34548)
  TeFit  = rprof_augmds.from_ped(conn, 'TeFit',  3.0, 5.9, \
                                 experiment='WLS', edition=1, verbose=1)
  

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

  args = parser.parse_args()

  if args.shot is None:
    print('No shot number given! Abort.')
    sys.exit()

  if args.experiment=='self' or args.experiment=='SELF':
    args.experiment=getpass.getuser()

  print('Shot: %d  Diagnostic: %s  Experiment: %s  Edition: %d' \
          % (args.shot, args.diag, args.experiment, args.edition))
  print('Requested time interval: %6.3f - %6.3f s' % (args.tbegin, args.tend))
    
  print('Connecting to MDSplus server: %s' % args.server)
  conn = MDSplus.Connection(args.server)

  print('Open shot: %s' % args.shot)
  conn.openTree(tree, args.shot)

  if args.diag=='PED':
    TeFit  = rprof_augmds.from_ped(conn, 'TeFit',  args.tbegin, args.tend, \
                                 experiment=args.experiment, edition=args.edition)
    TeData = rprof_augmds.from_ped(conn, 'TeData', args.tbegin, args.tend, \
                                 experiment=args.experiment, edition=args.edition)
  elif args.diag=='IDA':
    TeFit  = rprof_augmds.from_ida(conn, 'Te',     args.tbegin, args.tend, \
                                 experiment=args.experiment, edition=args.edition)
    TeData = None
    neFit  = rprof_augmds.from_ida(conn, 'ne',     args.tbegin, args.tend, \
                                 experiment=args.experiment, edition=args.edition)
  else:
    print('Unsupported diagnostic: %s -- Abort.' % argv.diag)
    sys.exit()

  if args.plot_profiles:
    TeData.plot_errors(1e3, color='grey', linewidth=0)
    TeFit.plot(1e3, color='red')
    plt.show()

