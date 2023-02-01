#!/usr/bin/env python3

def getHeadDataVmecNc(infile):
  import numpy as np
  import scipy.io.netcdf as ncdf

  if infile.endswith('nc'):
    f = ncdf.netcdf_file(infile)
    data = f.variables

    enfp =  np.copy(data['nfp'].data)
    a = np.copy(data['Aminor_p'].data)
    R0 = np.copy(data['Rmajor_p'].data)
    phipf = np.copy(data['phipf'].data)
    m = np.array(np.copy(data['xm'].data),int)
    n = np.array(-np.copy(data['xn'].data),int)
    psi_tor_a = phipf[0]
    empol = int(np.max(np.abs(m)))
    entor = int(np.max(np.abs(n)))
    m0b = 2*empol
    n0b = 2*entor

    del data
    f.close()

  else:
    enfp =  0
    a = 0.0
    R0 = 0.0
    psi_tor_a = 0.0
    m0b = 0
    n0b = 0

  return [enfp, psi_tor_a, a, R0, m0b, n0b]

if __name__ == "__main__":
  import sys

  infile = sys.argv[1]

  [enfp, psi_tor_a, a, R0, m0b, n0b] = getHeadDataVmecNc(infile)

  print('%3d' % enfp)
  print('%10.8f' % psi_tor_a)
  print('%10.8f' % a)
  print('%10.8f' % R0)
  print('%3d' % m0b)
  print('%3d' % n0b)
