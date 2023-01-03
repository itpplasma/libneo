#!/usr/bin/env python3
# -*- coding: utf-8 -*-

class eqdsk_file:

  def __init__ (self, filename:str):
    self.generate_from_file(filename)

  def generate_from_file(self, filename:str):
    """
    \note Might not be suitable for general format of eqdsk file.
    """
    from numpy import array

    with open(filename) as f:
      l = f.readlines()

    l = [ls.strip() for ls in l]

    self.filenamedateEQDSK = l[0][0:48]
    # l[0].split()[1] is dummy, not used
    # number r coordinates for the grid
    self.nrgr = int(l[0].split()[2]);
    # number z coordinates for the grid
    self.nzgr = int(l[0].split()[3]);
    # Read boxlength
    self.rboxlength = float(l[1].split()[0])
    self.zboxlength = float(l[1].split()[1])
    self.R0 = float(l[1].split()[2])
    self.rboxleft = float(l[1].split()[3])
    self.zbox_mid = float(l[1].split()[4])
    self.R = [self.rboxleft + self.rboxlength/(self.nrgr-1)*k for k in range(self.nrgr)]
    self.Z = [self.zbox_mid - self.zboxlength/2 + self.zboxlength/(self.nzgr-1)*k for k in range(self.nzgr)] # Needs to be transposed.
    # Read magnetic axis
    self.Rpsi0 = float(l[2].split()[0])
    self.Zpsi0 = float(l[2].split()[1])
    # Read psiaxis and psiedge
    self.PsiaxisVs = float(l[2].split()[2])
    self.PsiedgeVs = float(l[2].split()[3])
    # Read B0 and current
    self.Btor_at_R0=float(l[2].split()[4])
    self.Ip=float(l[3].split()[0])
    # Rest of line 3 (4 one-based) and the next line is ignored.

    startline = 5 # 0-based index
    # Read T(psi) - f (size nrgr)
    self.psieqprof = [1.0/(self.nrgr-1)*k for k in range(self.nrgr)]
    self.fprof = []
    line_number = startline
    while (len(self.fprof) < self.nrgr):
      self.fprof += [float(n) for n in l[line_number].split()]
      line_number += 1

    # Read p(psi) (size nrgr)
    self.ptotprof = []
    if self.nrgr % 5 == 0:
      line_number += 1
    while (len(self.ptotprof) < self.nrgr):
      self.ptotprof += [float(n) for n in l[line_number].split()]
      line_number += 1

    # Read TT'(psi) - ff' (size nrgr)
    self.fdfdpsiprof = []
    if self.nrgr % 5 == 0:
      line_number += 1
    while (len(self.fdfdpsiprof) < self.nrgr):
      self.fdfdpsiprof += [float(n) for n in l[line_number].split()]
      line_number += 1

    # Read p'  (size nrgr)
    self.dpressdpsiprof = []
    if self.nrgr % 5 == 0:
      line_number += 1
    while (len(self.dpressdpsiprof) < self.nrgr):
      self.dpressdpsiprof += [float(n) for n in l[line_number].split()]
      line_number += 1

    # Read psi (size nrgr x nzgr in file?, nzgr x nrgr in 'output')
    self.PsiVs = []
    if self.nrgr % 5 == 0:
      line_number += 1
    while (len(self.PsiVs) < self.nrgr*self.nzgr):
      self.PsiVs += [float(n) for n in l[line_number].split()]
      line_number += 1
    self.psinorm = [(PsiV-self.PsiaxisVs+self.PsiedgeVs)/(-self.PsiaxisVs+self.PsiedgeVs) for PsiV in self.PsiVs]
    self.PsiVs = array(self.PsiVs)
    self.PsiVs = self.PsiVs.reshape(self.nzgr,self.nrgr)

    # Read Q profile
    self.qprof = []
    if (self.nrgr*self.nzgr) % 5 == 0:
      line_number += 1
    while (len(self.qprof) < self.nrgr):
      self.qprof += [float(n) for n in l[line_number].split()]
      line_number += 1

    # Read npbound and nblimiter (in one line)
    if self.nrgr % 5 == 0:
      line_number += 1
    self.npbound = int(l[line_number].split()[0])
    self.nplimiter = int(l[line_number].split()[1])
    # Read boundary (LCFS)
    line_number += 1
    self.Lcfs = []
    while (len(self.Lcfs) < 2*self.npbound):
      self.Lcfs += [float(n) for n in l[line_number].split()]
      line_number += 1
    self.Lcfs = array(self.Lcfs)
    self.Lcfs = self.Lcfs.reshape(self.npbound, 2)
    self.Lcfs = self.Lcfs.transpose()

    # Read limiter (LCFS)
    if (2*self.npbound) % 5 == 0:
      line_number += 1
    self.Limiter = []
    while (len(self.Limiter) < 2*self.nplimiter):
      self.Limiter += [float(n) for n in l[line_number].split()]
      line_number += 1
    self.Limiter = array(self.Limiter)
    self.Limiter = self.Limiter.reshape(self.nplimiter, 2)
    self.Limiter = self.Limiter.transpose()

    if (2*self.nplimiter) % 5 == 0:
      line_number += 1
    self.ncoils = int(l[line_number])
    line_number += 1
    # Read PFcoilData
    self.PFcoilData = []
    while (len(self.PFcoilData) < 5*self.ncoils):
      self.PFcoilData += [float(n) for n in l[line_number].split()]
      line_number += 1
    self.PFcoilData = array(self.PFcoilData)
    self.PFcoilData = self.PFcoilData.reshape(self.ncoils, 5)

    self.Rcond = self.PFcoilData[:,0]
    self.Zcond = self.PFcoilData[:,1]
    self.DRcond = self.PFcoilData[:,2]
    self.DZcond = self.PFcoilData[:,3]
    self.Icond = self.PFcoilData[:,4]

    if (5*self.ncoils) % 5 == 0:
      line_number += 1
    # Read Brn
    self.Brn = []
    while (len(self.Brn) < self.nrgr*self.nzgr):
      self.Brn += [float(n) for n in l[line_number].split()]
      line_number += 1
    self.Brn = array(self.Brn)
    self.Brn = self.Brn.reshape(self.nzgr, self.nrgr)

    if (self.nrgr*self.nzgr) % 5 == 0:
      line_number += 1
    # Read Bzn
    self.Bzn = []
    while (len(self.Bzn) < self.nrgr*self.nzgr):
      self.Bzn += [float(n) for n in l[line_number].split()]
      line_number += 1
    self.Bzn = array(self.Bzn)
    self.Bzn = self.Bzn.reshape(self.nzgr, self.nrgr)

    self.rest = l[line_number:]

  def to_makegrid_input(self, filename):
    pass